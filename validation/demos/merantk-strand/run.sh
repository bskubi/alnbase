#!/usr/bin/env bash
# Run meRanGh and meRanGs, meRanTK's two genome aligners, on the four-strand fixture the other
# strand demos use, and check what queries/strand/merantk.toml makes of the BAMs they write.
#
# meRanTK is the only RNA tool in this series: it maps bisulfite-treated RNA for 5mC in mRNA,
# so it is a spliced aligner (hisat2 in meRanGh, STAR in meRanGs) rather than a DNA one. The
# fixture has no introns, which is deliberate -- what is under test is the tag, and the tag is
# written the same way either way.
#
# What this demo measures:
#
#   - Whether the rule's four-way table is right. YG names the converted genome index that won
#     and YR the conversion applied to the read, the same split as Bismark's XG and XR, so the
#     two together name the strand of origin per record. The run puts a record in all four
#     cells and checks each one.
#
#   - What the rule's comment got wrong. It said "in practice only OT and OB occur". That is
#     true of one mode out of three. In paired-end, YR is the mate number -- C2T on read 1 and
#     G2A on read 2, unconditionally (meRanGh.pl:3029-3060, meRanGs.pl:3056-3087) -- so half of
#     every run's records are CTOT or CTOB. In single-end, YR is a run-level constant set once
#     before the read loop (meRanGh.pl:1432-1433, meRanGs.pl:1483-1484): C2T for a -f run and
#     G2A for a -r one, so a -r run is CTOT and CTOB throughout and a -f run OT and OB. YR
#     never varies with the read in any mode.
#
#   - A bug in meRanGs that the rule's comment guessed at and this run confirms. meRanGh passes
#     the read direction into its single-end SEQ registration (meRanGh.pl:1988-2018) and stores
#     SEQ reference-forward in every mode. meRanGs' registration takes no such argument
#     (meRanGs.pl:2035-2064): it stores the read as sequenced on a C2T hit and reverse-
#     complemented on a G2A hit, which is right for a -f run and backwards for both halves of a
#     -r one. The records still carry NM:i:0 and MD:Z:60, computed against the converted genome,
#     so nothing in the BAM contradicts itself. alnbase's count of read bases differing from the
#     reference is what finds it.
#
#   MERANTK_SRC=/path/to/meRanTK/checkout MERANTK_ENV=/path/to/env \
#   PY=/path/to/python-with-pyarrow ALNBASE=/path/to/alnbase ./run.sh [workdir]
#
# MERANTK_SRC is a clone of https://github.com/icbi-lab/meRanTK (run against c06d726, v1.3.0);
# its scripts are in src/. MERANTK_ENV is the environment from ../../envs/merantk.yaml. samtools
# comes from PATH.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${MERANTK_SRC:?set MERANTK_SRC}" "${MERANTK_ENV:?set MERANTK_ENV}"
: "${PY:?set PY}" "${ALNBASE:?set ALNBASE}"
# perl-bio-samtools pins samtools 0.1.19, which is in the environment and whose `sort` takes a
# prefix rather than -o. Resolve the samtools this demo uses before the environment shadows it.
SAMTOOLS="${SAMTOOLS:-$(command -v samtools)}"
export PATH="$MERANTK_ENV/bin:$PATH"
PERL="$MERANTK_ENV/bin/perl"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }

# The reference is 1200 bases from a fixed recurrence with a CG planted every 17, so every
# read carries several cytosines in CG context and the four loci stay unique after either
# conversion. Four 200-base fragments, one per strand of origin; each pair is named for the
# strand its READ 1 was sequenced from. Identical to ../bsmap-strand/, ../bwameth-strand/,
# ../hisat-3n-strand/ and ../dnmtools-strand/.
say "The four pairs, each named for the strand its read 1 came from"
"$PY" - <<'PYEOF'
bases = "ACGT"
ref, x = [], 7
for _ in range(1200):
    x = (x * 31 + 17) % 1009
    ref.append(bases[x % 4])
for i in range(0, 1198, 17):
    ref[i], ref[i + 1] = "C", "G"
ref = "".join(ref)

with open("ref.fa", "w") as f:
    f.write(">chr1\n")
    for i in range(0, len(ref), 60):
        f.write(ref[i:i + 60] + "\n")

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s):
    return "".join(COMP[c] for c in reversed(s))

def bisulfite(s):
    """Unmethylated cytosines become T; a cytosine in CG stays C."""
    return "".join("T" if c == "C" and s[i + 1:i + 2] != "G" else c
                   for i, c in enumerate(s))

FRAG, LEN = 200, 60
MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}
pairs = [("OT", 100), ("CTOT", 350), ("OB", 600), ("CTOB", 850)]
f1, f2 = open("all_1.fq", "w"), open("all_2.fq", "w")
for name, start in pairs:
    frag = ref[start:start + FRAG]
    ot, ob = bisulfite(frag), bisulfite(rc(frag))
    strand = {"OT": ot, "CTOT": rc(ot), "OB": ob, "CTOB": rc(ob)}
    # Each read starts at its own strand's 5' end.
    r1, r2 = strand[name][:LEN], strand[MATE[name]][:LEN]
    f1.write("@%s\n%s\n+\n%s\n" % (name, r1, "I" * LEN))
    f2.write("@%s\n%s\n+\n%s\n" % (name, r2, "I" * LEN))
    print("%-5s fragment chr1 1-based %d-%d: read 1 from %s, read 2 from %s"
          % (name, start + 1, start + FRAG, name, MATE[name]))
f1.close()
f2.close()
PYEOF
"$SAMTOOLS" faidx ref.fa

# Both tools build a pair of converted indexes, C2T and G2A, and align to both. -f is the
# fastq of reads from the 5' end of the RNA and -r the fastq from the 3' end; giving both is a
# paired-end run, giving one alone is a single-end run in that direction. The three modes are
# run for each tool because YR says something different in each.
say "meRanGh (hisat2): index, then paired-end, -f alone and -r alone"
"$PERL" "$MERANTK_SRC/src/meRanGh.pl" mkbsidx -fa ref.fa -id ghidx -t 2 > ghidx.log 2>&1
gh() { "$PERL" "$MERANTK_SRC/src/meRanGh.pl" align -id ghidx -t 2 -ob "$@"; }
gh -f all_1.fq -r all_2.fq -o gh_pe -S gh_pe.sam > gh_pe.log 2>&1
gh -f all_1.fq              -o gh_f  -S gh_f.sam  > gh_f.log  2>&1
gh              -r all_2.fq -o gh_r  -S gh_r.sam  > gh_r.log  2>&1

say "meRanGs (STAR): the same three"
# The fixture's reference is 1200 bases, far below STAR's default suffix-array sizing.
"$PERL" "$MERANTK_SRC/src/meRanGs.pl" mkbsidx -fa ref.fa -id gsidx -t 2 \
  -star_genomeSAindexNbases 5 > gsidx.log 2>&1
gs() { "$PERL" "$MERANTK_SRC/src/meRanGs.pl" align -id gsidx -t 2 -ob "$@"; }
gs -f all_1.fq -r all_2.fq -o gs_pe -S gs_pe.sam > gs_pe.log 2>&1
gs -f all_1.fq              -o gs_f  -S gs_f.sam  > gs_f.log  2>&1
gs              -r all_2.fq -o gs_r  -S gs_r.sam  > gs_r.log  2>&1

RUNS="gh_pe gh_f gh_r gs_pe gs_f gs_r"

say "What each run wrote: the tags, the call the rule makes, and how SEQ is stored"
"$PY" - <<'PYEOF'
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
def rc(s):
    return "".join(COMP[c] for c in reversed(s))

def fastq(path):
    lines = open(path).read().split("\n")
    return dict((lines[i][1:].split()[0], lines[i + 1]) for i in range(0, len(lines) - 1, 4))
READS = (fastq("all_1.fq"), fastq("all_2.fq"))

# queries/strand/merantk.toml, written out so the run compares rather than asserts.
RULE = {("C2T", "C2T"): "OT", ("C2T", "G2A"): "CTOT",
        ("G2A", "C2T"): "OB",  ("G2A", "G2A"): "CTOB"}

# The strand each record was really sequenced from. Read 1 of a pair is the strand the pair is
# named for; read 2 is that strand's complement. A -f run is read 1 only, a -r run read 2 only.
MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}

for run in "gh_pe gh_f gh_r gs_pe gs_f gs_r".split():
    print("\n-- %s" % run)
    print("   pair  flag  pos   YG    YR    rule call  truth  SEQ stored")
    wrong_call = wrong_seq = 0
    for line in open("%s/%s.sam" % (run, run)):
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        pair, flag, pos, seq = f[0].split("/")[0], int(f[1]), int(f[3]), f[9]
        tags = dict(t.split(":Z:") for t in f[11:] if ":Z:" in t)
        yg, yr = tags.get("YG", "-"), tags.get("YR", "-")
        call = RULE.get((yg, yr), "unknown")
        # A paired run says which mate by the flag; a single-end run by which fastq it was given.
        mate = (0 if flag & 0x40 else 1) if flag & 0x1 else (1 if run.endswith("_r") else 0)
        truth = pair if mate == 0 else MATE[pair]
        asis = READS[mate][pair]
        stored = ("as sequenced" if seq == asis else
                  "reference-forward" if seq == rc(asis) else "neither")
        want = "reference-forward" if flag & 0x10 else "as sequenced"
        wrong_call += call != truth
        wrong_seq += stored != want
        print("   %-5s %-5d %-5d %-5s %-5s %-10s %-6s %-18s %s"
              % (pair, flag, pos, yg, yr, call, truth, stored,
                 "" if stored == want else "*** wrong way round"))
    print("   strand of origin wrong: %d;  SEQ stored the wrong way round: %d"
          % (wrong_call, wrong_seq))
PYEOF

say "alnbase on each, with queries/strand/merantk.toml"
"$ALNBASE" index ref.fa ref.aref > /dev/null
for run in $RUNS; do
  "$SAMTOOLS" sort -o "$run.bam" "$run/$run.sam" 2> /dev/null
  printf -- '-- %s\n' "$run"
  "$ALNBASE" query --parquet --only-hits \
    --query-file "$HERE/cg.toml" --query-file "$HERE/../../../queries/strand/merantk.toml" \
    -f qname,flags,strand,conv_strand,read_reverse \
    "$run.bam" ref.aref "$run.hits.parquet" 2>&1 \
    | grep -E 'records|differ' || true
done

# The two tools are meant to be interchangeable, so the hits are compared against each other
# as well as counted: same fixture, same tag, same rule, different aligner underneath.
say "CG hits: meRanGh against meRanGs, mode by mode"
"$PY" - <<'PYEOF'
import pyarrow.parquet as pq

def hits(run):
    t = pq.read_table("%s.hits_0_0.parquet" % run).to_pylist()
    return set((r["qname"].split("/")[0], r["refr_pos"], r["strand"]) for r in t), len(t)

for mode in ("pe", "f", "r"):
    gh, ngh = hits("gh_" + mode)
    gs, ngs = hits("gs_" + mode)
    print("   %-3s  meRanGh %2d hits, meRanGs %2d hits, %d in common, "
          "%d only in meRanGh, %d only in meRanGs"
          % (mode, ngh, ngs, len(gh & gs), len(gh - gs), len(gs - gh)))
PYEOF

say "meRanGs -r, the broken run, in full"
grep -v '^@' gs_r/gs_r.sam | cut -f1-9,12-
echo "   NM and MD are computed against the converted genome, so both records claim a perfect"
echo "   match. Only comparing SEQ with the real reference finds the problem."
