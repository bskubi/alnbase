#!/usr/bin/env bash
# Run abismal, dnmtools' own aligner, on the four-strand fixture the other strand demos use,
# and check what queries/strand/dnmtools.toml makes of the BAM it writes.
#
# This demo found the shipped rule wrong, and the round rewrote it. What it measures:
#
#   - Where CV comes from. The rule was written believing dnmtools translates CV in from
#     BSMAP's ZS or Bismark's XR. It does, but only inside `dnmtools format`
#     (bam_record_utils.cpp:833-913, called only from format-reads.cpp:339,348). Every
#     unformatted BAM carrying CV was therefore written by abismal, which assigns it
#     directly (abismal.cpp:480, :660, :688).
#
#   - What CV says. It is the conversion the READ shows as sequenced -- T-rich or A-rich --
#     not the reference strand that carried it. In the default directional mode it is the
#     mate number and nothing else: read 1 is always T, read 2 always A (abismal's own
#     README.md:112-117 says so). Only under -R, random PBAT, does it vary with the read.
#
#   - What the shipped rule made of that: CV:A:T -> '+' and CV:A:A -> '-', which is right on
#     a forward record and backwards on a reverse one. The conversion strand is the XOR of
#     the two, which is how abismal computes it itself (abismal.cpp:1269-1270). The rule
#     file now says that, and this run is the check.
#
#   - The larger problem, which no rule can fix. abismal stores SEQ as it came off the
#     sequencer on every record, reverse ones included (README.md:108-109), against the SAM
#     convention and against the invariant alnbase's walk relies on. The demo measures the
#     damage and then repairs the BAM with a five-line script, because the corrected rule is
#     only usable on a repaired file.
#
#   DNMTOOLS_ENV=/path/to/env PY=/path/to/python-with-pyarrow ALNBASE=/path/to/alnbase \
#     ./run.sh [workdir]
#
# DNMTOOLS_ENV is the environment from ../../envs/dnmtools.yaml. It provides abismal,
# dnmtools and samtools.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${DNMTOOLS_ENV:?set DNMTOOLS_ENV}" "${PY:?set PY}" "${ALNBASE:?set ALNBASE}"
export PATH="$DNMTOOLS_ENV/bin:$PATH"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }

# The reference is 1200 bases from a fixed recurrence with a CG planted every 17, so every
# read carries several cytosines in CG context and the four loci stay unique after either
# conversion. Four 200-base fragments, one per strand of origin; each pair is named for the
# strand its READ 1 was sequenced from. Identical to ../bsmap-strand/, ../bwameth-strand/
# and ../hisat-3n-strand/.
say "The four pairs, each named for the strand its read 1 came from"
python3 - <<'PYEOF'
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
samtools faidx ref.fa

# Default mode is directional: abismal maps read 1 T-rich and read 2 A-rich and nothing
# else. -R is its non-directional mode, where each end is tried both ways and CV records
# which way won. Both runs are kept, because the two say different things about CV.
say "abismal, paired-end: default (directional) and -R (random PBAT)"
abismal idx ref.fa ref.idx > idx.log 2>&1
abismal map -i ref.idx -o dir.sam  all_1.fq all_2.fq 2> dir.log
abismal map -i ref.idx -R -o rand.sam all_1.fq all_2.fq 2> rand.log
grep -c '^[^@]' dir.sam rand.sam || true

say "What abismal wrote, and what each rule makes of it"
"$PY" - <<'PYEOF'
ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s):
    return "".join(COMP[c] for c in reversed(s))

def fastq(path):
    lines = open(path).read().split("\n")
    return dict((lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 4))
READS = (fastq("all_1.fq"), fastq("all_2.fq"))

# Which strand each record was really sequenced from, and that strand's conversion strand:
# OT and CTOT are the converted top strand and its copy, so the informative reference base
# is a C and the walk runs '+'; OB and CTOB give '-'.
MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}
CONV = {"OT": "+", "CTOT": "+", "OB": "-", "CTOB": "-"}

def shipped(cv, rev):
    """queries/strand/dnmtools.toml as it stood before this round."""
    return {"T": "+", "A": "-"}.get(cv, "?")

def corrected(cv, rev):
    """The rule this round wrote: the conversion strand is CV XOR 0x10."""
    if cv not in ("T", "A"):
        return "?"
    return "+" if ((cv == "T") != rev) else "-"

def report(path, title):
    print("\n-- %s" % title)
    print("%-5s %-6s %-5s %-5s %-4s %-3s %-6s %-9s %-9s %s"
          % ("pair", "record", "from", "flag", "pos", "CV", "want", "shipped", "corrected",
             "SEQ stored"))
    bad_ship = bad_corr = bad_seq = 0
    rows = []
    for line in open(path):
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        qname, flag, pos, seq = f[0], int(f[1]), int(f[3]), f[9]
        tags = dict((t[:2], t[5:]) for t in f[11:])
        mate = 1 if flag & 0x40 else 2
        origin = qname if mate == 1 else MATE[qname]
        rev = bool(flag & 0x10)
        cv = tags.get("CV", "?")
        want = CONV[origin]
        s, c = shipped(cv, rev), corrected(cv, rev)
        bad_ship += s != want
        bad_corr += c != want
        # As sequenced, or reverse-complemented into the reference's orientation? The read
        # is in the FASTQ under its own name; compare against both.
        asis = READS[mate - 1][qname]
        stored = ("as sequenced" if seq == asis else
                  "reference-forward" if seq == rc(asis) else "neither")
        bad_seq += rev and stored == "as sequenced"
        rows.append((qname, mate, origin, flag, pos, cv, want, s, c, stored))
    for r in sorted(rows):
        print("%-5s read %d %-5s %-5d %-4d %-3s %-6s %-9s %-9s %s"
              % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9]))
    n = len(rows)
    print("   conversion strand wrong, shipped rule (CV alone):   %d of %d" % (bad_ship, n))
    print("   conversion strand wrong, corrected rule (CV ^ 0x10): %d of %d" % (bad_corr, n))
    print("   reverse records whose SEQ is NOT reference-forward:  %d of %d" % (bad_seq, n))
    return rows

report("dir.sam", "default mode: CV is the mate number and nothing else")
report("rand.sam", "-R: CV varies with the read, and the shipped rule halves")
PYEOF

# abismal's SEQ convention is documented and deliberate -- its own caller reads SEQ backwards
# to match -- but it is not SAM's, and every tool that assumes SAM's, alnbase included, will
# read those records wrong. Repairing it is mechanical, so the demo repairs it and carries on.
say "Repairing SEQ: reverse-complement the reverse records, as SAM requires"
samtools view -h rand.sam | "$PY" -c '
import sys
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
for line in sys.stdin:
    if not line.startswith("@"):
        f = line.rstrip("\n").split("\t")
        if int(f[1]) & 0x10:
            f[9] = "".join(COMP[c] for c in reversed(f[9]))
            f[10] = f[10][::-1]
            line = "\t".join(f) + "\n"
    sys.stdout.write(line)
' | samtools sort -o fixed.bam -
samtools index fixed.bam
samtools sort -o raw.bam rand.sam
samtools index raw.bam

say "alnbase on the raw BAM and on the repaired one, with each rule"
"$ALNBASE" index ref.fa ref.aref > /dev/null
cat > shipped.toml <<'EOF'
# queries/strand/dnmtools.toml as it stood before this round, kept here so the run can show
# what it did. Do not copy it; the file in queries/ is the corrected one.
[strand.conversion]
"+" = 'CV == "T"'
"-" = 'CV == "A"'
[strand.sequenced]
forward = "not is_reverse"
reverse = "is_reverse"
EOF
for bam in raw fixed; do
  for rule in shipped "$HERE/../../../queries/strand/dnmtools.toml"; do
    label="$([ "$rule" = shipped ] && echo shipped || echo corrected)"
    file="$([ "$rule" = shipped ] && echo shipped.toml || echo "$rule")"
    printf -- '-- %s BAM, %s rule\n' "$bam" "$label"
    "$ALNBASE" query --parquet --only-hits \
      --query-file "$HERE/cg.toml" --query-file "$file" \
      -f qname,flags,strand,conv_strand,read_reverse \
      "$bam.bam" ref.aref "$bam.$label.hits.parquet" 2>&1 \
      | grep -E 'records|differ' || true
  done
done

say "Where the hits land"
"$PY" - <<'PYEOF'
import pyarrow.parquet as pq

CONV = {"OT": "+", "CTOT": "+", "OB": "-", "CTOB": "-"}
MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}

def hits(bam, label):
    return pq.read_table("%s.%s.hits_0_0.parquet" % (bam, label)).to_pylist()

# The repaired BAM read with the corrected rule is the answer the other three are measured
# against: a hit is the same hit only if it is the same record, the same reference position
# and the same walk direction.
truth = set((r["qname"], r["flags"], r["refr_pos"], r["strand"])
            for r in hits("fixed", "corrected"))

for bam in ("raw", "fixed"):
    for label in ("shipped", "corrected"):
        t = hits(bam, label)
        off = 0
        for r in t:
            mate = 1 if r["flags"] & 0x40 else 2
            origin = r["qname"] if mate == 1 else MATE[r["qname"]]
            off += r["strand"] != CONV[origin]
        agree = sum(1 for r in t
                    if (r["qname"], r["flags"], r["refr_pos"], r["strand"]) in truth)
        print("%-6s BAM, %-9s rule: %3d CG hits over %d records, %3d walked along the "
              "wrong conversion strand, %2d of the %d right ones found"
              % (bam, label, len(t), len(set((r["qname"], r["flags"]) for r in t)), off,
                 agree, len(truth)))
PYEOF

# `dnmtools format` is the next step of the real pipeline, and it does more than normalise a
# tag: it merges the two mates of a fragment into ONE record whose CIGAR bridges the gap with
# an N, drops the paired bit, and writes CV:A:T on everything. No strand rule applies to that.
say "dnmtools format, for the record"
samtools sort -n -o rand.ns.bam rand.sam
dnmtools format -f abismal -B rand.ns.bam fmt.bam 2> fmt.log
samtools view fmt.bam | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$12"\t"$13}'
