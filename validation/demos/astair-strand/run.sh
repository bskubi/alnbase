#!/usr/bin/env bash
# Run asTair's caller on a TAPS fixture and check what queries/strand/astair.toml makes of the
# same BAM.
#
# asTair is the only non-bisulfite tool in this series. TAPS reads a MODIFIED cytosine as T and
# leaves an unmodified one as C -- the opposite way round from bisulfite -- and, because the
# chemistry does not convert the whole library, the reference is not converted either. So the
# aligner is an ordinary bwa mem rather than a three-letter one, and nothing at alignment time
# discards copy-strand reads. Every other strand demo here has had that filter for free.
#
# asTair also writes no conversion tag and reads none. Its caller takes the strand from a
# literal table of FLAG integers (caller.py:235-267), which is what astair.toml transcribes.
# What this demo measures:
#
#   - Whether that table is right, by mirroring it. taps.toml asks two questions at one anchor,
#     "a CG cytosine read as T" and "read as C", so alnbase's two hit counts are comparable
#     with asTair's MOD and UNMOD columns position by position rather than merely by label.
#
#   - What the rule's comment got wrong. It claimed asTair counts a record outside the six
#     FLAGs toward depth. It does not: pysam's pileup drops non-proper pairs before asTair sees
#     them, which the run shows by aligning the same fixture with and without a deliberate FF
#     pair and diffing the two .mods files.
#
#   - A bug in asTair. `astair call --ignore_orphans False` is a documented option
#     (caller.py:61) that cannot work: the branch it turns on calls list.extend() with two
#     arguments (caller.py:247-248, 259-260), which raises TypeError. Because clean_pileup
#     wraps its per-position body in `except Exception: continue` (caller.py:339), every
#     position is silently skipped, the run logs success and exits 0, and the user is handed an
#     empty .mods file.
#
#   ASTAIR_ENV=/path/to/env PY=/path/to/python-with-pyarrow ALNBASE=/path/to/alnbase \
#     ./run.sh [workdir]
#
# ASTAIR_ENV is the environment from ../../envs/astair.yaml, with asTair installed into it from
# a checkout; see that file. bwa and samtools come from the same environment.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${ASTAIR_ENV:?set ASTAIR_ENV}" "${PY:?set PY}" "${ALNBASE:?set ALNBASE}"
export PATH="$ASTAIR_ENV/bin:$PATH"
ASTAIR="$ASTAIR_ENV/bin/astair"
ASTAIR_PY="$ASTAIR_ENV/bin/python"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }

# The reference is the same 1200 bases every strand demo uses: a fixed recurrence with a CG
# planted every 17, so every read carries several cytosines in CG context and the four loci stay
# unique. What differs is the chemistry applied to the fragments.
#
# Every other planted CG is modified, alternating. Modifying all of them would leave no exact
# stretch longer than 16 bases, below bwa mem's default -k 19 seed, and nothing would map; the
# alternation also gives the fixture both modified and unmodified sites rather than a uniform
# modification level of 1.0.
say "The fixture: four probe pairs, 200 filler pairs, one FF pair"
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

MODIFIED = set(range(0, 1198, 34))   # every other planted CG, by reference offset

def taps(s, off, flip):
    """TAPS: a modified cytosine reads as T; an unmodified one stays C.

    off is the reference offset of s[0]; flip is True when s runs 3'->5' along the
    reference, as the bottom strand does, so that a base at reference position p is
    looked up in MODIFIED either way."""
    out = []
    for i, c in enumerate(s):
        p = off + (len(s) - 1 - i if flip else i)
        # The bottom strand's cytosine sits on the G of a planted CG, one to the right.
        out.append("T" if c == "C" and ((p in MODIFIED) if not flip else (p - 1 in MODIFIED))
                   else c)
    return "".join(out)

FRAG, LEN = 200, 60
MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}
f1, f2 = open("all_1.fq", "w"), open("all_2.fq", "w")

# The four probes, one pair per strand of origin, each named for the strand its read 1 came
# from. These are what the strand rule is judged on.
for name, start in [("OT", 100), ("CTOT", 350), ("OB", 600), ("CTOB", 850)]:
    frag = ref[start:start + FRAG]
    ot, ob = taps(frag, start, False), taps(rc(frag), start, True)
    strand = {"OT": ot, "CTOT": rc(ot), "OB": ob, "CTOB": rc(ob)}
    r1, r2 = strand[name][:LEN], strand[MATE[name]][:LEN]
    f1.write("@%s\n%s\n+\n%s\n" % (name, r1, "I" * LEN))
    f2.write("@%s\n%s\n+\n%s\n" % (name, r2, "I" * LEN))
    print("   %-5s chr1 1-based %d-%d: read 1 from %s, read 2 from %s"
          % (name, start + 1, start + FRAG, name, MATE[name]))

# Filler. bwa mem infers the insert-size distribution from the batch and will not set 0x2 on a
# handful of pairs, so with the four probes alone every record comes out 97/145/161/81 -- every
# one of them outside asTair's six-integer table, and asTair calls nothing at all. The filler
# exists only to give bwa a distribution to estimate.
n = 0
for start in range(0, 1000, 10):
    for which in ("OT", "OB"):
        frag = ref[start:start + FRAG]
        ot, ob = taps(frag, start, False), taps(rc(frag), start, True)
        r1, r2 = (ot[:LEN], rc(ot)[:LEN]) if which == "OT" else (ob[:LEN], rc(ob)[:LEN])
        n += 1
        f1.write("@fill%d_%s\n%s\n+\n%s\n" % (n, which, r1, "I" * LEN))
        f2.write("@fill%d_%s\n%s\n+\n%s\n" % (n, which, r2, "I" * LEN))
print("   filler pairs: %d" % n)

# One pair that is deliberately not a proper pair: both mates taken from the top strand, so
# they align FF and bwa leaves 0x2 clear. This is the record outside asTair's table, and it
# must be last in the file, because the comparison run below is these same fastqs with their
# final record dropped.
start = 100
ot = taps(ref[start:start + FRAG], start, False)
for f, r in ((f1, ot[:LEN]), (f2, ot[FRAG - LEN:])):
    f.write("@ff\n%s\n+\n%s\n" % (r, "I" * LEN))
f1.close()
f2.close()
print("   one FF pair, both mates from the top strand: not a proper pair")
PYEOF
samtools faidx ref.fa

# TAPS does not convert the reference, so this is a plain bwa mem against a plain reference.
say "bwa mem"
bwa index ref.fa 2> /dev/null
bwa mem -t 2 ref.fa all_1.fq all_2.fq 2> /dev/null | samtools sort -o pe.bam
samtools index pe.bam
# The same alignment without the FF pair, for the depth comparison further down.
head -n -4 all_1.fq > noff_1.fq; head -n -4 all_2.fq > noff_2.fq
bwa mem -t 2 ref.fa noff_1.fq noff_2.fq 2> /dev/null | samtools sort -o noff.bam
samtools index noff.bam

say "Flags: what bwa wrote, and where the four probes landed"
echo "   count flag"
samtools view pe.bam | awk '{print $2}' | sort -n | uniq -c | sed 's/^/  /'
echo "   probe pairs (qname, flag, pos):"
samtools view pe.bam | awk '$1 !~ /^fill/ {printf "   %-6s %-5s %s\n", $1, $2, $4}'

# asTair wants absolute paths: a relative -f is reported as "The genome reference fasta file
# does not exist", which is not what has gone wrong (simple_fasta_parser.py raises IndexError).
say "asTair call: the default run"
mkdir -p out
"$ASTAIR" call -i "$WORK/pe.bam" -f "$WORK/ref.fa" -d "$WORK/out" -co CpG -t 1 2>&1 | tail -3
echo "   -- stats"
sed 's/^/   /' out/pe_mCtoT_CpG.stats

say "The depth claim: the same fixture without the FF pair"
mkdir -p out_noff
"$ASTAIR" call -i "$WORK/noff.bam" -f "$WORK/ref.fa" -d "$WORK/out_noff" -co CpG -t 1 > /dev/null 2>&1
if diff -q <(sed 's/^noff/pe/' out_noff/noff_mCtoT_CpG.mods) out/pe_mCtoT_CpG.mods > /dev/null; then
  echo "   .mods files identical: the FF pair changed no position's MOD, UNMOD or TOTAL_DEPTH."
  echo "   pysam's pileup drops a non-proper pair before asTair sees it, so it never reaches"
  echo "   depth either -- contrary to what astair.toml's comment used to say."
else
  echo "   *** the .mods files differ; see out/ and out_noff/"
  diff <(sed 's/^noff/pe/' out_noff/noff_mCtoT_CpG.mods) out/pe_mCtoT_CpG.mods | head
fi

say "The bug: astair call --ignore_orphans False"
mkdir -p out_io
"$ASTAIR" call -i "$WORK/pe.bam" -f "$WORK/ref.fa" -d "$WORK/out_io" -co CpG -t 1 -io False 2>&1 | tail -2
echo "   exit status: $?"
echo "   -- the mods file it wrote, in full ($(wc -l < out_io/pe_mCtoT_CpG.mods) line)"
sed 's/^/   /' out_io/pe_mCtoT_CpG.mods
echo "   -- and the stats"
sed 's/^/   /' out_io/pe_mCtoT_CpG.stats
echo "   A successful-looking run, exit 0, and nothing in it. What actually happens, called"
echo "   directly so the exception is not swallowed:"
"$ASTAIR_PY" - <<'PYEOF' 2>&1 | sed 's/^/   /'
from astair.caller import flags_expectation
for io in (True, False):
    try:
        d, u = flags_expectation({}, 0, 'mCtoT', 'C', io, False, 'directional')
        print("ignore_orphans=%-5s -> %d desired, %d undesired tuples" % (io, len(d), len(u)))
    except Exception as e:
        print("ignore_orphans=%-5s -> %s: %s" % (io, type(e).__name__, e))
PYEOF
echo "   caller.py:247-248 and 259-260 call list.extend() with two arguments; caller.py:339"
echo "   catches every exception per position and continues, so the failure is silent."

say "alnbase with queries/strand/astair.toml"
"$ALNBASE" index ref.fa ref.aref > /dev/null
"$ALNBASE" query --parquet --only-hits \
  --query-file "$HERE/taps.toml" --query-file "$HERE/../../../queries/strand/astair.toml" \
  -f qname,flags,strand \
  pe.bam ref.aref hits.parquet 2>&1 | sed 's/^/   /'

say "The four probes, as the rule calls them"
"$PY" - <<'PYEOF'
import pyarrow.parquet as pq

# The strand each record was really sequenced from: read 1 of a pair is the strand the pair is
# named for, read 2 is that strand's complement. 0x40 is read 1, so 99 and 83 are read 1 and
# 147 and 163 are read 2.
TRUTH = {("OT", 99): "OT", ("OT", 147): "CTOT", ("CTOT", 83): "CTOT", ("CTOT", 163): "OT",
         ("OB", 83): "OB", ("OB", 163): "CTOB", ("CTOB", 99): "CTOB", ("CTOB", 147): "OB"}

seen = {}
for r in pq.read_table("hits_0_0.parquet").to_pylist():
    if not r["qname"].startswith("fill"):
        seen.setdefault((r["qname"], r["flags"]), r["strand"])

print("   pair  flag  sequenced from  informs about  rule says")
wrong = 0
for (q, f), called in sorted(seen.items()):
    origin = TRUTH[(q, f)]
    # A CTOT read is a PCR copy of an OT read, so it carries the TOP strand's conversion; a
    # CTOB read carries the bottom strand's. What a record informs about is its fragment's
    # original strand, not the strand it was itself sequenced from.
    informs = "+" if origin in ("OT", "CTOT") else "-"
    wrong += called != informs
    print("   %-5s %-5d %-15s %-14s %s%s"
          % (q, f, origin, informs, called, "" if called == informs else "   <- wrong"))
print("   conversion strand wrong: %d of %d" % (wrong, len(seen)))
print()
print("   CTOT and CTOB records are not the problem: read 2 of every ordinary fragment is one,")
print("   and asTair gets those right. 99 and 147 together are 'an FR pair whose read 1 is")
print("   forward', which is a top-strand fragment and both of its mates; 83 and 163 are the")
print("   same for a bottom-strand fragment. What breaks is a pair whose READ 1 came from the")
print("   copy strand, so that the pair's orientation is inverted relative to its fragment's")
print("   strand -- the CTOT and CTOB pairs above. A directional library never produces one,")
print("   and asTair is directional-only by design (--library directional|reverse; the docs")
print("   call non-directional 'under development'). This measures that documented limit.")
print("   What TAPS changes is only that nothing catches it earlier: the reference is not")
print("   converted, so bwa maps such a pair as happily as any other.")
PYEOF

say "alnbase against asTair, position by position"
"$PY" - <<'PYEOF'
import pyarrow.parquet as pq
from collections import defaultdict

mine = defaultdict(lambda: [0, 0])
for r in pq.read_table("hits_0_0.parquet").to_pylist():
    mine[r["refr_pos"]][0 if r["name"] == "CG_mod" else 1] += 1

theirs = {}
for line in open("out/pe_mCtoT_CpG.mods"):
    if line.startswith("#"):
        continue
    f = line.split("\t")
    theirs[int(f[1])] = (int(f[4]), int(f[5]))

both = set(mine) & set(theirs)
diff = [p for p in sorted(both) if tuple(mine[p]) != theirs[p]]
print("   asTair positions: %d   alnbase positions: %d   in common: %d   disagreeing: %d"
      % (len(theirs), len(mine), len(both), len(diff)))
print("   asTair  totals: mod %d unmod %d" % (sum(m for m, _ in theirs.values()),
                                              sum(u for _, u in theirs.values())))
print("   alnbase totals: mod %d unmod %d" % (sum(v[0] for v in mine.values()),
                                              sum(v[1] for v in mine.values())))
for p in diff:
    print("   disagree at %d: alnbase %s, asTair %s" % (p, tuple(mine[p]), theirs[p]))
for p in sorted(set(mine) - set(theirs)):
    print("   only alnbase reports %d: %s" % (p, tuple(mine[p])))
for p in sorted(set(theirs) - set(mine)):
    print("   only asTair reports %d: %s" % (p, theirs[p]))
print("   The positions asTair cannot reach are the ones whose trinucleotide SPECIFIC_CONTEXT")
print("   runs off the end of the contig; alnbase's CG query needs only the dinucleotide.")
PYEOF
