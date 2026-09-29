#!/usr/bin/env bash
# Run HISAT-3N on the same four-strand fixture the other strand demos use, and check what
# queries/strand/hisat-3n.toml makes of the BAM it writes. HISAT-3N's YZ:A is one character,
# + or -, the 3N index the read hit, so it reaches the conversion strand and no further; the
# file pairs it with FLAG 0x10 for the sequenced direction.
#
# What makes this rule different from the other two of its shape is where YZ comes from. In
# HISAT-3N's default non-directional mode it is INFERRED, by comparing how many C>T and how
# many G>A differences the alignment shows (alignment_3n.h:288-313), and when those counts
# are equal -- including when both are zero, which is what a fully methylated read looks like
# -- it falls back to the directional assumption (:299-308). The rule file says so; this demo
# measures it. The fixture is therefore aligned three times: with the usual partly converted
# reads, with fully methylated ones that carry no conversion at all, and with reads carrying
# one conversion of each kind, which tie at a non-zero count. The two ties are not equivalent
# downstream, because Yf:i:0 marks the first and nothing marks the second.
#
#   HISAT3N_ENV=/path/to/env PY=/path/to/python-with-pyarrow ALNBASE=/path/to/alnbase \
#     ./run.sh [workdir]
#
# HISAT3N_ENV is the environment from ../../envs/hisat-3n.yaml, built by
# ../../envs/build_hisat3n.sh. It provides hisat-3n, hisat-3n-build and samtools.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${HISAT3N_ENV:?set HISAT3N_ENV}" "${PY:?set PY}" "${ALNBASE:?set ALNBASE}"
export PATH="$HISAT3N_ENV/bin:$PATH"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }

# The reference is 1200 bases from a fixed recurrence with a CG planted every 17, so every
# read carries several cytosines in CG context and the four loci stay unique after either
# conversion. Four 200-base fragments, one per strand of origin; each pair is named for the
# strand its READ 1 was sequenced from. Identical to ../bsmap-strand/ and ../bwameth-strand/,
# except that each pair is written three times, at three amounts of conversion evidence.
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

def bisulfite(s, k=None):
    """Unmethylated cytosines become T; a cytosine in CG stays C. k caps how many are
    converted, so k=0 is a fully methylated read: every cytosine survives and the read is
    identical to the strand it came from."""
    sites = [i for i, c in enumerate(s) if c == "C" and s[i + 1:i + 2] != "G"]
    if k is not None:
        mid = len(s) // 2
        sites = sorted(sites, key=lambda i: abs(i - mid))[:k]
    out = list(s)
    for i in sites:
        out[i] = "T"
    return "".join(out)

def one_g_to_a(s):
    """Turn the G nearest the middle into an A, leaving the G of a CG alone. It is a plain
    substitution, but the aligner counts it as a conversion of the other kind, so a read
    carrying one of these and one converted cytosine ties with both counts at one rather
    than at zero."""
    mid = len(s) // 2
    for i in sorted(range(len(s)), key=lambda i: abs(i - mid)):
        if s[i] == "G" and s[i - 1:i] != "C":
            return s[:i] + "A" + s[i + 1:]
    return s

# One genomic fragment yields four sequences in two duplexes: OT with its PCR copy CTOT, and
# OB with its copy CTOB. A pair is one duplex read from both ends, so read 1 comes from one
# of its two strands and read 2 from the other.
#
# HISAT-3N's default mode is non-directional, so all four pairs are input it claims to
# handle. Each pair is written three times: `part`, the usual reads, where the conversion
# counts settle the strand; `full`, with nothing converted, so both counts are zero; and
# `tie`, with one converted cytosine and one G>A substitution, so both counts are one. The
# last two are the two ways the counts can fail to decide.
FRAG, LEN = 200, 60

def treat(s, k, subst):
    """Convert a strand of the fragment. A read is the first LEN bases of its own strand,
    so on the strand that was treated the two reads of a pair cover its first LEN and its
    last LEN bases; a small k has to be spent inside those windows, since converting the
    fragment as a whole would put it in the middle where no read would see it."""
    if k is None:
        return bisulfite(s, None)
    ends = [bisulfite(w, k) for w in (s[:LEN], s[-LEN:])]
    if subst:
        ends = [one_g_to_a(w) for w in ends]
    return ends[0] + s[LEN:-LEN] + ends[1]

MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}
pairs = [("OT", 100), ("CTOT", 350), ("OB", 600), ("CTOB", 850)]
for label, k, subst in (("part", None, False), ("full", 0, False), ("tie", 1, True)):
    f1 = open("%s_1.fq" % label, "w")
    f2 = open("%s_2.fq" % label, "w")
    for name, start in pairs:
        frag = ref[start:start + FRAG]
        ot = treat(frag, k, subst)
        ob = treat(rc(frag), k, subst)
        strand = {"OT": ot, "CTOT": rc(ot), "OB": ob, "CTOB": rc(ob)}
        # Each read starts at its own strand's 5' end.
        r1, r2 = strand[name][:LEN], strand[MATE[name]][:LEN]
        f1.write("@%s\n%s\n+\n%s\n" % (name, r1, "I" * LEN))
        f2.write("@%s\n%s\n+\n%s\n" % (name, r2, "I" * LEN))
        if label == "part":
            print("%-5s fragment chr1 1-based %d-%d: read 1 from %s, read 2 from %s"
                  % (name, start + 1, start + FRAG, name, MATE[name]))
    f1.close()
    f2.close()
print("written three times: `part`, CG methylated and CH converted; `full`, nothing "
      "converted; `tie`, one conversion and one G>A substitution in every read")
PYEOF
samtools faidx ref.fa

# hisat-3n-build writes a pair of three-letter indexes, the reference with C>T and with G>A;
# hisat-3n converts the read both ways, aligns to both, and reports one alignment. YZ names
# the index it came from. --no-spliced-alignment because HISAT is an RNA aligner by default
# and this fixture is genomic.
say "HISAT-3N, paired-end: the three fixtures"
hisat-3n-build --base-change C,T ref.fa ref3n > build.log 2>&1
for set in part full tie; do
  hisat-3n --index ref3n --base-change C,T --no-spliced-alignment \
    -1 "${set}_1.fq" -2 "${set}_2.fq" 2> "$set.log" | samtools sort -o "$set.bam" -
  samtools index "$set.bam"
done
# The same fully methylated reads under --directional-mapping, which is the library telling
# the aligner what the counts could not. Compared against the default run below.
hisat-3n --index ref3n --base-change C,T --no-spliced-alignment --directional-mapping \
  -1 full_1.fq -2 full_2.fq 2> full_dm.log | samtools sort -o full_dm.bam -

say "What HISAT-3N wrote for each record, and what the rule makes of it"
python3 - <<'PYEOF'
import subprocess

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s):
    return "".join(COMP[c] for c in reversed(s))

def fastq(path):
    lines = open(path).read().split("\n")
    return dict((lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 4))

ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))

# The rule in queries/strand/hisat-3n.toml, written out so the table below shows what it
# calls rather than asserting it. Two tables, because YZ reaches only the conversion strand:
# the sequenced direction has to come from the FLAG.
def rule(tags, flag):
    conv = tags.get("YZ")
    if conv not in ("+", "-"):
        return "unknown"
    return "%s:%s" % (conv, "reverse" if flag & 0x10 else "forward")

MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}
# The conversion strand each strand of origin sits on: a strand and its own PCR copy share
# one, which is exactly the split a one-character tag cannot cross.
ORIGIN_CONV = {"OT": "+", "CTOT": "+", "OB": "-", "CTOB": "-"}

def rows_of(bam, f1, f2):
    r1, r2 = fastq(f1), fastq(f2)
    rows = []
    sam = subprocess.run(["samtools", "view", bam], capture_output=True, text=True).stdout
    for line in sam.rstrip("\n").split("\n"):
        if not line:
            continue
        f = line.split("\t")
        qname, flag, pos, seq = f[0], int(f[1]), int(f[3]), f[9]
        tags = dict((t[:2], t[5:]) for t in f[11:])
        # Which fastq read this record holds, and whether SEQ was reverse-complemented on the
        # way in: both settled by comparing SEQ against the reads themselves, so that nothing
        # in the "really" columns comes from the FLAG or from YZ.
        as_read = r1[qname] if seq in (r1[qname], rc(r1[qname])) else r2[qname]
        is_r1 = as_read == r1[qname]
        revcomped = seq != as_read
        origin = qname if is_r1 else MATE[qname]
        window = ref[pos - 1:pos - 1 + len(seq)]
        mismatches = sum(1 for r, s in zip(window, seq)
                         if r != s and (r, s) not in (("C", "T"), ("G", "A")))
        want = "%s:%s" % (ORIGIN_CONV[origin], "reverse" if revcomped else "forward")
        rows.append(dict(qname=qname, is_r1=is_r1, flag=flag, pos=pos, origin=origin,
                         yz=tags.get("YZ", "-"), yf=tags.get("Yf", "-"),
                         revcomped=revcomped, mismatches=mismatches,
                         mapped=not flag & 0x4, call=rule(tags, flag), want=want))
    return rows

def report(label, bam, f1, f2):
    rows = rows_of(bam, f1, f2)
    print("\n-- %s" % label)
    print("%-5s %-9s %-6s %-5s %-3s %-3s %-24s %-6s %s"
          % ("pair", "really is", "origin", "FLAG", "YZ", "Yf", "rule says", "POS", "mism"))
    for r in sorted(rows, key=lambda r: (r["qname"], not r["is_r1"])):
        print("%-5s %-9s %-6s %-5d %-3s %-3s %-24s %-6d %d"
              % (r["qname"], "read 1" if r["is_r1"] else "read 2", r["origin"], r["flag"],
                 r["yz"], r["yf"], "%s (want %s)" % (r["call"], r["want"]), r["pos"],
                 r["mismatches"]))

    mapped = [r for r in rows if r["mapped"]]
    print("records: %d, of which unmapped: %d" % (len(rows), len(rows) - len(mapped)))
    bad = [r for r in mapped if r["call"] != r["want"]]
    print("mapped records the rule gets wrong (conversion strand or sequenced direction): "
          "%d of %d" % (len(bad), len(mapped)))
    # Where the disagreement is: the rule reports YZ faithfully, so a wrong call is a wrong
    # YZ, not a misreading of it.
    print("   of those, wrong because YZ names the other conversion strand: %d"
          % len([r for r in bad if r["yz"] != ORIGIN_CONV[r["origin"]]]))
    print("   of those, carrying Yf:i:0, the marker that no conversion was seen: %d"
          % len([r for r in bad if r["yf"] == "0"]))
    # The FLAG claim the rule depends on, since [strand.sequenced] reads nothing else.
    lying = [r for r in mapped if r["revcomped"] != bool(r["flag"] & 0x10)]
    print("mapped records whose 0x10 bit contradicts how SEQ was actually stored: %d"
          % len(lying))
    print("mapped records not stored reference-forward: %d"
          % len([r for r in mapped if r["mismatches"]]))
    # Both mates of a fragment sit on one conversion strand, so if YZ is the conversion
    # strand it must agree across a pair whenever both mates mapped.
    print("pairs with both mates mapped whose two YZ values differ: %d"
          % len(set(r["qname"] for r in mapped for o in mapped
                    if o["qname"] == r["qname"] and o["is_r1"] != r["is_r1"]
                    and o["yz"] != r["yz"])))
    return rows

report("partly converted reads, where the conversion counts decide", "part.bam",
       "part_1.fq", "part_2.fq")
full = report("fully methylated reads, where the counts are zero and zero",
              "full.bam", "full_1.fq", "full_2.fq")
report("one conversion of each kind, where the counts are one and one", "tie.bam",
       "tie_1.fq", "tie_2.fq")

# With no conversions to count, the default non-directional mode takes the same branch
# --directional-mapping does. If that is so, the two BAMs carry the same YZ on every record.
dm = rows_of("full_dm.bam", "full_1.fq", "full_2.fq")
key = lambda rows: sorted((r["qname"], r["flag"], r["yz"]) for r in rows)
print("\nfully methylated reads, default mode vs --directional-mapping: "
      "records whose YZ differs: %d of %d"
      % (len(set(key(full)) ^ set(key(dm))) // 2, len(full)))
PYEOF

### How much evidence the inference needs
#
# YZ is inferred from the counts of C>T and G>A differences, so the question this section
# answers is how many conversions it takes for the evidence to beat the assumption, and what
# the record looks like when it does not. The sweep is a 60 bp read from each of the four
# strands of origin every 10 bases along the reference, in four states, aligned single-end so
# that every read is read 1.
say "How much evidence the inference needs"
python3 - <<'PYEOF'
ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s):
    return "".join(COMP[c] for c in reversed(s))
def bisulfite(s, k=None):
    """As in the fixture above. A small k converts the cytosines nearest the middle of the
    read, so that the handful of differences this section turns on are counted rather than
    soft-clipped off an end."""
    sites = [i for i, c in enumerate(s) if c == "C" and s[i + 1:i + 2] != "G"]
    if k is not None:
        mid = len(s) // 2
        sites = sorted(sites, key=lambda i: abs(i - mid))[:k]
    out = list(s)
    for i in sites:
        out[i] = "T"
    return "".join(out)

def one_g_to_a(s):
    """Turn one G into an A: a plain substitution, which the aligner counts as a G>A
    conversion of the other kind. The G of a CG is left alone, since changing it would be a
    methylation difference rather than a mismatch, and the G nearest the middle is taken for
    the same reason the conversions are."""
    mid = len(s) // 2
    for i in sorted(range(len(s)), key=lambda i: abs(i - mid)):
        if s[i] == "G" and s[i - 1:i] != "C":
            return s[:i] + "A" + s[i + 1:]
    return s

LEN, STEP = 60, 10
# part: every CH converted, the usual case.  one: a single converted cytosine, the least
# evidence that can still decide.  full: nothing converted, so both counts are zero.
# tie: one converted cytosine and one G>A substitution, so both counts are one.
STATES = ("part", "one", "full", "tie")
with open("sweep.fq", "w") as out:
    for start in range(0, len(ref) - LEN, STEP):
        w = ref[start:start + LEN]
        for state in STATES:
            k = {"part": None, "one": 1, "full": 0, "tie": 1}[state]
            ot, ob = bisulfite(w, k), bisulfite(rc(w), k)
            if state == "tie":
                ot, ob = one_g_to_a(ot), one_g_to_a(ob)
            for name, seq in (("OT", ot), ("CTOT", rc(ot)), ("OB", ob), ("CTOB", rc(ob))):
                out.write("@%s_%s_%d\n%s\n+\n%s\n" % (state, name, start, seq, "I" * LEN))
PYEOF
hisat-3n --index ref3n --base-change C,T --no-spliced-alignment -U sweep.fq 2> sweep.log \
  | samtools sort -o sweep.bam -
python3 - <<'PYEOF'
import subprocess
from collections import Counter

ORIGIN_CONV = {"OT": "+", "CTOT": "+", "OB": "-", "CTOB": "-"}
STATES = ("part", "one", "full", "tie")

seen, mapped, wrong, marked, misplaced, clipped = (Counter() for _ in range(6))
mapq = {}
sam = subprocess.run(["samtools", "view", "sweep.bam"], capture_output=True,
                     text=True).stdout
for line in sam.rstrip("\n").split("\n"):
    f = line.split("\t")
    flag = int(f[1])
    if flag & 0x100 or flag & 0x800:      # one row per read, primary only
        continue
    state, name, start = f[0].split("_")
    seen[(state, name)] += 1
    if flag & 0x4:
        continue
    mapped[(state, name)] += 1
    tags = dict((t[:2], t[5:]) for t in f[11:])
    if tags.get("YZ") != ORIGIN_CONV[name]:
        wrong[(state, name)] += 1
        if tags.get("Yf") == "0":
            marked[(state, name)] += 1
    if int(f[3]) - 1 != int(start):
        misplaced[(state, name)] += 1
    if "S" in f[5] or "H" in f[5]:
        clipped[(state, name)] += 1
    mapq.setdefault((state, name), []).append(int(f[4]))

print("%-6s %-5s %-6s %-7s %-16s %-16s %-9s %-8s %s"
      % ("state", "from", "reads", "mapped", "conversion wrong", "of those, Yf:i:0",
         "misplaced", "clipped", "median MAPQ"))
for state in STATES:
    for n in ("OT", "CTOT", "OB", "CTOB"):
        k = (state, n)
        q = sorted(mapq.get(k, [0]))
        print("%-6s %-5s %-6d %-7d %-16d %-16d %-9d %-8d %d"
              % (state, n, seen[k], mapped[k], wrong[k], marked[k], misplaced[k],
                 clipped[k], q[len(q) // 2]))
tot = lambda c, s: sum(v for k, v in c.items() if k[0] == s)
print("\n%-6s %-8s %-8s %s" % ("state", "mapped", "wrong", "wrong records carrying Yf:i:0"))
for state in STATES:
    print("%-6s %-8d %-8d %d" % (state, tot(mapped, state), tot(wrong, state),
                                 tot(marked, state)))
PYEOF

say "alnbase, reading the rule and nothing else"
"$ALNBASE" index ref.fa ref.aref > /dev/null
for set in part full tie; do
  printf -- '-- %s\n' "$set"
  "$ALNBASE" query --parquet --only-hits \
    -f qname,flags,strand,conv_strand,read_reverse \
    --query-file "$HERE/cg.toml" \
    --query-file "$HERE/../../../queries/strand/hisat-3n.toml" \
    "$set.bam" ref.aref "$set.hits.parquet"
done

# The scan line above counts read bases that differ from the reference, allowing whichever
# conversion the declared strand permits. It is what surfaced BSMAP's wrong-strand records,
# and on `full` it cannot do the same job: a read with no conversion matches the reference
# whichever way it is walked, so the count stays at zero while half the records sit on the
# wrong strand. On `tie` the count is not zero -- but splitting the fixture into the records
# the rule got right and the ones it got wrong shows that what it is counting is the planted
# substitution, which both halves carry, and not the strand.
for half in right wrong; do
  case "$half" in
    right) keep='$1=="OT"||$1=="OB"' ;;
    wrong) keep='$1=="CTOT"||$1=="CTOB"' ;;
  esac
  samtools view -H tie.bam > "tie.$half.sam"
  samtools view tie.bam | awk "$keep" >> "tie.$half.sam"
  samtools view -b "tie.$half.sam" > "tie.$half.bam"
  printf -- '-- tie, the 4 records the rule got %s\n' "$half"
  "$ALNBASE" query --parquet --only-hits \
    --query-file "$HERE/cg.toml" \
    --query-file "$HERE/../../../queries/strand/hisat-3n.toml" \
    "tie.$half.bam" ref.aref "tie.$half.hits.parquet"
done

# Every hit's anchor is a cytosine in the read. Which reference base it sits on is decided by
# the conversion strand alnbase took from the rule: the C of a forward CG for '+', the G of
# the same CG for '-'. On a record whose YZ names the other strand the whole walk moves, so
# the anchors land on the other base of every CG -- which is what a wrong conversion strand
# costs, stated in the coordinates a caller would emit.
"$PY" - <<'PYEOF'
import glob
import pyarrow.parquet as pq

ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))
MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}
ORIGIN_CONV = {"OT": "+", "CTOT": "+", "OB": "-", "CTOB": "-"}

for name, label in (("part", "partly converted reads"),
                    ("full", "fully methylated reads"),
                    ("tie", "one conversion of each kind")):
    rows = pq.read_table(sorted(glob.glob("%s.hits_*.parquet" % name))).to_pylist()
    by_record = {}
    for r in rows:
        by_record.setdefault((r["qname"], r["flags"]), []).append(r)
    print("\n-- %s: %d CG hits over %d records" % (label, len(rows), len(by_record)))
    print("   %-5s %-5s %-6s %-11s %-13s %-4s %s"
          % ("pair", "FLAG", "strand", "conv_strand", "read_reverse", "hits",
             "reference base under each anchor"))
    for (qname, flag), rs in sorted(by_record.items()):
        base = set(ref[r["refr_pos"]] for r in rs)
        print("   %-5s %-5d %-6s %-11s %-13s %-4d %s"
              % (qname, flag, rs[0]["strand"],
                 "null" if rs[0]["conv_strand"] is None else rs[0]["conv_strand"],
                 str(rs[0]["read_reverse"]), len(rs), "".join(sorted(base))))
    want = {"+": "C", "-": "G"}
    print("   hits whose anchor is not the informative base for the strand walked: %d"
          % len([r for r in rows if ref[r["refr_pos"]] != want[r["strand"]]]))
    # The strand the record really came from is the pair's own or its mate's; both sit on one
    # conversion strand, so the fixture knows what every hit should have been walked as.
    off = [r for r in rows if r["strand"] != ORIGIN_CONV[r["qname"]]]
    print("   hits walked along the conversion strand the read did NOT come from: %d of %d"
          % (len(off), len(rows)))
PYEOF
