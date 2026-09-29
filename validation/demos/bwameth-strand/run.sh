#!/usr/bin/env bash
# Run bwa-meth on the same four-strand fixture the other strand demos use, and check what
# queries/strand/bwameth.toml makes of the BAM it writes. bwa-meth's YD:Z is one character,
# the converted contig the read hit, so it reaches the conversion strand and no further; the
# file pairs it with FLAG 0x10 for the sequenced direction. Three claims in that file's
# comment are measured here:
#
#   - 0x10 is conformant, because bwa-meth leaves bwa's FLAG alone and restores SEQ to the
#     reference's forward orientation (bwameth.py:500-503);
#   - YC:Z carries nothing the FLAG does not, being assigned from the mate number alone
#     (bwameth.py:205-209);
#   - directional.toml also applies to this BAM, and additionally recovers the four-way
#     strand of origin.
#
# The last is the one worth a demo: two rules of different shape, read from different
# evidence, walking the same records. The fourth section runs both and compares them hit for
# hit. The third section measures what this directional-only aligner does when copy-strand
# reads reach it anyway, which is the case the rule's own comment does not cover.
#
#   BWAMETH_ENV=/path/to/env PY=/path/to/python-with-pyarrow ALNBASE=/path/to/alnbase \
#     ./run.sh [workdir]
#
# BWAMETH_ENV is the environment from ../../envs/bwameth.yaml. It provides bwameth.py, bwa
# and samtools.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${BWAMETH_ENV:?set BWAMETH_ENV}" "${PY:?set PY}" "${ALNBASE:?set ALNBASE}"
export PATH="$BWAMETH_ENV/bin:$PATH"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }

# The reference is 1200 bases from a fixed recurrence with a CG planted every 17, so every
# read carries several cytosines in CG context and the four loci stay unique after either
# conversion. Four 200-base fragments, one per strand of origin; each pair is named for the
# strand its READ 1 was sequenced from. Identical to ../bsmap-strand/ and ../biscuit-strand/.
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

# One genomic fragment yields four sequences in two duplexes: OT with its PCR copy CTOT, and
# OB with its copy CTOB. A pair is one duplex read from both ends, so read 1 comes from one
# of its two strands and read 2 from the other.
#
# bwa-meth supports DIRECTIONAL libraries only (its own README.md:28), which is the OT and OB
# pairs. The nd set adds the two pairs whose read 1 came from the copy, so that what bwa-meth
# does with input it does not claim to handle can be measured rather than assumed.
FRAG, LEN = 200, 60
MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}
pairs = [("OT", 100), ("CTOT", 350), ("OB", 600), ("CTOB", 850)]
handles = {k: (open("%s_1.fq" % k, "w"), open("%s_2.fq" % k, "w"))
           for k in ("dir", "nd")}
for name, start in pairs:
    frag = ref[start:start + FRAG]
    ot = bisulfite(frag)
    ob = bisulfite(rc(frag))
    strand = {"OT": ot, "CTOT": rc(ot), "OB": ob, "CTOB": rc(ob)}
    # Each read starts at its own strand's 5' end.
    r1, r2 = strand[name][:LEN], strand[MATE[name]][:LEN]
    sets = ["nd"] + (["dir"] if name in ("OT", "OB") else [])
    for s in sets:
        f1, f2 = handles[s]
        f1.write("@%s\n%s\n+\n%s\n" % (name, r1, "I" * LEN))
        f2.write("@%s\n%s\n+\n%s\n" % (name, r2, "I" * LEN))
    print("%-5s fragment chr1 1-based %d-%d: read 1 from %s, read 2 from %s (%s)"
          % (name, start + 1, start + FRAG, name, MATE[name],
             "both runs" if "dir" in sets else "non-directional only"))
for f1, f2 in handles.values():
    f1.close()
    f2.close()
PYEOF
samtools faidx ref.fa

# bwameth.py converts the reference to a single C>T index whose contigs are the forward
# strand prefixed 'f' and the reverse strand prefixed 'r', converts read 1 C>T and read 2
# G>A, and hands the lot to bwa mem. YD is the prefix of the contig the read landed on.
say "bwa-meth, paired-end: the directional pairs, and all four"
bwameth.py index ref.fa > index.log 2>&1
bwameth.py --reference ref.fa dir_1.fq dir_2.fq 2> dir.log \
  | samtools sort -o dir.bam -
bwameth.py --reference ref.fa nd_1.fq  nd_2.fq  2> nd.log \
  | samtools sort -o nd.bam  -
samtools index dir.bam
samtools index nd.bam

say "What bwa-meth wrote for each record, and what the rule makes of it"
python3 - <<'PYEOF'
import subprocess

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s):
    return "".join(COMP[c] for c in reversed(s))

def fastq(path):
    lines = open(path).read().split("\n")
    return dict((lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 4))

ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))

# The rule in queries/strand/bwameth.toml, written out so the table below shows what it calls
# rather than asserting it. Two tables, because YD reaches only the conversion strand: the
# sequenced direction has to come from the FLAG.
CONV = {"f": "+", "r": "-"}
def rule(tags, flag):
    conv = CONV.get(tags.get("YD"))
    if conv is None:
        return "unknown"
    return "%s:%s" % (conv, "reverse" if flag & 0x10 else "forward")

MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}
# The conversion strand each strand of origin sits on: a strand and its own PCR copy share
# one, which is exactly the split a one-character tag cannot cross.
ORIGIN_CONV = {"OT": "+", "CTOT": "+", "OB": "-", "CTOB": "-"}

def report(label, bam, f1, f2):
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
        # in the "really" columns comes from the FLAG or from YD.
        as_read = r1[qname] if seq in (r1[qname], rc(r1[qname])) else r2[qname]
        is_r1 = as_read == r1[qname]
        revcomped = seq != as_read
        origin = qname if is_r1 else MATE[qname]
        window = ref[pos - 1:pos - 1 + len(seq)]
        mismatches = sum(1 for r, s in zip(window, seq)
                         if r != s and (r, s) not in (("C", "T"), ("G", "A")))
        want = "%s:%s" % (ORIGIN_CONV[origin], "reverse" if revcomped else "forward")
        rows.append(dict(qname=qname, is_r1=is_r1, flag=flag, pos=pos, origin=origin,
                         yd=tags.get("YD", "-"), yc=tags.get("YC", "-"),
                         revcomped=revcomped, mismatches=mismatches,
                         mapped=not flag & 0x4,
                         call=rule(tags, flag), want=want))

    print("\n-- %s" % label)
    print("%-5s %-9s %-6s %-5s %-3s %-3s %-24s %-6s %s"
          % ("pair", "really is", "origin", "FLAG", "YD", "YC", "rule says", "POS", "mism"))
    for r in sorted(rows, key=lambda r: (r["qname"], not r["is_r1"])):
        print("%-5s %-9s %-6s %-5d %-3s %-3s %-24s %-6d %d"
              % (r["qname"], "read 1" if r["is_r1"] else "read 2", r["origin"], r["flag"],
                 r["yd"], r["yc"], "%s (want %s)" % (r["call"], r["want"]), r["pos"],
                 r["mismatches"]))

    # bwa-meth tags YD only on mapped records (bwameth.py:469-473 returns early for the
    # rest), so an unmapped one is `unknown` to the rule and alnbase skips and counts it.
    # Everything below is therefore about the mapped records; 0x10, SEQ orientation and POS
    # all mean nothing on an unmapped one, whatever the file happens to hold there.
    mapped = [r for r in rows if r["mapped"]]
    print("records: %d, of which unmapped and so `unknown` to the rule: %d"
          % (len(rows), len(rows) - len(mapped)))
    bad = [r for r in mapped if r["call"] != r["want"]]
    print("mapped records the rule gets wrong (conversion strand or sequenced direction): "
          "%d of %d" % (len(bad), len(mapped)))
    # The FLAG claim the rule depends on, since [strand.sequenced] reads nothing else.
    lying = [r for r in mapped if r["revcomped"] != bool(r["flag"] & 0x10)]
    print("mapped records whose 0x10 bit contradicts how SEQ was actually stored: %d"
          % len(lying))
    print("mapped records not stored reference-forward: %d"
          % len([r for r in mapped if r["mismatches"]]))
    # The YC claim: assigned from the mate number alone, so it is CT on every first-in-pair
    # record and GA on every second-in-pair one, whatever strand the record came from.
    yc = set((bool(r["flag"] & 0x40), r["yc"]) for r in rows)
    print("distinct (first-in-pair, YC) combinations: %s"
          % ", ".join("%s->%s" % ("read1" if a else "read2", b) for a, b in sorted(yc)))
    # Both mates of a directional fragment sit on one conversion strand, so if YD is the
    # conversion strand it must agree across a pair whenever both mates mapped.
    print("pairs with both mates mapped whose two YD values differ: %d"
          % len(set(r["qname"] for r in mapped
                    for o in mapped
                    if o["qname"] == r["qname"] and o["is_r1"] != r["is_r1"]
                    and o["yd"] != r["yd"])))

report("the directional pairs, which is what bwa-meth supports", "dir.bam",
       "dir_1.fq", "dir_2.fq")
report("all four pairs, including the two bwa-meth does not claim to handle", "nd.bam",
       "nd_1.fq", "nd_2.fq")
PYEOF

### What a directional-only aligner does with copy-strand reads
#
# bwa-meth converts read 1 C>T and read 2 G>A (bwameth.py:205) on the assumption that read 1
# was sequenced from an original strand and read 2 from its copy. A CTOT or CTOB read
# arriving as read 1 gets the wrong conversion applied, and the question this section answers
# is what bwa mem then does with it: drop it, clip it, or place it somewhere. The sweep is a
# 60 bp read from each of the four strands of origin every 10 bases along the reference,
# aligned single-end, where every read is treated as read 1 and so converted C>T.
say "What bwa-meth does with copy-strand reads"
python3 - <<'PYEOF'
ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s):
    return "".join(COMP[c] for c in reversed(s))
def bisulfite(s):
    return "".join("T" if c == "C" and s[i + 1:i + 2] != "G" else c
                   for i, c in enumerate(s))

LEN, STEP = 60, 10
with open("sweep.fq", "w") as out:
    for start in range(0, len(ref) - LEN, STEP):
        w = ref[start:start + LEN]
        ot, ob = bisulfite(w), bisulfite(rc(w))
        for name, seq in (("OT", ot), ("CTOT", rc(ot)), ("OB", ob), ("CTOB", rc(ob))):
            out.write("@%s_%d\n%s\n+\n%s\n" % (name, start, seq, "I" * LEN))
PYEOF
bwameth.py --reference ref.fa sweep.fq 2> sweep.log | samtools sort -o sweep.bam -
python3 - <<'PYEOF'
import subprocess
from collections import Counter

ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))
CONV = {"f": "+", "r": "-"}
ORIGIN_CONV = {"OT": "+", "CTOT": "+", "OB": "-", "CTOB": "-"}

lines = open("sweep.fq").read().split("\n")
reads = dict((lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 4))

seen, mapped, wrong, misplaced, clipped, qcfail = (Counter() for _ in range(6))
mapq = {}
sam = subprocess.run(["samtools", "view", "sweep.bam"], capture_output=True,
                     text=True).stdout
for line in sam.rstrip("\n").split("\n"):
    f = line.split("\t")
    flag = int(f[1])
    if flag & 0x100 or flag & 0x800:      # one row per read, primary only
        continue
    name, start = f[0].split("_")
    start = int(start)
    seen[name] += 1
    if flag & 0x4:
        continue
    mapped[name] += 1
    tags = dict((t[:2], t[5:]) for t in f[11:])
    if CONV.get(tags.get("YD")) != ORIGIN_CONV[name]:
        wrong[name] += 1
    if int(f[3]) - 1 != start:
        misplaced[name] += 1
    if "S" in f[5] or "H" in f[5]:
        clipped[name] += 1
    if flag & 0x200:
        qcfail[name] += 1
    mapq.setdefault(name, []).append(int(f[4]))

print("%-5s %-6s %-7s %-16s %-9s %-8s %s"
      % ("from", "reads", "mapped", "conversion wrong", "misplaced", "clipped", "median MAPQ"))
for n in ("OT", "CTOT", "OB", "CTOB"):
    q = sorted(mapq.get(n, [0]))
    print("%-5s %-6d %-7d %-16d %-9d %-8d %d"
          % (n, seen[n], mapped[n], wrong[n], misplaced[n], clipped[n], q[len(q) // 2]))
print("records flagged QC-fail (0x200): %d" % sum(qcfail.values()))
PYEOF

say "alnbase, reading the rule and nothing else"
"$ALNBASE" index ref.fa ref.aref > /dev/null
# Two rules over the same BAM. bwameth.toml takes the conversion strand from YD and the
# sequenced direction from 0x10; directional.toml takes both from 0x10 and the mate bit, on
# the library's design rather than the aligner's evidence, and so also reaches conv_strand.
# The rule file's comment says both apply to a bwa-meth BAM. Whether they agree is measured.
for set in dir nd; do
  for rule in bwameth directional; do
    printf -- '-- %s, %s.toml\n' "$set" "$rule"
    "$ALNBASE" query --parquet --only-hits \
      -f qname,flags,strand,conv_strand,read_reverse \
      --query-file "$HERE/cg.toml" \
      --query-file "$HERE/../../../queries/strand/$rule.toml" \
      "$set.bam" ref.aref "$set.$rule.hits.parquet"
  done
done

# Every hit's anchor is a cytosine in the read. Which reference base it sits on is decided by
# the conversion strand alnbase took from the rule: the C of a forward CG for '+', the G of
# the same CG for '-'. Reading those bases back out of the reference is a check on the walk
# rather than on the label.
"$PY" - <<'PYEOF'
import glob
import pyarrow.parquet as pq

ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))
MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}

def load(name, rule):
    t = pq.read_table(sorted(glob.glob("%s.%s.hits_*.parquet" % (name, rule))))
    return t.to_pylist()

for name, label in (("dir", "the directional pairs"), ("nd", "all four pairs")):
    print("\n-- %s" % label)
    both = {}
    for rule in ("bwameth", "directional"):
        rows = both[rule] = load(name, rule)
        by_record = {}
        for r in rows:
            by_record.setdefault((r["qname"], r["flags"]), []).append(r)
        print("%s.toml: %d CG hits over %d records" % (rule, len(rows), len(by_record)))
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
        wrong = [r for r in rows if ref[r["refr_pos"]] != want[r["strand"]]]
        print("   hits whose anchor is not the informative base for the strand walked: %d"
              % len(wrong))

    # The claim in the rule file's comment: directional.toml also applies to a bwa-meth BAM,
    # and adds conv_strand. Two rules agree when they put the same anchors in the same
    # places; keyed on the record and the reference position, so a strand disagreement moves
    # the anchor from a C to the G beside it and shows up as a difference.
    key = lambda rows: sorted((r["qname"], r["flags"], r["refr_pos"], r["strand"])
                              for r in rows)
    a, b = key(both["bwameth"]), key(both["directional"])
    print("   hits where the two rules disagree: %d (of %d and %d)"
          % (len(set(a) ^ set(b)), len(a), len(b)))
    conv = both["directional"]
    unnamed = [r for r in conv if r["conv_strand"] not in (r["qname"], MATE[r["qname"]])]
    print("   directional.toml: hits with no strand of origin: %d; hits whose strand of "
          "origin is neither the pair's nor its mate's: %d"
          % (len([r for r in conv if r["conv_strand"] is None]), len(unnamed)))
PYEOF
