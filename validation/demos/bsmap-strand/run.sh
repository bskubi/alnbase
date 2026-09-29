#!/usr/bin/env bash
# Run BSMAP paired-end at the default -n 0 and at -n 1, and check what
# queries/strand/bsmap.toml makes of the SAM it writes. BSMAP's ZS is two characters, the
# conversion strand followed by the read's orientation within it, so unlike BISCUIT's YD it
# names all four strands from the tag alone. The claim this demo exists to measure is the
# one in that file's comment: that a paired-end run emits all four ZS values even at -n 0,
# the default, which is documented as mapping "only to 2 forward strands".
#
# Measuring it turned up something else. At -n 1, BSMAP puts a fraction of the copy-strand
# reads (CTOT and CTOB) on the opposite conversion strand, at the right position, reported
# as a unique alignment. The third section below measures the rate and tests the mechanism;
# ../../../docs/design/research/strand-determination.md writes it up.
#
#   BSMAP_ENV=/path/to/env PY=/path/to/python-with-pyarrow ALNBASE=/path/to/alnbase \
#     ./run.sh [workdir]
#
# BSMAP_ENV is the environment from ../../envs/bsmap.yaml. It provides bsmap and samtools.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${BSMAP_ENV:?set BSMAP_ENV}" "${PY:?set PY}" "${ALNBASE:?set ALNBASE}"
export PATH="$BSMAP_ENV/bin:$PATH"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }

# The reference is 1200 bases from a fixed recurrence with a CG planted every 17, so every
# read carries several cytosines in CG context and the four loci stay unique after either
# conversion. Four 200-base fragments, one per strand of origin; each pair is named for the
# strand its READ 1 was sequenced from.
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
# A DIRECTIONAL library holds the OT and OB pairs alone, and those two already show both
# conversion strands and both sequenced directions. A NON-DIRECTIONAL library adds the two
# pairs whose read 1 came from the copy instead.
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

# -p 1 and -S 1 so the run is reproducible; BSMAP otherwise seeds from the system clock.
say "BSMAP, paired-end, -n 0 (the default) and -n 1 (all four strands)"
bsmap -a dir_1.fq -b dir_2.fq -d ref.fa -o dir.sam -p 1 -S 1        > dir.log 2>&1
bsmap -a nd_1.fq  -b nd_2.fq  -d ref.fa -o nd.sam  -p 1 -S 1 -n 1   > nd.log  2>&1
samtools view -b -o dir.bam dir.sam
samtools view -b -o nd.bam  nd.sam

say "What BSMAP wrote for each record, and what the rule makes of it"
python3 - <<'PYEOF'
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s):
    return "".join(COMP[c] for c in reversed(s))

def fastq(path):
    lines = open(path).read().split("\n")
    return dict((lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 4))

ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))

# The rule in queries/strand/bsmap.toml, written out so the table below shows what it calls
# rather than asserting it. One table, four keys: ZS names the strand of origin outright.
RULE = {"++": "OT", "+-": "CTOT", "-+": "OB", "--": "CTOB"}

MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}

def report(label, sam, f1, f2):
    r1, r2 = fastq(f1), fastq(f2)
    rows = []
    for line in open(sam):
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        qname, flag, pos, seq = f[0], int(f[1]), int(f[3]), f[9]
        tags = dict((t[:2], t[5:]) for t in f[11:])
        # Which fastq read this record holds, and whether SEQ was reverse-complemented on
        # the way in: both settled by comparing SEQ against the reads themselves, so that
        # nothing in the "really" columns comes from the FLAG or from ZS.
        as_read = r1[qname] if seq in (r1[qname], rc(r1[qname])) else r2[qname]
        is_r1 = as_read == r1[qname]
        revcomped = seq != as_read
        origin = qname if is_r1 else MATE[qname]
        window = ref[pos - 1:pos - 1 + len(seq)]
        mismatches = sum(1 for r, s in zip(window, seq)
                         if r != s and (r, s) not in (("C", "T"), ("G", "A")))
        zs = tags.get("ZS", "-")
        rows.append(dict(qname=qname, is_r1=is_r1, flag=flag, pos=pos, zs=zs, origin=origin,
                         revcomped=revcomped, mismatches=mismatches,
                         call=RULE.get(zs, "unknown")))

    print("\n-- %s" % label)
    print("%-5s %-9s %-6s %-5s %-3s %-16s %-6s %s"
          % ("pair", "really is", "origin", "FLAG", "ZS", "rule says", "POS", "mismatches"))
    for r in sorted(rows, key=lambda r: (r["qname"], not r["is_r1"])):
        print("%-5s %-9s %-6s %-5d %-3s %-16s %-6d %d"
              % (r["qname"], "read 1" if r["is_r1"] else "read 2", r["origin"], r["flag"],
                 r["zs"], "%s (want %s)" % (r["call"], r["origin"]), r["pos"],
                 r["mismatches"]))

    # The claim in the rule's comment: -n 0 is documented as mapping "only to 2 forward
    # strands", yet a paired-end run emits all four ZS values, because read 2 is searched
    # against the complementary chain. Counted rather than asserted.
    print("distinct ZS values emitted: %s"
          % " ".join(sorted(set(r["zs"] for r in rows))))
    bad = [r for r in rows if r["call"] != r["origin"]]
    # BSMAP's FLAG is conformant, and the rule never reads it. Measured anyway, because a
    # four-way rule that agreed with a lying FLAG would be worth knowing about.
    lying = [r for r in rows if r["revcomped"] != bool(r["flag"] & 0x10)]
    zs_flag = [r for r in rows
               if (r["zs"][0] != r["zs"][1]) != bool(r["flag"] & 0x10)]
    print("records whose strand of origin the rule gets wrong: %d of %d" % (len(bad), len(rows)))
    print("records whose 0x10 bit contradicts how SEQ was actually stored: %d" % len(lying))
    print("records where 0x10 is not set exactly when ZS's two characters differ: %d"
          % len(zs_flag))
    print("records not stored reference-forward: %d"
          % len([r for r in rows if r["mismatches"]]))

report("-n 0, the default: the OT and OB pairs, which is all a directional library holds",
       "dir.sam", "dir_1.fq", "dir_2.fq")
report("-n 1: all four pairs", "nd.sam", "nd_1.fq", "nd_2.fq")
PYEOF

### How often BSMAP names the wrong conversion strand, and why
#
# One of the four pairs above comes back with its read 1 on the wrong conversion strand, so
# the question is how often that happens and to which reads. The sweep below is the same
# construction at scale: a 60 bp read from each of the four strands of origin every 10 bases
# along the reference, aligned single-end with -n 1 so that all four strands are searched.
#
# The mechanism is in BSMAP's source. A read is searched in two chains -- as sequenced
# (chain 0) and reverse-complemented (chain 1), align.cpp:92-93 -- against two converted
# reference chains, and the pair of indices becomes ZS (align.cpp:269). A hit is then
# recorded by AddHit (align.h:232-249), whose duplicate test is
#
#     if(!hitset[_ghit.chr>>1].insert(_ghit.loc).second) return 0; //hit already exist
#
# and `>>1` drops the reference-chain bit, so the key is the forward-reference position
# alone. Read chain 0 is searched first (align.cpp:211). So when a copy-strand read -- CTOT
# or CTOB, the reads whose correct alignment needs chain 1 -- happens to fit the *other*
# converted chain within the mismatch budget, that wrong-strand hit is registered first and
# the perfect chain-1 hit at the same position is then discarded as a duplicate. The read is
# placed correctly, and only its conversion strand and NM are wrong.
#
# The budget is 5 for a 60 bp read at the default -v 0.08 (align.cpp:501:
# (108-100)/100*60+0.5). The sweep tests that prediction against what BSMAP did.
say "How often BSMAP names the wrong conversion strand"
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
bsmap -a sweep.fq -d ref.fa -o sweep.sam -p 1 -S 1 -n 1 > sweep.log 2>&1
grep -E "aligned reads" sweep.log
python3 - <<'PYEOF'
from collections import Counter

ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s):
    return "".join(COMP[c] for c in reversed(s))
xt = lambda s: s.replace("C", "T")       # the C>T conversion BSMAP applies to both chains

RULE = {"++": "OT", "+-": "CTOT", "-+": "OB", "--": "CTOB"}
BUDGET = 5

lines = open("sweep.fq").read().split("\n")
reads = dict((lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 4))

total, wrong, called, misplaced = Counter(), Counter(), Counter(), 0
rows = []
for line in open("sweep.sam"):
    if line.startswith("@"):
        continue
    f = line.rstrip("\n").split("\t")
    name, start = f[0].split("_")
    start, pos = int(start), int(f[3])
    zs = dict((t[:2], t[5:]) for t in f[11:])["ZS"]
    call = RULE.get(zs, "unknown")
    total[name] += 1
    misplaced += pos - 1 != start
    if call != name:
        wrong[name] += 1
        called[(name, call)] += 1
    if name in ("CTOT", "CTOB"):
        # The candidate BSMAP looks at first: the read as sequenced, against the converted
        # chain its correct alignment does not use.
        w = ref[start:start + len(reads[f[0]])]
        other = rc(w) if name == "CTOT" else w
        mm = sum(1 for a, b in zip(xt(reads[f[0]]), xt(other)) if a != b)
        rows.append((mm <= BUDGET, call != name))

for n in ("OT", "CTOT", "OB", "CTOB"):
    print("%-4s reads aligned %3d, strand of origin wrong %2d" % (n, total[n], wrong[n]))
print("what the wrong ones were called: %s"
      % ", ".join("%s as %s (%d)" % (a, b, c) for (a, b), c in sorted(called.items())))
print("reads placed at the wrong position: %d" % misplaced)
print("copy-strand reads where 'the wrong-chain hit fits in %d mismatches' predicts what "
      "BSMAP did: %d of %d" % (BUDGET, sum(1 for p, a in rows if p == a), len(rows)))
PYEOF

say "alnbase, reading the rule and nothing else"
"$ALNBASE" index ref.fa ref.aref > /dev/null
for set in dir nd; do
  printf -- '-- %s\n' "$set"
  "$ALNBASE" query --parquet --only-hits \
    -f qname,flags,strand,conv_strand,read_reverse \
    --query-file "$HERE/cg.toml" \
    --query-file "$HERE/../../../queries/strand/bsmap.toml" \
    "$set.bam" ref.aref "$set.hits.parquet"
done

# Every hit's anchor is a cytosine in the read. Which reference base it sits on is decided by
# the conversion strand alnbase took from the rule: the C of a forward CG for '+', the G of
# the same CG for '-'. Reading those bases back out of the reference is a check on the walk
# rather than on the label. Unlike the BISCUIT demo, conv_strand is filled in here: a
# four-way tag reaches the strand of origin, so alnbase reports it.
"$PY" - <<'PYEOF'
import glob
import pyarrow.parquet as pq

ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))
MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}

for name, label in (("dir", "-n 0"), ("nd", "-n 1")):
    t = pq.read_table(sorted(glob.glob("%s.hits_*.parquet" % name)))
    rows = t.to_pylist()
    print("\n-- %s: %d CG hits over %d records"
          % (label, len(rows), len(set((r["qname"], r["flags"]) for r in rows))))
    print("%-5s %-5s %-6s %-11s %-13s %-4s %s"
          % ("pair", "FLAG", "strand", "conv_strand", "read_reverse", "hits",
             "reference base under each anchor"))
    by_record = {}
    for r in rows:
        by_record.setdefault((r["qname"], r["flags"]), []).append(r)
    for (qname, flag), rs in sorted(by_record.items()):
        base = set(ref[r["refr_pos"]] for r in rs)
        print("%-5s %-5d %-6s %-11s %-13s %-4d %s"
              % (qname, flag, rs[0]["strand"],
                 "null" if rs[0]["conv_strand"] is None else rs[0]["conv_strand"],
                 str(rs[0]["read_reverse"]), len(rs), "".join(sorted(base))))
    want = {"+": "C", "-": "G"}
    wrong = [r for r in rows if ref[r["refr_pos"]] != want[r["strand"]]]
    print("hits whose anchor is not the informative base for the strand alnbase walked: %d"
          % len(wrong))
    # Each pair's read 1 came from the strand the pair is named for and its read 2 from that
    # strand's complement, so conv_strand has one right answer per record.
    missing = [r for r in rows if r["conv_strand"] is None]
    unnamed = [r for r in rows
               if r["conv_strand"] not in (r["qname"], MATE[r["qname"]])]
    print("hits with no strand of origin: %d; hits whose strand of origin is neither the "
          "pair's nor its mate's: %d" % (len(missing), len(unnamed)))
PYEOF
