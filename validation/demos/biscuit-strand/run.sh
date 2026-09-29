#!/usr/bin/env bash
# Run BISCUIT paired-end, directional (-b 1) and non-directional (-b 0, the default), and
# check what queries/strand/biscuit.toml makes of the BAMs it writes. This is the first
# [strand.conversion] rule run against real output: a two-way tag for the conversion strand
# and the FLAG for the sequenced direction, rather than one tag that names all four strands.
#
#   BISCUIT_ENV=/path/to/env PY=/path/to/python-with-pyarrow ALNBASE=/path/to/alnbase \
#     ./run.sh [workdir]
#
# BISCUIT_ENV is the environment from ../../envs/biscuit.yaml with BISCUIT built into it by
# ../../envs/build_biscuit.sh (commit 0a5ceae, 1.10.3). It provides biscuit and samtools.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${BISCUIT_ENV:?set BISCUIT_ENV}" "${PY:?set PY}" "${ALNBASE:?set ALNBASE}"
export PATH="$BISCUIT_ENV/bin:$PATH"
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

say "BISCUIT, paired-end, -b 1 (directional) and -b 0 (non-directional, the default)"
biscuit index ref.fa > index.log 2>&1
biscuit align -b 1 ref.fa dir_1.fq dir_2.fq 2> dir.log | samtools view -b -o dir.bam
biscuit align       ref.fa nd_1.fq  nd_2.fq  2> nd.log  | samtools view -b -o nd.bam
samtools view -h -o dir.sam dir.bam
samtools view -h -o nd.sam nd.bam

say "What BISCUIT wrote for each record, and what the rule makes of it"
python3 - <<'PYEOF'
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s):
    return "".join(COMP[c] for c in reversed(s))

def fastq(path):
    lines = open(path).read().split("\n")
    return dict((lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 4))

ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))

# The rule in queries/strand/biscuit.toml, written out so the table below shows what it
# calls rather than asserting it. Unlike the four-way rules, it answers in two parts.
CONV = {"f": "+", "r": "-", "u": "unknown"}
SEQD = {False: "forward", True: "reverse"}          # [strand.sequenced], which reads 0x10

# What the reads were built from. Both columns follow from the strand of origin by biology
# (docs/design/strand-rules.md): OT and CTOT are the converted top strand and its copy, so
# the informative reference base is C; OT and CTOB carry top-strand sequence, so they read
# along the reference forward.
MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}
WANT_CONV = {"OT": "+", "CTOT": "+", "OB": "-", "CTOB": "-"}
WANT_SEQD = {"OT": "forward", "CTOB": "forward", "CTOT": "reverse", "OB": "reverse"}

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
        # nothing in the "really" columns comes from the FLAG the rule is trusting.
        as_read = r1[qname] if seq in (r1[qname], rc(r1[qname])) else r2[qname]
        is_r1 = as_read == r1[qname]
        revcomped = seq != as_read
        origin = qname if is_r1 else MATE[qname]
        window = ref[pos - 1:pos - 1 + len(seq)]
        mismatches = sum(1 for r, s in zip(window, seq)
                         if r != s and (r, s) not in (("C", "T"), ("G", "A")))
        rows.append(dict(qname=qname, is_r1=is_r1, flag=flag, pos=pos, yd=tags.get("YD", "-"),
                         origin=origin, revcomped=revcomped, mismatches=mismatches,
                         conv=CONV.get(tags.get("YD", "-"), "?"),
                         seqd=SEQD[bool(flag & 0x10)]))

    print("\n-- %s" % label)
    print("%-5s %-9s %-6s %-5s %-3s %-12s %-14s %-6s %s"
          % ("pair", "really is", "origin", "FLAG", "YD", "rule: conv", "rule: sequenced",
             "POS", "mismatches"))
    for r in sorted(rows, key=lambda r: (r["qname"], not r["is_r1"])):
        print("%-5s %-9s %-6s %-5d %-3s %-12s %-14s %-6d %d"
              % (r["qname"], "read 1" if r["is_r1"] else "read 2", r["origin"], r["flag"],
                 r["yd"], "%s (want %s)" % (r["conv"], WANT_CONV[r["origin"]]),
                 "%s (want %s)" % (r["seqd"], WANT_SEQD[r["origin"]]),
                 r["pos"], r["mismatches"]))

    # Which strands of origin BISCUIT put on each converted contig, counted rather than
    # assumed: this is the grouping that decides what a two-way tag can and cannot say.
    for yd in ("f", "r", "u"):
        got = sorted(set(r["origin"] for r in rows if r["yd"] == yd))
        if got:
            print("YD:A:%s records came from %s" % (yd, " and ".join(got)))
    bad_conv = [r for r in rows if r["conv"] != WANT_CONV[r["origin"]]]
    bad_seqd = [r for r in rows if r["seqd"] != WANT_SEQD[r["origin"]]]
    lying = [r for r in rows if r["revcomped"] != bool(r["flag"] & 0x10)]
    print("records whose conversion strand the rule gets wrong: %d of %d"
          % (len(bad_conv), len(rows)))
    print("records whose sequenced direction the rule gets wrong: %d of %d"
          % (len(bad_seqd), len(rows)))
    print("records whose 0x10 bit contradicts how SEQ was actually stored: %d"
          % len(lying))
    print("records not stored reference-forward: %d"
          % len([r for r in rows if r["mismatches"]]))

report("-b 1: the OT and OB pairs, which is all a directional library holds",
       "dir.sam", "dir_1.fq", "dir_2.fq")
report("-b 0 (the default): all four pairs", "nd.sam", "nd_1.fq", "nd_2.fq")
PYEOF

say "alnbase, reading the rule and nothing else"
"$ALNBASE" index ref.fa ref.aref > /dev/null
for set in dir nd; do
  printf -- '-- %s\n' "$set"
  "$ALNBASE" query --parquet --only-hits \
    -f qname,flags,strand,conv_strand,read_reverse \
    --query-file "$HERE/cg.toml" \
    --query-file "$HERE/../../../queries/strand/biscuit.toml" \
    "$set.bam" ref.aref "$set.hits.parquet"
done

# Every hit's anchor is a cytosine in the read. Which reference base it sits on is decided by
# the conversion strand alnbase took from the rule: the C of a forward CG for '+', the G of
# the same CG for '-'. Reading those bases back out of the reference is a check on the walk
# rather than on the label.
"$PY" - <<'PYEOF'
import glob
import pyarrow.parquet as pq

ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))

for name, label in (("dir", "-b 1"), ("nd", "-b 0")):
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
    print("records given a strand of origin by a rule that names only the conversion "
          "strand: %d" % len([r for r in rows if r["conv_strand"] is not None]))
PYEOF
