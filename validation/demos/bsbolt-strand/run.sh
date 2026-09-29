#!/usr/bin/env bash
# Run BSBolt paired-end, directional and -UN, and check what queries/strand/bsbolt.toml
# makes of the BAMs it writes -- in particular the claim that FLAG 0x10 is inverted on
# every *_G2A record while SEQ stays reference-forward.
#
#   BSBOLT_ENV=/path/to/env ALNBASE=/path/to/alnbase ./run.sh [workdir]
#
# BSBOLT_ENV is the environment from ../../envs/premethyst.yaml with BSBolt checked out
# into it by ../../envs/build_bsbolt.sh (commit ea4870e, v1.6.0). It provides bsbolt and
# samtools.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${BSBOLT_ENV:?set BSBOLT_ENV}" "${ALNBASE:?set ALNBASE}"
export PATH="$BSBOLT_ENV/bin:$PATH"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }

# The reference is 1200 bases from a fixed recurrence with a CG planted every 17, so every
# read carries several cytosines in CG context and the four loci stay unique after either
# conversion. Four 200-base fragments, one per strand of origin; each pair is named for the
# strand its READ 1 was sequenced from, which is the answer the rule has to produce for
# that record.
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

# Bisulfite converts each single strand on its own, and PCR then copies each converted
# strand, so one genomic fragment yields four sequences in two duplexes: OT with its copy
# CTOT, and OB with its copy CTOB. A pair is one duplex read from both ends, so read 1
# comes from one of its two strands and read 2 from the other.
#
# A DIRECTIONAL library sequences read 1 from the converted original only, so it holds the
# OT and OB pairs alone; those two already exercise all four of BSBolt's YS values, because
# read 2 of each is the complementary strand. A NON-DIRECTIONAL library (-UN) adds the two
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
             "directional and -UN" if "dir" in sets else "-UN only"))
for f1, f2 in handles.values():
    f1.close()
    f2.close()
PYEOF
samtools faidx ref.fa

say "BSBolt, paired-end, directional and -UN"
bsbolt Index -G ref.fa -DB db > index.log 2>&1
bsbolt Align -DB db -F1 dir_1.fq -F2 dir_2.fq -O dir > dir.log 2>&1
bsbolt Align -DB db -F1 nd_1.fq -F2 nd_2.fq -O nd -UN > nd.log 2>&1
samtools view -h -o dir.sam dir.bam
samtools view -h -o nd.sam nd.bam

say "What BSBolt wrote for each record, and which fastq read it really was"
python3 - <<'PYEOF'
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s):
    return "".join(COMP[c] for c in reversed(s))

def fastq(path):
    lines = open(path).read().split("\n")
    return dict((lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 4))

ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))

# The rule in queries/strand/bsbolt.toml, written out so the table below shows what it
# calls rather than asserting it.
RULE = {"W_C2T": "OT", "W_G2A": "CTOT", "C_C2T": "OB", "C_G2A": "CTOB", "WC": "unknown"}
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
        # the way in: both settled by comparing SEQ against the reads themselves, never
        # from the FLAG, since the FLAG is what is in question.
        as_read = r1[qname] if seq in (r1[qname], rc(r1[qname])) else r2[qname]
        is_r1 = as_read == r1[qname]
        revcomped = seq != as_read
        # Read 1 of a pair was sequenced from the strand the pair is named for; read 2 from
        # the other strand of the same duplex. That is the call the rule has to produce.
        want = qname if is_r1 else MATE[qname]
        window = ref[pos - 1:pos - 1 + len(seq)]
        mismatches = sum(1 for r, s in zip(window, seq)
                         if r != s and (r, s) not in (("C", "T"), ("G", "A")))
        rows.append((qname, is_r1, flag, tags["YS"], want, RULE[tags["YS"]],
                     revcomped, pos, mismatches))

    print("\n-- %s" % label)
    print("%-5s %-9s %-5s %-11s %-11s %-6s %-9s %-9s %-6s %s"
          % ("pair", "really is", "FLAG", "0x10 says", "SEQ really", "YS", "",
             "rule call", "POS", "mismatches"))
    for qname, is_r1, flag, ys, want, call, revcomped, pos, mm in sorted(
            rows, key=lambda r: (r[0], -r[1])):
        print("%-5s %-9s %-5d %-11s %-11s %-6s %-9s %-9s %-6d %d"
              % (qname, "read 1" if is_r1 else "read 2", flag,
                 "reversed" if flag & 0x10 else "forward",
                 "reversed" if revcomped else "forward",
                 ys, "want " + want, call, pos, mm))
    lying = [r for r in rows if r[6] != bool(r[2] & 0x10)]
    g2a = [r for r in rows if r[3].endswith("G2A")]
    wrong = [r for r in rows if r[4] != r[5]]
    print("records whose 0x10 bit contradicts how SEQ was actually stored: %d of %d"
          % (len(lying), len(rows)))
    print("  of which *_G2A records: %d of %d *_G2A records in all"
          % (len([r for r in lying if r[3].endswith("G2A")]), len(g2a)))
    print("records not stored reference-forward: %d"
          % len([r for r in rows if r[8]]))
    print("records the rule calls as something other than the strand they were built "
          "from: %d" % len(wrong))

report("directional: the OT and OB pairs, which is all a directional library holds",
       "dir.sam", "dir_1.fq", "dir_2.fq")
report("-UN: all four pairs", "nd.sam", "nd_1.fq", "nd_2.fq")
PYEOF

say "alnbase, reading the rule and nothing else"
"$ALNBASE" index ref.fa ref.aref > /dev/null
for set in dir nd; do
  printf -- '-- %s\n' "$set"
  "$ALNBASE" query --query-file "$HERE/origin.toml" \
    --query-file "$HERE/../../../queries/strand/bsbolt.toml" \
    "$set.bam" ref.aref "$set.tagged.bam"
  samtools view "$set.tagged.bam" \
    | awk '{yv="none"; for (i = 12; i <= NF; i++) if ($i ~ /^YV:Z:/) yv = substr($i, 6);
            printf "%-5s FLAG %-4d POS %-5d YV=%s\n", $1, $2, $4, yv}' \
    | sort
done
