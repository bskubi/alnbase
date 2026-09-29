#!/usr/bin/env bash
# Run Bismark --non_directional paired-end on one pair per strand of origin and check
# what queries/strand/bismark.toml makes of the BAM it writes.
#
#   BISMARK_ENV=/path/to/bismark-env ALNBASE=/path/to/alnbase ./run.sh [workdir]
#
# BISMARK_ENV is the environment from ../../envs/bismark.yaml: Bismark 0.25.1, bowtie2
# and samtools.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${BISMARK_ENV:?set BISMARK_ENV}" "${ALNBASE:?set ALNBASE}"
export PATH="$BISMARK_ENV/bin:$PATH"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }

# The reference is 1200 bases from a fixed recurrence with a CG planted every 17, so every
# read carries several cytosines in CG context and the four loci stay unique after either
# conversion. Four 200-base fragments, one per strand of origin; each pair is named for the
# strand its READ 1 was sequenced from, which is the answer the rule has to produce for
# that record.
say "The four pairs, each named for the strand its read 1 came from"
mkdir -p genome
python3 - <<'PYEOF'
bases = "ACGT"
ref, x = [], 7
for _ in range(1200):
    x = (x * 31 + 17) % 1009
    ref.append(bases[x % 4])
for i in range(0, 1198, 17):
    ref[i], ref[i + 1] = "C", "G"
ref = "".join(ref)

with open("genome/ref.fa", "w") as f:
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
# comes from one of its two strands and read 2 from the other. Which way round is exactly
# what a non-directional library leaves open, and what --non_directional is for; each pair
# here is named for the strand its read 1 was sequenced from.
FRAG, LEN = 200, 60
MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}
pairs = [("OT", 100), ("CTOT", 350), ("OB", 600), ("CTOB", 850)]
r1f = open("reads_1.fq", "w")
r2f = open("reads_2.fq", "w")
for name, start in pairs:
    frag = ref[start:start + FRAG]
    ot = bisulfite(frag)
    ob = bisulfite(rc(frag))
    strand = {"OT": ot, "CTOT": rc(ot), "OB": ob, "CTOB": rc(ob)}
    # Each read starts at its own strand's 5' end.
    r1, r2 = strand[name][:LEN], strand[MATE[name]][:LEN]
    r1f.write("@%s\n%s\n+\n%s\n" % (name, r1, "I" * LEN))
    r2f.write("@%s\n%s\n+\n%s\n" % (name, r2, "I" * LEN))
    print("%-5s fragment chr1 1-based %d-%d: read 1 from %s, read 2 from %s"
          % (name, start + 1, start + FRAG, name, MATE[name]))
r1f.close()
r2f.close()
PYEOF
samtools faidx genome/ref.fa

say "Bismark, paired-end, --non_directional so all four strands are aligned"
mkdir -p tmp
bismark_genome_preparation --bowtie2 genome > prep.log 2>&1
bismark --bowtie2 --non_directional --genome genome -1 reads_1.fq -2 reads_2.fq \
  -o . --temp_dir "$WORK/tmp" > align.log 2>&1
grep -E "Mapping efficiency|^(CT|GA)/(CT|GA)/(CT|GA):" reads_1_bismark_bt2_PE_report.txt

say "What Bismark wrote for each record, and which fastq read it really was"
samtools view -h -o out.sam reads_1_bismark_bt2_pe.bam
python3 - <<'PYEOF'
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s):
    return "".join(COMP[c] for c in reversed(s))

def fastq(path):
    lines = open(path).read().split("\n")
    return dict((lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 4))

r1, r2 = fastq("reads_1.fq"), fastq("reads_2.fq")
ref = "".join(l.strip() for l in open("genome/ref.fa") if not l.startswith(">"))

# The rule in queries/strand/bismark.toml, written out so the table below shows what it
# calls rather than asserting it.
RULE = {("CT", "CT"): "OT", ("GA", "CT"): "CTOT",
        ("CT", "GA"): "OB",  ("GA", "GA"): "CTOB"}
MATE = {"OT": "CTOT", "CTOT": "OT", "OB": "CTOB", "CTOB": "OB"}

rows = []
for line in open("out.sam"):
    if line.startswith("@"):
        continue
    f = line.rstrip("\n").split("\t")
    qname, flag, pos, seq = f[0], int(f[1]), int(f[3]), f[9]
    tags = dict((t[:2], t[5:]) for t in f[11:])
    # Which fastq read this record holds: SEQ is stored reference-forward, so the record
    # is either the read as sequenced or its reverse complement.
    is_r1 = seq in (r1[qname], rc(r1[qname]))
    # Read 1 of a pair was sequenced from the strand the pair is named for; read 2 from
    # the other strand of the same duplex. That is the call the rule has to produce.
    want = qname if is_r1 else MATE[qname]
    window = ref[pos - 1:pos - 1 + len(seq)]
    mismatches = sum(1 for r, s in zip(window, seq)
                     if r != s and (r, s) not in (("C", "T"), ("G", "A")))
    rows.append((qname, is_r1, flag, tags["XR"], tags["XG"], want,
                 RULE[(tags["XR"], tags["XG"])], pos, mismatches))

print("%-5s %-9s %-5s %-9s %-5s %-9s %-9s %-6s %s"
      % ("pair", "really is", "FLAG", "FLAG says", "XR/XG", "", "rule call",
         "POS", "mismatches"))
for qname, is_r1, flag, xr, xg, want, call, pos, mm in sorted(rows, key=lambda r: (r[0], -r[1])):
    print("%-5s %-9s %-5d %-9s %-5s %-9s %-9s %-6d %d"
          % (qname, "read 1" if is_r1 else "read 2", flag,
             "read 1" if flag & 0x40 else "read 2", xr + "/" + xg,
             "want " + want, call, pos, mm))
print()
swapped = [r for r in rows if r[1] != bool(r[2] & 0x40)]
wrong = [r for r in rows if r[5] != r[6]]
print("records whose 0x40/0x80 bit contradicts the fastq they came from: %d of %d (%s)"
      % (len(swapped), len(rows), ", ".join(sorted(set(r[0] for r in swapped))) or "none"))
print("records the rule calls as something other than the strand they were built from: %d"
      % len(wrong))
PYEOF

say "alnbase, reading the rule and nothing else"
"$ALNBASE" index genome/ref.fa ref.aref > /dev/null
"$ALNBASE" query --query-file "$HERE/origin.toml" \
  --query-file "$HERE/../../../queries/strand/bismark.toml" \
  reads_1_bismark_bt2_pe.bam ref.aref tagged.bam
printf -- '-- the YV tag each record came out with:\n'
samtools view tagged.bam \
  | awk '{yv="none"; for (i = 12; i <= NF; i++) if ($i ~ /^YV:Z:/) yv = substr($i, 6);
          printf "%-5s FLAG %-4d POS %-5d YV=%s\n", $1, $2, $4, yv}' \
  | sort
