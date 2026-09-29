#!/usr/bin/env bash
# Run BS-Seeker2 on one read per strand of origin and check what
# queries/strand/bs-seeker2-se.toml makes of the BAM it writes.
#
#   BSSEEKER2_ENV=/path/to/bs-seeker2-env ALNBASE=/path/to/alnbase ./run.sh [workdir]
#
# BSSEEKER2_ENV is the environment from ../../envs/bs-seeker2.yaml with BS-Seeker2 checked
# out into it by ../../envs/build_bsseeker2.sh (commit 0976fee, v2.1.8). It provides
# python 2.7, pysam < 0.9, bowtie2 and samtools.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${BSSEEKER2_ENV:?set BSSEEKER2_ENV}" "${ALNBASE:?set ALNBASE}"
export PATH="$BSSEEKER2_ENV/bin:$PATH"
BSS2="$BSSEEKER2_ENV/share/BSseeker2"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }

# The reference is 800 bases from a fixed recurrence with a CG planted every 17, so every
# read carries several cytosines in CG context and the four loci are unique after either
# conversion. The four reads are one per strand of origin, each named for the strand it was
# built from -- which is the answer the rule has to produce.
say "The four reads, each named for the strand it was built from"
python - <<'PYEOF'
bases = "ACGT"
ref, x = [], 7
for _ in range(800):
    x = (x * 31 + 17) % 1009
    ref.append(bases[x % 4])
for i in range(0, 798, 17):
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

# OT is the top strand as bisulfite leaves it; CTOT is the strand copied from it, and so
# the reverse complement of that sequence. OB and CTOB are the same pair on the bottom.
LEN = 60
loci = [("OT", 100), ("CTOT", 250), ("OB", 400), ("CTOB", 550)]
with open("reads.fq", "w") as f:
    for name, start in loci:
        top = ref[start:start + LEN]
        seq = {"OT":   bisulfite(top),
               "CTOT": rc(bisulfite(top)),
               "OB":   bisulfite(rc(top)),
               "CTOB": rc(bisulfite(rc(top)))}[name]
        f.write("@%s\n%s\n+\n%s\n" % (name, seq, "I" * LEN))
        print("%-5s from chr1 1-based %d-%d  %s" % (name, start + 1, start + LEN, seq))
PYEOF
samtools faidx ref.fa

say "BS-Seeker2, single-end, -t Y so all four strands are aligned"
mkdir -p db tmp
python "$BSS2/bs_seeker2-build.py" -f ref.fa --aligner=bowtie2 -p "$BSSEEKER2_ENV/bin/" \
  -d db > build.log 2>&1
python "$BSS2/bs_seeker2-align.py" -i reads.fq -g ref.fa -d db --aligner=bowtie2 \
  -p "$BSSEEKER2_ENV/bin/" -f sam -o out.sam -m 4 -t Y --temp_dir="$WORK/tmp" \
  > align.log 2>&1
grep -E "mapped to (Watson|Crick) strand$|Mappability" align.log | sed 's/^[^]]*] *//'

say "What BS-Seeker2 wrote, and what a conformant aligner would have written"
python - <<'PYEOF'
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s):
    return "".join(COMP[c] for c in reversed(s))

lines = open("reads.fq").read().split("\n")
read = dict((lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 4))
ref = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))

print("%-5s %-4s %-5s %-5s %-6s %-18s %s"
      % ("qname", "XO", "POS", "FLAG", "should", "SEQ is", "mismatches vs reference"))
for line in open("out.sam"):
    if line.startswith("@"):
        continue
    f = line.rstrip("\n").split("\t")
    qname, flag, pos, seq = f[0], int(f[1]), int(f[3]), f[9]
    xo = [t for t in f[11:] if t.startswith("XO:Z:")][0][5:]
    # A read is reverse relative to the reference exactly when the aligner had to
    # reverse-complement it to place it, which is what 0x10 is for.
    reversed_ = seq == rc(read[qname])
    # And SEQ is reference-forward if it matches the reference at POS once the two
    # substitutions bisulfite can make -- C>T on the top strand, G>A on the bottom, both
    # seen this way round because SEQ is stored reference-forward -- are allowed.
    window = ref[pos - 1:pos - 1 + len(seq)]
    mismatches = sum(1 for r, s in zip(window, seq)
                     if r != s and (r, s) not in (("C", "T"), ("G", "A")))
    print("%-5s %-4s %-5d %-5d %-6d %-18s %d"
          % (qname, xo, pos, flag, 16 if reversed_ else 0,
             "revcomp(read)" if reversed_ else "the read as given", mismatches))
PYEOF

say "alnbase, reading the rule and nothing else"
samtools view -b --no-PG -o aln.bam out.sam
"$ALNBASE" index ref.fa ref.aref > /dev/null
"$ALNBASE" query --query-file "$HERE/origin.toml" \
  --query-file "$HERE/../../../queries/strand/bs-seeker2-se.toml" aln.bam ref.aref tagged.bam
printf -- '-- the YV tag each record came out with:\n'
samtools view tagged.bam \
  | awk '{yv="none"; for (i = 12; i <= NF; i++) if ($i ~ /^YV:Z:/) yv = substr($i, 6);
          printf "%-5s FLAG %-3d POS %-4d YV=%s\n", $1, $2, $4, yv}' \
  | sort
