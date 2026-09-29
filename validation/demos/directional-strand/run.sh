#!/usr/bin/env bash
# Run a directional bisulfite library through Bismark and check what
# queries/strand/directional.toml makes of the BAM it writes.
#
# directional.toml is the odd one out among the strand rules: every other file describes
# one aligner's tag, and this one reads no tag at all. It is a claim about the LIBRARY --
# that read 1 was sequenced from a converted original and read 2 from that original's PCR
# copy -- so it applies to any conformant aligner's output and to none of the evidence any
# of them writes. That makes it hard to check, because there is nothing in the BAM to
# compare it against.
#
# Bismark is the way in. Its XR/XG pair names all four strands per record, from the read's
# own conversion and the index that won, and queries/strand/bismark.toml is already checked
# against real output. So on a directional Bismark BAM there are two independent four-way
# calls for every record -- one from the FLAG alone, one from the tags -- and a third
# answer, the strand each read was actually built from, that the fixture knows. This demo
# puts all three side by side.
#
# The single-end arm is here for a specific reason. bismark.toml's comment says that on
# single-end output 0x10 is set from the genome conversion rather than from SEQ's
# orientation, which is why that file has to read the tags. directional.toml reads 0x10.
# Whether it is therefore wrong on single-end Bismark output is a question the rule's own
# comment does not answer: it lists "Bismark 0.25.1 directional paired-end" and stops.
#
#   BISMARK_ENV=/path/to/env PY=/path/to/python-with-pyarrow ALNBASE=/path/to/alnbase \
#     ./run.sh [workdir]
#
# BISMARK_ENV is the environment from ../../envs/bismark.yaml: Bismark 0.25.1, bowtie2 and
# samtools.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${BISMARK_ENV:?set BISMARK_ENV}" "${PY:?set PY}" "${ALNBASE:?set ALNBASE}"
export PATH="$BISMARK_ENV/bin:$PATH"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }

# The reference is the 1200 bases every strand demo uses: a fixed recurrence with a CG
# planted every 17, so every read carries several cytosines in CG context and the loci stay
# unique after either conversion.
#
# A directional library is the two fragments' own strands and nothing else: each fragment
# gives an OT/CTOT duplex or an OB/CTOB duplex, read 1 always from the converted original.
# Four probe pairs make a table small enough to read; a sweep of 60 bp pairs every 10 bases
# along the reference, from both strands, turns the same question into a count.
say "A directional library: read 1 from the original strand, read 2 from its copy"
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

MATE = {"OT": "CTOT", "OB": "CTOB"}

def pair(name, start, frag, length):
    """The two reads of one directional fragment, each from its own strand's 5' end."""
    f = ref[start:start + frag]
    ot, ob = bisulfite(f), bisulfite(rc(f))
    strand = {"OT": ot, "CTOT": rc(ot), "OB": ob, "CTOB": rc(ob)}
    return strand[name][:length], strand[MATE[name]][:length]

FRAG, LEN = 200, 60
r1f, r2f = open("reads_1.fq", "w"), open("reads_2.fq", "w")
for name, start in (("OT", 100), ("OB", 600)):
    r1, r2 = pair(name, start, FRAG, LEN)
    r1f.write("@probe_%s\n%s\n+\n%s\n" % (name, r1, "I" * LEN))
    r2f.write("@probe_%s\n%s\n+\n%s\n" % (name, r2, "I" * LEN))
    print("probe_%-3s fragment chr1 1-based %d-%d: read 1 from %s, read 2 from %s"
          % (name, start + 1, start + FRAG, name, MATE[name]))
n = 0
for start in range(0, 1000, 10):
    for name in ("OT", "OB"):
        r1, r2 = pair(name, start, FRAG, LEN)
        r1f.write("@sweep_%s_%d\n%s\n+\n%s\n" % (name, start, r1, "I" * LEN))
        r2f.write("@sweep_%s_%d\n%s\n+\n%s\n" % (name, start, r2, "I" * LEN))
        n += 1
r1f.close()
r2f.close()
print("plus a sweep of %d pairs, 60 bp from each end of a 200 bp fragment every 10 bases"
      % n)

# The single-end arm is read 1 of every pair: in a directional library that is the
# converted original itself, OT or OB and never a copy.
with open("reads_1.fq") as src, open("se.fq", "w") as dst:
    dst.write(src.read())
PYEOF
samtools faidx genome/ref.fa

say "Bismark, directional (the default), paired-end and single-end"
mkdir -p tmp
bismark_genome_preparation --bowtie2 genome > prep.log 2>&1
bismark --bowtie2 --genome genome -1 reads_1.fq -2 reads_2.fq \
  -o . --temp_dir "$WORK/tmp" > align_pe.log 2>&1
bismark --bowtie2 --genome genome se.fq \
  -o . --temp_dir "$WORK/tmp" > align_se.log 2>&1
grep -E "Mapping efficiency" reads_1_bismark_bt2_PE_report.txt se_bismark_bt2_SE_report.txt

say "Three answers per record: the fixture's, the FLAG's, and the tags'"
python3 - <<'PYEOF'
import subprocess

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s):
    return "".join(COMP[c] for c in reversed(s))

def fastq(path):
    lines = open(path).read().split("\n")
    return dict((lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 4))

ref = "".join(l.strip() for l in open("genome/ref.fa") if not l.startswith(">"))

# The two rules, written out so the tables below show what each one calls rather than
# asserting it.
def directional(flag):
    """queries/strand/directional.toml: the mate bit picks original or copy, 0x10 the strand.
    An unpaired read has is_last_in_template false and so is treated as read 1."""
    rev, last = bool(flag & 0x10), bool(flag & 0x80)
    return {(False, False): "OT", (True, True): "CTOT",
            (True, False): "OB", (False, True): "CTOB"}[(rev, last)]

def bismark_rule(tags):
    """queries/strand/bismark.toml: XR is the read's own conversion, XG the index that won."""
    return {("CT", "CT"): "OT", ("GA", "CT"): "CTOT",
            ("CT", "GA"): "OB", ("GA", "GA"): "CTOB"}[(tags["XR"], tags["XG"])]

MATE = {"OT": "CTOT", "OB": "CTOB"}

def records(bam, reads):
    out = []
    sam = subprocess.run(["samtools", "view", bam], capture_output=True, text=True).stdout
    for line in sam.rstrip("\n").split("\n"):
        if not line:
            continue
        f = line.split("\t")
        qname, flag, pos, seq = f[0], int(f[1]), int(f[3]), f[9]
        tags = dict((t[:2], t[5:]) for t in f[11:])
        # Which fastq read this record holds, and whether SEQ was reverse-complemented on
        # the way in: both settled against the reads themselves, so nothing in the "really"
        # columns comes from the FLAG or from a tag.
        as_read = reads[0][qname] if seq in (reads[0][qname], rc(reads[0][qname])) \
            else reads[1][qname]
        is_r1 = as_read == reads[0][qname]
        # A pair is named for the strand its read 1 came from; read 2 is that strand's copy.
        pair_strand = qname.split("_")[1]
        truth = pair_strand if is_r1 else MATE[pair_strand]
        window = ref[pos - 1:pos - 1 + len(seq)]
        out.append(dict(
            qname=qname, flag=flag, pos=pos, is_r1=is_r1, truth=truth,
            xr=tags["XR"], xg=tags["XG"],
            revcomped=seq != as_read,
            flag_call=directional(flag), tag_call=bismark_rule(tags),
            mismatches=sum(1 for r, s in zip(window, seq)
                           if r != s and (r, s) not in (("C", "T"), ("G", "A")))))
    return out

def report(label, rows):
    print("\n-- %s: %d records" % (label, len(rows)))
    probes = [r for r in rows if r["qname"].startswith("probe_")]
    print("   %-9s %-9s %-5s %-5s %-6s %-11s %-9s %-6s %s"
          % ("read", "really is", "FLAG", "XR/XG", "truth", "directional", "bismark",
             "POS", "mism"))
    for r in sorted(probes, key=lambda r: (r["qname"], not r["is_r1"])):
        print("   %-9s %-9s %-5d %-5s %-6s %-11s %-9s %-6d %d"
              % (r["qname"], "read 1" if r["is_r1"] else "read 2", r["flag"],
                 r["xr"] + "/" + r["xg"], r["truth"], r["flag_call"], r["tag_call"],
                 r["pos"], r["mismatches"]))
    wrong_f = [r for r in rows if r["flag_call"] != r["truth"]]
    wrong_t = [r for r in rows if r["tag_call"] != r["truth"]]
    print("   directional.toml calls a strand the read was not built from: %d of %d"
          % (len(wrong_f), len(rows)))
    print("   bismark.toml calls a strand the read was not built from:     %d of %d"
          % (len(wrong_t), len(rows)))
    print("   records the two rules disagree about: %d"
          % len([r for r in rows if r["flag_call"] != r["tag_call"]]))
    # The FLAG claim directional.toml rests on, and the one bismark.toml's comment says
    # fails on single-end output: whether 0x10 tracks how SEQ was really stored.
    print("   records whose 0x10 contradicts how SEQ was actually stored: %d"
          % len([r for r in rows if r["revcomped"] != bool(r["flag"] & 0x10)]))
    print("   records not stored reference-forward: %d"
          % len([r for r in rows if r["mismatches"]]))

pe = fastq("reads_1.fq"), fastq("reads_2.fq")
se = fastq("se.fq"), {}
report("paired-end", records("reads_1_bismark_bt2_pe.bam", pe))
report("single-end", records("se_bismark_bt2.bam", se))
PYEOF

say "alnbase, reading each rule and nothing else"
"$ALNBASE" index genome/ref.fa ref.aref > /dev/null
for arm in pe se; do
  bam=$([ "$arm" = pe ] && echo reads_1_bismark_bt2_pe.bam || echo se_bismark_bt2.bam)
  for rule in directional bismark; do
    printf -- '-- %s, %s.toml\n' "$arm" "$rule"
    "$ALNBASE" query --parquet --only-hits \
      -f qname,flags,strand,conv_strand \
      --query-file "$HERE/cg.toml" \
      --query-file "$HERE/../../../queries/strand/$rule.toml" \
      "$bam" ref.aref "$arm.$rule.hits.parquet"
  done
done

# Every hit's anchor is a cytosine in the read; which reference base it sits on is decided
# by the conversion strand alnbase took from the rule -- the C of a forward CG for a top
# strand, the G of that same CG for a bottom one. Reading those bases back out of the
# reference checks the walk rather than the label.
"$PY" - <<'PYEOF'
import glob
import pyarrow.parquet as pq

ref = "".join(l.strip() for l in open("genome/ref.fa") if not l.startswith(">"))
TOP = {"OT", "CTOT"}

def load(arm, rule):
    return pq.read_table(sorted(glob.glob("%s.%s.hits_*.parquet" % (arm, rule)))).to_pylist()

for arm in ("pe", "se"):
    print("\n-- %s" % arm)
    both = {}
    for rule in ("directional", "bismark"):
        rows = both[rule] = load(arm, rule)
        wrong = [r for r in rows
                 if ref[r["refr_pos"]] != ("C" if r["conv_strand"] in TOP else "G")]
        print("   %-16s %6d CG hits over %4d records; anchors on the wrong base of the CG "
              "for the strand named: %d"
              % (rule + ".toml", len(rows),
                 len(set((r["qname"], r["flags"]) for r in rows)), len(wrong)))
    # Two rules agree when they put the same anchors in the same places, keyed on the record
    # and the reference position, so a strand disagreement moves the anchor from a C to the G
    # beside it and shows up as a difference.
    key = lambda rows: sorted((r["qname"], r["flags"], r["refr_pos"], r["conv_strand"])
                              for r in rows)
    a, b = key(both["directional"]), key(both["bismark"])
    print("   hits where the two rules disagree: %d (of %d and %d)"
          % (len(set(a) ^ set(b)), len(a), len(b)))
PYEOF
