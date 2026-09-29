#!/usr/bin/env bash
# Run a PBAT bisulfite library through Bismark and check what a PBAT strand rule makes of the
# BAM it writes.
#
# The rule under test is directional.toml with the two conversion strands exchanged: in a PBAT
# library read 1 comes off the copy strand and read 2 off the converted original, the other way
# round from a directional library. Like directional.toml it reads no tag, so it is a claim
# about the LIBRARY and there is nothing in the BAM to compare it against. Bismark's XR/XG
# names all four strands from evidence the FLAG rule never touches, and the fixture knows which
# strand every read was built from, so each record gets three answers.
#
# This round asked two things the directional round could not:
#
#   1. Nothing in the BAM says which library this was. The same records are run under
#      directional.toml as well, to measure what pointing at the wrong file costs.
#   2. Single-end output is where Bismark's non-conformant 0x10 finally shows. On a
#      directional library the bit Bismark sets from the genome conversion and the bit SAM
#      asks for coincide, which is why the directional round measured no contradiction. A
#      PBAT library is exactly the case where they part company.
#
# Both answers came back against the rule, and it is no longer shipped. Bismark's paired-end
# FLAG table swaps 0x40 and 0x80 on CTOT and CTOB pairs (bismark:8818-8866), and under --pbat
# every pair is one of those, so the output is in the directional arrangement: the PBAT rule is
# wrong on every record and directional.toml is right on every record. Single-end, 0x10
# contradicts how SEQ was stored on every record, and only bismark.toml gets the strand right.
# queries/strand/ therefore has no PBAT file; bismark.toml covers both --pbat layouts. The rule
# is written inline below so this comparison stays reproducible. See README.md.
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
# The reads are the directional demo's reads with the two mate files exchanged, which is
# literally what methylpy --pbat does to a PBAT run before aligning it (call_mc_pe.py:225-227):
# the same four strands exist, and which of each duplex is sequenced first is the whole
# difference between the two libraries.
say "A PBAT library: read 1 from the copy strand, read 2 from the converted original"
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

# A PBAT pair is named for the strand its read 1 came from -- a copy strand -- and read 2 is
# the converted original that copy was made from.
ORIGINAL = {"CTOT": "OT", "CTOB": "OB"}

def pair(name, start, frag, length):
    """The two reads of one PBAT fragment, each from its own strand's 5' end."""
    f = ref[start:start + frag]
    ot, ob = bisulfite(f), bisulfite(rc(f))
    strand = {"OT": ot, "CTOT": rc(ot), "OB": ob, "CTOB": rc(ob)}
    return strand[name][:length], strand[ORIGINAL[name]][:length]

FRAG, LEN = 200, 60
r1f, r2f = open("reads_1.fq", "w"), open("reads_2.fq", "w")
for name, start in (("CTOT", 100), ("CTOB", 600)):
    r1, r2 = pair(name, start, FRAG, LEN)
    r1f.write("@probe_%s\n%s\n+\n%s\n" % (name, r1, "I" * LEN))
    r2f.write("@probe_%s\n%s\n+\n%s\n" % (name, r2, "I" * LEN))
    print("probe_%-4s fragment chr1 1-based %d-%d: read 1 from %s, read 2 from %s"
          % (name, start + 1, start + FRAG, name, ORIGINAL[name]))
n = 0
for start in range(0, 1000, 10):
    for name in ("CTOT", "CTOB"):
        r1, r2 = pair(name, start, FRAG, LEN)
        r1f.write("@sweep_%s_%d\n%s\n+\n%s\n" % (name, start, r1, "I" * LEN))
        r2f.write("@sweep_%s_%d\n%s\n+\n%s\n" % (name, start, r2, "I" * LEN))
        n += 1
r1f.close()
r2f.close()
print("plus a sweep of %d pairs, 60 bp from each end of a 200 bp fragment every 10 bases"
      % n)

# The single-end arm is read 1 of every pair: in a PBAT library that is the copy strand
# itself, CTOT or CTOB and never an original.
with open("reads_1.fq") as src, open("se.fq", "w") as dst:
    dst.write(src.read())
PYEOF
samtools faidx genome/ref.fa

say "Bismark --pbat, paired-end and single-end"
mkdir -p tmp
bismark_genome_preparation --bowtie2 genome > prep.log 2>&1
bismark --bowtie2 --pbat --genome genome -1 reads_1.fq -2 reads_2.fq \
  -o . --temp_dir "$WORK/tmp" > align_pe.log 2>&1
bismark --bowtie2 --pbat --genome genome se.fq \
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

# The three rules, written out so the tables below show what each one calls rather than
# asserting it.
def pbat(flag):
    """The PBAT rule: read 1 is the copy strand, read 2 the converted original."""
    rev, last = bool(flag & 0x10), bool(flag & 0x80)
    return {(True, False): "CTOT", (False, False): "CTOB",
            (False, True): "OT", (True, True): "OB"}[(rev, last)]

def directional(flag):
    """queries/strand/directional.toml: the same rule with the conversion strands exchanged.
    It is here as the wrong file for this BAM, since nothing in the BAM tells them apart."""
    rev, last = bool(flag & 0x10), bool(flag & 0x80)
    return {(False, False): "OT", (True, True): "CTOT",
            (True, False): "OB", (False, True): "CTOB"}[(rev, last)]

def bismark_rule(tags):
    """queries/strand/bismark.toml: XR is the read's own conversion, XG the index that won."""
    return {("CT", "CT"): "OT", ("GA", "CT"): "CTOT",
            ("CT", "GA"): "OB", ("GA", "GA"): "CTOB"}[(tags["XR"], tags["XG"])]

ORIGINAL = {"CTOT": "OT", "CTOB": "OB"}

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
        # A PBAT pair is named for the copy strand its read 1 came from; read 2 is the
        # converted original.
        pair_strand = qname.split("_")[1]
        truth = pair_strand if is_r1 else ORIGINAL[pair_strand]
        window = ref[pos - 1:pos - 1 + len(seq)]
        out.append(dict(
            qname=qname, flag=flag, pos=pos, is_r1=is_r1, truth=truth,
            xr=tags["XR"], xg=tags["XG"],
            revcomped=seq != as_read,
            pbat_call=pbat(flag), dir_call=directional(flag), tag_call=bismark_rule(tags),
            mismatches=sum(1 for r, s in zip(window, seq)
                           if r != s and (r, s) not in (("C", "T"), ("G", "A")))))
    return out

def report(label, rows):
    print("\n-- %s: %d records" % (label, len(rows)))
    probes = [r for r in rows if r["qname"].startswith("probe_")]
    print("   %-10s %-9s %-5s %-5s %-6s %-6s %-11s %-8s %-6s %s"
          % ("read", "really is", "FLAG", "XR/XG", "truth", "pbat", "directional", "bismark",
             "POS", "mism"))
    for r in sorted(probes, key=lambda r: (r["qname"], not r["is_r1"])):
        print("   %-10s %-9s %-5d %-5s %-6s %-6s %-11s %-8s %-6d %d"
              % (r["qname"], "read 1" if r["is_r1"] else "read 2", r["flag"],
                 r["xr"] + "/" + r["xg"], r["truth"], r["pbat_call"], r["dir_call"],
                 r["tag_call"], r["pos"], r["mismatches"]))
    for name, key in (("pbat.toml       ", "pbat_call"),
                      ("directional.toml", "dir_call"),
                      ("bismark.toml    ", "tag_call")):
        print("   %s calls a strand the read was not built from: %d of %d"
              % (name, len([r for r in rows if r[key] != r["truth"]]), len(rows)))
    print("   records pbat.toml and bismark.toml disagree about: %d"
          % len([r for r in rows if r["pbat_call"] != r["tag_call"]]))
    # The FLAG claim pbat.toml rests on, and the one its comment says fails on single-end
    # output: whether 0x10 tracks how SEQ was really stored.
    print("   records whose 0x10 contradicts how SEQ was actually stored: %d"
          % len([r for r in rows if r["revcomped"] != bool(r["flag"] & 0x10)]))
    print("   records not stored reference-forward: %d"
          % len([r for r in rows if r["mismatches"]]))

pe = fastq("reads_1.fq"), fastq("reads_2.fq")
se = fastq("se.fq"), {}
report("paired-end", records("reads_1_bismark_bt2_pe.bam", pe))
report("single-end", records("se_bismark_bt2.bam", se))
PYEOF

# The PBAT rule itself is written here rather than taken from queries/strand/, because this
# run is why it is no longer shipped: it is the honest reading of a PBAT library under SAM's
# meaning of the flags, and no aligner produces a BAM it fits. Keeping the copy here keeps
# the comparison reproducible.
cat > pbat.toml <<'TOMLEOF'
[strand.origin]
CTOT = "is_reverse and not is_last_in_template"
CTOB = "not is_reverse and not is_last_in_template"
OT   = "not is_reverse and is_last_in_template"
OB   = "is_reverse and is_last_in_template"
TOMLEOF

say "alnbase, reading each rule and nothing else"
"$ALNBASE" index genome/ref.fa ref.aref > /dev/null
for arm in pe se; do
  bam=$([ "$arm" = pe ] && echo reads_1_bismark_bt2_pe.bam || echo se_bismark_bt2.bam)
  for rule in pbat directional bismark; do
    file=$([ "$rule" = pbat ] && echo pbat.toml || echo "$HERE/../../../queries/strand/$rule.toml")
    printf -- '-- %s, %s.toml\n' "$arm" "$rule"
    "$ALNBASE" query --parquet --only-hits \
      -f qname,flags,strand,conv_strand \
      --query-file "$HERE/cg.toml" \
      --query-file "$file" \
      "$bam" ref.aref "$arm.$rule.hits.parquet"
  done
done

# Every hit's anchor is a cytosine in the read; which reference base it sits on is decided
# by the conversion strand alnbase took from the rule -- the C of a forward CG for a top
# strand, the G of that same CG for a bottom one. Reading those bases back out of the
# reference checks the walk rather than the label, and it is where a rule that named the
# opposite conversion strand stops being a labelling mistake and becomes a wrong number.
"$PY" - <<'PYEOF'
import glob
import pyarrow.parquet as pq

ref = "".join(l.strip() for l in open("genome/ref.fa") if not l.startswith(">"))
TOP = {"OT", "CTOT"}

def load(arm, rule):
    return pq.read_table(sorted(glob.glob("%s.%s.hits_*.parquet" % (arm, rule)))).to_pylist()

key = lambda rows: sorted((r["qname"], r["flags"], r["refr_pos"], r["conv_strand"])
                          for r in rows)
place = lambda rows: sorted((r["qname"], r["flags"], r["refr_pos"]) for r in rows)

for arm in ("pe", "se"):
    print("\n-- %s" % arm)
    both = {}
    for rule in ("pbat", "directional", "bismark"):
        rows = both[rule] = load(arm, rule)
        wrong = [r for r in rows
                 if ref[r["refr_pos"]] != ("C" if r["conv_strand"] in TOP else "G")]
        # This is a consistency check on alnbase, not a check on the rule: it asks whether
        # the walk went the way the label says, so a rule that names the wrong strand and
        # then walks it consistently still scores 0. What catches a wrong rule is the
        # mismatch percentage above and the comparison with bismark.toml below.
        print("   %-16s %6d CG hits over %4d records; anchors on the base the strand named "
              "does not want: %d"
              % (rule + ".toml", len(rows),
                 len(set((r["qname"], r["flags"]) for r in rows)), len(wrong)))
    # Two kinds of disagreement, and they cost different things. A rule that names the other
    # conversion strand walks the read the other way, so its anchors move off the C onto the
    # G beside it -- a wrong number. A rule that gets the conversion strand right and the
    # strand of origin wrong puts every anchor in the same place and only mislabels it, which
    # still mirrors off_5p and off_3p but leaves the positions alone.
    for rule in ("pbat", "directional"):
        a, b = key(both[rule]), key(both["bismark"])
        pa, pb = place(both[rule]), place(both["bismark"])
        print("   %-16s vs bismark.toml: %5d hits anchored elsewhere, %5d more named a "
              "different strand in the same place (of %d and %d)"
              % (rule + ".toml", len(set(pa) ^ set(pb)),
                 len(set(a) ^ set(b)) - len(set(pa) ^ set(pb)), len(a), len(b)))
PYEOF
