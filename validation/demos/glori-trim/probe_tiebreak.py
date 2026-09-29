"""Which duplicate alignment's A-count does GLORI keep when two records share a read name?

GLORI's final BAM is a concatenation of a genome alignment and a transcriptome alignment
lifted back to genome coordinates (`run_GLORI.py:50, :71-73`), so one read can appear twice
at the same column. `pileup_genome_multiprocessing.py` deduplicates by read name and, when
the two records agree on the base, is meant to keep the larger per-read A-count -- the
conservative choice, because that count is what `--cutoff` filters on downstream.

It compares the stored count to itself instead (`:91` on the reverse strand, `:116` on the
forward), so the comparison is constant and the winner is decided by write order and strand
rather than by the counts. Each strand is probed twice here, once with the low count written
first and once with it written second; if the code were choosing on magnitude, both orders
would give the same answer.

A forward read's count is `query_sequence.count('A')` and a reverse read's is `.count('T')`,
which is the same base in the read as sequenced, so the reverse pair is built from T.

  GLORI=<checkout> python probe_tiebreak.py
"""
import os
import subprocess
import sys

import pysam

GLORI = os.environ.get("GLORI")
if not GLORI:
    sys.exit("set GLORI to a GLORI-tools checkout")

REF = "A" * 40
HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6", "SO": "coordinate"},
     "SQ": [{"SN": "chr1", "LN": len(REF)}]}
)

# Two records per read name, agreeing on the base at read offset 19 (reference position 20)
# but carrying 0 and 18 countable bases respectively.
LOW = {"fwd": "G" * 20, "rev": "G" * 20}
HIGH = {"fwd": "A" * 18 + "GG", "rev": "T" * 18 + "GG"}

cases = []
for strand, flag in (("fwd", 0), ("rev", 16)):
    cases.append((f"{strand}_low_first", flag, [LOW[strand], HIGH[strand]]))
    cases.append((f"{strand}_high_first", flag, [HIGH[strand], LOW[strand]]))

with pysam.AlignmentFile("tiebreak.bam", "wb", header=HEADER) as out:
    for name, flag, seqs in cases:
        for seq in seqs:
            assert seq[19] == "G"
            al = pysam.AlignedSegment(HEADER)
            al.query_name = name            # both records share one name
            al.query_sequence = seq
            al.query_qualities = pysam.qualitystring_to_array("I" * 20)
            al.reference_name = "chr1"
            al.reference_id = 0
            al.reference_start = 0
            al.mapping_quality = 60
            al.cigarstring = "20M"
            al.flag = flag
            out.write(al)
pysam.index("tiebreak.bam")

with open("ref40.fa", "w") as f:
    f.write(f">chr1\n{REF}\n")

subprocess.run([sys.executable, f"{GLORI}/pipelines/pileup_genome_multiprocessing.py",
                "-i", "tiebreak.bam", "-f", "ref40.fa", "-o", "tiebreak.pileup.txt"],
               check=True, capture_output=True)

names = {"+": ["fwd_low_first", "fwd_high_first"],
         "-": ["rev_low_first", "rev_high_first"]}
print("  each read name has two records with A-counts 0 and 18; the larger should win\n")
for line in open("tiebreak.pileup.txt"):
    chrom, pos, strand, ref, bases, acounts = line.rstrip("\n").split("\t")
    if pos != "20":
        continue
    for name, kept in zip(names[strand], acounts.split(",")):
        written = [0, 18] if name.endswith("low_first") else [18, 0]
        which = "first record" if int(kept) == written[0] else "second record"
        verdict = "ok" if int(kept) == 18 else "LOST THE LARGER COUNT"
        print(f"  {name:<16} written {written[0]:>2},{written[1]:>2}  "
              f"kept {kept:>2} ({which})  {verdict}")
