#!/usr/bin/env python3
"""Two 60 bp single-end BAMs for probing bam2bakR's mut_call.py filters.

The reference is the repeating tetramer ACGT, so a read starting at reference
index 20 has a T at every read offset j with j % 4 == 3: offsets 3, 7, ... 59.
That makes it trivial to place mutations at chosen offsets and to know exactly
how many T's the read contains.

  qual.bam  one T->C at Phred 10 and one at Phred 5, both well inside the read.
            Probes what --minQual actually compares against.
  ends.bam  EVERY T mutated to C, all at Phred 35. The true per-T conversion
            rate is 1.0 by construction, so any reported rate below 1.0 is the
            filter mismatch between the numerator and the denominator.
"""
import os
import sys

import pysam

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
L = 60
START = 20
REF = "".join("ACGT"[i % 4] for i in range(200))

os.makedirs(OUT, exist_ok=True)
with open(f"{OUT}/ref.fa", "w") as fh:
    fh.write(">chr1\n" + REF + "\n")

HEADER = {"HD": {"VN": "1.6", "SO": "queryname"},
          "SQ": [{"SN": "chr1", "LN": len(REF)}]}


def md_and_nm(seq):
    """MD/NM for an all-M alignment, so the demo needs no samtools calmd.

    mut_call.py calls pysam's get_aligned_pairs(with_seq=True), which raises
    ValueError unless the record carries an MD tag.
    """
    md, run, nm = "", 0, 0
    for j, base in enumerate(seq):
        ref = REF[START + j]
        if base == ref:
            run += 1
        else:
            md += f"{run}{ref}"
            run = 0
            nm += 1
    return md + str(run), nm


def write(path, seq, quals):
    md, nm = md_and_nm(seq)
    with pysam.AlignmentFile(path, "wb", header=HEADER) as out:
        a = pysam.AlignedSegment()
        a.query_name = "r1"
        a.flag = 0
        a.reference_id = 0
        a.reference_start = START
        a.mapping_quality = 60
        a.cigarstring = f"{L}M"
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array(
            "".join(chr(q + 33) for q in quals))
        a.set_tag("NM", nm, "i")
        a.set_tag("MD", md, "Z")
        out.write(a)


t_offsets = [j for j in range(L) if REF[START + j] == "T"]

# --- qual.bam: two mutations at different base qualities, both interior -----
seq = list(REF[START:START + L])
hi, lo = t_offsets[2], t_offsets[3]          # offsets 11 and 15
seq[hi] = seq[lo] = "C"
quals = [35] * L
quals[hi], quals[lo] = 10, 5
write(f"{OUT}/qual.bam", "".join(seq), quals)

# --- ends.bam: every T mutated, uniform high quality -----------------------
seq = list(REF[START:START + L])
for j in t_offsets:
    seq[j] = "C"
write(f"{OUT}/ends.bam", "".join(seq), [35] * L)

near_end = [j for j in t_offsets if j < 5 or j > L - 6]
print(f"read length {L}, T at offsets {t_offsets}")
print(f"qual.bam: T->C at offset {hi} (Phred 10) and offset {lo} (Phred 5)")
print(f"ends.bam: all {len(t_offsets)} T's mutated; "
      f"{len(near_end)} of them at offsets {near_end}, inside the "
      f"--minDist 5 end zones")
