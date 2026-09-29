"""What pysam does with the MD string HISAT-3N's repeat-index builder constructs.

`constructRepeatMD` (`alignment_3n.h:504-508`) appends `^` for a deletion without first
flushing the pending match run length, so the run before the deletion is dropped and merged
into the run after it. This probe does not run HISAT-3N: the repeat path only fires for
alignments to a repeat index, and `hisat-3n-build --repeat-index` will not build one for a
sequence with no repeats. It writes the MD string that branch constructs for one alignment
and asks what pysam and `samtools calmd` make of it.

  python probe_repeat.py

Needs pysam and samtools on PATH.
"""
import os
import subprocess
import tempfile

import pysam

REF = "ACGTACGTACGTACGTACGT"
HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6", "SO": "coordinate"},
     "SQ": [{"SN": "chr1", "LN": len(REF)}]}
)
WORK = tempfile.mkdtemp(prefix="hisat3n-repeat-md.")
with open(os.path.join(WORK, "ref.fa"), "w") as f:
    f.write(f">chr1\n{REF}\n")
subprocess.run(["samtools", "faidx", os.path.join(WORK, "ref.fa")], check=True)


def make_read(name, seq, cigar, md, nm):
    al = pysam.AlignedSegment(HEADER)
    al.query_name = name
    al.query_sequence = seq
    al.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    al.reference_name = "chr1"
    al.reference_id = 0
    al.reference_start = 0
    al.mapping_quality = 60
    al.cigarstring = cigar
    al.set_tags([("MD", md), ("NM", nm)])
    return al


def case(title, al):
    print(f"  --- {title}")
    print(f"    CIGAR {al.cigarstring}  MD {al.get_tag('MD')}")
    pairs = al.get_aligned_pairs(matches_only=True, with_seq=True)
    got = "".join(p[2] for p in pairs)
    want = "".join(REF[p[1]] for p in pairs)
    print(f"    pysam reads the reference as  {got}")
    print(f"    the reference actually is     {want}")
    print(f"    -> {'agrees' if got.upper() == want else 'MISREAD'}")
    bam = os.path.join(WORK, "one.bam")
    with pysam.AlignmentFile(bam, "wb", header=HEADER) as out:
        out.write(al)
    done = subprocess.run(["samtools", "calmd", "-b", bam, os.path.join(WORK, "ref.fa")],
                          check=True, capture_output=True)
    path = os.path.join(WORK, "fixed.bam")
    with open(path, "wb") as f:
        f.write(done.stdout)
    with pysam.AlignmentFile(path, "rb") as f:
        fixed = next(iter(f))
    print(f"    samtools calmd says MD {fixed.get_tag('MD')}")
    if done.stderr:
        print(f"    calmd stderr: {done.stderr.decode().strip()}")
    print()


# 5M2D13M with no mismatch: ref 5 and 6 (C, G) are deleted.
seq = REF[0:5] + REF[7:20]
case("correct MD, for comparison", make_read("a", seq, "5M2D13M", "5^CG13", 2))
case("what constructRepeatMD builds", make_read("b", seq, "5M2D13M", "^CG18", 2))
