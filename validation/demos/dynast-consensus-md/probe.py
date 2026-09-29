"""dynast's consensus caller writes MD and NM tags that do not describe its own CIGAR.

Runs dynast's real `call_consensus_from_reads` on hand-built reads -- no STAR, no
GTF, no real BAM needed -- and checks the tags it produces against `samtools calmd`,
which recomputes MD and NM from the alignment and the reference.

  DYNAST=<path to a dynast checkout> python probe.py

Needs pysam, numpy, ngs_tools, anndata (dynast's import chain) and samtools on PATH.
"""
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.environ.get("DYNAST", "."))

import pysam
from dynast.preprocessing.consensus import call_consensus_from_reads

REF = "ACGTACGTACGTACGTACGT"
HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6", "SO": "coordinate"},
     "SQ": [{"SN": "chr1", "LN": len(REF)}]}
)
WORK = tempfile.mkdtemp(prefix="dynast-consensus-md.")

with open(os.path.join(WORK, "ref.fa"), "w") as f:
    f.write(f">chr1\n{REF}\n")
subprocess.run(["samtools", "faidx", os.path.join(WORK, "ref.fa")], check=True)


def make_read(name, seq, cigar, md, nm, start=0):
    al = pysam.AlignedSegment(HEADER)
    al.query_name = name
    al.query_sequence = seq
    al.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    al.reference_name = "chr1"
    al.reference_id = 0
    al.reference_start = start
    al.mapping_quality = 60
    al.cigarstring = cigar
    al.set_tags([("MD", md), ("NM", nm)])
    return al


def calmd(al):
    """What samtools says MD and NM should be for this exact alignment."""
    bam = os.path.join(WORK, "one.bam")
    with pysam.AlignmentFile(bam, "wb", header=HEADER) as out:
        out.write(al)
    fixed = subprocess.run(
        ["samtools", "calmd", "-b", bam, os.path.join(WORK, "ref.fa")],
        check=True, capture_output=True,
    ).stdout
    path = os.path.join(WORK, "fixed.bam")
    with open(path, "wb") as f:
        f.write(fixed)
    with pysam.AlignmentFile(path, "rb") as f:
        r = next(iter(f))
    return r.get_tag("MD"), r.get_tag("NM")


def case(title, reads):
    print(f"--- {title}")
    for r in reads:
        print(f"  input      {r.reference_start:>3} {r.cigarstring:<10} "
              f"MD {str(r.get_tag('MD')):<12} NM {r.get_tag('NM')}")
    con = call_consensus_from_reads(reads, HEADER)
    md, nm = con.get_tag("MD"), con.get_tag("NM")
    want_md, want_nm = calmd(con)
    print(f"  dynast     {con.reference_start:>3} {con.cigarstring:<10} "
          f"MD {md:<12} NM {nm}")
    print(f"  samtools   {con.reference_start:>3} {con.cigarstring:<10} "
          f"MD {want_md:<12} NM {want_nm}")
    flags = []
    if md != want_md:
        flags.append("MD WRONG")
    if nm != want_nm:
        flags.append("NM WRONG")
    print(f"  -> {', '.join(flags) if flags else 'agrees'}\n")
    return con


print(f"reference chr1 = {REF}\n")

# Two identical reads per family: the consensus is unambiguous, so any
# disagreement with samtools is the tag writer's, not the base caller's.
def pair(seq, cigar, md, nm):
    return [make_read("a", seq, cigar, md, nm), make_read("b", seq, cigar, md, nm)]


# ref 5 deleted, ref 6 mismatched G->T.  Correct MD: 5^C0G13
case("deletion immediately followed by a mismatch",
     pair(REF[0:5] + "T" + REF[7:20], "5M1D14M", "5^C0G13", 2))

# ref 4 mismatched A->C, then ref 5 and ref 6 deleted.  Correct MD: 4A0^CG13
case("mismatch immediately followed by a 2 bp deletion",
     pair(REF[0:4] + "C" + REF[7:20], "5M2D13M", "4A0^CG13", 3))

# Control: a deletion with no adjacent mismatch.  Correct MD: 5^C14
case("isolated deletion (control)",
     pair(REF[0:5] + REF[6:20], "5M1D14M", "5^C14", 1))

# Control: a mismatch with no adjacent deletion.  Correct MD: 4A15
case("isolated mismatch (control)",
     pair(REF[0:4] + "C" + REF[5:20], "20M", "4A15", 1))

con = case("two non-overlapping reads in one UMI family",
           [make_read("a", REF[0:5], "5M", "5", 0),
            make_read("b", REF[12:20], "8M", "8", 0, start=12)])
print("  Read a covers ref 0-4, read b covers ref 12-19, and nothing covers the")
print("  7 bp between them. dynast emits that gap as an N operator -- a spliced")
print(f"  alignment. blocks = {con.get_blocks()}")
