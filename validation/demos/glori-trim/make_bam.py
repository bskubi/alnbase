"""Build a tiny BAM whose reads differ only in orientation, for GLORI's read-end trim.

Twenty 20 bp reads cover the same 20 bp reference exactly: ten on the forward strand and
ten on the reverse. Every read carries the same base at every position, so any difference
between the strands in GLORI's pileup output comes from the trim, not from the data.

The reference is all A except position 0 and position 19, which is where the 5' and 3'
trims are supposed to bite. Reads carry a G (a converted A) everywhere, so each column is
a clean count of "how many reads survived the trim here".
"""
import pysam

REF = "A" * 20
READ = "G" * 20

header = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6", "SO": "coordinate"},
     "SQ": [{"SN": "chr1", "LN": len(REF)}]}
)

with open("ref.fa", "w") as f:
    f.write(f">chr1\n{REF}\n")

with pysam.AlignmentFile("reads.bam", "wb", header=header) as out:
    for i in range(10):
        for strand in ("fwd", "rev"):
            al = pysam.AlignedSegment(header)
            al.query_name = f"{strand}{i}"
            al.query_sequence = READ
            al.query_qualities = pysam.qualitystring_to_array("I" * len(READ))
            al.reference_name = "chr1"
            al.reference_id = 0
            al.reference_start = 0
            al.mapping_quality = 60
            al.cigarstring = f"{len(READ)}M"
            al.flag = 16 if strand == "rev" else 0
            out.write(al)
pysam.index("reads.bam")
print("reads.bam: 10 forward + 10 reverse reads, all 20M at chr1:1-20")
