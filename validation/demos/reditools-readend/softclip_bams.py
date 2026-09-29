"""Same reads, byte-identical SEQ, two CIGARs. Does -mbp count from the read's
own 5' end, or from the first aligned base?

Reference is all A. Every read carries a G at read offset 2 (0-based), i.e. the
3rd base of the read as sequenced -- an A>G "edit". CLIP=2, so in the clipped
file that base is the FIRST aligned base.

  -mbp 0  : the edit is counted in both files.
  -mbp 2  : ignores the first 2 bases of each read. The edit sits at read
            offset 2, which is the 3rd base, so it SURVIVES in both files if
            the filter measures the read as sequenced.
"""
import pysam, pathlib, sys

REF_LEN, RLEN, NREADS, CLIP = 3000, 40, 300, 2
EDIT_OFF = 2                      # 0-based offset into the read as sequenced
START = 1000

d = pathlib.Path(sys.argv[1]); d.mkdir(parents=True, exist_ok=True)
ref = "A" * REF_LEN
(d / "ref.fa").write_text(">chr1\n" + "\n".join(ref[i:i+60] for i in range(0, REF_LEN, 60)) + "\n")
pysam.faidx(str(d / "ref.fa"))

seq = list("A" * RLEN)
seq[EDIT_OFF] = "G"
seq = "".join(seq)

hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": REF_LEN}]}
for name, cigar, pos in (("endtoend", [(0, RLEN)], START),
                         ("clipped", [(4, CLIP), (0, RLEN - CLIP)], START + CLIP)):
    with pysam.AlignmentFile(d / (name + ".bam"), "wb", header=hdr) as out:
        for n in range(NREADS):
            r = pysam.AlignedSegment()
            r.query_name, r.query_sequence, r.flag = "r%d" % n, seq, 0
            r.reference_id, r.reference_start, r.mapping_quality = 0, pos, 60
            r.cigar = cigar
            r.query_qualities = pysam.qualitystring_to_array("I" * RLEN)
            out.write(r)
    pysam.index(str(d / (name + ".bam")))

print("edit at read offset %d (0-based) = base %d of the read as sequenced" % (EDIT_OFF, EDIT_OFF + 1))
print("reference position of that base: %d (1-based)" % (START + EDIT_OFF + 1))
print("in clipped.bam it is the FIRST aligned base")
