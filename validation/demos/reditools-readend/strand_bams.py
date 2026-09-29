"""Reverse-strand read: is -mbp counted from the read's own 5' end (the
sequencing start, which for a FLAG-16 read is the RIGHT end in reference
coordinates) or from SEQ position 0 (the LEFT end)?

Reference all A. SEQ is 40 bp, stored reference-forward as SAM requires.
We place a T at SEQ index 2 and a T at SEQ index 37 -- as reference A>T
mismatches at reference positions 1003 and 1038.

For a FLAG-16 read the original read's 5' end corresponds to the RIGHT end of
SEQ. So:
  SEQ index 37 is read base 3 as sequenced   (near the 5' end)
  SEQ index  2 is read base 38 as sequenced  (near the 3' end)

-mbp 3 ("ignores the first 3 bases in each read") should therefore drop the
reference-1038 event and keep the reference-1003 one, if the filter really
means the read as sequenced.
"""
import pysam, pathlib, sys
REF_LEN, RLEN, NREADS, START = 3000, 40, 300, 1000
d = pathlib.Path(sys.argv[1]); d.mkdir(parents=True, exist_ok=True)
ref = "A" * REF_LEN
(d/"ref.fa").write_text(">chr1\n" + "\n".join(ref[i:i+60] for i in range(0,REF_LEN,60)) + "\n")
pysam.faidx(str(d/"ref.fa"))
seq = list("A"*RLEN); seq[2] = "T"; seq[37] = "T"; seq = "".join(seq)
hdr = {"HD":{"VN":"1.6","SO":"coordinate"},"SQ":[{"SN":"chr1","LN":REF_LEN}]}
for name, flag in (("fwd",0),("rev",16)):
    with pysam.AlignmentFile(d/(name+".bam"),"wb",header=hdr) as out:
        for n in range(NREADS):
            r = pysam.AlignedSegment()
            r.query_name, r.query_sequence, r.flag = "r%d"%n, seq, flag
            r.reference_id, r.reference_start, r.mapping_quality = 0, START, 60
            r.cigar = [(0,RLEN)]
            r.query_qualities = pysam.qualitystring_to_array("I"*RLEN)
            out.write(r)
    pysam.index(str(d/(name+".bam")))
print("ref 1003 = SEQ index 2 ; ref 1038 = SEQ index 37")
print("FLAG 16: ref1038 is read base 3 as sequenced, ref1003 is read base 38")
