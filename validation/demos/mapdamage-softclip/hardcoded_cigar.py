"""Case A: the same damaged reads, written twice with different CIGARs.

  endtoend.bam  CIGAR 40M, as a global aligner (bwa aln) emits
  clipped.bam   CIGAR 2S38M, as a local aligner emits when the two 5'-terminal
                bases mismatch

The SEQ strings are byte-identical between the two files, and the two 5'-terminal
bases are present in SEQ in both. Only the CIGAR and the reported alignment start
differ. Any tool that measures distance from the read's 5' end should give the
same answer for both.

Usage: hardcoded_cigar.py <outdir>
"""

import pathlib
import random
import sys

import pysam

REF_LEN, RLEN, NREADS = 5000, 40, 4000
DAMAGE = 0.40
CLIP = 2


def main(outdir):
    random.seed(7)
    d = pathlib.Path(outdir)
    d.mkdir(parents=True, exist_ok=True)

    ref = "".join(random.choice("ACGT") for _ in range(REF_LEN))
    (d / "ref.fa").write_text(
        ">chr1\n" + "\n".join(ref[i : i + 60] for i in range(0, REF_LEN, 60)) + "\n"
    )
    pysam.faidx(str(d / "ref.fa"))

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": REF_LEN}],
    }
    ete = pysam.AlignmentFile(d / "endtoend.bam", "wb", header=header)
    clp = pysam.AlignmentFile(d / "clipped.bam", "wb", header=header)

    truth = {i: [0, 0] for i in range(1, CLIP + 1)}
    for n, start in enumerate(
        sorted(random.randrange(0, REF_LEN - RLEN) for _ in range(NREADS))
    ):
        seq = list(ref[start : start + RLEN])
        for p in range(CLIP):
            if seq[p] == "C":
                truth[p + 1][0] += 1
                if random.random() < DAMAGE:
                    seq[p] = "T"
                    truth[p + 1][1] += 1
        seq = "".join(seq)

        for out, cigar, pos in (
            (ete, [(0, RLEN)], start),
            (clp, [(4, CLIP), (0, RLEN - CLIP)], start + CLIP),
        ):
            rec = pysam.AlignedSegment()
            rec.query_name = "r%d" % n
            rec.query_sequence = seq
            rec.flag = 0
            rec.reference_id = 0
            rec.reference_start = pos
            rec.mapping_quality = 60
            rec.cigar = cigar
            rec.query_qualities = pysam.qualitystring_to_array("I" * RLEN)
            out.write(rec)

    ete.close()
    clp.close()
    for name in ("endtoend.bam", "clipped.bam"):
        pysam.index(str(d / name))

    print("ground truth, against the read's own 5' end:")
    for p in range(1, CLIP + 1):
        c, t = truth[p]
        print("  position %d: %d reference C, %d deaminated -> %.4f" % (p, c, t, t / c))


if __name__ == "__main__":
    main(sys.argv[1])
