#!/usr/bin/env python3
"""Keep only read pairs whose mates do not overlap on the reference.

    select_nonoverlapping_pairs.py IN.bam OUT.bam

Every tool resolves overlapping mates its own way, and alnbase has its own
`overlap` command, so tool comparisons are made on fragments where the question
does not arise. A pair is kept when both mates are mapped to the same contig and
their aligned reference intervals do not intersect, or when they are on different
contigs. Unpaired and single-mapped records are kept as they are.

The input must have each pair's records adjacent (as Bismark writes them, or after
`samtools collate`). Writes the kept records in input order and prints counts.
"""

import sys

import pysam


def overlaps(a: pysam.AlignedSegment, b: pysam.AlignedSegment) -> bool:
    if a.is_unmapped or b.is_unmapped or a.reference_id != b.reference_id:
        return False
    return a.reference_start < b.reference_end and b.reference_start < a.reference_end


def main(src: str, dst: str) -> None:
    kept = dropped = 0
    with pysam.AlignmentFile(src) as bam, pysam.AlignmentFile(dst, "wb", template=bam) as out:
        pending = None
        for rec in bam:
            if not rec.is_paired:
                out.write(rec)
                kept += 1
                continue
            if pending is None:
                pending = rec
                continue
            if pending.query_name != rec.query_name:
                sys.exit(f"records of pair {pending.query_name} are not adjacent")
            if overlaps(pending, rec):
                dropped += 2
            else:
                out.write(pending)
                out.write(rec)
                kept += 2
            pending = None
        if pending is not None:
            sys.exit(f"pair {pending.query_name} has only one record")
    print(f"kept {kept} records, dropped {dropped} records in overlapping pairs", file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:3])
