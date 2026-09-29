01-basics: which records participate and what a plain resolution writes.

Templates (see the @CO lines in reads.sam for coordinates):
  pair_overlap    ordinary FR pair with a 20-base overlap -> resolved
  pair_apart      mates do not meet -> no overlap (oO:i:0)
  pair_overlap_3  3-base overlap, below --min-overlap 4 -> no overlap (oO:i:0)
  single_end      unpaired record -> single-ended, no overlap (oO:i:0)
  mate_unmapped   R2 unmapped -> excluded from planning, template single-ended
  with_secondary  secondary copy of R2 is written through, tagged, never clipped

QUAL characters: '?' = 30, 'I' = 40, '5' = 20, '+' = 10, '!' = 0.
reads.sam was produced by a small generator from the stated molecule layouts
(SEQ read off ref.fa; FLAG, RNEXT/PNEXT, TLEN, MC, SA, NM, MD, AS filled the way a
BWA-like aligner would) and checked by hand. Run: ./run.sh (needs samtools).
