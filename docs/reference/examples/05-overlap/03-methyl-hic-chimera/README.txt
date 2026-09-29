03-methyl-hic-chimera: mate overlap on a split (chimeric) read, the case reference-interval
tools cannot see.

hic_sup_overlap       R2 overlaps only R1's supplementary segment; the primaries are on
                      different contigs. alnbase recovers L=100 from the supplementary,
                      and clips R1's supplementary (R2 scores higher: MAPQ 60 vs 40).
hic_junction_overlap  the overlap spans the ligation junction; both reads are split.
                      R1 wins; R2's supplementary lies wholly inside the overlap and is
                      not written, R2's primary is clipped 30S20M and moves to chr2:121.
hic_inverted          the second locus entered against the reference; strands flip,
                      pairing still works.
colinear_control      a plain colinear pair.

The mpileup table shows depth over each overlap: on the chimeric templates samtools
mpileup counts the doubly observed bases twice with or without its overlap detection,
and after alnbase each base is counted once. On the colinear control mpileup's own
detection already removes the double count (depth 1 without -x, 2 with -x).

QUAL characters: '?' = 30, 'I' = 40.
