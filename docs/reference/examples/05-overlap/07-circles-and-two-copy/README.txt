07-circles-and-two-copy: molecules that visit a locus in an unexpected order or twice.

self_circle                 resolved normally (L=80): a sheared self-circle is a circular
                            permutation with every base present once.
sister_chromatid            resolved on the true group (L=80) although the primaries'
                            reference intervals coincide; the decoy group L=40 fails the
                            anchor test.
two_copy_past_end           refused past_end: the only group implies a 25 bp molecule,
                            but R2 has 30 aligned read positions.
two_copy_indistinguishable  resolved as a 50 bp colinear fragment: from alignments alone
                            nothing distinguishes it from one. A limit of any
                            alignment-based method.

The mpileup table shows depth; for sister_chromatid look at which positions each tool
counts once (see the reference text).
QUAL characters: '?' = 30, 'I' = 40.
