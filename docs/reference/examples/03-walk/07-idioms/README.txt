Idioms, each verified here:
1. Walk-first / walk-last read base via the flank: "fN" (mark ".+") and "Nf"
   (mark "+."), with the alias f = "{_:}", a pad or a clip. These are the
   first/last ALIGNED bases in WALK order,
   which for reads with 0x80 set are the sequencing 3'/5' ends. off_5p /
   off_3p are as sequenced for every FLAG, so for r2f/r2r "_N" fires on the
   base with off_3p = 0 and "N_" on the base with off_5p = 0.
   Beside a soft clip the flank column is a clip (':'), not a pad, so "_N"
   fires only where that end is unclipped and ":N" only where it is clipped
   (clip_only: clipped_first, not unclipped_first).
2. Reference context past the read end: "C_@CG" fires only when the G
   lies beyond the last aligned base.
3. Anchoring beside a deletion: "N.@NN" mark "+." (base 5' of the gap in
   walk order) and ".N@NN" mark ".+" (base 3' of it).
4. Reference-gap columns: with --insertions emit, every "N@." column is an
   inserted base, because soft-clipped bases never produce columns. A
   pattern bounded by aligned columns on both sides ("N.N", "N..N", ...)
   matches only insertions of that exact length; the one-sided "NN@N."
   matches the first base of an insertion of any length. The clipped
   records (clip_only, clip_3p) produce no N@. column at all: their clipped
   bases are clip columns in the flank, against the reference. Their
   walk-first / walk-last bases are aligned bases, and off_5p still counts
   the clipped bases.
