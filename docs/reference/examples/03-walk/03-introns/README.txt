CIGAR N. With context c, an intron of length n > 2c emits c reference
columns at each end (read ',') and one marker column (read ',' refr ',')
whose refr_pos is the first elided base in walk order. n <= 2c is emitted
whole with no marker. Every intron column has qual NULL and the off_5p of
the last read base before it. Run with the derived default (widest query
span = 3), then --splice-context 0 and 5.
