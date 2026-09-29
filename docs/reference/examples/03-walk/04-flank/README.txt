Pad (flank) columns: how many, what they carry, and contig edges.
By default the number of pad columns at each end is the span of the widest query
in the whole run, whatever query you care about; --end-context N sets it.
`col` reads "~" and accepts pads, so it fires on every pad column: its rows grow
with the widest query (see 02-query-language.md section 8.6, pads and the widest
query). The same BAM is run with a
span-1 query set (1 pad), span-3 (3 pads) and span-5 (5 pads).
Pads carry read_base '_', qual NULL, off_5p and off_3p NULL (a pad is not a
base of the read), and the REAL reference base when the
position is on the contig; off the contig refr_base is '_' and refr_pos is
emitted as-is (negative, or >= contig length).
