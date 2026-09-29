# A cytosine whose context runs off the end of the contig

`chrC` is 23 bp beginning with G and ending with C, with an ordinary CG at 11-12 (1-based)
as the control. The C at 23 has no base after it, and the G at 1 -- a cytosine on the
bottom strand -- has none before it. Their context columns are reference pads, so both
fall to `Cx` rather than being called CHH or silently dropped.

This is the boundary where off-by-one errors live. BISCUIT's `pileup` builds half-open
windows that end at the contig length, so it loses the last position of every contig and,
with `-g`, of every region; a run sharded into adjacent regions loses one position per
shard (BC-18, reproduced in `alnbase-validation/demos/biscuit-bugs/`).
