# A molecule base observed twice is one observation

`chrO` is `AAAACGAAAA` three times over (30 bp), so the top strand has a CG at 5-6, 15-16
and 25-26 (1-based). One 30 bp fragment is read 20 + 20 from both ends: R1 (FLAG 99)
aligns `20M` at 1 and R2 (FLAG 147) `20M` at 11, so the two mates both observe reference
11-20 — including the CG at 15-16. Every cytosine is methylated and every base is Q30.

`expected_raw.tsv` is what a caller sees with the duplicate left in place: **four** rows for
three cytosines, the one at 0-based 14 appearing twice, once from each mate at a different
`off_5p` (14 from R1, 15 from R2). Those two rows are not two molecules. They are one
molecule base counted twice, and counting it twice is what biases a methylation rate
towards whatever the overlapping region happens to say.

`expected.tsv` is the same file after `alnbase overlap`. R1 and R2 score equally, so R1
keeps the overlap and R2 gives up its read positions `e2 >= L - hi = 10`, which on a
reverse record is the CIGAR head: R2 becomes `10S10M` at 21. A soft-clipped base emits a
CLIP column, which only `~` matches, so the clipped copy makes no call and the duplicate
row is gone. The surviving row carries Q40, not Q30: the two agreeing observations were
summed and capped at `--qual-cap` (40).

## What this probes in other tools

`alnbase overlap` never reads `MC`, and these records deliberately carry none. It finds the
fragment length from reference positions both reads cover on opposite strands
(`L = e1 + e2 + 1`; `docs/reference/05-overlap.md` §4), so the estimate does not depend on a
tag the aligner may not have written. BISCUIT's BC-11 is the other behaviour: without `MC`
its overlap detection does not fire and the doubly observed base is counted twice, which is
exactly the four-row table above. The fixture holds both tables so the finding and its
absence are checked against each other rather than asserted.

The two mates walk the **same** strand here, not opposite ones: under
`queries/strand/directional.toml` FLAG 99 is OT and FLAG 147 is CTOT, and both are walked
along the top strand, so both call the reference C at 15 rather than the two cytosines of
the palindrome. That is why the rows collide in the first place.
