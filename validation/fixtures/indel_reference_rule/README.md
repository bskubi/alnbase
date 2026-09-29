# The same cytosines read with and without an indel beside them

`chrI` has a CAG at 5-7 and a CG at 18-19 (1-based). `plain` matches the reference;
`del` deletes the A at 6, putting the G directly after the C in read order; `ins` inserts
a base between the C at 18 and its G. All three give the same two calls, CHG and CG,
because `../context.toml` takes the context from the reference alone.

That invariance is the point, and it is a choice. Bismark projects the reference through
the CIGAR, so for `del` the deleted base is skipped and the call becomes CG, and for `ins`
the inserted base makes the context unknown (U). Running these same reads with
`--insertions emit` and `../../../presets/bismark/queries.toml` reproduces that answer instead.
Neither rule is wrong; they are different queries, and this is MD-1 against BM-1/BM-2
reduced to three reads. On `lambda-illumina` the disagreement is 76 calls in 1,066,039.
