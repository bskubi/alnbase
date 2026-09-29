# The same call at Phred 40 and at Phred 2

`good` and `tail` are the same 12 bp alignment over the same CG; `tail` carries Phred 2
over its last five bases, which includes the C. Both are called, and each call carries its
own base quality in the `qual` column.

No quality threshold is built in, so a caller-side filter is a predicate over these rows.
The reason to pin it is that a threshold is invisible when the input does not exercise it:
`lambda-stress` writes a flat Phred 40 on bases that are wrong 1 time in 200, which makes
MethylDackel's `-p` look free, and `lambda-illumina` writes a flat Phred 30, where it also
removes nothing. Two reads settle what neither simulation can.
