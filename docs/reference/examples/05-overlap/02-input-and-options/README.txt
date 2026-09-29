02-input-and-options: the grouping check, option validation, tag naming, and re-running.

check_grouping (src/overlap.rs) reads only the @HD line of the header: GO:query or
SO:queryname are accepted, anything else is refused. The record order itself is never
checked, which the interleaved run demonstrates: each run of equal QNAMEs is a template,
so an interleaved file silently becomes single-ended fragments.

QUAL '*' (absent qualities) stays '*': a base with no quality takes no part in the quality
arithmetic, so neither copy of it is edited. The overlap is still clipped.

--compression-level sets the BGZF level of the output: level 0 writes uncompressed blocks.

The last run feeds the output back in: the discarded copy is now soft-clipped, so no
overlap is found, and the previous run's tags are removed and replaced by oO:i:0.

QUAL characters: '?' = 30, 'I' = 40.
