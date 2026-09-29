# One record per FLAG a caller might filter on

The same 20 bp read over a single CG, written six times: plain, duplicate (0x400),
secondary (0x100), supplementary (0x800), QC-fail (0x200) and unmapped (0x4). The first
five are all walked and all produce the call, because filtering that samtools can do
belongs upstream of alnbase rather than inside it; choosing not to filter is a decision,
so it is pinned here rather than left to be discovered.

The unmapped record has no alignment to walk. It is skipped and counted rather than
guessed at or treated as fatal, and `expected_scan.txt` holds the count line alnbase
prints, which is the only place that behaviour is visible.
