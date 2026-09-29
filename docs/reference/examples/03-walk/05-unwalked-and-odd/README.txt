Records that are not walked (unmapped, off-reference), records whose SEQ or
QUAL is unusual, and a mapped record with SEQ '*' (skipped and counted
since 0.1.2; it used to crash the run).
Each record's purpose is in the @CO lines of reads.sam.
