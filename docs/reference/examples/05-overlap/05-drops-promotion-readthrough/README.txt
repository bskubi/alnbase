05-drops-promotion-readthrough: what happens to records that lose all their aligned bases.

- A supplementary with nothing left is not written at all.
- A primary with nothing left, when a supplementary of the same read survives: the
  5'-most survivor's alignment moves onto the primary record ("promotion"). Its hard
  clips become soft clips over the primary's full SEQ, NM/MD/AS come from the survivor,
  and the supplementary record is not written. The read keeps exactly one primary line
  (primary_inside_overlap: 20M60S at chr2:141, 80 bp SEQ).
- A primary with nothing left and no survivor is written unmapped: flag 0x4 set, 0x2 and
  0x10 cleared, MAPQ 0, CIGAR '*', NM/MD/AS/SA/XA removed, SEQ/QUAL reverse-complemented
  back to sequencing orientation if it was reverse, RNAME/POS set to the mate's.
- Bases past the end of the molecule (read-through) are clipped from both reads
  regardless of which end keeps the overlap, when they are aligned; when L <= read
  length, the losing read gives up everything.
- In a template where anything was clipped, dropped or promoted, mate fields
  (RNEXT/PNEXT/TLEN/0x2/0x8/0x20/MC) are rebuilt as samtools fixmate would, and SA is
  rebuilt for a read whose alignments changed. Templates left alone keep theirs.
- flagstat on the output counts 10 primary lines for 10 reads.
- --max-past-end alone cannot rescue read_through_aligned: the anchor test measures the
  same excess (see the reference text). With --max-anchor-slack raised too it resolves.
- adapter_indel: the terminal-indel cleanup inspects the ORIGINAL CIGAR, so an insertion
  inside the already-clipped adapter adds 3 more clipped bases (27M13S instead of 30M10S).

QUAL characters: '?' = 30, 'I' = 40.
