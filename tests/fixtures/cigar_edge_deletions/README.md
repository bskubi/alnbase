# A CIGAR that begins or ends with a deletion

`chrD` is 30 bp of A with a CG at 10-11 (1-based). All three reads are 10 bases of SEQ and
all three call the cytosine at 10; the deletions sit at the very edge of the CIGAR, where a
walk that assumes the first and last operation consume query bases goes wrong.

| read | CIGAR | POS | aligned to | deleted |
|---|---|---|---|---|
| `plain` | `10M` | 1 | 1-10 | — |
| `ends_del` | `10M2D` | 1 | 1-10 | 11-12 |
| `starts_del` | `2D10M` | 1 | 3-12 | 1-2 |

`samtools view` accepts both, which is worth knowing: neither is rejected on the way in, so
a caller meets them. In `ends_del` the deletion covers exactly the two reference bases the
cytosine's context needs, and the call is unchanged, because a deleted position still emits
a column carrying its reference base — the same reference rule `indel_reference_rule`
checks in the middle of a read, here at the edge where the read has run out. The trailing
`2D` also aligns nothing, so `off_3p` stays 0: the C is the last base sequenced.

`starts_del` is the mirror. Its first aligned base is reference 3 rather than POS, so a
caller that equates POS with the first aligned base is two positions out for the whole
read; the cytosine is nonetheless at reference 10 with `off_5p` 7.
