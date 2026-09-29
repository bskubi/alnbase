# A soft clip contributes no column but is still a cycle

`chrS` is 30 bp of A with a CG at 10-11 (1-based). All three reads call the same cytosine
at 10 and differ only in what surrounds it.

| read | CIGAR | POS | the point |
|---|---|---|---|
| `plain` | `10M` | 10 | the control; the C is the first base sequenced |
| `clip_before` | `5S10M` | 10 | five clipped bases precede the C, so it is the sixth cycle |
| `clip_hides_g` | `10M5S` | 1 | the C is the **last aligned base**; the G it needs lies under a soft-clipped `T` |

Two separate things are being pinned.

**Context is the reference's, not the read's.** In `clip_hides_g` the read's own base at
reference 11 is a `T`, and the reference's is a `G`. A soft-clipped base emits a flank
column that carries the real reference base and the clipped base's offsets but no quality
(`docs/reference/01-reference-and-codes.md` §1.1), and the `~` in the pattern's read row
admits it, so the reference context `CG` is read through the clip and the call lands. A
caller that took context from SEQ would see `CT` and call this a CHH, or drop it.

**Offsets count clipped bases.** `off_5p` is 0, 5 and 9 for the three reads, and
`off_5p + off_3p` is `len(SEQ) - 1` on every row: 9 for `plain` and 14 for the two clipped
reads. That is what makes `off_5p < n` mean "the first n sequencing cycles" for M-bias and
read-end trimming regardless of clipping, which is the same property `hard_clip_offsets`
checks for hard clips.
