# When the mates disagree, the surviving base can be the one that lost

The same 30 bp fragment and reference as [`../overlap_agreeing_mates/`](../overlap_agreeing_mates/),
with one base changed: at reference 15 R1 reads `C` (methylated) and R2 reads `T`
(converted). The mates are the same molecule, so one of them is wrong. Qualities are set so
that the two decisions `alnbase overlap` makes come apart:

- **Which end keeps the overlap** is decided once, for the whole overlap, by
  `mean(MAPQ) + --score-qual-weight x mean(base quality)` over each end's copy
  (§6.1). R1's ten overlap bases average Q33.5 (nine at Q35, the clash at Q20) and R2's
  average Q31 (nine at Q30, the clash at Q40), so **R1 keeps** and R2 is clipped to
  `10S10M`.
- **Who wins a disagreeing base** is decided per base, by `(quality, MAPQ, is_R2)` (§5).
  At reference 15 R2's Q40 beats R1's Q20, so **R2 wins**. Under the default
  `--mismatch-qual subtract` the winner keeps `40 - 20 = 20` and the loser is zeroed; under
  the default `--mismatch-base none` neither base is changed.

The two need not agree, and here they do not. R1 keeps the overlap, so the row that
survives is R1's `C` — the call that lost — reported as `CG_met` at **quality 0**. The
higher-quality observation said the cytosine was converted, and it was clipped away.

This is not a defect being hidden: quality 0 is the conflict, written where a caller can
act on it, and `expected.tsv` pins it so that it cannot change silently. But it does mean a
downstream `group by` that ignores `qual` will take the losing call at face value.
`--mismatch-base set-n` is the option that removes the row instead, since a read base of
`N` matches neither `read = "C~~"` nor `read = "T~~"`.

`expected_raw.tsv` shows what the same file gives with no overlap resolution: reference 15
yields both a `CG_met` and a `CG_unmet` row, from one molecule base. A tool that counts
both records a methylation level of 0.5 at a site where a single molecule was observed once.

One mismatch in ten compared positions is `mism_frac` 0.1, under the 0.15 default, so the
template still resolves; a second disagreement would push it to 0.2 and the whole template
would be passed through unmodified as `high_mismatch` (§4.5).
