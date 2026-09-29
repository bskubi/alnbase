# bam2bakR's mutation caller, probed (B2-1, B2-2), reproduced

bam2bakR converts aligned nucleotide-recoding RNA-seq (TimeLapse-seq, SLAM-seq, TUC-seq)
into the per-read tables that bakR fits. Its `_counts.csv` gives, for each read pair, the
number of each reference base observed (`nA`, `nC`, `nT`, `nG`) and the number of each of
the 21 substitution types. That is a genuinely molecule-level output and the same shape of
query alnbase performs, so the interest here is not that bam2bakR cannot do this — it can —
but in exactly where its counting rules sit.

Full entries: `docs/design/research/non-methylation-tool-findings.md` (B2-1..B2-5).

Tool: `simonlabcode/bam2bakR` at GitHub `f8b492e` (2025-09-16),
`workflow/scripts/mut_call.py`. Checked 2026-09-18.

```
PY=<python with pysam> MUTCALL=<.../workflow/scripts/mut_call.py> ./run.sh [workdir]
```

`expected.txt` is the output of that run. The reference is the repeating tetramer `ACGT`,
so a read starting at reference index 20 has a T at every read offset `j` with `j % 4 == 3`
and the number of T's in the read is known exactly. `make_bams.py` computes MD and NM
itself, so no `samtools calmd` step is needed.

## 1. `--minQual` is compared against Phred + 33 — bug

One read carrying a T→C at Phred 10 and another at Phred 5, both far from the read ends.

| `--minQual` | TC | nT | effective cutoff |
|---|---|---|---|
| 30 | 2 | 15 | Phred > −3 — no filter |
| 33 | 2 | 15 | Phred > 0 |
| 34 | 2 | 15 | Phred > 1 |
| 40 | 1 | 14 | Phred > 7 |
| 44 | 0 | 13 | Phred > 11 |
| 76 | 0 | 0 | Phred > 43 |

`mut_call.py` compares `b[2] + 33 > args.minQual`, where `b[2]` came from pysam's
`query_alignment_qualities` and is therefore already a decoded Phred score. The `+ 33`
re-encodes it as a Sanger ASCII code point before the comparison.

bam2bakR's `config/config.yaml:53` ships `minqual: 40` under the comment
*"Minimum base quality to call mutation"*. That is a Phred cutoff of **8**, not 40. A real
Phred-40 cutoff would need `minqual: 73`. The sharper consequence is at the other end: any
value of 33 or below silently disables the filter completely, so a user who moves from the
shipped 40 down to a seemingly more permissive 30 does not loosen the filter — they remove
it.

## 2. The mutation numerator and the trials denominator apply different filters — bug

A read in which **every** T is mutated to C, all at Phred 35. The true per-T conversion
rate is 1.000 by construction.

```
nT (trials) = 15   TC (mutations) = 12   reported rate = 0.800
```

`nT` (`mut_call.py:180-185`) filters on base quality alone. `TC` (`:214`) additionally
requires the base to be more than `--minDist` positions from either end of the aligned
block, and not to fall on a called SNP. The three T's at read offsets 3, 55 and 59 are
therefore trials that can never become mutations.

The bias is roughly `2 * minDist / readlength` — 20% here, and `--minDist` is never passed
by `workflow/scripts/mut_call.sh`, so it is 5 in every pipeline run. The author flags it in
a code comment at `:180` (*"I think this should be also filtered for closeness to read end
and presence of SNPs"*), but it appears nowhere in the pipeline's configuration or
documentation.

How much this matters downstream is worth stating honestly: bakR fits a mixture over reads,
and both mixture components are deflated by the same factor, so the estimated *fraction* of
new RNA is considerably more robust than the reported conversion *rate* is. The rate itself
is the quantity users read off directly to judge labelling efficiency, and it is wrong.

The `_cU.csv` output is affected more severely: its `trials` column (`:203-211`) applies
**no** filter at all, not even the base-quality one, while its `n` column (`:218-221`)
applies all three.
