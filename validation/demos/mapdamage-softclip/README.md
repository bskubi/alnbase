# mapDamage measures damage from the first aligned base (MD-1), reproduced

mapDamage2.0's misincorporation table is the ancient-DNA field's authentication figure:
the C>T frequency as a function of distance from the read's 5' end. That distance is
measured on the **aligned block**, not on the read as sequenced, so soft-clipped bases are
excluded from both the numerator and the denominator and every position index shifts by
the clip length. Because a local aligner clips a terminal base *precisely when it
mismatches*, and a damaged base is a mismatch, the excluded bases are the damaged ones.
The estimate is therefore biased downward, in one direction, by an amount set by the
aligner rather than by the sample.

Full entry: `docs/design/research/non-methylation-tool-findings.md` (MD-1).

Tool: mapDamage 2.2.x from GitHub `484bb30` (`main`), run with `--merge-libraries
--no-stats --no-plot`. Checked on 2026-09-18.

```
PY=<python with pysam> MAPDAMAGE=<mapDamage checkout> BWA_BIN=<dir with bwa and samtools> \
  ./run.sh [workdir]
```

`expected.txt` is the output of that run.

## Case A: the mechanism

4000 reads, 40 bp, with 40% deamination applied to the read's first two bases. The reads
are written to two BAMs with byte-identical SEQ, differing only in CIGAR: `40M` (what a
global aligner emits) and `2S38M` (what a local aligner emits when those two bases
mismatch). The soft-clipped bases are still present in SEQ in the second file.

| position | ground truth | `40M` | `2S38M` |
|---|---|---|---|
| 1 | 0.3788 | 0.3788 | 0.0000 |
| 2 | 0.4141 | 0.4141 | 0.0000 |

The end-to-end file recovers the truth exactly. The clipped file reports no damage
anywhere, and its position 1 holds what was position 3 of the read.

## Case B: what a real aligner costs

60 000 simulated ancient fragments, 35-90 bp, geometric damage decay from each terminus,
aligned with `bwa aln` (global, never clips) and `bwa mem` (local) from the same FASTQ.

High-damage library, 1748 of 59 948 mapped reads soft-clipped:

| position | ground truth | `bwa aln` | `bwa mem` | mem / truth |
|---|---|---|---|---|
| 1 | 0.4006 | 0.3988 | 0.3697 | 0.92 |
| 2 | 0.2194 | 0.2208 | 0.1929 | 0.88 |
| 3 | 0.1251 | 0.1204 | 0.1020 | 0.82 |
| 4 | 0.0639 | 0.0635 | 0.0538 | 0.84 |
| 5 | 0.0331 | 0.0349 | 0.0293 | 0.89 |

Low-damage library, 140 of 59 986 mapped reads soft-clipped:

| position | ground truth | `bwa aln` | `bwa mem` | mem / truth |
|---|---|---|---|---|
| 1 | 0.1062 | 0.1022 | 0.0991 | 0.93 |
| 2 | 0.0578 | 0.0562 | 0.0538 | 0.93 |
| 3 | 0.0305 | 0.0311 | 0.0296 | 0.97 |

Note that `bwa aln` is itself slightly low in the second table (0.1022 against 0.1062).
That part is the aligner's own ascertainment bias: a read carrying several deaminations
exceeds the permitted mismatch count and never maps at all. It is not mapDamage's doing,
and the mapDamage-attributable part of the deficit is the `mem` column minus the `aln`
column, about 3% relative at low damage and 7-15% relative at high damage.

## Why no quality-control step catches it

The misincorporation plot is the only check the field runs on this number. There is no
independent measurement of deamination to compare against, so a biased curve is confirmed
by the same code that produced it. mapDamage does record the clipped-read count in the `S`
column of `misincorporation.txt` and plots it, but `S` reports the *fraction of reads
clipped at that position*, which does not tell a reader that the substitution frequencies
printed beside it are diluted, and no documentation connects the two.
