# dynast's consensus caller, probed (DY-1, DY-2, DY-5), reproduced

`dynast consensus` collapses each UMI family into a single consensus alignment before
`dynast count` reads it back. That consensus alignment is synthesised from scratch —
sequence, CIGAR, MD and NM are all written by dynast — so it is worth checking that the
tags it writes describe the alignment it writes.

Full entries: `docs/design/research/non-methylation-tool-findings.md` (DY-1..DY-6).

Tool: `pachterlab/dynast` at GitHub `04bb0b8` (2025-08-11, v1.0.1),
`dynast/preprocessing/consensus.py`. Checked 2026-09-18.

```
DYNAST=<path to a dynast checkout> PY=<python with pysam, numpy, ngs_tools, anndata> \
    $PY probe.py
```

`expected.txt` is the output of that run. No aligner, GTF or real BAM is needed: the probe
calls `call_consensus_from_reads` directly on hand-built reads and compares its tags
against `samtools calmd`, which recomputes MD and NM from the alignment and the reference.
Every family is two byte-identical reads, so the consensus base calls are unambiguous and
any disagreement is the tag writer's.

## 1. The MD tag is not consistent with the CIGAR next to deletions — bug

| case | dynast | samtools |
|---|---|---|
| deletion then mismatch | `5^CG13` | `5^C0G13` |
| mismatch then 2 bp deletion | `4A0^C0G13` | `4A0^CG13` |
| isolated deletion (control) | `5^C14` | `5^C14` |
| isolated mismatch (control) | `4A15` | `4A15` |

Both failures come from one omission: the deletion branch (`consensus.py:122-131`) never
clears `md_zero`, which the mismatch branch sets at `:160`.

* After a deletion, `md_zero` is whatever the preceding match left it, so a mismatch that
  follows a deletion gets no `0` separator. `5^CG13` reads as a **two-base deletion**
  `^CG` rather than a one-base deletion followed by a mismatched G.
* After a mismatch, `md_zero` is still `True` when the next position is a deletion, so a
  `0` is emitted *between* the bases of a multi-base deletion. `4A0^C0G13` reads as a
  **one-base deletion** `^C` followed by a mismatched G.

Read literally, both tags describe a different alignment than the CIGAR does. pysam
happens to recover, because its MD parser is CIGAR-driven and takes the deletion length
from the `D` operator rather than from the MD string — which is why dynast's own
`count` step, which reads the reference base back with
`get_aligned_pairs(with_seq=True)`, is not affected. A consumer that sizes the deletion
from the MD string is not so lucky, and `samtools calmd` disagrees outright.

## 2. NM omits deleted bases entirely — bug

`nm += 1` fires only in the mismatch branch (`consensus.py:161`). Nothing counts the
deleted bases, and the `deletions` counter accumulated at `:88` is never read.

The control case is the cleanest statement of it: a consensus alignment with CIGAR
`5M1D14M` and no mismatch at all is written with **`NM:i:0`**. By the SAM specification NM
is the edit distance — mismatches plus inserted and deleted bases — so it should be 1.

This one is not rescued by a lenient parser. Any downstream step that filters or ranks on
edit distance sees a deletion-bearing consensus read as a perfect match.

## 3. A gap between reads of one UMI family becomes a splice junction — bug

Two reads of the same UMI, one covering ref 0-4 and one covering ref 12-19, with nothing
covering the 7 bp between:

```
dynast       0 5M7N8M     MD 13           NM 0
blocks = [(0, 5), (12, 20)]
```

`N` is the reference-skip operator: it asserts that the molecule was spliced across those
7 bases. The consensus builder emits it for every position where no read of the family
supplied a base (`consensus.py:117`), which conflates "intron" with "not sequenced". For a
UMI family whose reads do not tile the molecule contiguously — routine for 10x data, where
fragment starts vary within a UMI, and for paired reads with an unsequenced insert — the
output BAM carries junctions that are not in the data.
