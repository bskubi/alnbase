# GLORI's read-end trim and duplicate-alignment tie-break, reproduced

GLORI-tools calls m6A from an A-to-G converted library. Its pileup step
(`pipelines/pileup_genome_multiprocessing.py`) is built on `pysam`'s pileup rather than
`samtools mpileup` text, so it keeps a handle to each read and *can* express read-coordinate
predicates — and it uses that to offer a read-end trim and a per-read count of unconverted
A bases. Both of those are the kind of molecule-level predicate a pileup normally cannot
reach, so this is not an instance of the pileup structural signature. Two of them are wrong.

Full entries: `docs/design/research/non-methylation-tool-findings.md` (GL-1..GL-5).

Tool: `liucongcas/GLORI-tools` at GitHub `def608b`. Checked 2026-09-18.

```
GLORI=<GLORI-tools checkout> PY=<python with pysam and biopython> ./run.sh [workdir]
```

`expected.txt` is the output of that run.

## 1. The 3' trim is one base shorter on the forward strand — bug

`make_bam.py` writes twenty 20 bp reads that all align `20M` to the same 20 bp of reference:
ten forward and ten reverse, identical in every other respect. Any difference between the
strands in what survives the trim is the trim's doing.

```
=== --trim-head 2 --trim-tail 2 ===
  pos         1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
  forward     0  0 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10  0
  reverse     0  0 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10  0  0
```

The forward strand keeps reference position 19, which the reverse strand drops. The same
one-position gap appears at `--trim-tail 1`, `2` and `3`; at `--trim-tail 1` the forward
strand is not trimmed at all.

Three of the four branches at `:53-60` drop exactly N bases; the forward 3' branch uses `<`
where its counterparts use `<=` or a bare `<` against the other end, and drops N−1:

| branch | test | bases dropped |
|---|---|---|
| reverse 3' | `query_position < trim_tail` | `trim_tail` |
| reverse 5' | `query_length - trim_head <= query_position` | `trim_head` |
| forward 5' | `query_position < trim_head` | `trim_head` |
| forward 3' | `query_length - trim_tail < query_position` | `trim_tail - 1` |

Credit where it is due: GLORI *does* flip the trim axis for reverse-strand reads, so that
"tail" means the read's 3' end on both strands. That is the thing REDItools gets wrong
(RT-2). The bug here is an off-by-one inside an otherwise correct design.

Honest sizing: `run_GLORI.py:101` invokes the pileup with no `--trim-head` or `--trim-tail`,
so the packaged pipeline runs with both at 0 and is unaffected. The exposure is a user who
runs the pileup script directly, which is a documented way to use the toolkit, and who then
gets a strand-biased trim.

## 2. The duplicate-alignment tie-break compares a value to itself — bug

GLORI's final BAM is a concatenation of a genome alignment and a transcriptome alignment
lifted back to genome coordinates (`run_GLORI.py:50`, `:71-73`), so one read name can appear
twice in the same pileup column. The pileup deduplicates by read name; when the two records
agree on the base, it is meant to keep the larger per-read A-count, which is the
conservative choice because that count is what `--cutoff` filters on.

It reads the stored count into `lastAcount` and then compares `lastAcount` against the same
dictionary entry, which has not been updated yet:

```python
lastAcount = PF_positive_A.get(query_name)
...
elif lastRead == query_base:
    if lastAcount <= PF_positive_A.get(query_name):     # always True
```

The reverse-strand copy at `:91` uses `<` and is therefore always False. So the winner is
decided by strand and write order, never by the counts. `probe_tiebreak.py` writes each
strand twice, once with the low count first and once with it second:

```
  fwd_low_first    written  0,18  kept 18 (second record)  ok
  fwd_high_first   written 18, 0  kept  0 (second record)  LOST THE LARGER COUNT
  rev_low_first    written  0,18  kept  0 (first record)  LOST THE LARGER COUNT
  rev_high_first   written 18, 0  kept 18 (first record)  ok
```

Forward always keeps the second record, reverse always the first. Half the time that is the
read that looked better converted than it was.

Honest sizing: GLORI is single-end (there is no `-2`/`--fastq2` anywhere in `run_GLORI.py`
or `mapping_reads.py`), so this branch never sees an overlapping mate pair — only the
genome/transcriptome duplicate above. Those two records are the same read, so they usually
carry the same A-count and the wrong pick is harmless. They diverge when the two alignments
soft-clip differently, since the count is taken over `query_sequence`, which includes
soft-clipped bases.

## 3. Defaults disagree with their own help text — bug, trivial

`-q` defaults to 10 and its help says `default=30`; `-m/--max-depth` defaults to 10000000
and says `default=10000` (`:150`, `:152`). `run_GLORI.py:101` passes neither, so a packaged
run applies a base-quality floor of 10 while the help promises 30.

## Why this one is worth recording even though GLORI gets the architecture right

The per-read A-count is a genuine molecule-level predicate: how many unconverted A bases the
whole read still carries, used to throw out incompletely converted molecules. GLORI reaches
it by smuggling one integer per read through an extra text column and then, in
`m6A_pileup_formatter.py:174-190`, pre-binning coverage into thirteen fixed cumulative bins
(1..10, 15, 20, unlimited). `m6A_caller.py --cutoff` then picks one bin, defaulting to 3, so
the binomial test that produces the published calls runs on the binned columns and not the
raw ones.

That works, and it is more than most pileup tools attempt. What it cannot do is what the
bins cost: the cutoff must be one of the thirteen (`m6A_caller.py:211-213` rejects anything
else), the per-read counts are gone by the time the caller runs, and the A-count predicate
cannot be combined with a read-position predicate, because the trim was applied before the
counts were binned. Each additional predicate would need its own column and its own bins.
