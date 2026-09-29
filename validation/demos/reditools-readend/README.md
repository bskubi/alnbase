# REDItools' read-end trim filters, probed (RT-2, RT-3), reproduced

REDItools detects RNA editing by comparing each aligned base against the reference. Its
`-mbp`/`-Mbp` options exclude bases near the read ends, where random-hexamer priming
artefacts (5') and quality decay (3') concentrate — the same quantity alnbase calls
`off_5p`/`off_3p`. This demo probes what "near the read end" actually means, using pairs of
hand-built BAMs that differ in exactly one property.

Full entries: `docs/design/research/non-methylation-tool-findings.md` (RT-1..RT-6).

Tool: REDItools3 v3.7 from GitHub `8cfebea`, installed with `pip install .`. Checked
2026-09-18.

```
PY=<python with pysam and reditools installed> ./run.sh [workdir]
```

`expected.txt` is the output of that run.

## 1. Soft clipping does not shift the axis — clean

The same reads with byte-identical SEQ, written twice as `40M` and `2S38M`. The edit sits
at read base 3 in both. The `-mbp` threshold is identical across the two files (survives
`-mbp 2`, dropped at `-mbp 3`), so REDItools3 counts from the read's own first base rather
than from the first aligned base.

This is the right convention and worth stating plainly: it is exactly what mapDamage2.0
gets wrong (MD-1), and it matches alnbase's definition. Note that REDItools **v1** uses the
other convention — it measures from `query_alignment_start`/`query_alignment_end` — so the
same `-a 6-0` / `-mbp 6` intent filters differently in v1 than in v2/v3 on a soft-clipped
BAM. Neither convention is documented.

Hard clips are invisible to all three versions, since `query_length` excludes them.

## 2. The filter does not flip for reverse-strand reads — bug

One 40 bp read written twice, identical SEQ, differing only in FLAG (0 vs 16). SEQ is
stored reference-forward as SAM requires, so for the FLAG-16 record **SEQ index 37 is read
base 3 as sequenced** and **SEQ index 2 is read base 38**.

| file | `-mbp 3` → ref 1003 (SEQ idx 2) | → ref 1038 (SEQ idx 37) |
|---|---|---|
| `fwd` | DROP | kept |
| `rev` | DROP | kept |

The reverse file behaves identically to the forward one. `-mbp` therefore always trims from
SEQ position 0, which for a minus-strand read is the **3' end**. The help text says
*"Ignores the first `-mbp` bases in each read."* For roughly half the reads in any library
it ignores the last ones instead, and leaves the first ones — the random-priming-biased
ones the option exists to suppress — in place.

## 3. `-mbp` and `-Mbp` are not symmetric — bug

| | 3rd base from that end |
|---|---|
| `-mbp 2` → kept, `-mbp 3` → dropped | correct: ignores the first N |
| `-Mbp 3` → **kept**, `-Mbp 4` → dropped | ignores only the last N−1 |

Anyone setting `-mbp 6 -Mbp 6` for a symmetric trim gets 6 bases at one end and 5 at the
other. `-Mbp 1` removes nothing at all. The same inequality appears in REDItools2, and
REDItools3 ships a unit test asserting the off-by-one result.
