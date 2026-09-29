# An overlap the primaries cannot see

The molecule is a ligation product: chrA:11-30 joined to chrB:11-30, 40 bases long, read
30 from one end and 20 from the other. `chrA` is `AAAACGAAAA` four times over and `chrB`
is `TTTTCGTTTT` four times over (40 bp each), so each contig has a top-strand CG at
1-based 5-6, 15-16, 25-26 and 35-36, and the two sequences share nothing at the junction.

R1 is 30 bases and straddles the junction, so it aligns twice:

| record | FLAG | POS | MAPQ | CIGAR | molecule bases f |
|---|---|---|---|---|---|
| R1 primary | 97 | chrA:11 | 60 | `20M10S` | 0-19 |
| R1 supplementary | 2145 | chrB:11 | 40 | `20H10M` | 20-29 |
| R2 primary | 145 | chrB:11 | 60 | `20M` | 20-39 |

**The two primaries are on different contigs.** Their reference intervals do not intersect
and never can, so any method that looks for an overlap by intersecting the mates' primary
alignments finds none here. The shared molecule bases — f 20-29, chrB:11-20 — are held by
R1's *supplementary* and R2's primary.

## The arithmetic

alnbase indexes every base of every record by its position in the read, and uses the
reference only to find the fragment length (`docs/reference/05-overlap.md` §4). On the
supplementary the leading `20H` still counts, so its q 0-9 are end indices e1 20-29; R2 is
reverse, so its q 0-9 are e2 19-10. Each of the ten chrB:11-20 positions pairs a forward
R1 base with a reverse R2 base and gives `L = e1 + e2 + 1 = 40` — one group, n = 10.

`lo = 40 - 20 = 20`, `hi = min(30, 40) = 30`, predicted 10. R1's 3'-most aligned base is
e 29 and R2's is e 19, so reach1 = reach2 = 0 and past_end = 0; span_frac = 10/10 = 1 and
the ten bases agree, so mism_frac = 0. No refusal test fires (§4.5).

Score over f 20-29 is `mean(MAPQ) + mean(qual)`: R1's ten bases all sit on the MAPQ-40
supplementary, giving 40 + 30 = 70, against R2's 60 + 30 = 90. **R2 keeps the overlap**, so
R1 gives up every aligned base with e >= 20 — which is the whole supplementary. Its
`clip_back` of 10 covers all 10 aligned bases, so that record is marked drop and is not
written (§6.4, §7.1). The primary keeps all 20 of its aligned bases and is untouched: none
of them has e >= 20.

## What the two tables show

`expected_raw.tsv` has **five** rows for four cytosines. chrB 0-based 14 appears twice,
once from the supplementary at `off_5p` 24 (sequencing cycle 24 of a 30-base read, the
leading hard clip included) and once from R2 at `off_5p` 15. Those are one molecule base.
`expected.tsv` has four rows, the survivor carrying Q40 — the two agreeing Q30
observations summed and capped — and `expected_scan.txt` pins that alnbase now reads two
records rather than three, which is the dropped supplementary.

## What this probes in other tools

This is the BC-13 geometry. BISCUIT's overlap handling, like bamUtil's and fgbio's, works
from the mate pair's reference intervals; here they are on different contigs, so the
doubly observed base at chrB:15 is counted twice. samtools mpileup does not act on this
shape either (`docs/reference/05-overlap.md` §9.2 measures that on the equivalent
template). The failure is not a threshold anyone can loosen: the information needed is in
a record the interval comparison never looks at.

The records carry neither `MC` nor `SA`, on purpose. alnbase groups the template by QNAME
and derives every coordinate from FLAG and CIGAR, so it depends on no aligner-written tag
— the same argument as `overlap_agreeing_mates` makes for BC-11, extended to the chimeric
case. Methyl-HiC and other ligation-based bisulfite protocols produce this geometry
routinely, not as an edge case.
