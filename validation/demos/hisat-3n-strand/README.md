# HISAT-3N's `YZ`, checked against HISAT-3N — and what an inferred strand costs

[`queries/strand/hisat-3n.toml`](../../../queries/strand/hisat-3n.toml) takes the conversion
strand from `YZ:A` and the sequenced direction from FLAG 0x10. This runs HISAT-3N
2.2.1-3n-0.0.3, built from source at commit `f5dda37`, on the four-strand fixture the other
strand demos use, and checks the rule against what the aligner actually wrote.

`YZ` is unlike the other one-character conversion tags in this collection. BISCUIT's and
bwa-meth's name the converted contig the read was aligned to, a fact about the search.
HISAT-3N aligns to both three-letter indexes and then decides, by counting how many C→T and
how many G→A differences the alignment shows (`alignment_3n.h:288-313`). That is an
inference, and the demo is built around the case where it has nothing to infer from.

```
HISAT3N_ENV=<env> PY=<python with pyarrow> ALNBASE=<alnbase> ./run.sh [workdir]
```

The environment is [`../../envs/hisat-3n.yaml`](../../envs/hisat-3n.yaml), built by
[`../../envs/build_hisat3n.sh`](../../envs/build_hisat3n.sh). HISAT-3N is a branch of hisat2
rather than a separate project and is not packaged, so it is built from source.

## What is aligned

Four 200-base fragments, one per strand of origin, read from both ends — the same fixture as
[`../bsmap-strand/`](../bsmap-strand/) and [`../bwameth-strand/`](../bwameth-strand/), except
that each pair is written three times, at three amounts of conversion evidence:

| set | reads | C→T count | G→A count |
|---|---|---|---|
| `part` | CG methylated, every CH converted | 5–13 | 0 |
| `full` | nothing converted at all | 0 | 0 |
| `tie` | one converted cytosine and one G→A substitution | 1 | 1 |

`full` is what a fully methylated read looks like, and `tie` what a read with one conversion
and one sequencing error or SNP looks like. Both make the counts equal, which is the branch
that falls back to the directional assumption (`alignment_3n.h:299-308`).

## What HISAT-3N wrote

On `part` the rule is right on every record: `YZ` is `+` on the OT and CTOT pairs and `-` on
the OB and CTOB ones, identical across both mates of every pair as a conversion-strand tag
must be, 0x10 agreed with how SEQ was really stored on all eight, and SEQ was
reference-forward throughout. Nothing was unmapped, so nothing was skipped.

On `full`, four of the eight records name the other conversion strand:

```
pair  really is origin FLAG  YZ  Yf  rule says
CTOB  read 1    CTOB   99    +   0   +:forward (want -:forward)
CTOT  read 1    CTOT   83    -   0   -:reverse (want +:reverse)
OB    read 1    OB     83    -   0   -:reverse (want -:reverse)
OT    read 1    OT     99    +   0   +:forward (want +:forward)
```

The pattern is the fallback, not a mis-alignment: every record is at the right position with
MAPQ 60, and the wrong ones are exactly the two pairs whose read 1 came from a PCR copy. With
no conversions to count, HISAT-3N assigns the strand read 1 would have had in a directional
library — which is right for OT and OB and backwards for CTOT and CTOB. Aligning the same
reads with `--directional-mapping` gives **the same `YZ` on all eight records**: in the
absence of evidence the default non-directional mode *is* the directional mode, applied
without having been asked for.

`tie` produces the same four wrong records with `Yf:i:1` instead of `Yf:i:0`.

## How much evidence the inference needs

A 60 bp read from each of the four strands every 10 bases, aligned single-end, in each of the
four states. Every read mapped, at the right position, unclipped, MAPQ 60:

| state | reads | conversion strand wrong | of those, carrying `Yf:i:0` |
|---|---|---|---|
| `part` | 456 | 0 | 0 |
| `one` — a single converted cytosine | 456 | 0 | 0 |
| `full` — nothing converted | 456 | **228** | 228 |
| `tie` — one conversion of each kind | 456 | **228** | 0 |

One converted cytosine is enough: the inference is not a majority vote that needs depth, it
needs a single difference the other count does not match. When it has that it is right on
every read. When it does not, it is wrong on every copy-strand read — half the library, since
the two strands it gets right are the two a directional library would have contained.

**`Yf` marks one of the two ties and not the other.** `Yf` is the number of conversions of
the kind `YZ` names, so `Yf:i:0` says outright that no evidence was seen — a read whose strand
was assumed. That covers the fully methylated case exactly. It does not cover the equal-count
case, where `Yf:i:1` is indistinguishable from a genuine single conversion.

## alnbase, reading the rule and nothing else

- **A wrong conversion strand moves every call.** alnbase walks each read along the strand
  the rule declared, so on the 4 mis-tagged records of `full` all 27 of their CG hits anchor
  on the other base of the CG: a G under the CTOT records, which should have been walked as
  `+`, and a C under the CTOB ones. The hits are not dropped or flagged; they are emitted at
  coordinates one base away, on the opposite strand, with the same confidence as the correct
  ones. That is what an inferred strand costs, stated in the coordinates a caller would emit.

- **The check that caught BSMAP cannot catch this.** alnbase's scan line counts read bases
  that differ from the reference, allowing whichever conversion the declared strand permits.
  On `full` it reads `0 of 480` while half the records are on the wrong strand: a read with no
  conversion matches the reference whichever way it is walked, so there is nothing to count.
  On `tie` it reads `8 of 480` — but splitting the fixture into the records the rule got right
  and the ones it got wrong gives `4 of 240` for each half. What it is counting is the planted
  substitution, which both halves carry. Neither figure separates a wrong strand from a right
  one.

- Between them `Yf:i:0` and the scan line cover both ties, and neither covers both. The
  practical rule that follows is about the library rather than the record: pass
  `--directional-mapping` when the library is directional, so the assumption is one you made
  on purpose.

## A note on the source

Two comments in `alignment_3n.h` describe the opposite of what the code does. The header of
`makeYZ` says "if the conversion type 0 is less, the read is mapped to REF (+)" (`:285-286`)
and the declaration of `YZ` says `+` is for "conversionCount[0] is equal or smaller than
conversionCount[1]" (`:95-96`), while the code assigns `+` when `conversionCount[0] >=
conversionCount[1]` (`:309`) — that is, when the C→T count is the larger one. The code is
right, and this run confirms it: `part` reads, which carry only C→T differences, all come back
`YZ:A:+` on the OT side. Anyone writing a strand rule from those comments rather than from the
code would invert the tag.

## Caveats, and what is not exercised

- Only `--base-change C,T`. HISAT-3N is a general three-letter aligner and TAPS, SLAM-seq and
  other conversions go through the same `makeYZ`; the polarity of `YZ` for a different base
  change is not measured here.
- `--directional-mapping-reverse`, the undocumented PBAT mode, is not run.
- No spliced alignments: the fixture is genomic and `--no-spliced-alignment` is passed.
  HISAT-3N's RNA mode, which is its main use, would put introns in the walk.
- No gaps, no trimming, no repeats. The reference is 1200 unique bases, so nothing here
  measures what a low-MAPQ or multi-mapping read does to the inference.
- `hisat-3n-table`, HISAT-3N's own caller, is not run. This demo is about the tag, not about
  agreeing with the extractor that reads it.
