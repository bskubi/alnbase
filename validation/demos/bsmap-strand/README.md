# BSMAP's `ZS`, checked against BSMAP — and a conversion strand it gets wrong

`queries/strand/bsmap.toml` reads `ZS:Z`, whose two characters are the conversion strand and
the read's orientation within it, so unlike BISCUIT's `YD` it names all four strands of
origin from the tag alone. Two things are checked here: that the four values mean what the
file says, and the claim in its comment that a paired-end run emits all four of them even at
`-n 0`, the default, which the help text describes as mapping "only to 2 forward strands".

Both hold. A third thing turned up on the way: **at `-n 1`, BSMAP puts some copy-strand reads
on the opposite conversion strand**, at the right position, reported as a unique alignment.
That is written up below and in `../../../docs/design/research/strand-determination.md`.

Tool: BSMAP 2.90 from bioconda (`../../envs/bsmap.yaml`), whose source tarball is the one the
recipe pins, sha256 `8b5ae4ba…f007ac2`. Checked 2026-09-17 with alnbase 0.1.24.

```
BSMAP_ENV=... PY=python-with-pyarrow ALNBASE=... ./run.sh [workdir]
```

`expected.txt` is the output of that run.

## What is aligned

Bisulfite converts each single strand on its own and PCR then copies each converted strand,
so one fragment yields four sequences in two duplexes: `OT` with its copy `CTOT`, and `OB`
with its copy `CTOB`. A pair is one duplex read from both ends, so read 1 comes from one of
its two strands and read 2 from the other. Four 200 bp fragments are built on a 1200 bp
reference with a CG every 17 bases, one pair per strand of origin, each pair named for the
strand its **read 1** was sequenced from — the same fixture as `../biscuit-strand/`.

A directional library holds the `OT` and `OB` pairs alone, so the run aligns those two at the
default `-n 0` and all four at `-n 1`. Nothing in the "really is" columns comes from the FLAG
or from `ZS`: which fastq read a record holds, and whether SEQ was reverse-complemented on the
way in, are settled by comparing SEQ against the reads themselves.

## What BSMAP wrote

| pair | record | origin | FLAG | `ZS` | rule says |
|---|---|---|---|---|---|
| `OT` | read 1 | OT | 99 | `++` | OT |
| `OT` | read 2 | CTOT | 147 | `+-` | CTOT |
| `OB` | read 1 | OB | 83 | `-+` | OB |
| `OB` | read 2 | CTOB | 163 | `--` | CTOB |

**All four `ZS` values appear at `-n 0`**, from two pairs, because read 2 of each is searched
against the complementary chain. So the four-way rule is needed for a default BSMAP run and
not only for `-n 1`; a two-way reading of `ZS` would be wrong on half of a directional
library's records. 0 of 4 strands of origin wrong.

**The FLAG is conformant and the rule does not use it.** 0x10 agrees with how SEQ was really
stored on all twelve records, and is set exactly when `ZS`'s two characters differ, with no
exceptions. SEQ is reference-forward throughout, counting only differences bisulfite cannot
explain.

## The conversion strand BSMAP gets wrong

At `-n 1` one of the eight records is wrong: the `CTOT` pair's read 1, a CTOT-strand read,
comes back `ZS:Z:-+` — OB, the opposite conversion strand — at the correct position, with
`NM:i:5` and MAPQ 255. Its perfect alignment, `+-` with `NM:i:0`, exists and is never
reported. It is not a tie broken badly: `-r 2`, `-r 0`, `-v 10`, `-s 12` and every `-S` seed
including the clock give the same answer, single-end as well as paired.

Sweeping a 60 bp read from each strand of origin every 10 bases along the reference, 456 reads
aligned single-end at `-n 1`:

| strand of origin | aligned | strand of origin wrong |
|---|---|---|
| OT | 114 | 0 |
| CTOT | 114 | 5 |
| OB | 114 | 0 |
| CTOB | 114 | 5 |

Only the copy strands are affected, every one of them is called as the strand with **both**
`ZS` characters flipped (CTOT as OB, CTOB as OT), and none is misplaced. The mechanism is
that flip. BSMAP searches a read in two chains — as sequenced and reverse-complemented
(`align.cpp:92-93`) — against two converted reference chains, and the two indices become `ZS`
(`align.cpp:269`). Flipping both puts the read at the same forward position, so the two
candidates collide, and `AddHit` (`align.h:232-249`) keeps whichever arrived first:

```c
if(!hitset[_ghit.chr>>1].insert(_ghit.loc).second) return 0; //hit already exist
```

`>>1` drops the reference-chain bit, so the duplicate key is the forward position alone,
with no mismatch comparison. The as-sequenced chain is searched first (`align.cpp:211`), and
for a copy-strand read that is the wrong chain. So whenever a CTOT or CTOB read happens to
fit the opposite converted chain within the mismatch budget, the wrong-strand hit is
registered first and the perfect one is discarded as a duplicate. An OT or OB read is never
affected, because for it the chain searched first is the right one.

The budget is 5 mismatches for a 60 bp read at the default `-v 0.08` (`align.cpp:501`).
Predicting a misassignment from "the wrong-chain alignment fits in 5 mismatches" agrees with
what BSMAP did on 226 of the 228 copy-strand reads; the two exceptions sit exactly on the
boundary, where whether the candidate is found at all depends on the seed segments.

What this costs a caller is the whole read: every cytosine in it is read on the wrong strand.
Nothing in the BAM marks it — the position is right, the alignment is called unique, and only
the inflated `NM` hints at it. A directional library at `-n 0` is safe, since copy-strand
reads are outside the search space there; a non-directional or PBAT run at `-n 1` is not.

## alnbase, reading the rule and nothing else

The CG query in `cg.toml` anchors on a cytosine in the read whose reference context is CG.
Which reference base that anchor lands on is decided by the conversion strand alnbase took
from the rule: the C of a forward CG for a `+` record, the G of the same CG for a `−` one.
Reading those bases back out of the reference checks the walk rather than the label.

- 29 hits over 4 records (`-n 0`) and 56 over 8 (`-n 1`); **0 hits whose anchor is not the
  informative base for the strand alnbase walked**, in both runs — including the misassigned
  record, which alnbase walks exactly as `ZS` told it to.
- `conv_strand` is filled in on every row, which is what a four-way rule buys over
  `../biscuit-strand/`'s two-way one, and every value is the pair's strand or its mate's —
  except the 8 hits of the misassigned record, which is BSMAP's error faithfully reported.
- The scan line is the independent check, and it separates the two runs: `0 of 240 compared
  read bases differ from the reference` at `-n 0`, but `5 of 480` at `-n 1`. Those 5 are the
  misassigned read's bisulfite conversions, counted as mismatches because alnbase was told to
  walk the other strand. A run that reports a mismatch rate above zero on data this clean is
  telling its user to look.

## Caveats, and what is not exercised

- **RRBS mode (`-D`) is not exercised**, and it takes a different code path through
  `SnpAlign` with its own hit bookkeeping; the misassignment above was measured in WGBS mode
  only.
- **`methratio.py`, BSMAP's own caller, is not run here.** Whether it inherits the wrong
  strand or re-derives one is a separate question from what the aligner wrote.
- **BSMAPz and BSMAP 2.6 are not run.** `queries/strand/bsmap.toml` claims `ZS` is the same
  in all three from source; only 2.90 has been executed.
- Gapped alignment (`-g`), adapter trimming and quality trimming are off, so `ZS` has only
  been seen on ungapped full-length records.
