# The twelfth rule: a PBAT rule against Bismark `--pbat`, and why there is no PBAT file

`queries/strand/pbat.toml` was the last of the twelve strand rules to be run against real
output, and the only one its run did not confirm. **It is no longer shipped**; this round is
why. `run.sh` keeps a copy of its four conditions inline so the comparison stays reproducible.
Like [`directional.toml`](../../../queries/strand/directional.toml) it read no tag: it was a
claim about the **library** — that read 1 came off a copy strand and read 2 off the converted
original, the opposite way round from a directional library — so it applies to any conformant
aligner and rests on none of the evidence any of them writes.

The round was set up the way the [directional one](../directional-strand/) was, because the
same problem applies: with no tag to read there is nothing in the BAM to compare the rule
against. Bismark supplies the missing evidence, since `XR`/`XG` names all four strands from
the read's own conversion and the index that won, and the fixture knows which strand each read
was built from. Every record therefore gets three answers.

```
BISMARK_ENV=<env> PY=<python with pyarrow> ALNBASE=<alnbase> ./run.sh [workdir]
```

`BISMARK_ENV` is the environment from [`../../envs/bismark.yaml`](../../envs/bismark.yaml):
Bismark 0.25.1, bowtie2 and samtools.

The fixture is the 1200-base reference every strand demo uses, read as a PBAT library: the
directional demo's two mate files exchanged, which is literally what `methylpy --pbat` does to
a run. The same four strands exist either way; which member of each duplex is sequenced first
is the whole difference between the two libraries.

## What it found

**A PBAT rule is the wrong file for Bismark `--pbat` output, in both modes, and the shipped
file's own comment said otherwise** — it named "Bismark `--pbat` paired-end" as an input.
**Use [`bismark.toml`](../../../queries/strand/bismark.toml) for Bismark `--pbat` in either
layout**; `directional.toml` is also correct for the paired-end arm.

| arm | records | PBAT rule vs. truth | `directional.toml` vs. truth | `bismark.toml` vs. truth |
|---|---|---|---|---|
| paired-end | 402 | **402 wrong** | 0 wrong | 0 wrong |
| single-end | 201 | **201 wrong** | 201 wrong | 0 wrong |

### Paired-end: Bismark swaps the mate bits, so a PBAT run looks directional

Bismark's paired-end FLAG table (`bismark:8818-8866`) exchanges 0x40 and 0x80 on CTOT and CTOB
pairs, so that browsers do not discard them as discordant — the same swap already measured on
`--non_directional` output in [`../bismark-nondirectional-pe/`](../bismark-nondirectional-pe/),
where it affected four records of eight. Under `--pbat` **every** pair is CTOT or CTOB, so the
swap is universal: the record holding the FASTQ's read 1 is written with 0x80.

```
   read       really is FLAG  XR/XG truth  pbat   directional bismark  POS    mism
   probe_CTOT read 1    147   GA/CT CTOT   OB     CTOT        CTOT     241    0
   probe_CTOT read 2    99    CT/CT OT     CTOB   OT          OT       101    0
```

The mate bit is the only thing the PBAT rule has to tell a PBAT library from a directional one,
and Bismark has already inverted it, so the file is wrong on every record — and wrong in the
expensive direction, naming the opposite **conversion** strand, which reverses the walk:

```
-- pe, pbat.toml         4123 of 24120 compared read bases differ from the reference (17.09%)
-- pe, directional.toml     0 of 24120 compared read bases differ from the reference (0.00%)
```

alnbase's mismatch count is what surfaces it. On this clean fixture the right rule gives 0.00%
and the wrong one 17.09%, and the CG anchors move with it: 5228 hits differ between
the PBAT rule and `bismark.toml`, which is every hit of both, while `directional.toml` and
`bismark.toml` put 2601 anchors in identical places and name every one the same strand.

So **`directional.toml` is the file for a Bismark `--pbat` paired-end BAM**, despite the
library being PBAT. That is not a contradiction: a strand rule describes a BAM, and this BAM
has had its mate bits rewritten into the directional arrangement.

### Single-end: the lying 0x10 the directional round could not see

`bismark.toml`'s comment says Bismark sets 0x10 on single-end output from the **genome
conversion** rather than from SEQ's orientation. The [directional round](../directional-strand/)
looked for that and measured 0 contradictions over 202 records, because on a directional
library the two ways of setting the bit coincide. A PBAT library is exactly where they part
company, and here they do, on every record:

```
   records whose 0x10 contradicts how SEQ was actually stored: 201
```

That is the first direct measurement of this behaviour in the suite. Its consequence for the
two FLAG rules is not the same:

- The PBAT rule reads 0x10 as SAM defines it, so it takes a CTOT record (stored
  reverse-complemented, but flagged forward because it matched the CT index) for a CTOB one.
  Opposite conversion strand, reversed walk, 17.08% of read bases differing from the reference.
- `directional.toml` is wrong about the **origin** on all 201 records — it calls a CTOT read
  OT — but right about the conversion strand on all of them, because Bismark's non-conformant
  bit *is* the conversion strand. Its walk is correct and its anchors land exactly where
  `bismark.toml` puts them: 0 hits anchored elsewhere, 2576 of 1288 hits differing only in the
  strand name. The cost is the name and the `off_5p`/`off_3p` mirroring that follows from it,
  not the positions.

Neither file is right, so **`bismark.toml` is the file for Bismark `--pbat` single-end output**,
which is what the deleted file's comment already said for this arm. The arm is here because a
warning that has been measured is worth more than one that has been asserted.

## Why the file was deleted rather than kept with a narrower comment

Its comment named four inputs and all four turned out to be wrong.

The two Bismark modes went by measurement, above. **methylpy `--pbat` went the same way, by
reading methylpy 1.4.7's source** — the same standard the rules themselves were written to,
though not a run:

- paired-end, `call_mc_pe.py:225-227` exchanges `read1_files` and `read2_files` before
  anything else happens, so the aligner sees the converted originals as read 1;
- single-end, `call_mc_se.py:440` splits the input through
  `utilities.split_fastq_file_pbat` (`utilities.py:817`), which **reverse-complements every
  read** as it writes the chunks.

Either way what reaches the BAM is in the directional arrangement, so `directional.toml` is the
file — and this corrects a second line of our own documentation, which said methylpy's output
was "`directional.toml`, or `pbat.toml` for a `--pbat` run."

That leaves HISAT-3N `--directional-mapping-reverse`, which writes `YZ`, so
[`hisat-3n.toml`](../../../queries/strand/hisat-3n.toml) reads it, and bwa-meth, which accepts
directional libraries only and has no PBAT mode. No input remained.

The rule was not *wrong about PBAT*: it is the correct reading of a conformant BAM from a PBAT
library under SAM's meaning of the flags. What the round found is that a PBAT library does not
reach a BAM as one, because both pipelines built for PBAT normalise it into the directional
arrangement first. So the file's only remaining effect on a real user would have been to
attract the reasoning "my library was PBAT, so I want the PBAT rule" and hand them the opposite
conversion strand on every record. A shipped file that is right in principle and wrong on every
BAM anyone has is a trap, so it was deleted, and this page and
[`../../../queries/strand/README.md`](../../../queries/strand/README.md) say where it went.
Anyone who does turn up with a conformant PBAT BAM has the four conditions here and in
`run.sh`; forking `directional.toml` is a two-line change.

## What this round does not establish

- **No aligner's PBAT output has confirmed the rule.** See above. The round confirms the
  rule's logic against truth only in the negative sense that its call is exactly the
  conversion-strand inverse of the right one on this BAM.
- **methylpy was read, not run.** The two lines cited are unambiguous, but no methylpy BAM was
  produced here.
- **The fixture is not a real PBAT library.** It is the directional fixture's mates exchanged,
  which reproduces the strand assignment but not PBAT's random priming, its fragment-length
  distribution or its adapter chemistry.
- **One read of 202 is missing from the single-end BAM** and from 402 rather than 404 records
  paired-end. `sweep_CTOB_0` starts at reference position 1, and Bismark needs the two bases
  *before* the alignment start to make the call for a record on the bottom-strand index
  (`bismark:4315-4322`), so it returns early and writes nothing. The run report still counts
  the read as a unique alignment and reports 100% mapping efficiency, so report and BAM
  disagree by one at a chromosome edge. It is an edge case, not a defect in the strand logic,
  and the directional fixture never hit it because its bottom-strand reads start elsewhere.
- No indels, no soft or hard clipping, no repeats, no multimappers, no unmapped records.
