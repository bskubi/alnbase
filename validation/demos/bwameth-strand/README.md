# bwa-meth's `YD`, checked against bwa-meth — and what "directional only" costs

`queries/strand/bwameth.toml` reads `YD:Z`, one character naming the converted contig the
read landed on, so it reaches the conversion strand and no further; the file pairs it with
FLAG 0x10 for the sequenced direction. Three claims in its comment are measured here, and all
three hold on a directional library: 0x10 is conformant, `YC:Z` carries nothing the FLAG does
not, and `directional.toml` applies to the same BAM and agrees with it hit for hit while
additionally naming the strand of origin.

A fourth thing is measured because bwa-meth supports directional libraries only: what it does
when copy-strand reads arrive anyway. It drops 93% of them and places the rest **on the
opposite conversion strand**, at the right position, unclipped, at high MAPQ.

Tool: bwa-meth 0.2.10 from bioconda (`../../envs/bwameth.yaml`), aligning with `bwa mem`.
Checked 2026-09-17 with alnbase 0.1.24.

```
BWAMETH_ENV=... PY=python-with-pyarrow ALNBASE=... ./run.sh [workdir]
```

`expected.txt` is the output of that run.

## What is aligned

Bisulfite converts each single strand on its own and PCR then copies each converted strand,
so one fragment yields four sequences in two duplexes: `OT` with its copy `CTOT`, and `OB`
with its copy `CTOB`. A pair is one duplex read from both ends, so read 1 comes from one of
its two strands and read 2 from the other. Four 200 bp fragments on a 1200 bp reference with a
CG every 17 bases, one pair per strand of origin, each named for the strand its **read 1** was
sequenced from — the same fixture as `../biscuit-strand/` and `../bsmap-strand/`.

The `dir` run holds the `OT` and `OB` pairs, which is a directional library and all bwa-meth
claims to handle; the `nd` run adds the other two. Nothing in the "really is" columns comes
from the FLAG or from `YD`: which fastq read a record holds, and whether SEQ was
reverse-complemented on the way in, are settled by comparing SEQ against the reads themselves.

## What bwa-meth wrote

| pair | record | origin | FLAG | `YD` | `YC` | rule says |
|---|---|---|---|---|---|---|
| `OT` | read 1 | OT | 97 | `f` | `CT` | `+`, forward |
| `OT` | read 2 | CTOT | 145 | `f` | `GA` | `+`, reverse |
| `OB` | read 1 | OB | 81 | `r` | `CT` | `-`, reverse |
| `OB` | read 2 | CTOB | 161 | `r` | `GA` | `-`, forward |

0 of 4 wrong on either count. Three things in that table are the claims:

**`YD` is per fragment, not per record.** Both mates of each pair carry the same value, which
is what a conversion-strand tag should do: a strand and its own PCR copy hit the same
converted contig. That is also the split `YD` cannot cross, so the rule reports the strand of
origin as unknown and `conv_strand` comes out null.

**0x10 is conformant and is what the rule uses for direction.** bwa-meth leaves bwa's FLAG
untouched and restores SEQ from the original read, reverse-complemented when the record is
reverse (`bwameth.py:500-503`). On every mapped record 0x10 agreed with how SEQ was really
stored, and SEQ was reference-forward throughout, counting only differences bisulfite cannot
explain.

**`YC` carries nothing the FLAG does not.** It is `CT` on every first-in-pair record and `GA`
on every second-in-pair one, in both runs, because bwa-meth writes it from the mate number
while converting the fastq (`bwameth.py:205-209`) and never revisits it. It names the
conversion *applied to the read*, not the strand the read came from, so a rule that read it as
a strand would be right only on a directional library, and right there only by coincidence.

## What "directional only" costs

bwa-meth converts read 1 C→T and read 2 G→A on the assumption that read 1 came from an
original strand. A `CTOT` or `CTOB` read arriving as read 1 gets the wrong conversion applied.
Sweeping a 60 bp read from each strand of origin every 10 bases along the reference, aligned
single-end so every read is treated as read 1:

| strand of origin | reads | mapped | conversion strand wrong | misplaced | clipped | median MAPQ |
|---|---|---|---|---|---|---|
| OT | 114 | 114 | 0 | 0 | 0 | 60 |
| CTOT | 114 | 15 | **15** | 0 | 0 | 56 |
| OB | 114 | 114 | 0 | 0 | 0 | 60 |
| CTOB | 114 | 17 | **17** | 0 | 0 | 52 |

93% of the copy-strand reads are dropped, which is the honest outcome and what "directional
only" is meant to mean. The other 7% is the part worth knowing: **every copy-strand read that
maps at all is on the wrong conversion strand**, and nothing in its record says so. It is at
the right position, its CIGAR is a full-length match, its MAPQ is in the fifties, and no record
is flagged QC-fail. A read that survives does so precisely because the wrong conversion made it
look like a clean hit on the other contig.

This is not a bug in the sense of `../bsmap-strand/`, where BSMAP discarded a perfect
alignment it had already found. bwa-meth is documented as directional-only and the reads above
are outside what it accepts. The measurement is of what the documented limitation does in
practice: not a clean refusal but a 7% leak of silently wrong records, which is a reason to
filter a library by design rather than trust the aligner to drop what it cannot handle.

In the paired-end `nd` run the same thing shows at record scale: of the eight records, three
are unmapped and carry no `YD` at all — bwa-meth tags only mapped records
(`bwameth.py:469-473` returns early for the rest) — and one of the five that mapped is the
`CTOT` pair's read 1, called `-` where it should be `+`.

## alnbase, reading the rule and nothing else

The CG query in `cg.toml` anchors on a cytosine in the read whose reference context is CG.
Which reference base that anchor lands on is decided by the conversion strand alnbase took
from the rule: the C of a forward CG for a `+` record, the G of the same CG for a `−` one.
Reading those bases back out of the reference checks the walk rather than the label.

- 29 hits over 4 records (`dir`) and 37 over 5 (`nd`); **0 hits whose anchor is not the
  informative base for the strand alnbase walked**, under either rule, in either run.
- **The records the rule cannot classify are skipped and counted, not guessed at.** The `nd`
  run reads 8 records and reports `scanned 5, skipped 3 (given a null row)` — the three with
  no `YD`. The count is the signal that the rule does not cover this input, which is the point
  of declaring the rule rather than compiling it in.
- **`bwameth.toml` and `directional.toml` agree hit for hit**: 0 disagreements out of 29 and
  out of 37, keyed on record and reference position, so a strand disagreement would move an
  anchor from a C to the G beside it and show up. `directional.toml` additionally fills in
  `conv_strand` — OT, CTOT, OB, CTOB — which `YD` cannot reach. The agreement is worth reading
  precisely: both rules trace back to the contig bwa-meth chose, one through `YD` and one
  through the 0x10 that same choice set, so they agree on the misassigned `CTOT` record too.
  Two rules agreeing says the BAM is self-consistent, not that it is right.
- The scan line separates the two runs: `0 of 240 compared read bases differ from the
  reference` on `dir`, `5 of 300` on `nd`. Those 5 are the misassigned read's bisulfite
  conversions, counted as mismatches because alnbase was told to walk the other strand. Both
  rules give the same 5, for the same reason they agree everywhere else.

## Caveats, and what is not exercised

- **`--set-as-failed` is not used.** It marks one direction's records 0x200 for targeted
  libraries; whether that interacts with `YD` is untested here.
- **The chimera heuristic is not exercised.** bwa-meth fails a record QC and un-pairs it when
  its longest match is under 44% of the read (`bwameth.py:485-492`); no record here was near
  that, and `records flagged QC-fail (0x200): 0`.
- **bwa-mem2 is not run**, only `bwa mem`, bwameth.py's default. The tagging code is shared,
  but the alignments are not necessarily identical.
- **`MethylDackel` is not run here.** What bwa-meth's usual caller makes of a leaked
  copy-strand record is a separate question from what the aligner wrote.
- Gapped alignment and trimming are off, so `YD` has only been seen on ungapped full-length
  records.
