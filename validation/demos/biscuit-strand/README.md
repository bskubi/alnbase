# BISCUIT's conversion tag and its conformant FLAG, checked against BISCUIT

`queries/strand/biscuit.toml` is the first rule here that cannot name a read's strand of
origin. BISCUIT's `YD` distinguishes only two states, so the file declares the conversion
strand from the tag and the sequenced direction from FLAG 0x10 — two tables instead of one,
and `conv_strand` null on every row. That makes it a different claim from the four-way rules:
not "the tag names the strand" but "the tag names the conversion strand **and** this
aligner's 0x10 can be trusted for the rest". Both halves are measured here.

Tool: BISCUIT 1.10.3, GitHub commit `0a5ceae`, built by `../../envs/build_biscuit.sh` into
the environment from `../../envs/biscuit.yaml`. Checked 2026-09-17 with alnbase 0.1.24.

```
BISCUIT_ENV=... PY=python-with-pyarrow ALNBASE=... ./run.sh [workdir]
```

`expected.txt` is the output of that run.

## What is aligned

Bisulfite converts each single strand on its own and PCR then copies each converted strand,
so one fragment yields four sequences in two duplexes: `OT` with its copy `CTOT`, and `OB`
with its copy `CTOB`. A pair is one duplex read from both ends, so read 1 comes from one of
its two strands and read 2 from the other. Four 200 bp fragments are built on a 1200 bp
reference with a CG every 17 bases, one pair per strand of origin, each pair named for the
strand its **read 1** was sequenced from.

A directional library holds the `OT` and `OB` pairs alone, and those two already show both
conversion strands and both sequenced directions, because read 2 of each is the
complementary strand. So the run aligns those two with `-b 1` and all four with the default
`-b 0`, which is the mode `biscuit.toml` was written for.

Nothing in the "really is" columns comes from the FLAG: which fastq read a record holds, and
whether SEQ was reverse-complemented on the way in, are both settled by comparing SEQ against
the reads themselves, because the FLAG is part of what is under test.

## What BISCUIT wrote

| pair | record | origin | FLAG | `YD` | conversion | sequenced |
|---|---|---|---|---|---|---|
| `OT` | read 1 | OT | 97 | `f` | + | forward |
| `OT` | read 2 | CTOT | 145 | `f` | + | reverse |
| `OB` | read 1 | OB | 81 | `r` | − | reverse |
| `OB` | read 2 | CTOB | 161 | `r` | − | forward |

Every record's conversion strand and sequenced direction is the one its strand of origin
implies (`../../../docs/design/strand-rules.md`): 0 of 4 wrong with `-b 1`, 0 of 8 wrong
with `-b 0`. SEQ is reference-forward on all twelve records, counting only differences that
bisulfite cannot explain.

**The FLAG really is conformant.** 0x10 agrees with how SEQ was actually stored on every
record in both runs — which is the claim that lets this file read 0x10 at all, and is exactly
what Bismark single-end, BSBolt and BS-Seeker2 violate (`../bsbolt-strand/`,
`../bs-seeker2-strand/`).

**The contig grouping was counted, not assumed.** `YD:A:f` held `OT` and `CTOT`; `YD:A:r`
held `OB` and `CTOB`. A strand and its PCR copy hit the same converted contig, which is
precisely why a two-way tag cannot reach the four-way answer. The file's comment had this
wrong — it said `OT` and `CTOB` share the `f` contig, which is the grouping by *sequenced
direction*, not by contig — and the demo corrected it.

## alnbase, reading the rule and nothing else

The CG query in `cg.toml` anchors on a cytosine in the read whose reference context is CG.
Which reference base that anchor lands on is decided by the conversion strand alnbase took
from the rule: the C of a forward CG for a `+` record, the G of the same CG for a `−` one.
Reading those bases back out of the reference checks the walk rather than the label.

- 29 hits over 4 records (`-b 1`) and 56 over 8 (`-b 0`); **0 hits whose anchor is not the
  informative base for the strand alnbase walked**.
- `conv_strand` is null on every row, and **0 records were given a strand of origin** — the
  rule names two strands, so alnbase declines to invent the other split rather than reporting
  a guess as a measurement.
- `read_reverse` follows 0x10 here, as the rule's `[strand.sequenced]` table says it should.
- The scan line — `0 of 240` and `0 of 480 compared read bases differ from the reference` —
  says the records were placed and walked correctly rather than merely labelled consistently.

## Caveats, and what is not exercised

- **FLAG 0x2 is absent from every record, and that is not a BISCUIT behaviour.** BISCUIT's
  aligner is a bwa derivative that sets the proper-pair bit from an insert-size distribution
  it estimates per batch, and with two pairs it has none: `dir.log` says `[M:mem_pestat]
  There are not enough pairs for insert size inference`. A real library carries the bit
  normally. `queries/strand/records/biscuit.sam` therefore keeps its conformant
  99/147/83/163; apart from 0x2 those are the FLAGs BISCUIT wrote here.
- **`YD:A:u` is not exercised.** It is emitted when no conversion event was seen anywhere in
  the alignment, and every read here carries many. Those records are the rule's `unknown`
  key, which alnbase counts and skips.
- Single-end BISCUIT (`-b 1`/`-b 3`/`-b 0` for SE), `biscuit pileup`, and the `-f`
  strand-restriction options are not exercised. `YD` is per record and the rule never looks
  at the pairing bits, so single-end is covered by the same reasoning as paired-end, but it
  has not been run.
