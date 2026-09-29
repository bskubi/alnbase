# Bismark's `--non_directional` paired-end BAM, checked against the rule

`queries/strand/bismark.toml` was confirmed against real Bismark BAMs in both layouts, but
it carried one explicit exception: **`--non_directional` paired-end was not claimed.** In
that mode Bismark deliberately swaps the 0x40 and 0x80 bits on CTOT and CTOB pairs, so that
genome browsers do not throw them out as discordant (`bismark:8818-8867`). Reading the
source said `XR` and `XG` should survive the swap, because they are per-record — but a tag
that "should" survive a deliberate flag rewrite is exactly the kind of claim only the
aligner's own output can settle. This directory runs the aligner.

Tool: Bismark 0.25.1 with bowtie2 (`../../envs/bismark.yaml`). Checked on 2026-09-17 with
alnbase 0.1.22.

```
BISMARK_ENV=... ALNBASE=... ./run.sh [workdir]
```

`expected.txt` is the output of that run. `origin.toml` is a one-tag query file that writes
the strand alnbase settled on into `YV`, so each record's call can be read off the BAM.

## The case

Bisulfite converts each single strand on its own and PCR then copies each converted strand,
so one genomic fragment yields four sequences in two duplexes: OT with its copy CTOT, and
OB with its copy CTOB. A pair is one duplex read from both ends, so read 1 comes from one
of its two strands and read 2 from the other. Which way round is precisely what a
non-directional library leaves open.

Four fragments, one per strand of origin, each 200 bases of a 1200-base reference, and each
pair **named for the strand its read 1 was sequenced from**. Read 2 then comes from the
other strand of the same duplex, so the answer the rule owes each record is fixed in
advance by construction.

| pair | read 1 from | read 2 from |
|---|---|---|
| OT | OT | CTOT |
| CTOT | CTOT | OT |
| OB | OB | CTOB |
| CTOB | CTOB | OB |

All four pairs aligned, one to each of Bismark's four indexes, and **all eight records were
called as the strand they were sequenced from**. The `XR`/`XG` table in the rule holds
per record under `--non_directional`, which is what was in question.

## The swap, which is the reason to read the tags

| pair | record | FLAG | what 0x40/0x80 says | `XR`/`XG` | rule's call |
|---|---|---|---|---|---|
| OT | read 1 | 99 | read 1 | CT/CT | OT |
| OT | read 2 | 147 | read 2 | GA/CT | CTOT |
| CTOT | read 1 | 147 | **read 2** | GA/CT | CTOT |
| CTOT | read 2 | 99 | **read 1** | CT/CT | OT |
| OB | read 1 | 83 | read 1 | CT/GA | OB |
| OB | read 2 | 163 | read 2 | GA/GA | CTOB |
| CTOB | read 1 | 163 | **read 2** | GA/GA | CTOB |
| CTOB | read 2 | 83 | **read 1** | CT/GA | OB |

Four of the eight records — every record of the CTOT and CTOB pairs — carry a mate bit
naming the other mate. The demo counts this rather than asserting it: it works out which
fastq read each record actually holds by comparing `SEQ` against both reads of the pair, and
reports four of eight.

Two things follow. First, an OT pair and a CTOT pair produce **the same two BAM records**,
distinguished only by which fastq each came from; a `--non_directional` paired-end BAM is
therefore indistinguishable, record by record, from a directional one, which is why
`queries/strand/records/bismark.sam` needs no new rows for this mode. Second, anything that
reads 0x40 to decide which strand a mate came from — "read 1 is OT, read 2 is CTOT" — has
those four records exactly backwards. The tags do not move: `XR_tag_1` is printed on the
same line as `flag_1` (`bismark:9192-9214`), so each record keeps the conversion state of
the read it actually holds.

`SEQ` is reference-forward in all eight records, checked directly by counting how many of
each record's 60 bases differ from the reference at `POS` once the two substitutions
bisulfite can make are allowed, and getting zero every time.

alnbase's own scan line is the independent check that the calls are right rather than
merely self-consistent: **0 of 480 compared read bases differ from the reference**. A record
walked in the wrong direction, or placed as the wrong strand, would disagree with the
reference at most of its bases rather than none of them.

## What this closes, and what it does not

`bismark.toml` now claims all of Bismark's layouts: single-end and paired-end, directional
and `--non_directional`. Nothing here exercises `--pbat`, which is a single-end matter and
covered by the same argument as the other single-end modes, or `--strandID`, whose `YS:Z`
names the strand directly and would want its own rule file.
