# BSBolt's strand tag and its inverted FLAG, checked against BSBolt

`queries/strand/bsbolt.toml` makes an unusually strong claim for a file written by reading
source: that BSBolt's FLAG 0x10 is *wrong* — set from the converted reference contig the
read matched rather than from the read's own orientation — on every `*_G2A` record, which
in a directional paired-end library is read 2 of every pair. Two other things in the repo
rest on that claim: the file exists in the shape it does because of it, and
`../premethyst-bugs/README.md` cites it. It had only ever been read, never measured. This
directory runs the aligner.

Tool: BSBolt 1.6.0, GitHub master `ea4870e` (`../../envs/premethyst.yaml` with
`../../envs/build_bsbolt.sh`), with samtools. Checked on 2026-09-17 with alnbase 0.1.23.

```
BSBOLT_ENV=... ALNBASE=... ./run.sh [workdir]
```

`expected.txt` is the output of that run. `origin.toml` is a one-tag query file that writes
the strand alnbase settled on into `YV`, so each record's call can be read off the BAM.

## The case

Four fragments, one per strand of origin, each 200 bases of a 1200-base reference, and each
pair **named for the strand its read 1 was sequenced from**. Read 2 then comes from the
other strand of the same duplex, so the answer the rule owes each record is fixed in advance
by construction.

A directional library sequences read 1 from the converted original only, so it holds the OT
and OB pairs alone — and those two already exercise all four of BSBolt's `YS` values,
because read 2 of each is the complementary strand. The `-UN` run adds the two pairs whose
read 1 came from the copy. Both runs were made, because `bsbolt.toml` claims the rule holds
in both.

All records aligned, and **every one was called as the strand it was sequenced from**, in
both runs.

## The FLAG, which is the reason to read the tag

| pair | record | FLAG | what 0x10 says | how SEQ is stored | `YS` | rule's call |
|---|---|---|---|---|---|---|
| OT | read 1 | 65 | forward | forward | `W_C2T` | OT |
| OT | read 2 | 129 | forward | **reverse-complemented** | `W_G2A` | CTOT |
| OB | read 1 | 113 | reversed | reverse-complemented | `C_C2T` | OB |
| OB | read 2 | 177 | reversed | **forward** | `C_G2A` | CTOB |

The `-UN` run adds the CTOT and CTOB pairs, which produce the same four FLAG values against
the mirrored `YS` values, for eight records in all.

The bit is inverted on **every `*_G2A` record and no other**: 2 of 2 in the directional run,
4 of 4 under `-UN`. The demo measures this rather than asserting it — it works out how SEQ
was actually stored by comparing it against the fastq read itself, never from the FLAG,
since the FLAG is what is in question.

`SEQ` is reference-forward in all of them, checked by counting how many of each record's 60
bases differ from the reference at `POS` once the two substitutions bisulfite can make are
allowed, and getting zero every time. That is the whole reason a `YS`-only rule works: the
tag names the strand, the record is already placed and oriented, and 0x10 need never be
consulted. alnbase's own scan line is the independent check that the calls are right rather
than merely self-consistent: **0 of 240 compared read bases differ from the reference** in
the directional run and **0 of 480** under `-UN`. A record walked in the wrong direction, or
placed as the wrong strand, would disagree with the reference at most of its bases rather
than none of them.

0x20 is not set on either Watson record, matching the source (`bwamem.c:854,868` sets the
mate's `is_rev` from the *aligned* contig too, and leaves it 0 on the Watson branch). Nothing
in the rule depends on it.

## What this closes, and one thing it does not show

The claim `bsbolt.toml` was built around is now measured on BSBolt's own output, in both the
directional and the `-UN` layout, and `../premethyst-bugs/README.md`'s note that BSBolt
paired-end BAMs need their read-2 FLAGs repaired before alnbase reads them is superseded:
they do not, as long as `bsbolt.toml` is the strand rule.

FLAG 0x2 is absent from every record here, but that is **not** a BSBolt behaviour and should
not be read as one. bwa sets 0x2 from an insert-size distribution it estimates per batch, and
with two pairs it has none to estimate from — `dir.log` says so directly ("skip orientation
FR as there are not enough pairs"). A real library would carry the bit normally.

Not covered: unmapped records, whose `YS` is `WC` and which the rule sends to `unknown`; and
BSBolt's own `CallMethylation`, which is a separate comparison.
