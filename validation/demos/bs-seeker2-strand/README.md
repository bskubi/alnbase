# BS-Seeker2's strand tag, checked against BS-Seeker2

`queries/strand/bs-seeker2-se.toml` was written by reading BS-Seeker2's source, and until
now nothing had run a BS-Seeker2 BAM through it. That is the one thing the records files
in `queries/strand/records/` cannot do: they are written from the same reading as the rule,
so they catch a rule that stops agreeing with its own documentation, not an aligner that
writes something other than what its source appeared to say. This directory runs the
aligner.

Tool: BS-Seeker2 v2.1.8, GitHub master `0976fee` (`../../envs/bs-seeker2.yaml` and
`../../envs/build_bsseeker2.sh`), with bowtie2 and python 2.7 + pysam 0.8.4. Checked on
2026-09-17 with alnbase 0.1.21.

```
BSSEEKER2_ENV=... ALNBASE=... ./run.sh [workdir]
```

`expected.txt` is the output of that run. `origin.toml` is a one-tag query file that writes
the strand alnbase settled on into `YV`, so each record's call can be read off the BAM.

## The case

Four reads, one per strand of origin, each built from a known 60-base window of an 800-base
reference and **named for the strand it was built from**, which is the answer the rule has
to produce. They are aligned single-end with `-t Y`, BS-Seeker2's non-directional mode:
a directional run reports only the two `FW` classes, so the other two `XO` values would
never appear and half the rule would go unexercised.

| read | built as | `XO` written | the rule's call | alnbase's `YV` |
|---|---|---|---|---|
| OT | the top strand, bisulfite-converted | `+FW` | OT | OT |
| CTOT | the reverse complement of that | `+RC` | CTOT | CTOT |
| OB | the bottom strand, bisulfite-converted | `-FW` | OB | OB |
| CTOB | the reverse complement of that | `-RC` | CTOB | CTOB |

All four aligned at the window they came from, and each got the call its name demands. The
four-way reading of `XO` — first character the genome strand, second the read's direction
within it — is what BS-Seeker2 actually writes.

## The FLAG, which is why the rule reads the tag

`bs_align/output.py:54-61` sets `0x10` from the **genome strand** rather than from the
record's orientation, and the run confirms it:

| read | `XO` | FLAG written | FLAG a conformant aligner writes |
|---|---|---|---|
| OT | `+FW` | 0 | 0 |
| CTOT | `+RC` | **0** | 16 |
| OB | `-FW` | 16 | 16 |
| CTOB | `-RC` | **16** | 0 |

The two `RC` classes are inverted, exactly as reading the source predicted. `SEQ` is stored
reference-forward in all four cases — the demo checks that directly, by counting how many
of each record's 60 bases differ from the reference at `POS` once the two substitutions
bisulfite can make are allowed, and getting zero every time. So the records are placed and
oriented correctly; it is only the flag that lies about how the read was sequenced.

A rule built on `0x10` would therefore give the wrong sequenced direction for every `RC`
read, which is to say for every read from a complementary strand. BS-Seeker2's own caller
reads `is_reverse` and nothing else (`bs_seeker2-call_methylation.py:335-340`). Reading
`XO` avoids the question: it names the strand of origin, from which the walk direction and
the sequenced direction both follow, so alnbase never consults `0x10` on this input.

alnbase's own scan line is the independent check that this is right rather than merely
self-consistent: **0 of 240 compared read bases differ from the reference**. A record
walked in the wrong direction, or placed as the wrong strand, would disagree with the
reference at most of its bases rather than none of them.

## What this closes, and what it does not

`bs-seeker2-se.toml` is no longer marked `verify`, and its records file now carries the
FLAGs BS-Seeker2 writes rather than the conformant ones it was first given — a detail that
did not affect any call, since the rule reads only `XO`, but did misrepresent the aligner
in a file whose whole purpose is to say what the aligner writes.

Paired-end BS-Seeker2 still has no rule. Its `XO` values are `+FR`/`-FR`/`+RF`/`-RF`, one
per pair rather than per record (`bs_align/bs_pair_end.py:531-545,607-612`), so the mate bit
has to be combined with the tag; nothing here tests that.
