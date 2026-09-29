# Strand rules

One file per aligner, each declaring how to recover a read's original strand from that
aligner's flags and tags. The syntax is specified in
[`../../docs/design/strand-rules.md`](../../docs/design/strand-rules.md).

**Status: in use.** Since 0.1.18 a run gets its strand from one of these files and from
nowhere else: pass the one that matches the BAM, like any other query file, and a run that
reads records without one is refused. `directional.toml` is the plain directional rule, for
an aligner that writes no strand tag; it is the rule alnbase used to apply on its own, and
every row and tag it produces is byte-for-byte what that rule produced.

```
alnbase query --query-file my-queries.toml --query-file queries/strand/bismark.toml in.bam ref.aref out.bam
```

Every file below was written by reading the aligner's source. A file marked `verify` has
not yet been run against that aligner's own output. None is so marked at present.

Every file here ships a `records/<name>.sam` holding one record per strand it names, each
record named for the call the rule must make. A test runs them on every build, so a rule
that stops agreeing with the table in its own comment fails rather than goes on producing
calls on the wrong strand. See [`records/README.md`](records/README.md), which also says
what that test cannot catch.

## Why files

Supporting a new aligner is writing one of these, not recompiling alnbase, and adapting one
is forking it — copy the file, change the conditions, pass your copy instead. Nothing about
any aligner is compiled into the binary, so there is no rule alnbase can apply that you
cannot read here, and none it can express that you cannot write.

A rule is a property of the BAM, not of the analysis, so pointing at the wrong file is a
real mistake: every file below says which aligner, version and mode it was read from, and
what it would get wrong if used elsewhere.

| File | Input | Evidence | Gives |
|---|---|---|---|
| `directional.toml` | any conformant aligner on a directional library | FLAG 0x10 and the mate bit | all four strands |
| `unconverted.toml` | any assay that converts nothing: variants, motifs, RNA modifications read as mismatches, damage | FLAG 0x10, for the read's direction only | neither — everything walks `+` |
| `bismark.toml` | Bismark 0.25.1, single-end or paired-end, directional or `--non_directional` | `XR` + `XG` | all four strands |
| `bsbolt.toml` | BSBolt 1.6.0, single-end or paired-end, directional or `-UN` | `YS` | all four strands |
| `bsmap.toml` | BSMAP 2.90 (2.6 and BSMAPz from source), including `-n 1` | `ZS` | all four strands |
| `bs-seeker2-se.toml` | BS-Seeker2 2.1.8, single-end only | `XO` | all four strands |
| `merantk.toml` | meRanTK 1.3.0 (meRanGh and meRanGs), single-end or paired-end; SEQ repaired for a meRanGs `-r` run | `YG` + `YR` | all four strands |
| `biscuit.toml` | BISCUIT 1.10.3, single-end or paired-end, any `-b` | `YD` + FLAG 0x10 | conversion strand |
| `bwameth.toml` | bwa-meth 0.2.10, directional only | `YD` + FLAG 0x10 | conversion strand |
| `hisat-3n.toml` | HISAT-3N 2.2.1-3n-0.0.3, any mode | `YZ` + FLAG 0x10 | conversion strand |
| `dnmtools.toml` | abismal 3.3.0 (in dnmtools 1.5.0), before `dnmtools format`, SEQ repaired | `CV` + FLAG 0x10 | conversion strand |
| `astair.toml` | asTair 3.3.3 (TAPS), paired-end or single-end, to mirror its own caller | whole FLAG | conversion strand |

methylpy's output BAM is conformant and carries no conversion tag, so it is
`directional.toml` — including for a `--pbat` run, which reaches the aligner already in the
directional arrangement: paired-end methylpy exchanges the two mate files
(`call_mc_pe.py:225-227`), single-end it reverse-complements every read
(`utilities.py:817`). Read from methylpy 1.4.7 in
`alnbase-validation/demos/pbat-strand/`, not run.

"All four strands" means the rule yields OT, CTOT, OB or CTOB, from which the walk
direction and the read's sequenced direction both follow by biology; "conversion strand"
means the input only distinguishes two of the four, so the file declares the walk direction
and the sequenced direction separately and reports the strand of origin as unknown.

The split those inputs cannot cross is a strand and its own PCR copy, which hit the same
converted contig — measured on BISCUIT's output in
`alnbase-validation/demos/biscuit-strand/`, where `YD:A:f` held OT and CTOT and `YD:A:r` held
OB and CTOB. Such a file has to read FLAG 0x10 for the sequenced direction, which is a claim
about the aligner and not a default; on BISCUIT the bit agreed with how SEQ was really stored
on every record of both a `-b 1` and a `-b 0` run.

Two rules can apply to one BAM. bwa-meth is directional-only, so `bwameth.toml` and
`directional.toml` both describe its output — one from the aligner's own evidence, the other
from the library's design — and on a real bwa-meth BAM they put the same anchors at the same
reference positions, 0 disagreements over 29 and 37 hits in
`alnbase-validation/demos/bwameth-strand/`. `directional.toml` additionally names the strand of
origin. That agreement is a consistency check, not a correctness one: both rules trace back to
the contig bwa-meth chose, so they agree even on a record it placed on the wrong strand.

`directional.toml` is the one file with no tag to read, so a real run has nothing in the BAM to
compare it against; Bismark supplies one. Its `XR`/`XG` names all four strands from the read's
own conversion and the index that won, which is evidence the FLAG rule never touches, so on a
directional Bismark BAM the two calls are independent and the fixture holds a third answer
besides. Over 402 paired-end and 202 single-end records in
`alnbase-validation/demos/directional-strand/`, both rules called every record the strand it was
built from and disagreed about none. That run is also why the file now names Bismark
single-end: Bismark sets 0x10 there from the genome conversion rather than from SEQ's
orientation, but on a directional library the two coincide, because matching the bottom-strand
index is the same thing as having been sequenced from the bottom strand.

**There is no PBAT file, and that is a result rather than an omission.** One was written and
then deleted, because a PBAT library does not reach a BAM as one: both pipelines built for PBAT
put the reads back into the directional arrangement first. Bismark `--pbat` swaps 0x40 and 0x80
on CTOT and CTOB pairs and under `--pbat` every pair is one of those; methylpy `--pbat`
exchanges the mate files (`call_mc_pe.py:225-227`) or reverse-complements every read
(`utilities.py:817`). So **Bismark `--pbat` output is `bismark.toml` in both layouts**, and
`directional.toml` is also correct for the paired-end arm. On real `--pbat` output in
`alnbase-validation/demos/pbat-strand/` the deleted PBAT rule named the opposite conversion
strand on 402 of 402 paired-end records — 17.09% of read bases then differ from the reference,
against 0.00% under `directional.toml`. The library does not determine the rule; the BAM does,
and that demo is the worked example of the difference.

That demo is also where Bismark's single-end 0x10 was finally caught contradicting how SEQ was
stored, on 201 of 201 records, which is why single-end `--pbat` needs `bismark.toml` and no
FLAG rule can serve. The directional round could not see it, because on a directional library
the genome conversion and the read's orientation coincide.

`unconverted.toml` is the one file here that is not about a conversion chemistry at all. A
strand rule is mandatory for every run that reads records, and an assay that converts nothing
still needs one; this is it. It walks every record along the reference forward, so a pattern
written against the reference means the same thing on every record, and it keeps the read's own
orientation in `read_reverse` so `off_5p` and `off_3p` still count from the sequencer's 5' end.
Reaching for `directional.toml` instead makes the walk follow each fragment's orientation and
reverse-complements roughly half the data.

A tag can also be an inference rather than a record of what the aligner did. BISCUIT's `YD`
and bwa-meth's name the converted contig the read was aligned to; HISAT-3N aligns to both
three-letter indexes and then decides by counting the read's C→T and G→A differences, so its
`YZ` is a conclusion drawn from the same alignment the caller is about to read. It is right
whenever there is one difference to count — measured over 456 reads in
`alnbase-validation/demos/hisat-3n-strand/` — and when the counts tie, which a fully methylated
read makes them do, it falls back to the strand a directional library would have had. That is
wrong for every read from a PCR copy, half the sweep, at the right position and MAPQ 60. A
rule cannot see the difference, because the tag looks the same either way.

A tag can also say less than its name suggests. dnmtools' `CV` is the conversion the *read*
shows as it came off the sequencer, not the reference strand that carried it, so the
conversion strand is `CV` combined with 0x10 — which is how abismal, the aligner that writes
the tag, computes it itself. `dnmtools.toml` originally read `CV` alone and was backwards on
every reverse record, 4 of 8 in `alnbase-validation/demos/dnmtools-strand/`; it now declares the
combination and gets 0 of 8 wrong. Outside abismal's PBAT modes `CV` is the mate number and
nothing more, the same emptiness bwa-meth's `YC` has.

A rule need not read a tag at all. asTair writes none and reads none: its caller takes the
strand from a literal table of FLAG integers (`caller.py:235-267`), and `astair.toml` is that
table written out, so an alnbase run mirrors asTair's extraction rather than describing an
aligner's output. Mirroring it reproduces asTair's own `.mods` file — 261 positions in common,
none disagreeing, mod 710 against 711 and unmod 1954 against 1954 in
`alnbase-validation/demos/astair-strand/`. A rule of this shape has to be total in a way a tag
rule does not: any FLAG outside the six is neither strand, so the file declares `unknown` for
them and the run skips and counts them, two records of 410 in that demo.

A four-way tag can still be wrong. BSMAP's `ZS` names all four strands, and at `-n 1` it
names the opposite conversion strand for some copy-strand reads, measured at 10 of 228 in
`alnbase-validation/demos/bsmap-strand/`. A rule reports what the tag says, which is the
honest thing for it to do; what surfaces such a record is alnbase's count of read bases that
differ from the reference, which stops being zero on data that is otherwise clean.

## Reading the tag rather than the FLAG

Three of these aligners set FLAG 0x10 from the reference strand the read matched rather than
from the read's own orientation: Bismark single-end, BSBolt on every `*_G2A` record (read 2
of every directional pair, measured on its output in
`alnbase-validation/demos/bsbolt-strand/`), and BS-Seeker2's `RC` classes (measured in
`alnbase-validation/demos/bs-seeker2-strand/`). All of them still store SEQ in
the reference's forward orientation, so a rule that names the strand of origin places and
orients the record correctly without consulting 0x10 at all — which is why the file for each
of those three is a plain declaration and not a special case. `bsbolt.toml` in particular
removes the need to repair read-2 FLAGs before alnbase reads a BSBolt BAM
(`alnbase-validation/demos/premethyst-bugs/README.md`).

The mate bit can lie in the same way. Under `--non_directional`, Bismark swaps 0x40 and
0x80 on CTOT and CTOB pairs so that browsers do not discard them as discordant, so on those
pairs every record's mate bit names the other mate — measured on its output in
`alnbase-validation/demos/bismark-nondirectional-pe/`, four of eight records. `XR` and `XG`
stay with the record they describe, so `bismark.toml` is unaffected; a rule that inferred
the strand from "read 1 or read 2" would have those records backwards. Under `--pbat` every
pair is CTOT or CTOB, so the swap applies to all of them, which is what makes a PBAT-shaped
rule wrong on that output and `directional.toml` right.

## Not covered

- **BS-Seeker2 paired-end**, whose `XO` is per pair rather than per record.
- **abismal output with SEQ as written.** abismal stores SEQ exactly as it appears in the
  input FASTQ, on reverse records too, which its own caller compensates for by reading SEQ
  backwards without complementing. A strand rule cannot express that, because it declares
  which strand a record belongs to and not how its bases are stored. Reverse-complement SEQ
  and reverse QUAL on every 0x10 record first; alnbase's count of read bases differing from
  the reference is near a third of all bases until you do.
- **dnmtools after `format`**, which merges the two mates of a fragment into one N-gapped
  record, drops the paired bit and writes `CV:A:T` on everything.
- **Copy-strand reads in a bwa-meth BAM.** bwa-meth accepts directional libraries only. A
  CTOT or CTOB read that reaches it is usually dropped, but the fraction that maps — 32 of 228
  in `alnbase-validation/demos/bwameth-strand/` — is placed on the opposite conversion strand at
  the right position and high MAPQ, with nothing in the record to mark it. No rule can recover
  that; the library has to be filtered before alignment.
- **A non-directional pair in a TAPS BAM.** asTair is directional-only by design, so this is
  a limit rather than a defect. Its six-integer table is a claim about fragment orientation,
  which is why it handles the CTOT and CTOB *records* that are read 2 of every ordinary
  fragment; a pair whose read 1 came off the copy strand inverts that orientation and both
  mates are counted on the wrong side, measured at 4 of 8 probe records in
  `alnbase-validation/demos/astair-strand/`. Nothing upstream catches it, since TAPS leaves the
  reference unconverted, and only 7 of 24480 read bases differ from the reference, so
  alnbase's mismatch count does not mark it either. The fix is upstream of alignment.
- **Non-directional input with no evidence at all** — Bismark `--non_directional` under
  `--old_flag`. The honest rule is `unknown` for every record, which means the input cannot be
  used; no file here guesses.
- **A HISAT-3N record whose conversion counts tied.** HISAT-3N does not leave such a record
  unmarked: it writes the `YZ` a directional library would have implied, which is wrong for
  every copy-strand read. `Yf:i:0` marks the fully methylated case and nothing marks the
  equal-counts one, so no rule over `YZ` can decline them. Pass `--directional-mapping` when
  the library is directional, so the assumption is one you made on purpose.
