# alnbase command-line reference

Every option, what it does, and when you would reach for it. For the concepts
behind them — columns, patterns, queries, tags — read
[`alnbase.md`](alnbase.md) first; this is the reference you come back to.

`alnbase <command> --help` prints the same information in short form, and
`-h` gives a one-line summary per option. Where the two disagree, the binary
is right and this file is stale.

- [Common shapes](#common-shapes)
- [`index`](#index)
- [`info`](#info)
- [`query`](#query)
  - [Dry runs: understanding a query](#dry-runs-understanding-a-query)
  - [The walk](#the-walk)
  - [Threads](#threads)
  - [Data you did not expect](#data-you-did-not-expect)
  - [Parquet output](#parquet-output)
- [`extract`](#extract)
- [`dump-query`](#dump-query)
- [`overlap`](#overlap)
- [Exit status and diagnostics](#exit-status-and-diagnostics)

---

## Common shapes

**`-` means standard input or output** everywhere a BAM is read or written.
A parquet output is not a stream: its `OUT` is a path prefix (see
[Parquet output](#parquet-output)).

**Repeatable options accumulate.** `--query-file` and `--field` can be given
more than once and combine in order.

**Nothing is filtered for you.** alnbase does not skip duplicate, secondary,
supplementary or QC-fail reads, and has no MAPQ or base-quality threshold.
`samtools view -F 0x900 -q 20` upstream does that, and the parquet output
carries `qual` and `mapq` so downstream filtering is a `WHERE` clause. Two
filters that disagree are worse than one.

---

## `index`

```
alnbase index <INPUT_FASTA> <OUTPUT_INDEX>
```

Convert a reference FASTA into alnbase's own format, once per reference. The
result is memory-mapped and `Seq`-encoded, so every worker reads sequence
without locking or decoding.

Accepts plain or bgzip-compressed FASTA. Soft-masking is not preserved: the
index records which base is at a position, not whether a masker lowercased it.
A FASTA with two contigs of the same name is refused, since a BAM could not say
which one it means. So is any sequence character other than an IUPAC base
(`A C G T U R Y S W K M B D H V N`, either case): `*`, `-`, `.`, `X` or a digit
means a protein, alignment or damaged FASTA, and the error names the line, the
contig and the position.

The index is written to `OUTPUT_INDEX.partial` and renamed when complete, so a
failed run leaves no file behind and does not touch an existing index.

The index also stores each contig's MD5, the digest SAM puts in `@SQ M5`. This
is index format version 3: an index built by alnbase 0.1.9 or earlier is
refused with a message saying so. Rebuild it with `alnbase index`.

**Use the FASTA the reads were aligned against.** An index built from a
different file produces calls that are self-consistent and wrong, which is
much harder to notice than a crash. `query` catches it only as far as the BAM
header allows: a contig whose length differs stops the run, and so does one
whose `@SQ M5` differs from the index's MD5. Many aligners write no `M5`
(`samtools dict` does), and without it a contig with the same name and length
but different sequence cannot be detected.

---

## `info`

```
alnbase info [--all] [--seq CHR:FROM-TO] <INDEX>
```

| | |
|---|---|
| *(no options)* | format version, contig count, total bases, file size, then the contigs with their lengths and MD5s (the `@SQ M5` values `samtools dict` writes) |
| `--all` | list every contig instead of the first and last five |
| `--seq CHR:FROM-TO` | print the sequence over a region |

`--seq` uses `samtools faidx` coordinates — 1-based, inclusive, 60-column
lines — so the two outputs diff directly. The output is uppercase, so compare
case-insensitively:

```bash
diff <(alnbase info --seq chr19:3034740-3034780 mm10.aref | tr a-z A-Z) \
     <(samtools faidx mm10.fa chr19:3034740-3034780   | tr a-z A-Z)
```

This is the tool for "are these two references the same file", which is worth
checking before believing any disagreement between alnbase and another caller.

---

## `query`

```
alnbase query [OPTIONS] [BAM] [REFR] [OUT]
```

Walk every record and write it back out, in input order, with the tags the
query files declare. With `--parquet`, write hit rows instead.

The three positional arguments are the input BAM, the reference index, and the
output. They are not required for a dry run (`--explain`, `--trace`,
`--list-codes`).

### Queries

| | |
|---|---|
| `--query-file FILE` | TOML file of queries and tags. Repeatable; files combine, and a tag may be declared in only one of them |

The long `--help` for this option prints the whole file format and a worked
example.

### Dry runs: understanding a query

None of these opens a BAM, except `--trace-records`.

| | |
|---|---|
| `--explain` | each query rendered column-aligned, as a tree, and in English, plus the normal form its predicate compiled to. With more than one query, ends with a one-line summary per query. A query whose read row places a pad (`_`, or a group such as `{._}`) gets a note that for read 2 the pattern's left and right run against the sequenced 3' and 5' ends |
| `--trace READ@REFR` | run the queries over a synthetic column pair, e.g. `--trace 'TTTCGTTT@TTTCGTTT'`, showing how far each pattern got and which query fired |
| `--trace-records N` | the same report over the first N records of a real BAM, so the columns include real CIGARs, indels, clips and read ends. Each record ends with its reference concordance, and the traced records are checked against `--max-discordance` like a scan (see below) |
| `--trace-grid <bool>` | show the per-column progress grid in a trace (default true; turn it off for very wide records) |
| `--list-codes` | every code a pattern row may hold, and which letters are still free for aliases |

`--explain` before a long run; `--trace-records` when a specific read does
something you did not expect. The second is the one to reach for when
comparing against another caller: find the position, trace the read, and see
which pattern broke where.

### The walk

These change which columns exist, so they change what a pattern can match.

| | default | |
|---|---|---|
| `--insertions skip\|emit` | `skip` | `skip` keeps the two reference bases flanking an insertion adjacent, so a pattern can span it. `emit` is needed to write patterns *about* insertions, at the cost of breaking patterns that span one |
| `--end-context N` | widest query's span | flank columns shown past each end of the aligned read: reference side the real base; read side a clip (`:`) for each soft-clipped base nearest the aligned part, then pads (`_`). The default is far enough for every query, the widest included, to fit a whole window into the pads. **Edge case:** with the default, adding a wider query adds pads, so a query that can match windows made only of pads (read row `~` or `_` throughout) gains hits. Set this explicitly to fix the pads whatever else is queried, or exclude all-pad windows in the query with `and not` ([details](reference/02-query-language.md#86-pads-flank-contig-ends-and-record-boundaries)) |
| `--splice-context N` | widest query's span | reference bases shown at each end of an intron. The default makes a junction edge behave like the end of a read. Set it explicitly when matching a splice motif or the `,@,` marker (0 puts it between the exons), so the layout does not move with the widest query in the run. **Edge case:** as with `--end-context`, a query that can match windows made only of intron columns gains hits when a wider query raises the default ([details](reference/02-query-language.md#86-pads-flank-contig-ends-and-record-boundaries)) |

Soft-clipped bases are never aligned columns: they were sequenced but not
aligned, and appear only as clip columns in the flank, which a `bases` tag never
marks. (The `--soft-clips` option was removed in 0.1.2.)

`off_5p` and `off_3p` count skipped insertions, soft-clipped bases and (since
0.1.4) hard-clipped bases, and are measured from the ends of the read as
sequenced, for read 2 as for read 1: `off_5p = 0` is the first base the
sequencer read, and a supplementary `60H20M` reports `off_5p = 60` for its first
base, as its primary does. Patterns, by contrast, run
along the conversion strand, so for read 2 a pad before a base (`_N`) marks its
sequenced 3' end.

A query fires on any window of these columns it matches, including windows made
only of pads, intron columns or the `,@,` marker: `_N` fires on a read's end base,
and a one-column `_` on every pad. Set `--end-context` and `--splice-context`
explicitly and every query's hits are independent of the other queries in the
run; with the defaults, a query that can match all-pad or all-intron windows
depends on the widest one ([pads and the widest query](reference/02-query-language.md#86-pads-flank-contig-ends-and-record-boundaries)).

### Threads

| | default | |
|---|---|---|
| `-@, --threads N` | 1 | tagging threads |
| `--reader-threads N` | `0`, meaning min(4, `--threads`) | BGZF decompression. One reader feeds every worker, so this is usually the ceiling: a single BGZF thread tops out around 150–250 MB/s |
| `--writer-threads N` | `0`, meaning min(4, `--threads`) | BGZF compression for the output BAM. One writer takes every record in input order, so without threads of its own it serialises the taggers |

Raise the reader and writer threads until throughput stops improving before
raising `--threads`; they compete for the same cores.

### Data you did not expect

| | |
|---|---|
| `--permissive` | warn and count instead of stopping |
| `--max-discordance FRACTION` | stop when more than this fraction of compared read bases disagree with the reference. Default 0.25; 0 to 1; `1` turns the check off |
| `--require-m5` | stop at the first record on a contig whose `@SQ` line has no `M5`, even with `--permissive` |
| `--overwrite-tags` | replace a tag an input record already carries |

By default a record on a contig the reference lacks stops the run at the first
one, and so does a record whose aux block cannot be read. `--permissive`
turns both into counted warnings: off-reference records are written untagged,
and a damaged record is still tagged, with the new tags written at the
*front* of its aux block where readers reach them before the damage. Neither
changes what is written for healthy records.

`--max-discordance` catches a wrong reference, a BAM aligned to a different
assembly, or shifted coordinates. Every walked record counts read bases that
disagree with the reference, over columns where both sides show one unambiguous
base; a reference C read as T is not counted, so conversion does not raise the
rate. Correct alignments disagree at a few percent at most. The count is shared by
every thread: once the run has compared 100,000 bases, a rate above the threshold
stops it, removing its partial output as for any failure, and a run that ends
sooner is judged at the end on everything it compared (from 1,000 bases). Before
0.1.6 each worker was judged alone, so an input spread thinly over many threads
could finish with a high rate. The run summary always reports the rate, as
`N of M compared read bases differ from the reference (x.xx%; a reference C read
as T is not counted)`. The check applies to the tagged BAM, to `--parquet` and to
`--trace-records`, not to `extract`, `--explain` or `--trace`. A trace prints each
record's rate, stops with the same error when the traced records compared at least
1,000 bases above the threshold, and prints a warning when they are above it on
fewer.

The reference is always matched to the BAM header by contig name and length, and
by MD5 wherever the header's `@SQ` line carries `M5` (a disagreement stops the run).
Most aligners write no `M5`, and a masked analysis set or a patched assembly has the
same names and lengths as the original, so without one nothing tells them apart.
`--require-m5` makes the MD5 check mandatory: a record on a contig without `M5`
stops the run. Contigs no record is on are not checked. `samtools dict REF.fa`
prints the `M5` values for a FASTA, `alnbase info` prints the index's, and
`samtools reheader` puts an edited header on a BAM.

`--overwrite-tags` is off by default because a record that already carries,
say, `XM` was probably tagged by another tool, and quietly replacing its calls
is a mistake that surfaces much later. Giving your tag a different name is
often better: a tag name containing a lowercase letter (`xm`) is reserved for
local use and cannot collide.

### Parquet output

`--parquet` writes hit rows instead of a tagged BAM, sharded into
`--threads × --shards-per-worker` files named `{stem}_{worker}_{shard}.{ext}`
beside `OUT`: `calls.parquet` becomes `calls_0_0.parquet`, `calls_0_1.parquet`,
and so on, and an `OUT` with no extension just gains the suffix. The tag tables
in the query file are not used.
Everything below applies only with `--parquet`, and is refused without it.

**Columns**

| | |
|---|---|
| `-f, --field FIELD` | record columns to write, replacing the default (`core`) |
| `-F, --add-field FIELD` | add to the default instead of replacing it |
| `--only-hits` | drop the rows for records that matched nothing |

Both are comma-separated and repeatable, and order is preserved, so they set
column order too. `-F NM:i` is `core` plus that tag. Aux tags are given as
`TAG:type` using SAM's type codes: `i c C s S I` → Int64, `f d` → Float64,
`Z A H B` → Utf8. The long `--help` for `--field` lists the named columns and
groups. Among them, `conv_strand` gives the strand of origin (`OT`, `OB`,
`CTOT`, `CTOB`), which splits `strand` (`+` for OT and
CTOT, `-` for OB and CTOB) by mate, and is null when the run's strand rule
names only two strands; `read_reverse` gives the direction the read was
sequenced in, which is what `off_5p` and `off_3p` count from and is not the
same question as the FLAG's `is_reverse`; `hard_clip_5p` and `hard_clip_3p` give the number of
hard-clipped bases at the read's sequenced 5' and 3' ends (0 when there are
none), which `off_5p` and `off_3p` already include; `soft_clip_5p` and
`soft_clip_3p` do the same for soft clips, so `off_5p - soft_clip_5p -
hard_clip_5p` counts from the first aligned base; and `qual_phred` is null
when QUAL is `*` (as is the hit column `qual`).

Every record read gets a row, so `record_id` is its position in the input and
the table lines up one-to-one with the BAM. A record that matched nothing, or
could not be walked, gets a row with the query columns null; `--only-hits`
drops those.

**Sharding**

| | default | |
|---|---|---|
| `--partition-by FIELD` | `qname` | which record attributes decide a record's file |
| `--shards-per-worker N` | 1 | files per worker; total files is `--threads × this` |

Within a run, records with the same key always share a file. Across runs, the
same key is only guaranteed the same file when `--threads` and
`--shards-per-worker` are both unchanged, since the file is picked by the hash
modulo the total file count.

Partitioning by anything other than `qname` can put a fragment's two ends in
different files, which matters if a later step resolves mate overlap file by
file. alnbase warns when the key allows it.

**Format and size**

| | default | |
|---|---|---|
| `--output-format parquet\|ipc` | `parquet` | IPC is append-only, so it can be written to a FIFO and read batch by batch as the scan runs |
| `--batch-rows N` | 500000 | rows buffered in memory before they are written, in total across all output files: each file holds an even share. Memory scales with this, whatever the file count. At least 1 |
| `--row-group-rows N` | 50000 | rows per parquet row group. At least 1 |
| `--output-profile fast\|small` | `fast` | `fast` is LZ4 with chunk statistics; `small` is zstd with page statistics |
| `--compression none\|lz4\|snappy\|zstd` | from profile | overrides the profile |
| `--compression-level N` | 3 | zstd only |
| `--statistics none\|chunk\|page` | from profile | what a reader can prune without decompressing |
| `--column-encoding COL=ENC` | per column | overrides the encoding for one column |

`--compression none` is usually the *slowest* option, not the fastest: the
file is roughly 17× larger and writing those bytes costs more than LZ4 costs
to avoid them. It is worth it only when the output goes to a pipe or a tmpfs.

**Manifest and metadata**

Every coordinate in a hit table (`refr_pos`, `pos`, `end_pos`, `mate_pos`, ...)
is 0-based. Each file says so: its parquet footer carries the key-value entries
`format_version` (currently `7`), `coordinate_base` (`0`) and `alnbase_manifest`,
a JSON description of the run (an IPC stream carries the same keys as schema
metadata). DuckDB reads them with `parquet_kv_metadata('calls_*.parquet')`.

The manifest records the alnbase version, a `run_id` shared by the run's files,
the command line, the input and reference paths, the library, the output settings
(format, partition key, worker and shard counts, fields, capture layout, and the
file names in slot order, so `files[shard]` is the file holding that `shard`
value), the text of every query file, the walk with its defaults resolved
(`end_context`, `splice_context`, whether each was given, the widest span,
`--max-discordance`, `--require-m5`), and every `@SQ` contig in header order with
its length, MD5 and whether the reference has it.

When the run succeeds, the same JSON plus a `result` object (record counts, the
concordance totals, and the rows written to each file) is written last, to
`{stem}.manifest.json` beside `OUT` (`calls.parquet` → `calls.manifest.json`).
It is the completeness marker: a run removes any manifest already at that path
before it writes anything, and a failed run leaves none. Read the file list from
the manifest rather than globbing: files an earlier run with more shards left
under the same names are not removed, and a run that finds any says so in a
warning.

---

## `extract`

```
alnbase extract [OPTIONS] <BAM> <OUT>
```

Turn a `bases` tag back into hit rows, with no reference and no walking. The
definition comes from the query file alnbase stored in the BAM's header.

| | |
|---|---|
| `--tag TAG` | which tag to extract; needed only when more than one `bases` tag is defined |
| `--query-file FILE` | read the tag with these files instead of the stored ones. Repeatable |
| `--permissive` | as for `query` |
| `-@, --threads N`, `--reader-threads N` | as for `query` |
| *Parquet output* | the same options as `query --parquet`, and they apply without needing that flag; the output gets the same metadata and manifest, with the tag and its definition in place of the walk |

When the BAM has been through alnbase more than once, each tag resolves to the
most recent run that defines it — tag by tag, not run by run, so a tag only
the older run declared still resolves to that run's definition.

`--query-file` is also how you read a tag written by another tool: a file
whose `bases` table maps that tool's characters to queries is its definition.

Rows carry `name`, offsets from both ends of the read as sequenced (hard clips
included), the reference position, the base quality (null when QUAL is `*`) and
the read base. `refr_base` is filled in when the query
determines it — `T~~@CHH` can only fire over a C — and null otherwise.

---

## `dump-query`

```
alnbase dump-query [--pg ID] [--file N] <BAM>
```

Print the query file stored in a BAM's header, with its `@PG` line as a
comment on top. The output below the comments is byte-identical to the file
the run read, and is itself a valid query file.

| | |
|---|---|
| `--pg ID` | which run, by `@PG` ID. Defaults to the most recent |
| `--file N` | which of that run's files, counting from 0. Needed only when a run read more than one |

---

## `overlap`

```
alnbase overlap [OPTIONS] <INPUT> <OUTPUT>
```

The two reads of a fragment overlap, and the overlap is one molecule observed
twice. This works out each fragment's length, keeps one copy of the overlap,
and merges the qualities. Run it before `query`.

**Input must be grouped by name** — `samtools sort -n` or the cheaper
`samtools collate -O`. Anything else is rejected rather than warned about: an
ungrouped input would resolve the few templates that happen to be adjacent and
miss the rest, which looks like success.

Unlike the incumbents, it keys on the molecule rather than on reference
intervals, because a ligation product is not colinear with the reference.

**The output stays a valid, mate-consistent SAM.** A supplementary alignment
with nothing left after clipping is not written. A primary alignment with
nothing left hands its place to the read's 5'-most surviving supplementary,
whose alignment moves onto the primary record (hard clips become soft clips, so
the primary still carries the whole read). If nothing of the read survives, it
is written unmapped: flag 0x4, 0x2 and 0x10 cleared, MAPQ 0, CIGAR `*`, SEQ and
QUAL back in sequencing orientation, placed at its mate's RNAME/POS. Each read
keeps exactly one primary line. In every template where an alignment was
clipped, dropped or promoted, the mate fields (RNEXT, PNEXT, TLEN, flags 0x2,
0x8, 0x20, `MC`, and `MQ` where present) are rebuilt the way
`samtools fixmate` builds them (TLEN is 5' end to 5' end), so running
`fixmate` afterwards changes none of them. The `SA` tags of a read whose
alignments changed are rebuilt to list the surviving alignments. Templates
that were not changed are written with their mate fields as they were. A
missing quality string (QUAL `*`) is left as `*` and takes no part in quality
arithmetic.

**Basics**

| | default | |
|---|---|---|
| `--refr FILE` | — | needed only by `--stale-tags recompute` |
| `-@, --threads N` | 1 | compression threads for reading and writing |
| `--compression-level N` | 3 | BGZF level of the output, 0 (uncompressed) to 9 |

**Establishing the fragment length**

| | default | |
|---|---|---|
| `--length-tolerance N` | 0 | permitted disagreement between one read's records |
| `--min-overlap N` | 4 | shorter overlaps are left alone |
| `--min-support N` | 3 | reference-matched pairs needed to believe a length |
| `--min-span-frac F` | 0.7 | fraction of the implied overlap the agreeing pairs must span |
| `--max-anchor-slack N` | 5 | how far the supporting pairs may stop short of either read's last aligned base. This is the test that does most of the work |
| `--max-past-end N` | 0 | bases a read may still have aligned past the end of the implied molecule |
| `--max-mismatch-frac F` | 0.15 | above this, nothing is merged |
| `--ambiguity-ratio F` | 0.8 | a rival length this well supported blocks the merge |
| `--indel-window N` | 10 | a rival this close in length is reported as an indel disagreement rather than an ambiguous length |
| `--on-unresolved pass\|drop` | `pass` | `pass` writes the template unmodified and tagged, so coverage is kept and the reads stay filterable |

**Consensus, where the two copies agree or disagree**

| | default | |
|---|---|---|
| `--mismatch-qual none\|zero-both\|zero-loser\|subtract` | `subtract` | what to do to two disagreeing bases' qualities |
| `--mismatch-base none\|set-n` | `none` | whether to replace disagreeing calls with N |
| `--match-qual none\|max\|sum-capped` | `sum-capped` | what to do when they agree |
| `--qual-cap N` | 40 | ceiling on a combined quality |

The cap is not a formality: two mates are the same molecule, so a
pre-amplification polymerase error appears in both and the observations are
not independent. Summing to the SAM maximum would assert a confidence the
evidence does not support.

**Which copy survives**

| | default | |
|---|---|---|
| `--keep r1\|r2\|score` | `score` | `score` picks whichever end scores higher on MAPQ and base quality |
| `--score-qual-weight F` | 1 | weight on mean base quality relative to mean MAPQ |
| `--clip-mode soft\|hard` | `soft` | soft clipping is enough to keep bases out of a call and, unlike hard, can be undone |
| `--keep-terminal-indels` | off | do not clip terminal indels off the kept 3' end |
| `--terminal-indel-window N` | 3 | matched bases within which an indel counts as terminal |

**Tags and reporting**

| | default | |
|---|---|---|
| `--stale-tags recompute\|strip` | `strip` | clipping invalidates NM, MD and AS; `recompute` recalculates NM and MD (needs `--refr`) and drops AS, `strip` drops all three on the grounds that absent beats wrong. XA is always dropped from a clipped record |
| `--length-tag TAG` | `XL` | where the fragment length is recorded |
| `--tag-prefix C` | `o` | first character of the two-character tags this tool writes. Lower case deliberately: `XM` and `XG` are Bismark's call string and genome-conversion tag, and writing there would corrupt the thing downstream is about to read |
| `--no-tag` | off | do not write the length tag or the prefixed tags (`oO oN oS oX oK oV ov`). Tags that describe an alignment this step changed are still kept correct: NM/MD/AS/XA per `--stale-tags`, SA, MC and MQ |
| `--stats FILE` | | write the run's statistics |
| `--report-unresolved` | off | list the templates whose length could not be established |

---

## Exit status and diagnostics

An error exits 1, except a command-line syntax error (an unknown or removed
option, a bad value such as `--batch-rows 0`), which exits 2. There are no
finer codes, so do not write a pipeline that reads more meaning into the
number.

A failed `query` run writing a BAM removes its own partial output (unless it
was writing to standard output) rather than leaving a file that looks complete — htslib writes an end-of-file marker however a BAM is closed,
so a truncated one is not otherwise detectable. A failed `--parquet` or
`extract` run removes every file it was writing, including the ones that
were complete, and leaves no `{stem}.manifest.json`. Files from an earlier run
with more shards are not touched; a successful run warns about them.

When the system refuses to start a thread, usually because of the per-user
process/thread limit (`ulimit -u`), the run stops with an error that says so
rather than crashing; lower `--threads` or raise the limit.

Counts that are printed only when non-zero, at the end of a run:

| | |
|---|---|
| records skipped | unmapped, on a contig the reference lacks, or mapped with no SEQ (`*`) |
| unplaced hits | a hit anchored past either end of the read, which a per-base tag has nowhere to record |
| damaged aux blocks | records whose aux could not be read in full |
| non-UTF-8 aux text | values read as null rather than lossily converted |

A zero count prints nothing, so a quiet run means none of these happened.
A `query` run that compared any read bases also ends with the reference
concordance line, `N of M compared read bases differ from the reference (x.xx%;
a reference C read as T is not counted)`, whatever the rate (see
`--max-discordance`).
