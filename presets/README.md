# Presets

One directory per tool. Each holds that tool's extraction, written out in full:

| file | what it is | what runs it |
| --- | --- | --- |
| `queries.toml` | what counts as a methylation call | `alnbase query -q` |
| `run.sql` | how calls are grouped, and the file written | `alnbase-agg run` |
| `README.md` | each choice the tool makes, and where in its source it is made |

Nothing here is compiled into either program. A preset is a pair of ordinary input
files, so reproducing a tool you do not see here means writing two files, not
patching and re-releasing anything — and forking one of these is `cp -r`. The same
goes for the SQL a run reads: `alnbase-agg lib` prints where it lives, and a copy
next to `run.sql` is the one that gets read.

## Why they exist

A preset is a tool's implicit extraction query made explicit. Every caller in this
ecosystem makes a series of choices — which reference character decides a context,
what happens at a contig's last base, whether a deleted base is skipped or used,
which reads are silently dropped — and almost none of them are written down in one
place. Reading a preset is reading that list.

That makes them useful twice over:

- **As a comparison.** `run.sql` *is* the query: there is nothing it compiles to
  and nothing filled in elsewhere. Two presets side by side are a diff between two
  tools' definitions of the same measurement, in a language that does not hide
  anything. `alnbase-agg show run.sql` prints it with the library files it reads
  substituted in, if you want the whole thing in one piece.
- **As a starting point.** "Bismark's calls, but with my own read filter" is an
  edit to one file.

A preset is a specification, not a test harness. It is not run against the tool it
names and its outputs are not diffed against that tool's: reproducing another
caller's numbers is not a goal here, and where a tool has a correctness bug,
agreeing with it would mean reproducing the bug. Bugs are demonstrated separately,
by minimal reproductions under `alnbase-validation/demos/`, which isolate a defect better
than a dataset-wide diff that mixes every other difference in with it.

## Running one

Each recipe's header has the two commands, because both halves need flags that
match the tool: Bismark sees inserted bases and MethylDackel does not, so one runs
with `--insertions emit` and the other with `--insertions skip`. A preset run with
the wrong walk options is not that tool.

```sh
cd presets/methyldackel
alnbase query -q queries.toml --insertions skip \
    -F qname,is_first_in_template,ref_name,strand,mapq,flags,is_paired,\
is_proper_pair,is_mate_unmapped,qual \
    reads.bam calls
alnbase-agg show run.sql                  # the whole query, before running it
alnbase-agg run run.sql -s hits='calls_*_*.parquet' -s sample=sample
```

## What a preset does not cover

**Mate overlap.** Tools differ in how they handle the region where two mates of a
fragment overlap, and those rules are not reproduced here: alnbase has its own
`overlap` command, and reimplementing each tool's version would be reimplementing
the thing under comparison. A recipe's `[unit] by = "fragment"` collapses a
position seen by both mates into one call and drops it if they disagree, which is
its own rule and not any tool's. Compare on fragments without mate overlap, and
read each preset's README for what its tool does there.

**Configurable filters.** Where a tool's defaults are options rather than
hardcoded behaviour — MAPQ, base quality, flag masks — the preset carries the
defaults in `[filter.hit]` and says so, so that changing them on one side is one
edit on the other. Filters a tool hardcodes belong in `queries.toml`, because they
are part of what it means by a call.

**Statistical calls.** Nothing here writes a column a model produces rather than a
count. ALLC's significance column is the example. Where an established converter
turns a count table into such a format, the route is to write the count table —
`formats/bismark_cov.sql`, or `formats/sites_parquet.sql` when nothing standard
fits — and run the converter yourself. A preset does not run it, because a preset
is a specification rather than a pipeline.
