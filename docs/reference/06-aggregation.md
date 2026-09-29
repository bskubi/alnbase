# Aggregation

`alnbase-agg` turns hit rows into the products other tools expect: per-position
tables, windows, per-cell summaries, per-read summaries. It is a Python package in
`python/`, not part of the Rust binary, and it aggregates nothing itself — it
runs a SQL script that you write and lets DuckDB do the work.

```
alnbase-agg explain python/examples/cg_ch.toml --hits 'calls_*_*.parquet'
alnbase-agg run     python/examples/cg_ch.toml --hits 'calls_*_*.parquet' -t 8 -m 16GB
```

## Why it is not in the binary

Two reasons, and neither is tidiness.

DuckDB already does the group-by, the spilling, the parallelism and the predicate
and projection pushdown into parquet. Writing that again in Rust to save a
read-and-write cycle would be reimplementing a query engine for a constant
factor.

And the long tail of output formats is where the release pressure is. A new
format must never mean recompiling and re-releasing alnbase, the same reason
strand rules ship as query files rather than compiled-in presets. So writers live
behind a registry (`alnbase_agg.writers.register`), and the package that owns a
format registers it: amethyst-facet writes Amethyst H5 and premethyst cellInfo,
because it also converts, windows and deletes those files.

A plugin declares an entry point in the `alnbase_agg.writers` group, and importing
the module it names is the registration:

```toml
[tool.poetry.plugins."alnbase_agg.writers"]
amethyst = "amethyst_facet.alnbase.writers"
```

So installing amethyst-facet is all it takes for a recipe naming `amethyst-h5` to
work; `alnbase-agg formats` lists what arrived, and names the formats a recipe may
legally ask for that nothing installed can write. A plugin that fails to import
warns rather than stopping a run that only needed a TSV.

## The shape of a recipe

| table | what it declares |
| --- | --- |
| `[context.NAME]` | `c` and `t`: which hit rows mean methylated, which unmethylated |
| `[filter.hit]` | `where`, a predicate over hit-row columns |
| `[unit]` | `by = "fragment"` or `"read"`: what a read filter judges |
| `[cell]` | `barcode`, a format string over parquet columns |
| `[read_filter.CTX]` | `rules`, a ladder over the read's observed site count |
| `[output.NAME]` | `by`, the group key; `contexts`; `where`; `units`; `positions`; `derived` |
| `[output.NAME.FORMAT]` | one writer over that aggregation |

Contexts are declared, never chosen from a menu. `c` and `t` take a list of
alnbase query names or a SQL predicate, so NOMe's `GCH`, a single `CAC`, or a
T>C conversion measure over the same hit rows are all ordinary contexts.

## The group key is the only thing that changes

`by` is a list of dimensions, and every product level is a different one:

| `by` | one row per | serves |
| --- | --- | --- |
| `["site"]` | contig, position | bedGraph, `.cov`, bedMethyl |
| `["cell", "site"]` | cell × position | Amethyst base-resolution datasets |
| `["cell", "window"]` | cell × fixed bin | window sums, with `window = 500` |
| `["cell"]` | cell | cellInfo, per-cell QC |
| `["read"]` | fragment or record | per-read conversion, read-level QC |
| `["site", "strand"]` | cytosine | stranded bedMethyl, unmerged CG calls |
| `["site", "strand", "context"]` | cytosine × context | CGmap |

Windows are half-open, so a position belongs to exactly one window and no
boundary is double counted. Sliding windows are refused and pointed at
amethyst-facet's `facet agg`, which computes them from a written
base-resolution file.

Counts arrive as one column pair per context — `CG_c`, `CG_t` — and a predicate
or derived expression refers to them as `CG.c`, `CG.t`, `CG.sites` (their sum) and
`CG.level` (`c / sites`, NULL when there is nothing to divide by). `m` and `u`
are accepted spellings of `c` and `t`.

### `strand` and `context` are dimensions, not options

`strand` puts alnbase's conversion strand in the key — the strand the cytosine
itself is on, not FLAG 0x10. It costs nothing: one reference base has one strand,
so grouping by it splits no group that was not already distinct, and a per-site
file written from a stranded product still has one line per cytosine. What it buys
is a strand column that is a measurement rather than a label, which is what a
CGmap needs and what lets a bedMethyl show the two strands apart. It needs the
hit table's `strand` column, which is one of alnbase's defaults.

`context` turns the counts from wide to long. Without it a product carries one
column pair per context; with it the context is a key column named `ctx` and the
counts are plain `c` and `t`, one row per context:

```toml
[output.cytosines]
by = ["site", "strand", "context"]   # chr, pos, strand, ctx, c, t
```

That is the only shape in which a per-trinucleotide file over a genome is
tractable. The sixteen `C**` contexts wide would be 32 columns across every
position — some 1.2 × 10⁹ rows for a human genome, nearly all of them zero — and
long it is one `ctx` column with a row only where something was seen. In a long
product a predicate names the counts without a context (`where = "sites > 1"`,
`derived = { class = "substr(ctx, 1, 2)" }`); `CG.sites` is refused, because every
row is already one context.

Where the trinucleotide comes from is alnbase's side of the seam: patterns are
walk-oriented, so a query `refr = "CGA"` matches that trinucleotide read 5'→3' on
the cytosine's own strand, and its name lands in the hit table. A query file
enumerating the sixteen is an ordinary query file — nothing here or in the binary
knows what a trinucleotide is.

## Two things a product can also count

`units = true` on a cell-keyed product adds `units`, `pairs`, `singles`, `passed`,
`rejected` and `unclassified`: how many reads the cell had, and how the read
filter judged them. The first three are counted before the filter chooses, so they
are totals a QC table can be trusted to add up, and the last three add up to
`units` because a ladder that cannot classify a read says so rather than guessing.

That is also how a read threshold is expressed. Other tools have a `min_reads`
option somewhere near whichever writer needed it; here it is

```toml
[output.cg_sites]
by = ["cell", "site"]
units = true
where = "passed >= 1000"
```

— a predicate, like every other filter, which means it composes with the counts
(`where = "passed >= 1000 AND CG.sites >= 1"`) and shows up in the generated SQL
where you can see it. A cell whose every read was rejected still gets a row in a
cell-keyed product, with zero counts: that cell is a fact about the run, and
dropping it would mean re-aggregating to ask about it.

`positions = true` adds `{CONTEXT}_positions`, the distinct positions the group
covered. At one position that is `sites` by definition, so it is refused on a
product keyed by `site`; at any grain that spans several positions it is a
different question, because a cell covering one CG twenty times has twenty calls
and one position covered. premethyst's cellInfo means the second by "coverage",
and Amethyst's window datasets carry it alongside the call counts — which is why
it is a column a product asks for rather than something a writer recomputes from
sums it cannot recover a `DISTINCT` from.

## The formats it writes

A writer is `[output.NAME.FORMAT]` with a `path`, and one aggregation can feed
several of them. These need nothing installed, because a product is already a
table and DuckDB writes it:

| format | columns | coordinates | grain |
| --- | --- | --- | --- |
| `bismark-cov` | chr, start, end, percent, **c**, **t** | 1-based, start == end | `site` |
| `bedgraph` | chr, start, end, value | 0-based half-open | `site` or `window` |
| `methyldackel-bedgraph` | chr, start, end, percent, **c**, **t** | 0-based half-open | `site` |
| `bedmethyl` | ENCODE's eleven | 0-based half-open | `site` |
| `cgmap` | chr, C/G, pos, class, dinucleotide, level, mC, total | 1-based | `["site", "strand", "context"]` |
| `csv`, `parquet` | the product as it stands | 1-based | any |

Bismark's coverage file puts the **methylated** count before the unmethylated
one, and its percentages are Perl's fifteen significant digits
(`66.6666666666667`); both are reproduced exactly, so a diff against Bismark's
own output is a diff about methylation rather than about column order or
floating point. `decimals` overrides the formatting.

`methyldackel-bedgraph` holds the same six fields and is a separate format rather
than an option, because the two differ in exactly the two ways that are invisible
when wrong: MethylDackel's coordinates are 0-based half-open where Bismark's are
1-based with `start == end`, and its percentage is a C cast that **truncates**, so
two thirds methylated is `66` where Bismark writes `66.6666666666667`. A
coordinate base chosen by a flag is how a file comes out wrong at every position
and right-looking at any one of them.

The bed-like three report one number per line, so a product counting more than
one context has to say which with `context = "CG"` — or declare one
`[output.*.FORMAT]` per context. Either count shape feeds them: a writer does not
know whether it was handed `CG_c` or a `ctx` column. A position where that context was not observed
is left out rather than written as zero: an uncovered cytosine is not an
unmethylated cytosine, and every one of these formats is read that way.

`bedgraph` also takes `value`, naming any other column of the product — a
coverage count, or something a `derived` expression computed — because a bedGraph
is a browser's idea of a number per interval and nothing about it is specific to
methylation. It is the one format here a window product can write.

`bedmethyl` takes `strand`, defaulting to `"."` — the honest value for a product
whose contexts do not separate the two cytosines of a CG, and what strand-merged
bedMethyl carries. A product grouped by `strand` writes each cytosine's own strand
instead, and naming the option as well is refused rather than silently preferred.

`cgmap` is the one format here that is not one context per file: its fourth and
fifth columns name the context on the line, which is why it wants the long shape.
Both are read off the cytosine's trinucleotide — the class is CG, CHG or CHH and
the dinucleotide is the cytosine plus its 3′ neighbour — so they come from the
context's name, and a context named `CHG` is refused rather than guessed at,
because `CHG` does not say which base H was. `class` and `dinucleotide` name
columns to use instead. Column 2 is the base on the Watson strand, so a minus-strand
cytosine is written `G`, at its own coordinate. BS-Seeker2 writes the level as
Python's `str(round(level, 2))`, which is neither `%.2f` nor a bare float: it
rounds ties to the even digit (5/8 is `0.62`, not `0.63`) and keeps a trailing
`.0` on a whole number. Both are reproduced, so a diff against BS-Seeker2's own
output is a diff about methylation.

`amethyst-h5` and `premethyst-cellinfo` come from amethyst-facet, which owns those
formats; `alnbase-agg formats` lists what is actually installed.

**ALLC** is not written here, and will not be. Its seventh column is a per-site
significance call — methylpy computes it with a binomial test against the
non-conversion rate — and this tool counts rather than infers; writing a
placeholder there would put a number in a column that means something it did not
measure. Nothing is lost: ALLCools converts a plain count table, deriving the
strand and the context from a reference FASTA, so a `.cov` file is one command
away from an ALLC.

## Conversions: the command, not the run

A recipe can say which converter follows an output, and `alnbase-agg convert`
prints the commands:

```toml
[output.CpG.bismark-cov]
path = "out/sample.CG.cov"

[output.CpG.bismark-cov.convert.allc]
reference = "genome.fa"
prefix = "out/sample.CG"
```

```sh
alnbase-agg convert recipe.toml > convert.sh   # needs no hit table
alnbase-agg run recipe.toml 'calls_*_*.parquet'
sh convert.sh
```

The commands are **printed, not run**. Three reasons, in order of how much they
matter:

- A recipe stays an ordinary data file. A preset is a pair of inputs you fork
  with `cp -r`, which stops being true the moment a recipe someone sent you can
  execute a command.
- The failure surface stays here. A converter missing from PATH or one that
  renamed a flag should not arrive as "the recipe failed" after an hour of
  aggregation that in fact succeeded. Ordering steps is what a shell script and
  Nextflow are for.
- The value stays the converter's. ALLC's seventh column is ALLCools'
  `np.round(mc / cov)` — a 0 or a 1, not methylpy's binomial test — and a command
  you run keeps whose number it is visible.

What the recipe is actually for is the part a person gets wrong. `table-to-allc`
is configured with **0-based column indices**, and the recipe is what chose the
column layout:

| source format | indices |
| --- | --- |
| `bismark-cov` | `--chrom 0 --pos 1 --mc 4 --uc 5` |
| `methyldackel-bedgraph` | `--chrom 0 --pos 2 --mc 4 --uc 5` |
| `cgmap` | `--chrom 0 --pos 2 --mc 6 --cov 7` |

`--pos` must name a **1-based** column, because the converter fetches
`[pos - 1, pos)` from the FASTA to decide which strand the cytosine is on. That
is why the two six-column formats differ: Bismark's coverage file is 1-based, so
its start is the position, while MethylDackel's bedGraph is 0-based half-open, so
its start is one too low and its *end* is the position. Given the wrong one the
converter raises nothing — on a test contig, `--pos 1` over a bedGraph turned
`chr1 11 + CGN` into `chr1 10 - NNN`: wrong position, wrong strand, wrong
context, exit status 0.

A `gzip = true` writer whose path does not end in `.gz` is refused, because
these converters infer compression from the extension and would read the
compressed bytes as text.

> **`allcools table-to-allc` is currently broken upstream**, in release 1.1.1 and
> on master. Its `@doc_params` decorator passes `input_path_doc=` while the
> docstring interpolates `{table_to_allc_input_path}`, so importing the module
> raises `KeyError` and every invocation dies before reading a byte
> (`ALLCools/table_to_allc.py:227,274`). Other subcommands are unaffected. The
> indices above were verified by calling the converter with that decorator
> neutralised: `bismark-cov` and `methyldackel-bedgraph` then produce
> byte-identical ALLC files.

## Filters are CTEs, not passes

A filter conditioned on an aggregate — premethyst's read rule, ScaleMethyl's
`CH.level <= 0.5` — is not a second program. It is a CTE and a join in the same
query, and DuckDB decides what to materialise:

```
src → hit → call → resolved → unit_records → unit_counts → judged → counted
```

`src` applies `[filter.hit]`. `call` collapses mates, so a position both mates
cover is one observation; a fragment that contradicts itself there has that
position dropped, because keeping either call would be a coin flip and
`alnbase overlap` is where mates are reconciled properly. `judged` applies the
read-filter ladders.

A ladder must tile the site counts exactly once — no gap, no overlap — so no read
falls through to an unwritten default. A read observing none of a context's sites
is kept by a rung the tool inserts if you leave it out; `explain` reports it,
because discarding those reads is what premethyst does and a recipe should have
to say which it wants. A rung whose predicate evaluates to NULL leaves the read
unclassified rather than guessing, the same way alnbase skips and counts a record
its strand rule cannot name.

## What the manifest decides

Four things are read from the run rather than assumed:

- **the coordinate base**, so the shift to 1-based output is derived and appears
  once in the generated SQL instead of as a literal `+ 1`;
- **the query names**, so a context naming a query the run never declared is
  refused while you are still looking at the recipe, rather than aggregating to
  zero everywhere;
- **the contigs and their lengths**, for a reference check;
- **the alnbase version and run id**, so a product can say what produced it.

The parquet footer carries the run's settings and is always there; the
`{stem}.manifest.json` beside the output is written only on success, so its
absence means the run may have been interrupted. That is a warning, not an error
— a hit table can be perfectly usable without a finished run's blessing.

## The reference check

A coordinate convention is the one thing here that can be wrong by exactly one
and produce output that looks perfect: every file writes, every count is
plausible, every plot has the right shape, and every cytosine is reported one
base to the right of where it is. Nothing downstream can catch it, because
nothing downstream knows where the cytosines were supposed to be.

The reference can. Hit rows carry `refr_base`, the reference base at the column
the row is anchored on, so the check needs no notion of what a context means:

```
alnbase-agg run recipe.toml --hits 'calls_*_*.parquet' -r genome.fa
```

Two checks, because two different things go wrong and they want different
instruments.

**Is this the reference the run used?** Exact and free — the manifest records
every `@SQ` contig with its length, so comparing that against the FASTA index
catches a wrong assembly, and catches `chr1` versus `1` completely rather than
statistically. No base is read. A contig the FASTA lacks is a warning (a run
against a full assembly is often checked against a primary-only FASTA, and the
shared chromosomes are still worth checking); a shared contig of a different
length is fatal, because it means a different assembly.

**Are the positions shifted?** Sampled, and scored at several candidate shifts
rather than only at the one intended. A report saying *"0 of 2000 match as
written, 2000 of 2000 match one base to the left"* names the bug; one saying only
"0 of 2000 match" leaves you to guess between a shift, the wrong assembly and a
bug in the checker. `--verify-sites` sizes the sample (0 checks contigs only) and
`--verify-shards` bounds how many files it is drawn from, since a shift is
systematic and more shards buy accuracy nobody needs.

The positions checked are the ones the query reports: the offset comes from the
same `Manifest.shift` the generated SQL uses, so the check cannot pass a run whose
SQL shifts differently.

One wrinkle, because getting it backwards would make the check pass on shifted
data: alnbase walks a read 5'→3' along its conversion strand, so on a
minus-strand read `refr_base` is the *complement* of what the FASTA holds there.
The comparison complements those rows, identifying them from the hit table's
`strand` column — the walk's own direction, not the FLAG's `is_reverse`, which
differs from it for exactly one mate of every pair. Without a `strand` column the
check accepts either orientation and says so, which matters most exactly where it
is weakest: a CG's C and the G opposite it are complements, so a one-base shift
inside a CG cannot be told from a correct call. Write the hits with
`-F strand` for the strict check.

A failure refuses the run before anything is written, because an hour of
aggregation that produces shifted files is worse than a refusal;
`--allow-reference-mismatch` downgrades it for a run whose reference is a near
relation of the one the reads were aligned to. Passing no reference at all is
allowed and says so on stderr — silence would leave you to notice the absence of
a note. Writers with somewhere to record provenance record the one-line result:
`amethyst-h5` puts it in `/metadata/reference_check`, so a file carries the
evidence its coordinates were checked rather than leaving it in a terminal that
has since been closed.

## Sharding does not enter the plan

`read_parquet` takes a glob and is correct whatever `--partition-by` the run
used, so nothing has to reason about whether a grouping happens to be
shard-local. What that costs is memory for the group-by, which DuckDB spills;
what it buys is that a hit table cannot be aggregated wrongly because it was
partitioned by the wrong key.

One connection also means one thread pool and one memory budget. `-t/--threads`
and `-m/--memory-limit` are **totals for the process**, so a caller running N
aggregations at once divides by N — several DuckDBs each sizing a pool from the
machine's core count is how a node with a handful of concurrent jobs ends up
oversubscribed by an order of magnitude.

## `explain` prints the query, not a description of it

A recipe cannot express everything, so the way out of the grammar is to take the
SQL and edit it. `explain` prints the statement that runs, indented and
commented; a test asserts that pasting it back gives the same rows, because if
that ever stops being true the escape hatch has closed.

This is also what makes presets useful beyond convenience: write one recipe per
comparator tool and **diffing two presets' generated SQL is a precise statement of
how the two tools' implicit extraction queries differ**, in a language a reader
can check.

## Python API

Products come back as DuckDB relations, so an analysis need not write a file and
read it back:

```python
import alnbase_agg as aa

recipe = aa.Recipe.load("recipe.toml")
agg = aa.Aggregation.open(recipe, hits="calls_*_*.parquet", threads=8)
print(agg.explain())
df = agg.relation("cg_sites").df()      # or .pl(), .arrow(), .fetchall()
```

With one output the relation is the standalone statement `explain` printed. With
several, the shared prelude is created as temp tables under the same names the
CTEs have, so the work happens once and each product's final `SELECT` is the same
text either way.
