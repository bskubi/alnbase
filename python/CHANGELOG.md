# Changelog — alnbase-agg

This package versions separately from the alnbase binary, which has its own
`CHANGELOG.md` at the repository root.

## 0.8.1 — 2026-09-29

- **Licence declared** (`pyproject.toml`). The package is MIT-licensed, matching
  `LICENSE` at the repository root. The SPDX `license` string needs setuptools 77
  or later, so the build requirement moved up from 68.

## 0.8.0 — 2026-09-23

**Every filter the target formats need is now expressible without raw SQL.**
Going through what Bismark, MethylDackel, premethyst, ScaleMethyl, Amethyst,
scbs and the hyper-editing pipeline actually filter on turned up five things
the syntax could not say. All five are closed.

**Blocks can be judged on reads and on positions.** `mCG.reads.lt(10000)` and
`mCG.sites.lt(100)` join the quantity set at block level, over a new
`block_extent` view. Neither is derivable from `num` and `den` — a block with
2000 hits may be one read seen 2000 times or 2000 reads seen once — and
premethyst's `unique_reads >= 10000` and its `calls-filter -G` minimum covered
CG *positions* are the two commonest cell cuts in this ecosystem.

**A BED file works in both directions.** `in_regions(path)` is a blacklist and
`outside_regions(path)` restricts a run to a capture panel, which is what
premethyst's `context-extract` does. The second is not reachable by inverting
the first without complementing the BED by hand. `regions()` is renamed to
`in_regions()` with no alias left behind.

**A position can be vetoed on pooled evidence.** `at_position(test)` judges a
different relation at the same position, summing across blocks and strands, and
rejects the position in every block:

```python
run.drop_sites(at_position(snp.relate("G", "A").ratio.gt(0.5)))
```

That is the SNP veto named as a canonical site stage in the design notes: a C
read as a T may be conversion or a C-to-T variant, and the evidence is on the
other strand. Pooling matters — testing block by block would veto a position
because one cell saw one variant read.

**`both()` is the conjunction.** Filters compose by union, which is right: each
is an independent reason to throw something away. `both` covers the cases where
one reason is itself a conjunction — scbs drops a site that is neither wholly
methylated nor wholly unmethylated, a condition on both sides of the relation
at once that no single comparison can state.

**A ladder needs no mechanism of its own.** premethyst's read cut is a
threshold conditional on how much evidence there is, and it falls out of `both`
plus what already existed:

```python
run.drop_reads(
    both(mCH.total.ge(3), mCH.total.le(5), mCH.main.gt(1)),
    both(mCH.total.gt(5), mCH.ratio.ge(0.4)),
)
```

What is *not* recovered: the deleted recipe language required a ladder's rungs
to tile the range exactly once, no gap and no overlap. Here an overlap is
harmless, because the parts are unioned, and a gap simply means those reads are
kept — and nothing will tell you that you left one.

**`quantity_value` names every branch.** The `ELSE` used to fall through to the
ratio, so adding a fifth quantity to the `ENUM` and forgetting a branch would
have computed a ratio silently. It raises now.

**A toy dataset ships.** `examples/toy/` holds a 48-base reference and 30 hits,
written out one per line with no randomness, plus `tour.py`, which walks seven
concepts printing what it computed beside what it should be. Four assays share
the one table, deliberately: the hardest thing to believe is that a measure is
only a set of query names. `tests/test_toy.py` pins that regenerating gives
identical bytes, that every claim the tour makes holds, and that the fixture
stays under 40 hits and 100 bases so it cannot quietly stop being checkable.
The tour gained two steps for the new filters, so it now makes eleven claims.

**The chained example writes its ladder with `both`.**
`examples/single_cell_methylation.py` spelled premethyst's read cut as a block
of raw SQL passed to `fails(...)`, because nothing else could say it. It is now
two `both(...)` calls, and the hand-written `.sql` twin keeps its `CASE` so the
pair still shows the same run written two ways. The outputs remain byte-identical.

**A header's second `REQUIRES` line is no longer ignored.** `_requires`
returned at the first one, so a file splitting a long list over two lines would
have had half of it honoured and nothing said about the rest -- the include
would simply be absent when the format asked for the view. Every `REQUIRES`
line in the header now counts, and the header still ends at the first line that
is not a comment. `formats/amethyst_h5.sql` had its `.python` line indented
directly under `REQUIRES`, where it read as a continuation of the declaration;
it is documentation, the Python half being paired by name, so it moved up into
the usage block with the other directives.

**The prose caught up with four rounds of renaming.** A sweep of every comment,
docstring and example turned up references to `site_ratio`, `fails_above_rate`,
`fails_listed`, `fails_ratio`, `merge_cg`, `positive`/`negative`, `num`/`rest`,
`depth`, `count`, `unit`, `cell` as the partition and `recipe` as a file format,
plus a handful of wrong counts — "29 hits" for 30, "four measures" for six,
"five classes" for six, "80 bases" for 79. All are corrected, in the package,
the examples, the test names and both published explainers.

## 0.7.1 — 2026-09-23

**A relation's sides are `main` and `other`.** The principle is that the first
side is the one in focus, and the names now say so:

```sql
CREATE OR REPLACE VIEW relation (relation, measure, main, other) AS VALUES
  ('mCG', 'mCG', ['methylated'], ['unmethylated']);
```

```python
mCH.main.gt(2)      mCH.other.gt(5)      mCG.total.lt(3)      mCH.ratio.gt(0.4)
```

`num` was worse than it looked: "numerator" is ratio vocabulary, and the point
of a relation is that a ratio is one *reading* of it rather than what it is.
`rest` was worse still — true when the second side is the complement, false
when you write `relate("G", "A")` and side two is a named category rather than
the rest of anything.

`signal`/`background` was considered and rejected: both sides of a relation are
real observations. An unmethylated CG is a genuine negative call, not
background, and in `relate("G", "A")` neither side is noise.

The change of vocabulary between layers is now deliberate and documented. A
relation is declared and filtered in relation words — `main`, `other` — while a
file format writes a ratio, where `num` and `den` are the right words. Each
layer speaks its own, and `quantity_value` is where they meet.

`Relation.main_categories` and `.other_categories` hold the category lists.

**The package README was rewritten.** It still described `unit`, `cell`,
`fails_rate`, `positive`/`negative` and nine separate filter macros, none of
which exist.

## 0.7.0 — 2026-09-23

**The quantities are named for the relation's own two sides.** `count` and
`depth` are gone; a relation now yields `num`, `rest`, `total` and `ratio`:

```python
mCH.num.gt(2)       # more than two methylated CH calls, whatever the ratio
mCH.rest.gt(5)      # more than five unmethylated ones
mCG.total.lt(3)     # saw fewer than three CG positions at all
mCH.ratio.gt(0.4)   # more than 40% methylated
```

You declare `num=[...]` and `rest=[...]`, so `num` and `rest` are how many hits
fell on each — there is no second vocabulary to learn. `total` is the two
together, which is the `den` a file format writes.

`depth` in particular was a trap: per read it meant "how much of the relation
this read saw", and at site level the same word reads as coverage — how many
reads hit a position. Same word, two groupings, two numbers. `total` says the
arithmetic instead of implying a biological quantity.

`Relation.num` and `.rest` are now the quantities, so the category lists that
define them are `num_categories` and `rest_categories`.

## 0.6.2 — 2026-09-23

**A depth threshold no row can meet is refused.** The counting views are
sparse: a read, block or site gets a row because it had a hit, so its depth is
never below 1. `depth.lt(1)` therefore binds, runs and matches nothing — which
is the failure this library works hardest to prevent, and which the `ENUM`
guards on measure and relation names exist to stop in the neighbouring case.

`mCG.depth.lt(1)`, `.le(0)` and `.eq(0)` now raise, and the message says what
the call meant: a read that saw none of the measure has no row in the counting
views at all, so the filter for it is `missing()`, which reads `src` instead.

The guard found three filters in this repository's own tests and one in a check
harness that could never have fired.

## 0.6.1 — 2026-09-23

**The comparison methods are the standard two-letter names.** `gt`, `ge`, `lt`,
`le`, `eq`, `ne` replace `above`, `at_least`, `below`, `at_most`, `equals` and
`not_equals`:

```python
run.drop_reads(mCH.ratio.gt(0.4), mCH.consecutive.gt(3))
run.drop_blocks(mCG.depth.lt(1000))
```

They are the names everyone already knows from `operator`, from shell test, and
from every comparison API in the language, so there is nothing to learn and no
guessing at whether a boundary is included — `ge` is unambiguous where
`at_least` had to be read twice.

The genotype measure in the examples and tests is renamed from `gt` to `snp`,
because a measure called `gt` beside a method called `gt` is a sentence nobody
should have to parse.

## 0.6.0 — 2026-09-23

**A filter is three choices, and each is written where it belongs.** Eight
threshold macros with the direction baked into their names became two:

```sql
read_fails (relation, quantity, comparison, threshold)
block_fails(relation, quantity, comparison, threshold)
site_fails (relation, quantity, comparison, threshold)
```

`quantity` is `'count'`, `'depth'` or `'ratio'`; `comparison` is one of
`< <= > >= = <>`, guarded by an `ENUM` so a typo fails at bind time. So "blocks
whose depth is below 1000" is `block_fails('mCG', 'depth', '<', 1000)` and needs
no macro of its own. What the operators mean is one `compare` macro, so no
filter can disagree with another about a boundary.

This is more expressive as well as shorter. `> 0.4` and `>= 0.4` are different
filters and premethyst's read cut uses the second; with the direction in a macro
name there was no way to ask for it, which is why that cut had to be written out
as raw SQL.

In Python the three choices are three parts of one expression — the quantity is
a property, the comparison a method, and the grouping comes from the verb:

```python
run.drop_reads(mCH.ratio.above(0.4), mCH.consecutive.above(3))
run.drop_blocks(mCG.depth.below(1000))
run.drop_sites(mCG.depth.below(5))
```

The same `mCG.depth.below(5)` renders at whichever level it is handed to. A
quantity that is not a plain aggregate over one group's own totals keeps a macro
of its own and says where it applies: `consecutive` is a window over positions
and exists only for reads, `blocks` is a count across blocks and exists only for
sites. Using one at the wrong level is refused by name rather than binding and
quietly dropping nothing.

**Position filters.** A run can now drop sites, which nothing before could
express:

```python
run.drop_sites(mCG.depth.below(5))        # too little coverage here
run.drop_sites(mCG.blocks.below(3))       # only one cell saw this position
run.drop_sites(regions("blacklist.bed"))  # ENCODE blacklist, SNPs, off-target
```

A site is one position in one block, so these judge a position per block; a
filter meant to apply everywhere, such as a BED file, rejects the position in
every block. `lib/count_sites.sql` grew a second indirection — `site_relation`
over `site_all` — so a run narrows sites by redefining a view, exactly as
`merge_strands` folds coordinates by redefining another.

This is not a format's `min_depth`, which stays: that is what one file shows and
is set per output. Site filters are what the run keeps, are recorded, and appear
in `report()` beside every other filter.

**A relation, not a ratio.** What a run declares is a relation between two
disjoint sets of categories; ratio, count and depth are three ways of reading
one:

```sql
CREATE OR REPLACE VIEW relation (relation, measure, num, rest) AS VALUES
  ('mCG', 'mCG', ['methylated'], ['unmethylated']);
```

`rest` no longer repeats `num` the way `den` did, matching how the categories it
is built from do not overlap, and `count` and `depth` no longer have to be keyed
on a measure because they are not about ratios. `Measure.ratio()` is
`Measure.relate()`, and the session variable a format reads is `relation`.

## 0.5.2 — 2026-09-23

**`block` takes several columns.** A run whose hits carry a condition as well
as a barcode makes one block out of both:

```python
Run(hits, block=["condition", "XB"])      # ctrl.AACGTT
Run(hits, block=["condition", "XB"], block_sep="_")
```

They are joined rather than made into a struct, which is the opposite of what
`read` does with several columns, and deliberately: nothing ever prints a
`read`, it is only compared and counted, while a block is written out — it is
the barcode column of a cellInfo file and the name of a dataset inside an
Amethyst H5, and a struct there would be unreadable.

`concat_ws` skips a NULL rather than poisoning the whole key, so a row missing
one part lands in the block named by the rest. A run that would rather see that
fail writes the expression itself, which `block` has always accepted.

## 0.5.0 — 2026-09-23

**A measure holds any number of categories; a ratio relates two of them.** The
`Call` and `Count` classes are gone, and with them the reserved category names
`positive` and `negative`. Both were imprecise: a mathematician would call a
ratio a count too, and a class named for its output sat oddly beside argument
names describing its input.

What replaces them is one class and one relation. A `Measure` sorts alnbase
query names into categories you name yourself. A `Ratio` names one set of
categories against another, and is what every proportion filter tests and every
file format writes:

```python
mCG = Measure("mCG", methylated="CG", unmethylated="TG")
gt  = Measure("gt", A="A>A", C="A>C", G="A>G", T="A>T")

mCH.rate(0.4)                 # two categories: the only ratio, unnamed
gt.ratio("G", "A").rate(0.9)  # G among G and A
gt.ratio("G")                 # G among every category
gt.ratio(["G", "A"], "C")     # either side may gather several
```

A measure with exactly two categories has one possible ratio, declares it under
the measure's own name, and can be handed straight to a format, so the common
case never mentions the word. A measure with more has to say which pair,
because there is no default worth guessing.

In SQL this is a second declaration table beside `measure`:

```sql
CREATE OR REPLACE VIEW ratio (ratio, measure, num, den) AS VALUES
  ('mCG', 'mCG', ['methylated'], ['methylated', 'unmethylated']);
```

`num` and `den` are lists, and `den` is written out in full rather than implied.
Filters and formats take a ratio name, so every macro is back to two arguments:
`fails_rate('mCH', 0.4)`. `fails_shallow` and `fails_missing` still take a
measure, being about depth rather than proportion. `read_category` and
`sample_category` hold the per-category counts; `read_ratio`, `sample_ratio`,
`site_ratio` and `sample_summary` hold `num` and `den` per declared ratio.

**A query name belongs to at most one category *of a given measure*, and to no
category of the measures it has nothing to do with.** The previous wording said
"exactly one", which was wrong: a NOMe-seq run declares an accessibility measure
over its GC queries and a methylation measure over its CG queries, sharing no
query names, and a run that declares only mCG has CH query names in no category
at all. Only the within-measure partition is required, and only because it lets
a denominator be a plain sum.

**`unit` and `cell` are now `read` and `block`.** `unit` said nothing about
what it was a unit of. `cell` was too narrow: the same column is a cell in a
single-cell run, a spot in a spatial one, and one library otherwise. `block`
names what the column actually is — a partition over reads — and carries no
assay with it. `drop_cells` is `drop_blocks`, `cell_fails_*` is
`block_fails_*`, `cell_measure` and `cell_summary` are `block_ratio` and
`block_summary`, and `Run(hits, cell=...)` is `Run(hits, block=...)`.

`sample` was tried first and rejected: it reads as the biological specimen, so
a run of three donors by a thousand cells would have to call each cell a
sample. DuckDB also quotes `sample` in query plans, which broke the streaming
sort check until it was taught to strip quotes; `block` is plain there. A
session variable called `sample` still means what it always did — the output
filename prefix — and is untouched.

**`ref_name` is `contig` inside the library.** `ref_name` beside `refr_pos` was
a spelling difference with no meaning behind it. `src` now writes `ref_name AS
contig` and the library reads `contig` throughout, so the run works against
today's hit tables and keeps working when alnbase renames its own column.

**`Run.report()`** replaces `Run.counts()`, which collided with the old `Count`
class and said less about what it does.

`check_sorted` now tolerates a quoted identifier in the query plan, which
DuckDB emits for a column named `sample`.

## 0.4.0 — 2026-09-23

**A run is a SQL script.** The TOML recipe format, the compiler that turned it
into SQL, and the writer registry are gone. What replaces them is a runner for
`.sql` files and a library of SQL to include.

The reason is that nothing in an aggregation needed a configuration language.
alnbase's query files earn TOML because the patterns they declare are spatial —
a layout over aligned read and reference columns has no natural expression in a
query language. A filter threshold and a measure definition have no such
excuse: they are ordinary parameters to an ordinary query. The technical
argument for generating SQL turned out to be false as well. DuckDB refuses a
bound parameter inside `CREATE MACRO`, from which it does not follow that
thresholds must be written into the SQL text: `SET VARIABLE` / `getvariable()`
is late bound and works in macro bodies, in view bodies, in `read_parquet` and
as a `COPY` target, so a shipped SQL file is parameterisable with no templating
and no builder at all.

### What a run is now

- **`alnbase_agg.script`** executes a `.sql` file. It adds two directives to
  what `con.execute()` already does: `.read` (DuckDB's own CLI spelling for
  "execute that file here") and `.python`, for the formats SQL cannot finish
  alone. Both resolve against the including file's directory first and the
  shipped library second, so forking a format file is a copy and an edit.
- **Values arrive as session variables**, set with `-s name=value` or written
  into the script. A value that reads as a number becomes one; a barcode like
  `01` stays a string. `getvariable` on an unset name is NULL, so a format file
  gives an option a default with `coalesce` and there is no required-flag list
  anywhere.
- **Statements run one at a time**, so a failure is quoted with its file, its
  line and DuckDB's own message. A statement that returns rows has its result
  printed, which is what makes a QC query in a run script worth writing.
- **`alnbase-agg show`** prints a script with every `.read` substituted: the
  query a run will issue, in one piece, before it runs.

### The SQL library

- **`lib/filters.sql`** — `unit_measure` and `cell_measure`, then nine filters.
  Each returns what *fails* it, so they compose by union and a run reads as one
  `WHERE ... NOT IN (...)`: `fails_rate`, `fails_below_rate`, `fails_count`,
  `fails_shallow`, `fails_missing`, `fails_consecutive`, `fails_listed`,
  `cell_fails_shallow`, `cell_fails_rate`. `fails_consecutive` finds a run of
  adjacent positive calls along a read, which a rate provably cannot: four in a
  row and four scattered give the same rate, and only the first looks like a
  patch of failed conversion.
- **A mistyped measure name is now an error rather than a filter that rejects
  nothing.** `CREATE TYPE ... AS ENUM` is derived from the run's own `measure`
  declarations and filters cast against it. This cannot be done with an
  assertion inside a query: DuckDB's optimizer prunes `error()` out of a
  predicate, a `UNION ALL` branch and an unused projection alike.
- **`lib/count_sites.sql`** — `site_stranded` (one row per cytosine per cell),
  `site_measure`, `cell_summary`, and the `pos1`/`pos0` macros that hold the
  coordinate convention. No format file writes `+ 1`.
- **`lib/merge_strands.sql`** — collapses a CG's two cytosines onto the plus-strand
  one. A separate include because it is a choice, not a fact.
- **`lib/load_contigs.sql`** — contig order from a `.fai`, because `ORDER BY chr` is
  lexical and puts chr10 between chr1 and chr2.

### Formats

`bismark_cov`, `bedgraph`, `methyldackel_bedgraph`, `bedmethyl`, `cgmap`,
`cellinfo`, `amethyst_h5` (SQL plus a Python half) and `sites_parquet`. Each
reproduces its tool's conventions rather than a generic one: Bismark's
fifteen-significant-digit percentage, MethylDackel's truncating C cast,
BS-Seeker2's banker's rounding with a trailing `.0`, Amethyst v2's `(chr S10,
pos, c, t)` dtype with no `pct`.

The Amethyst writer checks three things that would otherwise produce a file
that opens cleanly and is wrong: that the stream is sorted on the cut key, that
no contig name exceeds the `S10` the dtype stores, and that no contig in the
data is missing from the reference index.

### The reference check

`alnbase-agg check --hits ... --reference genome.fa` is its own subcommand
rather than a flag on `run`. The question it answers — are these coordinates
where they say they are — is about a hit table and a FASTA, and is asked once
per alnbase run rather than once per output.

### Removed

`Recipe`, `compile_recipe`, `Aggregation`, the writer registry and its entry
point group, and `alnbase-agg explain` / `convert`. The presets ship `run.sql`
in place of `recipe.toml`.
