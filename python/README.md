# alnbase-agg

Aggregate [alnbase](../README.md) hit tables with DuckDB. **A run is a SQL
script.** There is no recipe format, no expression vocabulary and no query
builder: you write the query, and this runs it.

```sh
alnbase-agg run run.sql -s hits='calls_*_*.parquet' -s out=pbmc.CG.cov
```

## What a run looks like

A few declarations and an include. This is a complete, working run:

```sql
-- 1. The hits, and the hit-level filter.
CREATE OR REPLACE VIEW src AS
SELECT *,
       ref_name AS contig,
       {'qname': qname, 'mate': is_first_in_template} AS read,
       'bulk' AS block
FROM read_parquet(getvariable('hits'))
WHERE qual >= 20;

-- 2. What the alnbase query names mean. The only place biology appears.
CREATE OR REPLACE VIEW measure (measure, name, category) AS VALUES
  ('mCG', 'CG', 'methylated'), ('mCG', 'TG', 'unmethylated');

-- 3. What is set against what. `main` is the side in focus.
CREATE OR REPLACE VIEW relation (relation, measure, main, other) AS VALUES
  ('mCG', 'mCG', ['methylated'], ['unmethylated']);

.read lib/filters.sql

-- 4. Which observations survive.
CREATE OR REPLACE VIEW resolved AS SELECT * FROM src
WHERE read NOT IN (SELECT read FROM read_fails('mCG', 'ratio', '>', 0.9));

-- 5. Count the hits into sites, and write.
.read lib/count_sites.sql
SET VARIABLE relation = 'mCG';
.read formats/bismark_cov.sql
```

`examples/bulk_bismark.sql` is a runnable version of it that also merges the
two strands of each CG. `examples/single_cell_methylation.sql` is the same
shape with a block key, four read filters, a block filter and three output
formats in it. `examples/toy/` is smaller than either: a 48-base reference and
30 hits you can check by hand, and `tour.py` beside it prints every concept
once with the answer next to it.

## The two directives

A line starting with `.read` or `.python` is a directive rather than SQL.
`.read` is DuckDB's own CLI spelling for "execute that file here". `.python`
is its counterpart for the one or two formats SQL cannot finish alone.

Both resolve against the including file's directory first, then against the
shipped library. So forking a format is a copy and an edit:

```sh
cp "$(alnbase-agg lib)/formats/bismark_cov.sql" ./formats/
$EDITOR formats/bismark_cov.sql        # this copy is now the one that is read
```

## Values

DuckDB session variables, which are late bound, so one format file written once
produces every output it is asked for:

```sql
SET VARIABLE relation = 'mCG';
SET VARIABLE out = 'sample.CG.cov';
.read formats/bismark_cov.sql

SET VARIABLE relation = 'mCH';
SET VARIABLE out = 'sample.CH.cov';
.read formats/bismark_cov.sql
```

`-s name=value` on the command line sets one before the script starts. A value
that reads as a number becomes one, so `-s min_depth=5` is an integer; a
barcode like `01` stays a string.

## What ships

`alnbase-agg lib` prints the directory. In it:

| file | what it defines |
| --- | --- |
| `lib/filters.sql` | `read_category`, `read_relation`, `block_*`, and the filters |
| `lib/count_sites.sql` | `site_stranded`, `site_all`, `site_relation`, `block_summary`, `pos1`/`pos0`, the site filters |
| `lib/merge_strands.sql` | collapse a CG's two cytosines onto one coordinate |
| `lib/load_contigs.sql` | contig order from a `.fai`, for anything that sorts |
| `formats/bismark_cov.sql` | Bismark `.cov`: 1-based, six columns |
| `formats/bedgraph.sql` | Bismark bedGraph: 0-based, four columns |
| `formats/methyldackel_bedgraph.sql` | MethylDackel's: 0-based, truncated percent |
| `formats/bedmethyl.sql` | ENCODE bedMethyl: eleven columns, keeps strand |
| `formats/cgmap.sql` | BS-Seeker2 CGmap: every context in one file |
| `formats/cellinfo.sql` | premethyst/Amethyst cellInfo: one line per block |
| `formats/amethyst_h5.sql` + `.py` | Amethyst HDF5 v2: one dataset per block |
| `formats/sites_parquet.sql` | the site table itself, losing nothing |

Every filter returns what *fails* it, so they compose by union and a run reads
as one `WHERE ... NOT IN (...)`. A filter is three choices — what it groups on,
what it computes, and when it fails — so three macros cover the whole grid:

```sql
read_fails (relation, quantity, comparison, threshold)
block_fails(relation, quantity, comparison, threshold)
site_fails (relation, quantity, comparison, threshold)
```

`quantity` is `'main'`, `'other'`, `'total'` or `'ratio'`; `comparison` is one
of `< <= > >= = <>`. So the bisulfite non-conversion cut is
`read_fails('mCH', 'ratio', '>', 0.4)` and the empty-droplet cut is
`block_fails('mCG', 'total', '<', 1000)`.

A quantity that is not a plain aggregate over one group's own totals keeps a
macro of its own:

| filter | rejects |
| --- | --- |
| `read_fails_consecutive(rel, op, n)` | a read with too long a run of adjacent main-side calls |
| `read_fails_missing(m)` | a read with no hits of the measure at all |
| `read_fails_listed(path)` | a read whose name appears in a file |
| `block_fails_reads(rel, op, n)` | a block that contributed too few reads |
| `block_fails_sites(rel, op, n)` | a block that covered too few positions |
| `site_fails_blocks(rel, op, n)` | a position too few blocks covered |
| `site_fails_in_regions(path)` | a position inside any region of a BED file |
| `site_fails_outside_regions(path)` | a position no region of a BED file names |
| `site_fails_position(rel, q, op, x)` | a position another relation vetoes, on evidence pooled across blocks |

Filters compose by union, so `both(...)` in the Python front end is the other
connective: a rejection that is a conjunction rather than one comparison. It is
what says "neither wholly methylated nor wholly unmethylated", and it is how a
threshold that depends on how much evidence there is — premethyst's ladder — is
written, as one call per rung. In a hand-written script a ladder is a `CASE`.
Nothing checks that the rungs tile the range: an overlap is harmless, and a gap
quietly keeps the reads it should have judged.

## The other commands

```sh
alnbase-agg show run.sql               # the script with every .read substituted
alnbase-agg lib --list                 # what the library holds
alnbase-agg check --hits 'calls_*.parquet' --reference genome.fa
```

`check` is the coordinate check, and it involves no script: it looks the
reference base up in the FASTA and scores the hits at several candidate shifts.
An off-by-one is the one error here that writes cleanly, counts plausibly and is
wrong at every base.

## Installing

```sh
pip install -e .                       # duckdb, click, pyfaidx
pip install -e '.[amethyst]'           # and h5py, numpy, pyarrow for the H5 writer
```
