# alnbase

Query aligned read/reference columns.

alnbase compares a read against the reference it aligned to, column by column,
and answers questions you write down: *where does this read show a T over a
reference CG?* The answers come back either as tags on the reads, or as a table
of one row per hit.

It was built for bisulfite and enzymatic methylation data, where the question
is some variation on "which cytosines were converted, in which sequence
context", but nothing in it is specific to methylation. A query is a pattern
over aligned columns, and the patterns say what you want them to say.

```toml
# calls.toml -- a protected CpG cytosine, and a converted one
[query.Z]
read = "C~"
refr = "CG"

[query.z]
read = "T~"
refr = "CG"

[tag.XM.bases]
fill = "."
Z = "Z"
z = "z"
```

```bash
alnbase index mm10.fa mm10.aref
alnbase query --query-file calls.toml --query-file queries/strand/directional.toml \
  -@ 8 in.bam mm10.aref out.bam
alnbase extract out.bam calls.parquet
```

The second file is the **strand rule**: how this aligner records the strand a
read was converted on. `queries/strand/` ships one per aligner, every run
declares one, and nothing about any aligner is built into the binary.

Every record comes out in input order carrying an `XM` string, one character
per base. `extract` turns those calls back into a table without needing the
reference, because the query file travels in the BAM's header.

## Why

A methylation caller has the rules baked into it: what counts as a CpG, what to
do when a deletion falls inside the context, whether to look past the end of a
read. Those decisions are correct for the common case and invisible when they
are not.

alnbase makes them something you write in a file, and that file travels with
the output. `alnbase dump-query out.bam` prints the definitions a BAM was
tagged with, six months later.

Three consequences worth knowing:

- **Reference context is not limited to the read.** A cytosine at the last base
  still has its CG context.
- **Indels are explicit.** Columns are CIGAR-aligned, so what happens when a
  deletion lands inside a context window is decided by your pattern.
- **No filters that other tools already have.** alnbase does not skip
  duplicates, secondaries or low-MAPQ reads: `samtools view -F/-q` upstream
  does that, and a second implementation would be one more thing to keep in
  step.

## Installing

Rust 1.85 or newer, a C compiler, and the headers htslib needs.

```bash
cargo build --release
```

For a faster binary, `RUSTFLAGS="-C target-cpu=x86-64-v3"` enables the SIMD
CRAM decoders; the result needs a Haswell-era or newer CPU. For profiling,
`cargo build --profile profiling` is the same build with debug info.

## Commands

| | |
|---|---|
| `alnbase index FASTA INDEX` | convert a reference FASTA into alnbase's format; rebuild indexes made by 0.1.1 or earlier |
| `alnbase info INDEX` | what an index holds: version, contigs, total bases, each contig's MD5 |
| `alnbase query [OPTS] BAM REFR OUT` | tag a BAM, in input order; `--parquet` writes hit rows instead |
| `alnbase extract [OPTS] BAM OUT` | turn a tag back into hit rows, no reference needed |
| `alnbase dump-query BAM` | print the query file stored in a BAM's header |
| `alnbase overlap [OPTS] IN OUT` | resolve mate overlap before calling |

`-` works for standard input and output, so the tools pipe:

```bash
samtools collate -O in.bam | alnbase overlap - - \
  | alnbase query --query-file calls.toml --query-file queries/strand/directional.toml \
      - mm10.aref - > out.bam
```

## Documentation

- [`docs/cli-reference.md`](docs/cli-reference.md) — every command-line
  option, its default, and when to reach for it.
- [`docs/alnbase.md`](docs/alnbase.md) — the full guide: the column model, the
  query file format, every command, and how to inspect a query before running
  it.
- [`docs/reference/06-aggregation.md`](docs/reference/06-aggregation.md) — turning
  hit rows into per-position, window, per-cell and per-read products with
  `alnbase-agg`, a Python package in [`python/`](python/) that compiles a
  declarative recipe into one DuckDB query.
  [`python/examples/cg_ch.toml`](python/examples/cg_ch.toml) is a commented recipe
  to copy.
- [`docs/bismark-xm.toml`](docs/bismark-xm.toml) — Bismark's eight `XM` codes,
  including `u`/`U` for cytosines whose context the reference does not
  determine. A working file to copy.

## Before trusting the output

alnbase has not yet been compared against Bismark or MethylDackel on real data.

## Status

Prototype. The query file format and the BAM tag layout are settled enough to
build on; command-line details may still move.

## Licence

Not yet chosen — see the TODO in `Cargo.toml`. Until then, all rights reserved
by the authors.
