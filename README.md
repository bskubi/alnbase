# alnbase

Compare aligned reads to the reference genome to find mismatches, base conversions, indels, and splice junctions in specific sequence contexts.

## Example: methylation calling

`CG_methylation.toml`
```
[query.methylated]
read = "C~"
refr = "CG"

[query.unmethylated]
read = "T~"
refr = "CG"
```

Run `alnbase query`
```
# Build reference index to generate mm10.aref
alnbase index mm10.fa mm10.aref

# Run alnbase query
alnbase query --parquet --query-file CG_methylation.toml --query-file query/strand/bwameth.toml input.bam mm10.aref output.parquet
```

Example results, in `output_0_0.parquet` (one row per hit)
```
ref_name  strand  name          off_5p  off_3p  refr_pos  qual  read_base  refr_base
chr1      +       methylated         2      11        12    40  C          C
chr1      +       unmethylated       6       7        16    40  T          C
chr1      +       methylated        13       0        23    40  C          C
```

`read_5p` and `read_3p` are the offset of the hit with respect to the 5' and 3' ends of the read respectively. `refr_pos` is the position of the hit with respect to the reference genome. `qual` is the PHRED score. `strand` is not based solely on the`is_reverse` flag, but is a reconstruction of the original native DNA strand that the read reports on. `name` is the name of the query that produced the hit. `read_base` and `refr_base` are the bases at the hit position. Note that any other BAM attributes, flags or tags can be written out as well.
## Why use `alnbase`?
+ Short, declarative queries
+ Many aligners work through a strand file
+ Per-read, single-cell-compatible output, not just pileups
+ Explicit edge case handling (indels, introns, read ends, clipped bases)

## Install

Requires [rust](https://rust-lang.org/tools/install/)
```
git clone https://github.com/bskubi/alnbase
cd alnbase
cargo test
cargo build --release
./target/release/alnbase --help
```

Note: alnbase builds htslib from source. This requires a C compiler, cmake, libclang, and the zlib, bzip2, xz, curl and OpenSSL development headers. On Debian or Ubuntu: `sudo apt-get install build-essential cmake libclang-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev`.
## Documentation

See `./docs/`

## Commands

| command | what it does |
| ------- | ------------ |
| `alnbase index` | build the reference index (`.aref`) from a FASTA |
| `alnbase info` | report what is in a reference index |
| `alnbase query` | run query files over a BAM, writing BAM tags or parquet hit rows |
| `alnbase extract` | write parquet hit rows from a tag in a BAM tagged by `alnbase query` |
| `alnbase dump-query` | print the query files stored in a tagged BAM's header |
| `alnbase overlap` | resolve bases observed twice in one fragment, where the mates overlap |

## Repository layout

```
alnbase/
├── src/               the alnbase program (Rust), with its unit tests
├── query/
│   ├── strand/        one strand file per aligner, and records/ that test them
│   └── extract/       query files that reproduce other tools' calls (Bismark, MethylDackel)
├── python/            alnbase-agg: turns hit tables into bedGraph, .cov, CGmap, Amethyst H5, ...
├── docs/              user documentation
├── Cargo.toml
└── LICENSE
```

## Status

Pre-release. Formats may change.

## License

MIT licence; see `LICENSE`.
