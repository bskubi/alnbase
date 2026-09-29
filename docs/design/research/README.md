# Phase-0 research notes

Source-level surveys written on 2026-09-16 to inform
[`../aggregation-and-validation.md`](../aggregation-and-validation.md).

- `methylation-tools.md`: bulk methylation extractors and their implicit queries.
- `singlecell-formats.md`: Amethyst H5, premethyst cellInfo, and other single-cell formats.
- `write_amethyst_v2_minimal.py`: minimal Amethyst v2.0.0 writer, run once with h5py
  3.6.0 and checked with `h5ls`; not yet read back by the R package.
- `other-assays.md`: non-methylation read-vs-reference assays and their comparators.
- `strand-determination.md`: how eleven aligners say which strand a read came from, written
  2026-09-17 to scope declarative strand handling.

Citations of the form `path:line` refer to shallow clones of each tool's repository,
made in a temporary `clones/` directory that is not part of this repo. The
commits and versions are listed at the top of each file, so any citation can be
resolved against the upstream repository at that commit.

Comparator behaviours here were found by reading source unless marked "verified by
execution". Treat suspected bugs as hypotheses until a validator reproduces them.
