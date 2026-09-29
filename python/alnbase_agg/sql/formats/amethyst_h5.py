"""Amethyst HDF5 v2, one dataset per block. The Python half.

Run by `.python formats/amethyst_h5.py`, after `formats/amethyst_h5.sql` has
defined the sorted `amethyst` view. Two names are already in scope: `con`, the
connection carrying every view the run has defined, and `var`, which reads a
session variable.

All this does is walk a sorted stream and call `create_dataset` per block, which
is the only thing HDF5 needs that SQL cannot do. Nothing larger than one block
is ever in memory, so the peak cost is the biggest block rather than the genome.

Two things are checked before a byte is written, because both produce a file
that opens cleanly and is wrong:

* **the sort**, from the query plan, since a stream not cut on `block` would
  write a block's rows into several datasets and truncate all but the last;
* **contigs the reference index did not name**, which the view's inner join
  drops -- a hit table and a FASTA from different assemblies.

A third hazard, a contig name too long for the `chr` field, is handled rather
than refused: the field is widened to fit, which is what facet does.
"""

import h5py
import numpy as np
import pyarrow.compute as pc

from alnbase_agg.stream import check_sorted, groups

# `amethyst_facet/h5/dataset.py:14-17`, field for field and in that order.
#
# `c` before `t` is deliberate and is worth knowing, because premethyst's own
# files have them the other way round. premethyst writes v1, whose dataset is
# `(chr, pos, pct, t, c)` (`calls2h5.py:34-41`), matching the `.cov` files it
# reads. This is v2, which drops `pct` and is defined by amethyst-facet, whose
# `observations_v2_dtype` puts `c` first. Both readers that matter take fields
# by name -- facet's `Dataset.convert_dtype` copies "by field name, never by
# position", and Amethyst's R side turns the dataset into a data.table with the
# file's own column names -- so the order is a matter of matching the format's
# definition rather than of correctness.
OBSERVATIONS_V2 = [("chr", "S10"), ("pos", "<i8"), ("c", "<i8"), ("t", "<i8")]
VERSION = "amethyst2.0.0"

#: The narrowest `chr` field the format uses. Wider is written when a contig
#: name needs it: S10 fits a primary assembly, but truncating is not a safe
#: fallback, because chr1_KI270706v1_random and chr1_KI270707v1_random both
#: become b"chr1_KI270" and their observations are then summed together. facet
#: widens for the same reason (`h5/dataset.py:20-37`), so a file written here
#: and a file written there hold the same widths.
CHR_MIN_WIDTH = 10

# The dataset is `/<context>/<barcode>/1`; `1` is the conventional name for
# base-level observations, with windowed aggregations alongside it under their
# own names. amethyst-facet computes those, so nothing here writes them.
BASE_LEVEL = "1"

SQL = "SELECT block, chr, pos, c, t FROM amethyst"

out = var("out")
context = var("context") or var("relation")
if out is None or context is None:
    raise ValueError(
        "formats/amethyst_h5.py needs `out` and one of `context` or `relation`: "
        "SET VARIABLE out = 'sample.h5'; SET VARIABLE context = 'CG';"
    )

check_sorted(con, SQL, ["block"])

# What the view's inner join against `contigs` would drop. Counted from
# `site_relation` rather than from `amethyst`, which is the join's output and so
# cannot show what is missing from it.
orphans = con.execute(
    "SELECT DISTINCT contig FROM site_relation "
    "WHERE relation = getvariable('relation') "
    "  AND contig NOT IN (SELECT name FROM contigs) "
    "ORDER BY contig LIMIT 5"
).fetchall()
if orphans:
    names = ", ".join(r[0] for r in orphans)
    raise ValueError(
        f"the hit tables cover contig(s) the reference index does not name: "
        f"{names}. That is usually a FASTA and a BAM from different assemblies, "
        f"or a `.fai` for one chromosome. Writing the file would drop those "
        f"positions without saying so."
    )

written = 0
# Append, because an Amethyst file holds a group per context and a run writing
# both CG and CH reads this file twice. The context group is replaced rather
# than added to, so re-running a script rebuilds it rather than doubling it --
# but HDF5 does not reclaim the space a deleted group used, so a file rebuilt
# many times grows. Delete it first for a clean build.
with h5py.File(out, "a") as h5:
    if context in h5:
        del h5[context]
    if "metadata/version" not in h5:
        # h5py stores a Python str as a variable-length UTF-8 scalar, which is
        # what `h5/handles.py:43-56` reads back and decodes.
        h5.create_dataset("metadata/version", data=VERSION)

    for (block,), table in groups(con.execute(SQL).to_arrow_reader(8192), ["block"]):
        longest = pc.max(pc.utf8_length(table.column("chr"))).as_py() or 0
        width = max(CHR_MIN_WIDTH, longest)
        dtype = [(f, f"S{width}" if f == "chr" else k) for f, k in OBSERVATIONS_V2]

        rows = np.empty(table.num_rows, dtype=dtype)
        for column, kind in dtype:
            values = table.column(column).to_numpy(zero_copy_only=False)
            rows[column] = values.astype(kind) if column == "chr" else values

        # compression_opts=9 is what premethyst's calls2h5 uses, so a file from
        # here is the same size as a file from there.
        h5.create_dataset(
            f"{context}/{block}/{BASE_LEVEL}",
            data=rows,
            compression="gzip",
            compression_opts=9,
        )
        written += 1

print(f"{out}: {written} block(s) under /{context}")
