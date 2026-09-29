"""alnbase-agg: aggregate alnbase hit tables with DuckDB.

There is no recipe format to learn. A run is a SQL script, and this package is
what runs one: a directive for including another file, a way to pass values in,
a library of SQL to include, and one helper for the formats that are not a
single table.

    alnbase-agg run block.sql -s hits='calls_*_*.parquet'

The pieces are importable, so an analysis can stop at a relation rather than a
file:

    import duckdb, alnbase_agg as aa

    con = duckdb.connect()
    aa.run_script(con, "run.sql", values={"hits": "calls_*_*.parquet"})
    df = con.execute("SELECT * FROM site_relation").df()

`run_script` is the whole entry point. `groups` is what a format uses when its
output is not one file -- an HDF5 with a dataset per block, a folder with a file
per block -- and `verify_reference` is the coordinate check, which needs no
script at all.
"""

from .manifest import Contig, Manifest, ManifestError  # noqa: F401
from .script import SQL_DIR, ScriptError, expand, parse_value  # noqa: F401
from .script import run as run_script  # noqa: F401
from .stream import NotGrouped, check_sorted, groups  # noqa: F401

__all__ = [
    "Contig",
    "Manifest",
    "ManifestError",
    "NotGrouped",
    "SQL_DIR",
    "ScriptError",
    "check_sorted",
    "expand",
    "groups",
    "parse_value",
    "run_script",
]


def __getattr__(name):
    """Expose the reference check without importing pyfaidx on every import.

    `verify_reference` is the only thing here with a dependency a user might
    not have installed, and the common case -- running a script that writes a
    `.cov` file -- never touches it.
    """
    if name in ("verify_reference", "VerifyError"):
        from . import verify

        return getattr(verify, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
