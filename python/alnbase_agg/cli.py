"""The `alnbase-agg` command.

Four subcommands, and the first one is most of the tool:

    alnbase-agg run   script.sql -s hits='calls_*.parquet'   # do it
    alnbase-agg show  script.sql                             # the SQL, expanded
    alnbase-agg lib                                          # where the library is
    alnbase-agg check --hits 'calls_*.parquet' --reference genome.fa

`run` executes a script. `show` prints the same script with every `.read`
substituted, which is the query a run will actually issue -- something to read
before running, to paste into a DuckDB shell, or to attach to a methods
section. `lib` prints the directory the shipped SQL lives in, because forking a
format is `cp` and then editing, and you have to know what to copy. `check`
does not involve a script at all: it reads hit tables and a FASTA and reports
whether the coordinates are where they claim to be.
"""

from __future__ import annotations

import sys
from pathlib import Path

import click
import duckdb

from . import script as scriptlib
from .script import SQL_DIR, ScriptError


def _values(settings) -> dict:
    """`-s name=value` pairs as a dict, with each value typed."""
    out = {}
    for item in settings:
        name, _, text = item.partition("=")
        if not _:
            raise click.BadParameter(
                f"{item!r} is not name=value. Example: -s min_depth=5", param_hint="-s"
            )
        out[name.strip()] = scriptlib.parse_value(text)
    return out


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(package_name="alnbase-agg")
def main() -> None:
    """Aggregate alnbase hit tables with DuckDB. A run is a SQL script."""


@main.command()
@click.argument("script", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("-s", "--set", "settings", multiple=True, metavar="NAME=VALUE",
              help="Set a session variable the script reads with getvariable().")
@click.option("-I", "--include", multiple=True, metavar="DIR",
              type=click.Path(exists=True, file_okay=False, path_type=Path),
              help="Another directory to resolve `.read` against. Repeatable.")
@click.option("--db", metavar="PATH",
              help="Run against a persistent database instead of memory, so the "
                   "views survive the run and can be queried afterwards.")
@click.option("-v", "--verbose", is_flag=True,
              help="Print every directive followed and every statement executed.")
@click.option("-q", "--quiet", is_flag=True,
              help="Do not print the result of statements that return rows.")
def run(script: Path, settings, include, db, verbose, quiet) -> None:
    """Execute a run script."""
    con = duckdb.connect(db or ":memory:")
    # The trace goes to stderr and results to stdout, so a script whose QC
    # queries are the point can be redirected without the trace in the file.
    trace = (lambda line: click.echo(line, err=True)) if verbose else None
    show = None if quiet else (lambda text: click.echo(text + "\n"))
    try:
        scriptlib.run(
            con,
            script,
            include=(*include, SQL_DIR),
            values=_values(settings),
            trace=trace,
            show=show,
        )
    except ScriptError as exc:
        raise SystemExit(f"alnbase-agg: {exc}")


@main.command()
@click.argument("script", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("-I", "--include", multiple=True, metavar="DIR",
              type=click.Path(exists=True, file_okay=False, path_type=Path))
def show(script: Path, include) -> None:
    """Print the script with every `.read` substituted."""
    try:
        click.echo(scriptlib.expand(script, include=(*include, SQL_DIR)))
    except ScriptError as exc:
        raise SystemExit(f"alnbase-agg: {exc}")


@main.command(name="lib")
@click.option("--list", "listing", is_flag=True, help="List the files, not the path.")
def lib_path(listing) -> None:
    """Print where the shipped SQL library lives.

    Copy a file out of it, edit the copy, and put it next to your script: a
    `.read` resolves against the script's own directory first, so your version
    wins with no flag and no registry.
    """
    if not listing:
        click.echo(SQL_DIR)
        return
    for path in sorted(SQL_DIR.rglob("*")):
        if path.is_file():
            click.echo(path.relative_to(SQL_DIR))


@main.command()
@click.option("--hits", required=True, metavar="GLOB",
              help="The hit tables, e.g. 'calls_*_*.parquet'.")
@click.option("--reference", required=True, metavar="FASTA",
              type=click.Path(exists=True, dir_okay=False, path_type=Path),
              help="The FASTA the alignment was made against.")
@click.option("--sites", default=2000, show_default=True,
              help="How many positions to sample.")
@click.option("--strict/--no-strict", default=True, show_default=True,
              help="Exit non-zero when the check fails, rather than warning.")
def check(hits, reference, sites, strict) -> None:
    """Check that hit-table positions match the reference genome.

    A coordinate off by one is the one error here that produces output which
    looks perfect and is wrong at every base: every file writes, every count is
    plausible, and nothing downstream knows where the cytosines were supposed
    to be. This looks the reference base up in the FASTA and scores the hits at
    several candidate shifts, so the report says which shift, not merely that
    something is wrong.
    """
    import logging

    from .manifest import Manifest, ManifestError
    from .verify import VerifyError, verify_reference

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    con = duckdb.connect()
    try:
        # The manifest is read first and passed in, so the shift the check
        # applies is the run's own declared coordinate base rather than a
        # number this command worked out for itself.
        manifest = Manifest.read(con, hits)
        report = verify_reference(
            con, manifest, hits, reference, sites=sites, strict=strict
        )
    except (ManifestError, VerifyError) as exc:
        raise SystemExit(f"alnbase-agg: {exc}")
    click.echo(report)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
