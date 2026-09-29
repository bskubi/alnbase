"""What the hit table says about itself.

An aggregation must not guess at four things, and the alnbase run recorded all
four, so this module reads them rather than assuming:

- **the coordinate base**, so the shift to a 1-based output is derived instead of
  written as a literal `+ 1` somewhere;
- **the query names**, so a script naming a query the run never declared fails
  while the user is still looking at the script;
- **the contigs and their lengths**, for a reference check;
- **the alnbase version and run id**, so an output can say what produced it.

There are two copies of the manifest and they are not interchangeable. Every
parquet file carries one in its footer, written before the first row: it
describes the run's *settings*. A `{stem}.manifest.json` beside the output is
written last and only on success: its presence is the *completeness* marker.

So the footer is what we read (it is always there, and a glob gives us a file to
read it from) and the sidecar is what we check for (its absence means the run may
have been interrupted, which is worth saying). Absence is a warning rather than
an error: a hit table can be perfectly usable without a finished run's blessing,
and refusing to look at one would be a worse failure than an honest note.
"""

from __future__ import annotations

import json
import tomllib
import warnings
from dataclasses import dataclass, field
from pathlib import Path

# The footer key alnbase writes its manifest under; `src/manifest.rs`.
MANIFEST_KEY = "alnbase_manifest"


class ManifestError(Exception):
    """The hit table could not be described, phrased for the person who ran it."""


@dataclass(frozen=True)
class Contig:
    """One `@SQ` line from the BAM header, as the run recorded it.

    `length` is the header's length, not a FASTA's. That is what lets
    `alnbase_agg.verify.check_contigs` find a mismatched assembly without reading
    a base. It
    is 0 when the header declared no length, and a comparison against 0 proves
    nothing, so callers test for it before using it.

    The other two fields describe the reference the *run itself* read, and both
    are None for a run that read none. `in_reference` says whether that reference
    had the contig. `md5` is the reference's digest for it, falling back to the
    header's own `M5` when there was no reference or it lacked the contig -- so
    two `md5` values are comparable only when both came from the same source.
    Nothing in this package reads `md5` yet; it is carried because the run
    recorded it and an error message may want to name it.
    """

    name: str
    length: int
    md5: str | None = None
    in_reference: bool | None = None


@dataclass(frozen=True)
class Manifest:
    """One alnbase run's settings, as its output records them."""

    alnbase_version: str
    format_version: str
    coordinate_base: int
    run_id: str
    command: str
    partition_by: tuple[str, ...]
    fields: tuple[str, ...]
    contigs: tuple[Contig, ...]
    query_names: frozenset[str]
    # The reference index the run itself read, as given on its command line; None
    # for a run that read none. Reported when a reference check fails, because
    # "the run read mm10" is usually the whole explanation.
    reference: str | None = None
    # The file the footer was read from, for error messages that have to name one.
    source: Path | None = None
    complete: bool = True
    raw: dict = field(default_factory=dict, repr=False)

    # -- construction -------------------------------------------------------

    @classmethod
    def from_json(cls, doc: dict, source: Path | None = None, complete: bool = True) -> "Manifest":
        """Build a manifest from either copy of the JSON, footer or sidecar.

        The two copies are not the same document -- a finished run's sidecar adds
        a `result` block counting the records it read -- but they agree on every
        field read here, so one reader serves both. `source` is the file to name
        in errors, and `complete` records whether a sidecar was the source.

        A missing or malformed required field raises `ManifestError` rather than
        producing a manifest with a hole in it. That is the point of the class:
        an aggregation that cannot read the coordinate base must stop, because
        the alternative is to assume one and emit an output that is silently a
        base out. Fields the format really does leave optional -- `command`,
        `contigs`, `query_files`, `reference` -- are read with defaults, and the
        methods that use them say what they do when they are absent.
        """
        try:
            output = doc["output"]
            return cls(
                alnbase_version=doc["alnbase_version"],
                format_version=doc["format_version"],
                coordinate_base=int(doc["coordinate_base"]),
                run_id=doc["run_id"],
                command=doc.get("command", "query"),
                partition_by=tuple(output.get("partition_by") or ()),
                fields=tuple(output.get("fields") or ()),
                contigs=tuple(
                    Contig(
                        name=c["name"],
                        length=int(c["length"]),
                        md5=c.get("md5"),
                        in_reference=c.get("in_reference"),
                    )
                    for c in doc.get("contigs") or ()
                ),
                query_names=_query_names(doc.get("query_files") or ()),
                reference=doc.get("reference"),
                source=source,
                complete=complete,
                raw=doc,
            )
        except (KeyError, TypeError, ValueError) as exc:
            where = f"{source}: " if source else ""
            raise ManifestError(
                f"{where}this does not look like an alnbase manifest ({exc})"
            ) from None

    @classmethod
    def read(cls, con, hits: str) -> "Manifest":
        """Read the manifest for a hit table given as a path or a glob.

        `con` is a DuckDB connection, used because it already knows how to expand
        a glob and read a parquet footer -- which keeps this package's
        dependencies to DuckDB alone rather than adding pyarrow for one lookup.
        """
        files = expand(con, hits)
        footer = _footer_manifest(con, files[0])
        sidecar = _sidecar(files)
        complete = sidecar is not None
        if not complete:
            warnings.warn(
                f"no manifest file beside {files[0].name}; alnbase writes "
                f"{{stem}}.manifest.json only when a run finishes, so this hit "
                f"table may be from an interrupted run. Reading the settings from "
                f"the parquet footer instead.",
                stacklevel=2,
            )
        doc = sidecar if sidecar is not None else footer
        return cls.from_json(doc, source=files[0], complete=complete)

    # -- what callers ask it ------------------------------------------------

    @property
    def contig_lengths(self) -> dict[str, int]:
        """Contig name to header length, for lookup by name.

        The reference check walks this against a FASTA index. A length of 0 means
        the header declared none; see `Contig`.
        """
        return {c.name: c.length for c in self.contigs}

    def shift(self, base: int = 1) -> int:
        """How far a position moves when reported in `base`.

        The number itself, for a caller that needs it rather than the SQL: the
        reference check looks up the position an aggregation would report, and
        working the offset out a second time is how the check and the query come
        to disagree about which one of them is right.
        """
        return base - self.coordinate_base

    def position_sql(self, column: str = "refr_pos", base: int = 1) -> str:
        """`refr_pos` shifted from the run's coordinate base to `base`.

        The whole reason this is a method: the shift appears once, derived from
        what the run declared, so no output can be one base out because someone
        wrote the constant by hand. A shift of zero is emitted as the bare column
        so that generated SQL does not carry a pointless `+ 0`.
        """
        delta = self.shift(base)
        if delta == 0:
            return column
        return f"{column} + {delta}" if delta > 0 else f"{column} - {-delta}"

    def unknown_queries(self, names) -> list[str]:
        """Which of `names` the run never declared; empty when all are known.

        Returns nothing when the manifest listed no query files at all, because
        then we cannot tell an unknown name from an unrecorded one, and inventing
        an error out of missing evidence is worse than letting DuckDB count zero.
        """
        if not self.query_names:
            return []
        return sorted(set(names) - self.query_names)


def expand(con, hits: str) -> list[Path]:
    """The files a path or glob names, in sorted order.

    The order is load-bearing, which is why it is sorted here rather than left to
    the directory listing. `Manifest.read` reads the footer from the first file,
    and the reference check reads its column list from the first file and then
    samples shards at even intervals along the list. An unstable order would make
    both of those depend on which files the filesystem happened to return first,
    so two runs over one glob could check different shards and disagree.

    No match is an error, not an empty list. A caller given an empty list would
    aggregate nothing and report a successful run of zero rows, which looks the
    same as a real result and is the failure a mistyped glob actually produces.
    """
    rows = con.execute(
        "SELECT file FROM glob(?) ORDER BY file", [hits]
    ).fetchall()
    files = [Path(r[0]) for r in rows]
    if not files:
        raise ManifestError(
            f"no files match {hits!r}. A sharded alnbase run writes "
            f"{{stem}}_{{worker}}_{{shard}}.parquet, so the glob usually looks "
            f"like 'calls_*_*.parquet'"
        )
    return files


def _footer_manifest(con, path: Path) -> dict:
    """The manifest JSON stored in one parquet file's footer metadata.

    alnbase writes this before the first row, so a file carries it even when the
    run was interrupted. That is what makes it the copy this package reads: the
    sidecar is the fuller document, but it exists only after a run succeeds.

    `decode(key)` is in the query because DuckDB hands back both the key and the
    value of parquet key-value metadata as BLOB. Comparing the key to a string
    without decoding it matches no rows, and the failure would look exactly like
    a file written by some other program.
    """
    rows = con.execute(
        "SELECT value FROM parquet_kv_metadata(?) WHERE decode(key) = ?",
        [str(path), MANIFEST_KEY],
    ).fetchall()
    if not rows:
        raise ManifestError(
            f"{path} carries no {MANIFEST_KEY!r} in its footer, so it was not "
            f"written by alnbase -- or was written by a version predating the "
            f"manifest. There is nothing to say what its coordinates mean, which "
            f"is the one thing this pipeline cannot guess."
        )
    value = rows[0][0]
    text = value.decode("utf-8") if isinstance(value, (bytes, bytearray)) else str(value)
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise ManifestError(f"{path}: the footer manifest is not valid JSON ({exc})") from None


def _sidecar(files: list[Path]) -> dict | None:
    """The `{stem}.manifest.json` beside a run's files, if it is there.

    A sharded run names its files `{stem}_{worker}_{shard}.parquet`, so the stem
    is the name with those two numeric indices taken off, and an unsharded file
    keeps its own stem. The first candidate that exists wins; the files are tried
    in order only because a glob may reach across directories.

    Nothing here checks that the sidecar belongs to the run whose footer was
    read. Both documents carry a `run_id` and comparing them is the check that
    would settle it, but neither this function nor `Manifest.read` does so. A
    sidecar left behind by an earlier run over the same output prefix would
    therefore be read in preference to the footer, and its settings believed.
    Treat that as a known gap rather than a deliberate choice.
    """
    for path in files:
        stem = path.name.split(".")[0]
        parts = stem.rsplit("_", 2)
        if len(parts) == 3 and all(p.isdigit() for p in parts[1:]):
            stem = parts[0]
        candidate = path.parent / f"{stem}.manifest.json"
        if candidate.is_file():
            try:
                return json.loads(candidate.read_text())
            except json.JSONDecodeError as exc:
                raise ManifestError(f"{candidate}: not valid JSON ({exc})") from None
    return None


def _query_names(query_files) -> frozenset[str]:
    """The `[query.NAME]` keys declared across the run's query files.

    The manifest stores each query file's whole text rather than a list of
    names, so the names come from parsing it -- the same TOML alnbase parsed. A
    file that will not parse is skipped rather than fatal: we are collecting
    names to give a better error message, and failing to do so must not become an
    error of its own.
    """
    names: set[str] = set()
    for entry in query_files:
        text = entry.get("text") if isinstance(entry, dict) else None
        if not text:
            continue
        try:
            doc = tomllib.loads(text)
        except tomllib.TOMLDecodeError:
            continue
        names.update((doc.get("query") or {}).keys())
    return frozenset(names)
