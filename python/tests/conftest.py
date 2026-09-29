"""A hand-written hit table, small enough to reason about row by row.

Every test below asserts against numbers you can get by reading this file, so
the data is spelled out rather than generated. The columns are the ones an
alnbase run selects with

    -F qname,is_first_in_template,ref_name,strand,mapq,qual,off_5p,off_3p

plus `XB`, a cell barcode moved into a tag upstream.

The six reads are chosen so that each filter in `lib/filters.sql` rejects
exactly one of them, which is what makes a failure point at the filter that
broke rather than at the fixture:

    rA  cellA  ordinary
    rB  cellA  two of three CH methylated -- a ratio threshold
    rC  cellB  four CH in a row, but a low overall ratio -- read_fails_consecutive
    rD  cellB  ordinary
    rE  cellC  clean, but named in the exclusion file -- read_fails_listed
    rF  cellC  ordinary, and cellC is the thin block -- a block total threshold

chr1 carries the CH calls, chr2 the CG ones. The CG site at chr2:200/201 is
covered on both strands, which is what `lib/merge_strands.sql` collapses.
"""

from __future__ import annotations

import pathlib
import sys

import duckdb
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

COLUMNS = (
    "qname VARCHAR, is_first_in_template BOOLEAN, XB VARCHAR, name VARCHAR, "
    "ref_name VARCHAR, refr_pos BIGINT, strand VARCHAR, qual UTINYINT, "
    "mapq UTINYINT, off_5p BIGINT, off_3p BIGINT"
)

# qname, mate, block, query name, contig, 0-based position, strand
ROWS = [
    # -- cellA ------------------------------------------------------------
    ("rA", True, "cellA", "THH", "chr1", 100, "+"),
    ("rA", True, "cellA", "THG", "chr1", 110, "+"),
    ("rA", True, "cellA", "CHH", "chr1", 120, "+"),
    ("rA", True, "cellA", "CG",  "chr2", 200, "+"),
    ("rA", True, "cellA", "TG",  "chr2", 210, "+"),

    # Two of three CH methylated: a read that did not convert.
    ("rB", True, "cellA", "CHH", "chr1", 100, "+"),
    ("rB", True, "cellA", "CHG", "chr1", 110, "+"),
    ("rB", True, "cellA", "THH", "chr1", 120, "+"),
    ("rB", True, "cellA", "CG",  "chr2", 200, "+"),

    # -- cellB ------------------------------------------------------------
    # Rate 4/10, below the 0.5 threshold, but four adjacent: the case a rate
    # provably cannot express.
    ("rC", True, "cellB", "CHH", "chr1", 300, "+"),
    ("rC", True, "cellB", "CHH", "chr1", 310, "+"),
    ("rC", True, "cellB", "CHG", "chr1", 320, "+"),
    ("rC", True, "cellB", "CHH", "chr1", 330, "+"),
    ("rC", True, "cellB", "THH", "chr1", 340, "+"),
    ("rC", True, "cellB", "THH", "chr1", 350, "+"),
    ("rC", True, "cellB", "THG", "chr1", 360, "+"),
    ("rC", True, "cellB", "THH", "chr1", 370, "+"),
    ("rC", True, "cellB", "THH", "chr1", 380, "+"),
    ("rC", True, "cellB", "THH", "chr1", 390, "+"),

    ("rD", True, "cellB", "THH", "chr1", 300, "+"),
    ("rD", True, "cellB", "CHH", "chr1", 310, "+"),
    ("rD", True, "cellB", "CG",  "chr2", 200, "+"),
    # The minus-strand cytosine of the same CG, one base to the right.
    ("rD", True, "cellB", "TG",  "chr2", 201, "-"),

    # -- cellC ------------------------------------------------------------
    ("rE", True, "cellC", "TG",  "chr2", 200, "+"),
    ("rE", True, "cellC", "THH", "chr1", 400, "+"),

    ("rF", True, "cellC", "CG",  "chr2", 201, "-"),
    ("rF", True, "cellC", "THH", "chr1", 400, "+"),
]

#: What every row shares. Varying them is a per-test concern, so they are
#: constant here and the hit-level filter in a run script sees clean rows.
QUAL, MAPQ, OFF_5P, OFF_3P = 35, 40, 20, 20

#: The declarations a run makes. Repeated in the test scripts rather than
#: shared, because a run script is meant to be readable on its own.
MEASURES = """
CREATE OR REPLACE VIEW measure (measure, name, category) AS VALUES
  ('mCG', 'CG',  'methylated'), ('mCG', 'TG',  'unmethylated'),
  ('mCH', 'CHG', 'methylated'), ('mCH', 'CHH', 'methylated'),
  ('mCH', 'THG', 'unmethylated'), ('mCH', 'THH', 'unmethylated');

CREATE OR REPLACE VIEW relation (relation, measure, main, other) AS VALUES
  ('mCG', 'mCG', ['methylated'], ['unmethylated']),
  ('mCH', 'mCH', ['methylated'], ['unmethylated']);
"""

#: The adapter every test script starts with.
SRC = """
CREATE OR REPLACE VIEW src AS
SELECT *,
       ref_name AS contig,
       {'qname': qname, 'mate': is_first_in_template} AS read,
       XB AS block
FROM read_parquet(getvariable('hits'));
"""


@pytest.fixture(scope="session")
def hits(tmp_path_factory) -> str:
    """The hit table, as a parquet file. Returns its path."""
    path = tmp_path_factory.mktemp("hits") / "calls_0_0.parquet"
    con = duckdb.connect()
    con.execute(f"CREATE TABLE h ({COLUMNS})")
    con.executemany(
        "INSERT INTO h VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [(*row, QUAL, MAPQ, OFF_5P, OFF_3P) for row in ROWS],
    )
    con.execute("COPY h TO ? (FORMAT parquet)", [str(path)])
    return str(path)


@pytest.fixture(scope="session")
def exclude(tmp_path_factory) -> str:
    """A one-line exclusion list naming rE."""
    path = tmp_path_factory.mktemp("lists") / "exclude.txt"
    path.write_text("rE\n")
    return str(path)


@pytest.fixture(scope="session")
def fai(tmp_path_factory) -> str:
    """A two-contig FASTA index, in reference order rather than alphabetical."""
    path = tmp_path_factory.mktemp("ref") / "genome.fa.fai"
    path.write_text("chr1\t1000\t6\t60\t61\nchr2\t1000\t1024\t60\t61\n")
    return str(path)


@pytest.fixture
def con():
    """A fresh in-memory DuckDB connection."""
    connection = duckdb.connect()
    yield connection
    connection.close()


@pytest.fixture
def write(tmp_path):
    """Write a script into the test's own directory and return its path."""

    def _write(text: str, name: str = "run.sql") -> pathlib.Path:
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    return _write
