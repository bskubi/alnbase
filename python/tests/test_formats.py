"""The output formats, checked as the bytes they actually write.

Each test asserts the file's content rather than a relation, because the
things these formats get wrong -- a coordinate off by one, a percentage
rounded the wrong way, a column in the wrong order -- are only visible in the
file. The expected values are worked out by hand from the fixture in conftest.

After the read filters in `resolved` below, two reads survive: rA and rD. The
CG site at chr2:200/201 is covered by rA (plus strand) and rD (both strands),
so merged it is 2 methylated of 3 calls.
"""

from __future__ import annotations

import pytest

from alnbase_agg import script as scriptlib
from conftest import MEASURES, SRC

# The run every format test shares: source, measures, filters, one read filter
# and one block filter, then the site counts.
PIPELINE = (
    SRC
    + MEASURES
    + """
.read lib/filters.sql

CREATE OR REPLACE VIEW resolved AS SELECT * FROM src
WHERE read NOT IN (
        SELECT read FROM read_fails('mCH', 'ratio', '>', 0.5)
  UNION ALL SELECT read FROM read_fails_consecutive('mCH', '>', 3))
  AND block NOT IN (SELECT block FROM block_fails('mCH', 'total', '<', 3));

.read lib/count_sites.sql
"""
)

MERGE = """
SET VARIABLE merge = 'mCG';
.read lib/merge_strands.sql
"""


def run(con, write, hits, body: str, **values) -> None:
    path = write(PIPELINE + body)
    scriptlib.run(con, path, values={"hits": hits, **values})


def lines(path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text().splitlines()]


# --- the site table -------------------------------------------------------


def test_sites_are_per_cytosine_per_block_before_merging(con, write, hits):
    run(con, write, hits, "")
    rows = con.execute(
        "SELECT block, contig, refr_pos, strand, num, den FROM site_relation "
        "WHERE relation = 'mCG' ORDER BY block, refr_pos"
    ).fetchall()
    assert rows == [
        ("cellA", "chr2", 200, "+", 1, 1),
        ("cellA", "chr2", 210, "+", 0, 1),
        ("cellB", "chr2", 200, "+", 1, 1),
        ("cellB", "chr2", 201, "-", 0, 1),
    ]


def test_merge_strands_moves_the_minus_strand_cytosine_onto_its_partner(con, write, hits):
    run(con, write, hits, MERGE)
    rows = con.execute(
        "SELECT block, refr_pos, strand, num, den FROM site_relation "
        "WHERE relation = 'mCG' ORDER BY block, refr_pos"
    ).fetchall()
    assert rows == [
        ("cellA", 200, "+", 1, 1),
        ("cellA", 210, "+", 0, 1),
        ("cellB", 200, "+", 1, 2),
    ]


def test_merge_leaves_a_non_palindromic_measure_alone(con, write, hits):
    """Adding two unrelated cytosines together is exactly what must not happen."""
    run(con, write, hits, MERGE)
    before = con.execute(
        "SELECT count(*) FROM site_stranded WHERE measure = 'mCH'"
    ).fetchone()
    after = con.execute(
        "SELECT count(*) FROM site_relation WHERE relation = 'mCH'"
    ).fetchone()
    assert before == after


# --- bed-like -------------------------------------------------------------


def test_bismark_cov(con, write, hits, tmp_path):
    out = tmp_path / "s.cov"
    run(
        con, write, hits,
        MERGE + "SET VARIABLE relation = 'mCG';\n.read formats/bismark_cov.sql\n",
        out=str(out),
    )
    assert lines(out) == [
        # 1-based, start == end, percent as Perl stringifies it.
        ["chr2", "201", "201", "66.6666666666667", "2", "1"],
        ["chr2", "211", "211", "0", "0", "1"],
    ]


def test_bismark_cov_min_depth(con, write, hits, tmp_path):
    out = tmp_path / "s.cov"
    run(
        con, write, hits,
        MERGE + "SET VARIABLE relation = 'mCG';\n.read formats/bismark_cov.sql\n",
        out=str(out), min_depth=2,
    )
    assert [row[1] for row in lines(out)] == ["201"]


def test_bedgraph_is_zero_based_and_half_open(con, write, hits, tmp_path):
    out = tmp_path / "s.bedGraph"
    run(
        con, write, hits,
        MERGE + "SET VARIABLE relation = 'mCG';\n.read formats/bedgraph.sql\n",
        out=str(out),
    )
    assert lines(out) == [
        ["chr2", "200", "201", "66.6666666666667"],
        ["chr2", "210", "211", "0"],
    ]


def test_methyldackel_truncates_the_percentage(con, write, hits, tmp_path):
    """Two thirds is 66, not 67: MethylDackel's percentage is a C cast."""
    out = tmp_path / "s.bedGraph"
    run(
        con, write, hits,
        MERGE + "SET VARIABLE relation = 'mCG';\n"
                ".read formats/methyldackel_bedgraph.sql\n",
        out=str(out),
    )
    assert lines(out) == [
        ["chr2", "200", "201", "66", "2", "1"],
        ["chr2", "210", "211", "0", "0", "1"],
    ]


def test_bedmethyl_has_eleven_columns_and_keeps_the_strand(con, write, hits, tmp_path):
    out = tmp_path / "s.bedmethyl"
    run(
        con, write, hits,
        "SET VARIABLE relation = 'mCG';\n.read formats/bedmethyl.sql\n",
        out=str(out),
    )
    rows = lines(out)
    assert all(len(row) == 11 for row in rows)
    assert rows[0] == ["chr2", "200", "201", "mCG", "2", "+",
                       "200", "201", "0,0,0", "2", "100.00"]
    # Not merged, so the minus-strand cytosine is its own line.
    assert [row[5] for row in rows] == ["+", "-", "+"]


# --- CGmap ----------------------------------------------------------------

TRINUC_MEASURES = """
CREATE OR REPLACE VIEW measure (measure, name, category) AS VALUES
  ('CGA', 'CG', 'methylated'), ('CGA', 'TG', 'unmethylated');

CREATE OR REPLACE VIEW relation (relation, measure, main, other) AS VALUES
  ('CGA', 'CGA', ['methylated'], ['unmethylated']);
"""


def test_cgmap_reads_its_context_columns_off_the_relation_name(con, write, hits, tmp_path):
    out = tmp_path / "s.CGmap"
    path = write(
        SRC + TRINUC_MEASURES
        + "\n.read lib/filters.sql\n"
        "CREATE OR REPLACE VIEW resolved AS SELECT * FROM src;\n"
        ".read lib/count_sites.sql\n"
        ".read formats/cgmap.sql\n"
    )
    scriptlib.run(con, path, values={"hits": hits, "out": str(out)})
    # Not merged: each cytosine is its own line, at its own position, and
    # column 2 is the base on the Watson strand rather than the read's.
    assert lines(out) == [
        ["chr2", "C", "201", "CG", "CG", "0.75", "3", "4"],
        ["chr2", "G", "202", "CG", "CG", "0.5", "1", "2"],
        ["chr2", "C", "211", "CG", "CG", "0.0", "0", "1"],
    ]


def test_cgmap_refuses_a_relation_that_is_not_a_trinucleotide(con, write, hits):
    """An empty CGmap is the failure that looks like success."""
    import duckdb

    path = write(
        SRC + MEASURES
        + "\n.read lib/filters.sql\n"
        "CREATE OR REPLACE VIEW resolved AS SELECT * FROM src;\n"
        ".read lib/count_sites.sql\n"
        ".read formats/cgmap.sql\n"
    )
    with pytest.raises(Exception) as exc:
        scriptlib.run(con, path, values={"hits": hits, "out": "/dev/null"})
    assert "trinucleotide" in str(exc.value)


def test_cgmap_level_keeps_pythons_trailing_zero(con, write, hits, tmp_path):
    """`str(round(1.0, 2))` is '1.0', not '1'."""
    out = tmp_path / "s.CGmap"
    path = write(
        SRC + TRINUC_MEASURES
        + "\n.read lib/filters.sql\n"
        "CREATE OR REPLACE VIEW resolved AS SELECT * FROM src WHERE name = 'CG';\n"
        ".read lib/count_sites.sql\n"
        ".read formats/cgmap.sql\n"
    )
    scriptlib.run(con, path, values={"hits": hits, "out": str(out)})
    assert {row[5] for row in lines(out)} == {"1.0"}


# --- per block -------------------------------------------------------------


def test_cellinfo(con, write, hits, tmp_path):
    out = tmp_path / "s.cellInfo.txt"
    run(
        con, write, hits,
        MERGE + "SET VARIABLE cg_relation = 'mCG';\n"
                "SET VARIABLE ch_relation = 'mCH';\n"
                ".read formats/cellinfo.sql\n",
        out=str(out),
    )
    assert lines(out) == [
        # block, all calls, CG sites, CG %, CH sites, CH %
        ["cellA", "5", "2", "50.00", "3", "33.33"],
        ["cellB", "4", "1", "50.00", "2", "50.00"],
    ]


def test_block_summary_counts_positions_not_calls(con, write, hits):
    run(con, write, hits, MERGE)
    assert con.execute(
        "SELECT den, sites FROM block_summary "
        "WHERE block = 'cellB' AND relation = 'mCG'"
    ).fetchone() == (2, 1)


# --- parquet --------------------------------------------------------------


def test_sites_parquet_round_trips(con, write, hits, tmp_path):
    out = tmp_path / "s.parquet"
    run(con, write, hits, MERGE + ".read formats/sites_parquet.sql\n", out=str(out))
    columns = [
        r[0] for r in con.execute(f"DESCRIBE SELECT * FROM '{out}'").fetchall()
    ]
    assert columns == ["block", "contig", "refr_pos", "strand", "relation",
                       "num", "den"]


# --- Amethyst -------------------------------------------------------------


@pytest.fixture
def h5py_():
    return pytest.importorskip("h5py")


def test_amethyst_h5(con, write, hits, fai, tmp_path, h5py_):
    pytest.importorskip("pyarrow")
    out = tmp_path / "s.h5"
    run(
        con, write, hits,
        MERGE + ".read lib/load_contigs.sql\n"
                "SET VARIABLE relation = 'mCG';\n"
                "SET VARIABLE context = 'CG';\n"
                ".read formats/amethyst_h5.sql\n"
                ".python formats/amethyst_h5.py\n",
        out=str(out), fai=fai,
    )
    with h5py_.File(out, "r") as h5:
        assert h5["metadata/version"][()].decode() == "amethyst2.0.0"
        assert sorted(h5["CG"].keys()) == ["cellA", "cellB"]
        rows = h5["CG/cellB/1"][:]
        assert list(rows.dtype.names) == ["chr", "pos", "c", "t"]
        # cellB has one merged CG: 1 methylated of 2 calls, at 1-based 201.
        assert rows.tolist() == [(b"chr2", 201, 1, 1)]


def test_amethyst_h5_holds_a_group_per_context(con, write, hits, fai, tmp_path, h5py_):
    pytest.importorskip("pyarrow")
    out = tmp_path / "s.h5"
    run(
        con, write, hits,
        MERGE + ".read lib/load_contigs.sql\n"
                "SET VARIABLE context = 'CG';\n"
                "SET VARIABLE relation = 'mCG';\n"
                ".read formats/amethyst_h5.sql\n"
                ".python formats/amethyst_h5.py\n"
                "SET VARIABLE context = 'CH';\n"
                "SET VARIABLE relation = 'mCH';\n"
                ".read formats/amethyst_h5.sql\n"
                ".python formats/amethyst_h5.py\n",
        out=str(out), fai=fai,
    )
    with h5py_.File(out, "r") as h5:
        assert sorted(k for k in h5.keys() if k != "metadata") == ["CG", "CH"]


def test_amethyst_refuses_a_contig_the_reference_index_lacks(
    con, write, hits, tmp_path, h5py_
):
    """A FASTA and a BAM from different assemblies, caught before writing."""
    pytest.importorskip("pyarrow")
    partial = tmp_path / "one.fai"
    partial.write_text("chr1\t1000\t6\t60\t61\n")
    with pytest.raises(Exception) as exc:
        run(
            con, write, hits,
            MERGE + ".read lib/load_contigs.sql\n"
                    "SET VARIABLE relation = 'mCG';\n"
                    "SET VARIABLE context = 'CG';\n"
                    ".read formats/amethyst_h5.sql\n"
                    ".python formats/amethyst_h5.py\n",
            out=str(tmp_path / "s.h5"), fai=str(partial),
        )
    assert "chr2" in str(exc.value)
