"""The read- and block-level filters in `lib/filters.sql`, against the six
reads in conftest. The site-level ones are exercised in `test_pipeline.py`.

Every assertion below is a read name, because the fixture is built so that
each filter rejects exactly one thing. A test failing with the wrong read name
says which filter broke; a test failing with an empty set usually means the
measure declarations changed.
"""

from __future__ import annotations

import pytest

from alnbase_agg import script as scriptlib
from alnbase_agg.script import ScriptError
from conftest import MEASURES, SRC


@pytest.fixture
def loaded(con, hits, write):
    """A connection with `src`, `measure` and the filter library in place."""
    path = write(SRC + MEASURES + "\n.read lib/filters.sql\n")
    scriptlib.run(con, path, values={"hits": hits})
    return con


def failing(con, call: str) -> set:
    """The read names a filter rejects."""
    rows = con.execute(f"SELECT DISTINCT read.qname FROM {call}").fetchall()
    return {r[0] for r in rows}


def failing_blocks(con, call: str) -> set:
    return {r[0] for r in con.execute(f"SELECT DISTINCT block FROM {call}").fetchall()}


# --- the aggregates -------------------------------------------------------


def test_a_numerator_category_counts_toward_the_denominator_too(loaded):
    """rB is two methylated CH of three CH seen. The numerator's category is
    named in `den` as well, so a call raises both counts and a ratio is a
    proportion rather than an odds ratio."""
    row = loaded.execute(
        "SELECT num, den FROM read_relation "
        "WHERE read.qname = 'rB' AND relation = 'mCH'"
    ).fetchone()
    assert row == (2, 3)


def test_a_category_with_no_hits_counts_as_zero(loaded):
    """rE covers one CH position and methylates none of it: a ratio of 0, not
    no ratio. A sum over an empty FILTER is NULL without the coalesce."""
    assert loaded.execute(
        "SELECT num, den FROM read_relation "
        "WHERE read.qname = 'rE' AND relation = 'mCH'"
    ).fetchone() == (0, 1)


def test_block_relation_reads_src_not_resolved(loaded):
    """A block is judged on everything sequenced from it, filters or no."""
    assert loaded.execute(
        "SELECT den FROM block_relation WHERE block = 'cellB' AND relation = 'mCH'"
    ).fetchone() == (12,)


# --- read-level -----------------------------------------------------------


def test_read_fails_on_a_ratio(loaded):
    assert failing(loaded, "read_fails('mCH', 'ratio', '>', 0.5)") == {"rB"}


def test_a_greater_than_is_strict(loaded):
    """rD is exactly 1/2, and a threshold that rejects it would be a surprise."""
    assert "rD" not in failing(loaded, "read_fails('mCH', 'ratio', '>', 0.5)")


def test_the_other_direction_selects_rather_than_rejects(loaded):
    """The hyper-editing shape: everything below the threshold is discarded."""
    assert failing(loaded, "read_fails('mCH', 'ratio', '<', 0.4)") == {"rA", "rE", "rF"}


def test_main_counts_one_side(loaded):
    assert failing(loaded, "read_fails('mCH', 'main', '>', 1)") == {"rB", "rC"}


def test_read_fails_on_a_total(loaded):
    assert failing(loaded, "read_fails('mCH', 'total', '<', 2)") == {"rE", "rF"}


def test_fails_missing_catches_what_no_threshold_can_reach(loaded):
    """rC covers no CG at all, so no CG threshold would ever see it."""
    assert failing(loaded, "read_fails_missing('mCG')") == {"rC"}


def test_fails_consecutive_finds_what_a_ratio_cannot(loaded):
    """rC's overall ratio is 0.4 and four of its calls are adjacent."""
    assert failing(loaded, "read_fails_consecutive('mCH', '>', 3)") == {"rC"}
    assert "rC" not in failing(loaded, "read_fails('mCH', 'ratio', '>', 0.5)")


def test_fails_consecutive_is_not_fooled_by_a_scattered_read(loaded):
    """rB is 2 of 3 methylated but its longest run is 2."""
    assert "rB" not in failing(loaded, "read_fails_consecutive('mCH', '>', 3)")


def test_read_fails_listed(loaded, exclude):
    assert failing(loaded, f"read_fails_listed('{exclude}')") == {"rE"}


# --- block-level -----------------------------------------------------------


def test_block_fails_on_a_total(loaded):
    assert failing_blocks(loaded, "block_fails('mCH', 'total', '<', 3)") == {"cellC"}


def test_block_fails_on_a_ratio(loaded):
    assert failing_blocks(loaded, "block_fails('mCH', 'ratio', '>', 0.45)") == {"cellA"}


# --- the guard ------------------------------------------------------------


def test_a_mistyped_measure_is_an_error_not_an_empty_result(loaded):
    """The silent failure this design is most exposed to, made loud."""
    with pytest.raises(duckdb_error()) as exc:
        loaded.execute("SELECT * FROM read_fails('mCG_typo', 'ratio', '>', 0.5)")
    assert "mCG_typo" in str(exc.value)


def test_the_guard_is_derived_from_the_run_s_own_declarations(loaded):
    """No list of valid names is written anywhere; the ENUM is the measures."""
    values = loaded.execute(
        "SELECT unnest(enum_range(NULL::measure_name))"
    ).fetchall()
    assert {v[0] for v in values} == {"mCG", "mCH"}


def duckdb_error():
    import duckdb

    return duckdb.Error


# --- composition ----------------------------------------------------------


def test_every_level_composes_in_one_where(con, hits, exclude, write):
    """The whole point of filters returning failures rather than passes."""
    path = write(
        SRC + MEASURES + "\n.read lib/filters.sql\n"
        "CREATE OR REPLACE VIEW resolved AS SELECT * FROM src\n"
        "WHERE read NOT IN (\n"
        "        SELECT read FROM read_fails('mCH', 'ratio', '>', 0.5)\n"
        "  UNION ALL SELECT read FROM read_fails_consecutive('mCH', '>', 3)\n"
        f"  UNION ALL SELECT read FROM read_fails_listed('{exclude}'))\n"
        "  AND block NOT IN (SELECT block FROM block_fails('mCH', 'total', '<', 3));\n"
    )
    scriptlib.run(con, path, values={"hits": hits})
    survivors = {
        r[0] for r in con.execute("SELECT DISTINCT read.qname FROM resolved").fetchall()
    }
    # rB rate, rC run, rE listed, rF in the shallow block.
    assert survivors == {"rA", "rD"}


def test_a_filter_naming_a_view_that_does_not_exist_fails_at_load(con, write, hits):
    """A table macro binds its tables when created, so the order is checked."""
    path = write(SRC + "\n.read lib/filters.sql\n")
    with pytest.raises(ScriptError):
        scriptlib.run(con, path, values={"hits": hits})
