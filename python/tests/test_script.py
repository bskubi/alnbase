"""The runner: include resolution, variables, and what a failure says."""

from __future__ import annotations

import duckdb
import pytest

from alnbase_agg import script as scriptlib
from alnbase_agg.script import SQL_DIR, ScriptError, parse_value


# --- values ---------------------------------------------------------------


@pytest.mark.parametrize(
    "text, expected",
    [
        ("5", 5),
        ("-2", -2),
        ("0.4", 0.4),
        ("true", True),
        ("false", False),
        ("calls_*.parquet", "calls_*.parquet"),
        # A barcode that looks like a number keeps its leading zero, and
        # `1e3` keeps its spelling: converting and printing back has to give
        # the original text, or the text wins.
        ("01", "01"),
        ("1e3", "1e3"),
        ("", ""),
    ],
)
def test_parse_value(text, expected):
    assert parse_value(text) == expected
    assert type(parse_value(text)) is type(expected)


def test_variables_are_typed_where_it_matters(con):
    """An integer must arrive as an integer: a threshold is compared, not read."""
    scriptlib.set_variables(con, {"min_depth": parse_value("5")})
    assert con.execute("SELECT typeof(getvariable('min_depth'))").fetchone()[0] \
        in ("INTEGER", "BIGINT")
    assert con.execute("SELECT 6 >= getvariable('min_depth')").fetchone()[0] is True


def test_unset_variable_is_null_not_an_error(con):
    """Which is what lets a format file give an option a default."""
    assert con.execute("SELECT getvariable('never_set')").fetchone()[0] is None
    assert con.execute(
        "SELECT coalesce(getvariable('never_set'), 1)"
    ).fetchone()[0] == 1


def test_refuses_a_name_that_is_not_an_identifier(con):
    with pytest.raises(ScriptError, match="not usable as a variable name"):
        scriptlib.set_variables(con, {"min depth": 1})


# --- includes -------------------------------------------------------------


def test_read_resolves_against_the_including_file(con, write, tmp_path):
    (tmp_path / "lib").mkdir()
    (tmp_path / "lib" / "thing.sql").write_text("CREATE VIEW t AS SELECT 1 AS x;")
    path = write(".read lib/thing.sql\n")
    scriptlib.run(con, path)
    assert con.execute("SELECT x FROM t").fetchone() == (1,)


def test_a_local_copy_shadows_the_shipped_one(con, write, tmp_path, hits):
    """Forking a library file is a copy and an edit, with no flag to remember."""
    (tmp_path / "lib").mkdir()
    (tmp_path / "lib" / "count_sites.sql").write_text(
        "CREATE VIEW site_relation AS SELECT 'mine' AS marker;"
    )
    path = write(".read lib/count_sites.sql\n")
    scriptlib.run(con, path, include=(SQL_DIR,))
    assert con.execute("SELECT marker FROM site_relation").fetchone() == ("mine",)


def test_missing_include_names_everywhere_it_looked(con, write):
    path = write(".read lib/nope.sql\n")
    with pytest.raises(ScriptError) as exc:
        scriptlib.run(con, path)
    assert "cannot find 'lib/nope.sql'" in str(exc.value)
    assert str(SQL_DIR) in str(exc.value)


def test_a_cycle_is_refused_rather_than_recursing(con, write, tmp_path):
    (tmp_path / "a.sql").write_text(".read b.sql\n")
    (tmp_path / "b.sql").write_text(".read a.sql\n")
    with pytest.raises(ScriptError, match="reads itself"):
        scriptlib.run(con, tmp_path / "a.sql")


def test_python_runs_against_the_same_connection(con, write, tmp_path):
    (tmp_path / "w.py").write_text(
        "rows = con.execute('SELECT x FROM t').fetchall()\n"
        "open(var('out'), 'w').write(str(rows))\n"
    )
    out = tmp_path / "out.txt"
    path = write(
        f"CREATE VIEW t AS SELECT 7 AS x;\n"
        f"SET VARIABLE out = '{out}';\n"
        f".python w.py\n"
    )
    scriptlib.run(con, path)
    assert out.read_text() == "[(7,)]"


# --- errors ---------------------------------------------------------------


def test_a_failing_statement_is_quoted_with_its_file_and_line(con, write):
    path = write("CREATE VIEW a AS SELECT 1;\n\nSELECT * FROM missing_view;\n")
    with pytest.raises(ScriptError) as exc:
        scriptlib.run(con, path)
    message = str(exc.value)
    assert "run.sql" in message
    assert "SELECT * FROM missing_view" in message
    assert "missing_view" in message


def test_load_order_is_enforced_by_duckdb_not_by_convention(con, write, hits):
    """lib/filters.sql reads `src`, so reading it first names `src`."""
    path = write(".read lib/filters.sql\n")
    with pytest.raises(ScriptError) as exc:
        scriptlib.run(con, path, values={"hits": hits})
    assert "measure" in str(exc.value) or "src" in str(exc.value)


# --- show -----------------------------------------------------------------


def test_expand_substitutes_every_read(write, tmp_path):
    (tmp_path / "lib").mkdir()
    (tmp_path / "lib" / "x.sql").write_text("SELECT 'inner';")
    path = write("SELECT 'outer';\n.read lib/x.sql\n")
    text = scriptlib.expand(path)
    assert "'outer'" in text and "'inner'" in text and ".read" not in text


def test_expand_leaves_python_as_a_line(write, tmp_path):
    (tmp_path / "w.py").write_text("pass\n")
    path = write(".python w.py\n")
    text = scriptlib.expand(path)
    assert ".python w.py" in text and "w.py" in text


def test_results_are_shown_and_ddl_is_not(con, write):
    shown = []
    path = write("CREATE VIEW v AS SELECT 1 AS n;\nSELECT n FROM v;\n")
    scriptlib.run(con, path, show=shown.append)
    assert len(shown) == 1
    assert "n" in shown[0] and "1" in shown[0]


def test_trace_records_what_ran(con, write, tmp_path):
    (tmp_path / "lib").mkdir()
    (tmp_path / "lib" / "x.sql").write_text("SELECT 1;")
    lines = []
    path = write(".read lib/x.sql\n")
    scriptlib.run(con, path, trace=lines.append)
    assert any(line.startswith(".read ") for line in lines)
    assert any("SELECT 1" in line for line in lines)
