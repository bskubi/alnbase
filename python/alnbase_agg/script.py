"""Running a `.sql` file, which is the whole of an aggregation.

There is no recipe format. A run is a SQL script: it names the hit tables, says
what the query names mean, says which observations survive, and writes files.
Everything this module adds to `duckdb.connect().execute()` is two include
directives and a way to pass values in, because those are the only two things a
single `execute()` call cannot do.

# The two directives

A line whose first non-blank characters are `.read` or `.python` is a directive
rather than SQL. `.read` is DuckDB's own CLI spelling for "execute that file
here", reused rather than reinvented; `.python` is its counterpart for the few
formats that SQL cannot write on its own.

    .read lib/filters.sql          -- more SQL, executed at this point
    .python formats/amethyst_h5.py -- Python, run against this connection

Both take one path, resolved in this order:

1. against the directory of the file doing the reading, so a script's own
   library sits beside it and a shipped format file can read its neighbours;
2. against each directory on the include path, which is how the library this
   package ships is found without anybody writing an absolute path.

The order is what makes forking work. `alnbase-agg lib` prints the shipped
directory; copy `formats/bismark_cov.sql` next to your script, and your copy is
the one `.read formats/bismark_cov.sql` finds, with no flag and no registry.

# Values

A directive takes no arguments, and it does not need to, because DuckDB has
session variables and they are late bound:

    SET VARIABLE relation = 'mCG';
    SET VARIABLE out = 'sample.CG.cov';
    .read formats/bismark_cov.sql

    SET VARIABLE relation = 'mCH';
    SET VARIABLE out = 'sample.CH.cov';
    .read formats/bismark_cov.sql

The format file reads `getvariable('relation')`. It is read when the query runs,
not when the file was written, so one file written once produces both outputs.
This is also how the command line reaches the script: `-s hits=calls_*.parquet`
is a `SET VARIABLE` executed before the first line of the script.

Two consequences worth knowing. `getvariable` on a name that was never set
returns NULL rather than raising, so a file gives an optional value a default
with `coalesce(getvariable('min_depth'), 1)` and there is no such thing as a
required-flag list to maintain. And a variable holds a typed value, so `-s
min_depth=5` must arrive as the integer 5 and not the string `'5'`; see
`parse_value`.

# Order of definition

DuckDB binds a view's *macros* late but its *tables and views* eagerly:
`CREATE VIEW v AS SELECT * FROM not_yet_defined` fails at once. That is a
feature here, because it turns the library's load order into an error message
instead of a convention. `lib/filters.sql` reads `src` and `measure`, so it
must be read after the script defines them; `lib/count_sites.sql` reads `resolved`,
so it comes after that. Reading them in the wrong order names the missing view.
"""

from __future__ import annotations

import re
from pathlib import Path

#: The directory holding the SQL this package ships: `lib/` and `formats/`.
#: On the include path of every run, and the thing `alnbase-agg lib` prints.
SQL_DIR = Path(__file__).parent / "sql"

#: A directive line: a leading dot, a verb, and one path. Anything else is SQL.
#: The path runs to the end of the line and is stripped, so it may contain
#: spaces; it may not contain a comment, which is why `--` is not special here.
_DIRECTIVE = re.compile(r"^\s*\.(read|python)\s+(\S.*?)\s*$")


class ScriptError(Exception):
    """A script could not be run, phrased for the person who wrote it."""


def parse_value(text: str):
    """A command-line `-s name=value` string as the value it obviously is.

    A session variable is typed, and the type shows: `HAVING sum(den) >=
    getvariable('min_depth')` compares a count against whatever was set, so a
    string `'5'` and an integer 5 are not the same argument.

    The rule is deliberately narrow. A value is converted only when converting
    it and printing it back gives the original text, so `5` becomes an integer
    and `0.4` a float, while `01` and `1e3` stay strings -- a barcode that
    happens to look like a number keeps its leading zero. `true` and `false`
    become booleans. Everything else, including every path and glob, is a
    string, which is what makes the common case need no thought at all.
    """
    if text in ("true", "false"):
        return text == "true"
    for cast in (int, float):
        try:
            value = cast(text)
        except ValueError:
            continue
        if str(value) == text:
            return value
    return text


def set_variables(con, values: dict) -> None:
    """`SET VARIABLE` for each name, with the value bound rather than spliced.

    Bound, so a value containing a quote is a value and not a syntax error.
    DuckDB refuses a bound parameter inside `CREATE MACRO`, which is what makes
    the variable route necessary in the first place, but `SET VARIABLE` itself
    takes one happily.
    """
    for name, value in values.items():
        if not name.isidentifier():
            raise ScriptError(
                f"{name!r} is not usable as a variable name. A name is what you "
                f"would write inside getvariable('...'): letters, digits and "
                f"underscores, not starting with a digit."
            )
        con.execute(f"SET VARIABLE {name} = $value", {"value": value})


def resolve(path: str, *, relative_to: Path, include=()) -> Path:
    """Where a directive's path points, or an error naming everywhere we looked.

    `relative_to` is the directory of the file containing the directive, and it
    wins, so a copied-and-edited file shadows the shipped one of the same name.
    """
    candidate = Path(path)
    if candidate.is_absolute():
        if candidate.exists():
            return candidate
        raise ScriptError(f"no such file: {candidate}")

    tried = [relative_to / candidate] + [Path(d) / candidate for d in include]
    for option in tried:
        if option.exists():
            return option
    looked = "\n".join(f"    {p}" for p in tried)
    raise ScriptError(f"cannot find {path!r}. Looked for:\n{looked}")


#: How many rows of a statement's result are shown before the rest are
#: summarised. A run script's SELECTs are meant to be small -- how many reads
#: each filter rejected, how many blocks survived -- and a script that means to
#: emit a large table writes a file.
RESULT_ROWS = 40


def run(con, path, *, include=(SQL_DIR,), values=None, trace=None, show=None) -> None:
    """Execute one script against `con`, following its directives.

    `values` are set before the first line runs, so a script may use them in
    `read_parquet` and in a `CREATE VIEW` body.

    `trace`, if given, is called with a one-line description of every directive
    followed and every statement executed: a run's own account of what it did,
    which is the honest replacement for asking a `.sql` file what filters it
    applied.

    `show`, if given, is called with the formatted result of every statement
    that returns rows. That is what makes a QC query in a run script worth
    writing -- `SELECT count(*) FROM read_fails('mCH', 'ratio', '>', 0.5)` reports what the
    filter threw away -- and it is also how `COPY` reports how many rows it
    wrote to each file.
    """
    if values:
        set_variables(con, values)
    _run_file(
        con, Path(path).resolve(), tuple(Path(d) for d in include), trace, show, ()
    )


def _run_file(con, path: Path, include, trace, show, stack) -> None:
    """One file, statement by statement, with the include stack for errors."""
    if path in stack:
        chain = " -> ".join(p.name for p in (*stack, path))
        raise ScriptError(f"a file reads itself, directly or not: {chain}")

    text = _read(path)
    buffered, first_line = [], 1

    def flush():
        """Execute what has accumulated, attributing a failure to its statement."""
        if buffered and "".join(buffered).strip():
            _execute(con, "\n".join(buffered), path, first_line, trace, show)
        buffered.clear()

    for number, line in enumerate(text.splitlines(), start=1):
        directive = _DIRECTIVE.match(line)
        if directive is None:
            if not buffered:
                first_line = number
            buffered.append(line)
            continue
        flush()
        first_line = number + 1
        verb, argument = directive.groups()
        target = resolve(argument, relative_to=path.parent, include=include)
        if trace:
            trace(f".{verb} {target}")
        if verb == "read":
            _run_file(con, target, include, trace, show, (*stack, path))
        else:
            _run_python(con, target)
    flush()


def _read(path: Path) -> str:
    try:
        return path.read_text()
    except OSError as exc:
        raise ScriptError(f"cannot read {path}: {exc}") from None


def _execute(con, sql: str, path: Path, first_line: int, trace, show=None) -> None:
    """Run a block of SQL one statement at a time.

    One at a time for the error message. DuckDB will accept the whole block,
    but then a failure names neither which of the statements failed nor where
    it came from, and a script that reads four library files is otherwise very
    hard to place an error in.
    """
    try:
        statements = con.extract_statements(sql)
    except Exception as exc:  # a syntax error, which has its own position
        raise ScriptError(f"{path}, from line {first_line}:\n    {exc}") from None

    for statement in statements:
        query = statement.query.strip()
        if trace:
            trace(_summarise(query))
        try:
            con.execute(query)
            result = _format_result(con) if show else None
        except Exception as exc:
            raise ScriptError(
                f"{path}, in the block starting at line {first_line}:\n"
                f"{_indent(query)}\n"
                f"  {type(exc).__name__}: {exc}"
            ) from None
        if result:
            show(result)


def _format_result(con) -> str | None:
    """A statement's result as a small aligned table, or None if it has none.

    DDL returns no result set at all, so most of a script prints nothing. What
    does print is the two things a run wants to say out loud: what a QC query
    counted, and how many rows each `COPY` wrote.
    """
    if con.description is None:
        return None
    names = [d[0] for d in con.description]
    rows = con.fetchmany(RESULT_ROWS)
    if not rows:
        return None
    extra = con.fetchone() is not None
    cells = [names] + [["" if v is None else str(v) for v in row] for row in rows]
    widths = [max(len(row[i]) for row in cells) for i in range(len(names))]
    lines = ["  ".join(v.ljust(w) for v, w in zip(row, widths)).rstrip() for row in cells]
    lines.insert(1, "  ".join("-" * w for w in widths))
    if extra:
        lines.append(f"... ({RESULT_ROWS} rows shown)")
    return "\n".join(lines)


def _run_python(con, path: Path) -> None:
    """Execute a format's Python against this connection.

    The file is run, not imported: it is part of the run in the same way the
    SQL files are, and giving it a module identity would only invite somebody
    to import it for a side effect. It is handed `con`, already carrying every
    view the script has defined, and `var`, which reads a session variable --
    so a Python writer takes its output path the same way a SQL one does.
    """
    def var(name, default=None):
        value = con.execute("SELECT getvariable($n)", {"n": name}).fetchone()[0]
        return default if value is None else value

    namespace = {
        "__name__": "__main__",
        "__file__": str(path),
        "con": con,
        "var": var,
    }
    try:
        exec(compile(path.read_text(), str(path), "exec"), namespace)
    except Exception as exc:
        raise ScriptError(f"{path}: {type(exc).__name__}: {exc}") from exc


def expand(path, *, include=(SQL_DIR,), depth=0) -> str:
    """The script with every `.read` substituted, as one block of text.

    What `alnbase-agg show` prints. A run is a query, so the query is something
    a person can be shown before it runs, paste into a DuckDB shell, or attach
    to a methods section. `.python` cannot be substituted and is left as the
    line it is, marked with where it resolved to.
    """
    path = Path(path).resolve()
    include = tuple(Path(d) for d in include)
    out = [f"-- {'-' * 8} {path} {'-' * 8}"]
    for line in _read(path).splitlines():
        directive = _DIRECTIVE.match(line)
        if directive is None:
            out.append(line)
            continue
        verb, argument = directive.groups()
        target = resolve(argument, relative_to=path.parent, include=include)
        if verb == "read":
            out.append(expand(target, include=include, depth=depth + 1))
            out.append(f"-- {'-' * 8} back in {path.name} {'-' * 8}")
        else:
            out.append(f"{line}    -- {target}")
    return "\n".join(out)


def _indent(text: str, prefix: str = "    ") -> str:
    return "\n".join(prefix + line for line in text.splitlines())


def _summarise(query: str, width: int = 72) -> str:
    """One line describing a statement, for the trace."""
    flat = " ".join(query.split())
    return flat if len(flat) <= width else flat[: width - 1] + "…"
