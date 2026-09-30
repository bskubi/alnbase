"""Building a run by chained method calls, which emits the same SQL script.

A run is a SQL script; `alnbase_agg.script` executes one. This module is a
second way to write that script, for people who would rather call a method than
write a `SET VARIABLE` and a `.read` beside it.

Nothing here is a new language. The rule is one sentence:

    a method is a file in the library, and a keyword argument is one of that
    file's session variables.

`run.bismark_cov(mCG, out="s.cov")` emits `SET VARIABLE relation`, `SET VARIABLE
out` and `.read formats/bismark_cov.sql`, and that is the whole of what it does.
So the library documents the API: `alnbase-agg lib --list` names every method,
and the VARIABLES block at the top of each file names its arguments. Adding a
format file adds a method with no Python written anywhere.

The one exception is `Measure`, below, which mirrors a table the library
already requires rather than inventing anything. It is the only concept in this
module that is not a spelling of a file.

# The script is the output, not an implementation detail

`Run.sql()` returns the script, and it is an ordinary one: `alnbase-agg run`
executes it, a DuckDB shell accepts most of it, and a methods section can quote
it. So the builder is removable. Generate the script once, save it, delete the
Python, and the run still works -- which is the property that keeps this an
alternative front end rather than a fork of the toolchain.

# What the chain adds over writing the script by hand

**A format call becomes one call.** In the script a format takes its arguments
by setting session variables beside the include, so two statements and a
directive say what `bismark_cov(mCG, out=...)` says.

**A variable cannot leak into the next format.** Session variables are global
and stay set, so a script that passes `min_depth` to one format has quietly
passed it to every later one as well. `Run` remembers what each call set and
emits `RESET VARIABLE` for anything the next call does not set.

**The filters are named once.** A hand-written script lists its filters in the
`resolved` view and then lists them again in the QC query that reports what each
threw away, with the thresholds written twice and nothing keeping them equal.
Here `drop_reads` records the list and `report()` renders it from that.

**The load order is guaranteed.** `lib/filters.sql` must be read after `src` and
`measure` exist and before `resolved` is defined, because DuckDB binds a macro's
tables when the macro is created; reading it early yields filters that silently
match nothing. The file has to warn about that in prose. Here the declarations
are recorded and emitted in the right order by whichever call first needs them,
so the hazard is not reachable.

**A file's dependencies come from the file.** Each one names what it needs on a
`REQUIRES` line in its header -- `formats/bismark_cov.sql` needs
`lib/count_sites.sql`, and `formats/amethyst_h5.sql` needs contig order as well
-- and the builder reads those first and once. So a chain says nothing about
counting sites or loading a FASTA index: a run that writes no HDF5 never
mentions a `.fai`, and a run that writes one cannot forget it.
"""

from __future__ import annotations

from pathlib import Path

from alnbase_agg.script import SQL_DIR, ScriptError

#: What identifies one read, for a hit table that selected these fields. It is
#: the read and not the alignment record: a chimeric Methyl-HiC read is several
#: records, and judging them apart would filter half a molecule. A run whose
#: hit table lacks these columns fails naming the missing one, which is the
#: thing it needed to be told.
READ = {"qname": "qname", "mate": "is_first_in_template"}

#: How the library's column names are reached from a hit table's. `src` writes
#: both, so a run works today and keeps working when alnbase renames its own.
ALIASES = {"contig": "ref_name"}


class Var:
    """A value read from a session variable when the script runs, not now.

    The common case needs none of this: a path known to the Python that builds
    the script is written into the script. `Var` is for the other case, where
    the generated script is meant to be re-run against a different block with
    `alnbase-agg run -s block=other`, and so must not have this block's name
    baked into it.

        out=var("block") + ".CG.cov"     -->  (getvariable('block') || '.CG.cov')
    """

    def __init__(self, sql: str):
        self.sql = sql

    def __add__(self, other) -> "Var":
        return Var(f"({self.sql} || {literal(other)})")

    def __radd__(self, other) -> "Var":
        return Var(f"({literal(other)} || {self.sql})")

    def __repr__(self) -> str:
        return f"var({self.sql})"


def var(name: str) -> Var:
    """`getvariable('name')`, for a value the run supplies rather than this code."""
    return Var(f"getvariable({literal(name)})")


def literal(value) -> str:
    """A Python value as the SQL literal for it.

    Quoting is doubling a quote, which is the whole of SQL string escaping and
    is worth doing by hand here because the result has to be readable: this
    text is what `Run.sql()` shows somebody, so a parameter marker would make
    the script unreadable and unrunnable on its own.
    """
    if isinstance(value, Var):
        return value.sql
    if value is None:
        return "NULL"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value)
    text = str(value)
    return "'" + text.replace("'", "''") + "'"


class Fails:
    """One filter: a table-returning macro call, and a label for the report.

    `level` says what the filter judges, so passing a block filter to
    `drop_reads` is caught here rather than producing a view that compares a
    block key against a read identifier and drops nothing.
    """

    def __init__(self, sql: str, label: str, *, level: str = "read",
                 measure=None, relation=None):
        self.sql = sql
        self.label = label
        self.level = level
        self.measure = measure
        self.relation = relation

    def __repr__(self) -> str:
        return f"Fails({self.sql!r})"


def fails(sql: str, label: str | None = None, *, level: str = "read") -> Fails:
    """A filter written as SQL: any table expression with a `read` or `block`
    column.

    The escape hatch, and the way a `CASE` ladder is written, since spelling
    one out says more than a name for it would.
    """
    return Fails(sql, label or " ".join(sql.split())[:40], level=level)


def in_regions(path) -> "Test":
    """Sites inside any region of a BED file: a blacklist, or known SNPs.

    Three columns are read and the rest of the line ignored, and the
    coordinates are BED's own: 0-based and half open. The file is read at query
    time, so regenerating it and re-running changes the answer without touching
    the script.
    """
    return Test("in_regions", (literal(path),), levels=("site",),
                label="in regions")


def outside_regions(path) -> "Test":
    """Sites outside every region of a BED file: a capture panel.

    The other direction, and not reachable by inverting the first: restricting
    a run to the positions a site file names is what premethyst's
    `context-extract` does, and doing it without this meant complementing the
    BED by hand.
    """
    return Test("outside_regions", (literal(path),), levels=("site",),
                label="outside regions")


def at_position(test: "Test") -> "Test":
    """Judge a test at the position and reject it in every block.

        run.drop_sites(at_position(snp.relate("G", "A").ratio.gt(0.5)))

    The SNP veto. A C read as a T may be bisulfite conversion or a C-to-T
    variant, and on that strand the two look identical; the evidence is a
    different relation at the same position on the other strand, in whichever
    blocks happened to cover it. A filter keyed to one block's own row cannot
    reach that, so this one pools the evidence across blocks and strands and
    drops the position everywhere.
    """
    if test.macro is not None:
        raise ScriptError(
            f"at_position() takes a plain quantity test, not {test.label!r}; "
            "the quantities with a macro of their own are already about "
            "something other than one group's totals"
        )
    return Test("position", test.args, relation=test.relation,
                levels=("site",), label=f"at a position where {test.label}")


def both(*tests: "Test") -> "Test":
    """A test that fails only where every one of its parts fails.

        run.drop_sites(both(mCG.main.gt(0), mCG.other.gt(0)))

    Filters compose by union, which is right: each is an independent reason to
    throw something away. `both` is the other connective, for the cases where
    one reason is a conjunction -- scbs drops a site that is neither wholly
    methylated nor wholly unmethylated, which is a condition on both sides of
    the relation at once and no single comparison can say it.

    It is also what makes a ladder expressible. premethyst's read cut is a
    threshold that depends on how much evidence there is, and it is two
    conjunctions unioned:

        run.drop_reads(
            both(mCH.total.ge(3), mCH.total.le(5), mCH.main.gt(1)),
            both(mCH.total.gt(5), mCH.ratio.ge(0.4)),
        )

    Note what is not checked: the recipe language this replaced required a
    ladder's rungs to tile the range exactly once, with no gap and no overlap.
    Here an overlap is harmless, because the parts are unioned, and a gap
    simply means those reads are kept. Nothing will tell you that you left one.
    """
    if len(tests) < 2:
        raise ScriptError("both() needs at least two tests to combine")
    for test in tests:
        if not isinstance(test, Test):
            raise ScriptError(
                f"both() combines tests, not {type(test).__name__}"
            )
    return _Both(tests)


def listed(path) -> "Test":
    """Reads whose name appears in a file: a blacklist, a duplicate set, or
    the output of another tool entirely.

    One of the few filters that is about no measure, so it is a function here
    rather than a method. The file is read at query time, so regenerating it
    and re-running changes the answer without touching the script.
    """
    return Test("listed", (literal(path),), levels=("read",), label="listed")


# ---------------------------------------------------------------------------
# Measures, relations, quantities
# ---------------------------------------------------------------------------


class Measure:
    """alnbase query names, sorted into named buckets called categories.

        mCG = Measure("mCG", methylated="CG", unmethylated="TG")
        gt  = Measure("gt", A="A>A", C="A>C", G="A>G", T="A>T")

    This is the `measure` table of `lib/filters.sql` as a Python object, and
    the only place biology appears in a run. You choose the category names; no
    name is reserved.

    A query name belongs to at most one category *of this measure*. It may
    belong to no category at all, and to a category of some other measure: a
    NOMe-seq run declares one measure over its accessibility queries and
    another over its methylation queries, sharing nothing.

    A count on its own is rarely the question. `relate()` puts two sets of
    categories against each other, and a measure with exactly two categories
    has one possible relation, reachable without naming it.
    """

    def __init__(self, name: str, categories: dict | None = None, **by_name):
        table = {}
        for category, entries in {**(categories or {}), **by_name}.items():
            table[category] = [entries] if isinstance(entries, str) else list(entries)
        if not name:
            raise ScriptError("a measure needs a name")
        if len(table) < 2:
            raise ScriptError(
                f"measure {name!r} has {len(table)} categor"
                f"{'y' if len(table) == 1 else 'ies'}; it needs at least two, "
                "because a count with nothing to compare against is just a count"
            )
        seen = {}
        for category, entries in table.items():
            for entry in entries:
                if entry in seen:
                    raise ScriptError(
                        f"measure {name!r} puts query name {entry!r} in both "
                        f"{seen[entry]!r} and {category!r}. Within one measure a "
                        "query name belongs to at most one category, so that a "
                        "denominator can be a plain sum. Two independent readouts "
                        "are two measures, not two categories."
                    )
                seen[entry] = category
        self.name = name
        self.categories = table

    def __repr__(self) -> str:
        return f"Measure({self.name!r}, {list(self.categories)})"

    def rows(self) -> list[str]:
        """The `(measure, name, category)` rows this measure contributes."""
        return [
            f"({literal(self.name)}, {literal(entry)}, {literal(category)})"
            for category, entries in self.categories.items()
            for entry in entries
        ]

    def relate(self, main=None, other=None, *, name: str | None = None) -> "Relation":
        """Put one set of categories against another.

            mCG.relate()                   the only pair, for two categories
            snp.relate("G", "A")           G against A
            snp.relate("G")                G against everything else
            snp.relate(["G", "A"], "C")    either side may gather several

        `main` is the side in focus and `other` is what it is set against. The
        two are disjoint and neither repeats the other, exactly as the
        categories they are built from do not overlap. `main` and `other` are
        also the two counts you read back off it; `total` is both, and `ratio`
        is the side in focus over that total.
        """
        if main is None:
            if len(self.categories) != 2:
                raise ScriptError(
                    f"measure {self.name!r} has {len(self.categories)} categories "
                    f"({', '.join(self.categories)}), so there is no single "
                    "relation to mean. Say which: m.relate('a', 'b')."
                )
            main, other = list(self.categories)[0], list(self.categories)[1]
        main = self._check(main)
        if other is None:
            other = [c for c in self.categories if c not in main]
        else:
            other = self._check(other)
        overlap = set(main) & set(other)
        if overlap:
            raise ScriptError(
                f"relation on {self.name!r} puts {sorted(overlap)} on both sides; "
                "the two sides are disjoint, and the denominator is their sum"
            )
        return Relation(self, name or _relation_name(self, main, other), main, other)

    def _check(self, categories) -> list[str]:
        if isinstance(categories, str):
            categories = [categories]
        for category in categories:
            if category not in self.categories:
                known = ", ".join(repr(c) for c in self.categories)
                raise ScriptError(
                    f"measure {self.name!r} has no category {category!r}; "
                    f"it has {known}"
                )
        return list(categories)

    # -- the only relation, for a measure that has one --------------------

    @property
    def main(self) -> "Quantity":
        """How many hits fell on the side in focus."""
        return self.relate().main

    @property
    def other(self) -> "Quantity":
        """How many hits fell on the side it is set against."""
        return self.relate().other

    @property
    def total(self) -> "Quantity":
        """How many hits fell on the relation at all, both sides together."""
        return self.relate().total

    @property
    def ratio(self) -> "Quantity":
        """The first side over both sides."""
        return self.relate().ratio

    @property
    def consecutive(self) -> "Quantity":
        """The longest run of adjacent main-side calls along the read."""
        return self.relate().consecutive

    @property
    def reads(self) -> "Quantity":
        """How many reads the block contributed. Block level only."""
        return self.relate().reads

    @property
    def sites(self) -> "Quantity":
        """How many distinct positions the block covered. Block level only."""
        return self.relate().sites

    @property
    def blocks(self) -> "Quantity":
        """How many blocks covered the position. Site level only."""
        return self.relate().blocks

    def __getattr__(self, name: str):
        """`m.whatever(...)` for `<level>_fails_whatever('m', ...)`.

        Reached only for a macro with no quantity above, which means a filter
        added to a forked library. Such a call gets no completion and no
        signature, but it works, so forking does not mean leaving this front
        end. The level still comes from the `drop_*` call, as it does for
        everything else.
        """
        if name.startswith("_"):
            raise AttributeError(name)
        return self.relate().__getattr__(name)

    def missing(self) -> "Test":
        """Reads with no hits of this measure at all, which no threshold can
        reach. premethyst discards exactly these; whether that is right is a
        real choice, so it is spelled rather than assumed."""
        return Test("missing", (literal(self.name),), measure=self, levels=("read",))


def _relation_name(measure: "Measure", main, other) -> str:
    """A readable name for a relation nobody named.

    It appears in the script, in the report and in a DuckDB error, so it is
    built out of the category names rather than being a number.
    """
    if len(measure.categories) == 2 and len(main) == 1 and len(other) == 1:
        return measure.name          # the only relation: it *is* the measure
    return f"{measure.name}.{'+'.join(main)}_vs_{'+'.join(other)}"


class Relation:
    """Two sets of a measure's categories, put against each other.

    Everything a filter tests and every file format writes is an
    interpretation of one of these, reached as a property:

        mCH.ratio.gt(0.4)          the non-conversion cut
        mCH.ratio.ge(0.4)       premethyst's, whose boundary differs
        mCG.total.lt(3)            too little evidence to judge
        mCH.main.gt(2)             two methylated calls is two too many

    A measure with exactly two categories has one possible relation and
    declares it under the measure's own name, which is why the common case
    above never says the word.
    """

    def __init__(self, measure: Measure, name: str, main: list, other: list):
        self.measure = measure
        self.name = name
        #: The category names on each side. `main` and `other` themselves are
        #: the quantities -- how many hits fell on each -- so the lists that
        #: define them are spelled out.
        self.main_categories = main
        self.other_categories = other

    def __repr__(self) -> str:
        return (f"Relation({self.name!r}, {self.main_categories} vs "
                f"{self.other_categories})")

    def row(self) -> str:
        """The `(relation, measure, main, other)` row this relation contributes."""
        sides = ", ".join(
            "[" + ", ".join(literal(c) for c in side) + "]"
            for side in (self.main_categories, self.other_categories)
        )
        return f"({literal(self.name)}, {literal(self.measure.name)}, {sides})"

    def _quantity(self, name, macro=None, levels=("read", "block", "site")):
        return Quantity(self, name, macro=macro, levels=levels)

    @property
    def main(self) -> "Quantity":
        """How many hits fell on the side in focus -- the categories you named
        `main`. Two methylated CH calls is two too many however few sites the
        read covered, which is the thing a ratio cannot say."""
        return self._quantity("main")

    @property
    def other(self) -> "Quantity":
        """How many hits fell on the side it is set against."""
        return self._quantity("other")

    @property
    def total(self) -> "Quantity":
        """Both sides together: how much of the relation was seen at all. The
        companion to every ratio filter, since a ratio over two sites is not
        evidence. This is the `den` a file format writes."""
        return self._quantity("total")

    @property
    def ratio(self) -> "Quantity":
        """The side in focus, over the total."""
        return self._quantity("ratio")

    @property
    def consecutive(self) -> "Quantity":
        """The longest run of adjacent first-side calls along the read.

        The test a ratio provably cannot express: four methylated CH in a row
        and four scattered give the same ratio, and only the first looks like a
        patch of failed conversion. A window over positions rather than an
        aggregate, so it has a macro of its own and exists only for reads.
        """
        return self._quantity("consecutive", macro="consecutive", levels=("read",))

    @property
    def reads(self) -> "Quantity":
        """How many reads the block contributed.

        Not a hit count: a block with 2000 hits may be one read seen 2000 times
        or 2000 reads seen once, and `unique_reads >= 10000` is the commonest
        block cut in this ecosystem. Block level only, because at read level it
        is 1 and at site level it is `total`.
        """
        return self._quantity("reads", macro="reads", levels=("block",))

    @property
    def sites(self) -> "Quantity":
        """How many distinct positions the block covered.

        The other commonest block cut -- premethyst's `calls-filter -G` is a
        minimum number of covered CG *positions*, not of calls. Block level
        only, for the same reason.
        """
        return self._quantity("sites", macro="sites", levels=("block",))

    @property
    def blocks(self) -> "Quantity":
        """How many blocks covered the position at all.

        A count across blocks rather than over one site's own rows, so it has a
        macro of its own and exists only for sites. A position one block saw is
        not a position.
        """
        return self._quantity("blocks", macro="blocks", levels=("site",))

    def __getattr__(self, name: str):
        """A macro in a forked library, called by its own name."""
        if name.startswith("_"):
            raise AttributeError(name)

        def call(*args) -> "Test":
            rendered = (literal(self.name), *(literal(a) for a in args))
            return Test(
                name, rendered, relation=self,
                label=f"{name} {' '.join(str(a) for a in args)}".strip(),
            )

        return call


class Quantity:
    """A number a filter can be about, before anything is compared to it.

    A filter is three choices -- what it groups on, what it computes, and when
    it fails -- and this is the middle one. The comparison methods below supply
    the third; the `drop_reads` / `drop_blocks` / `drop_sites` call supplies the
    first.
    """

    def __init__(self, relation: Relation, name: str, *, macro=None, levels=()):
        self.relation = relation
        self.name = name
        self.macro = macro
        self.levels = levels

    def __repr__(self) -> str:
        return f"Quantity({self.relation.name!r}.{self.name})"

    def _test(self, op: str, threshold) -> "Test":
        _refuse_unsatisfiable(self, op, threshold)
        if self.macro:
            args = (literal(self.relation.name), literal(op), literal(threshold))
            return Test(
                self.macro, args, relation=self.relation, levels=self.levels,
                label=f"{self.name} {op} {threshold}",
            )
        args = (
            literal(self.relation.name),
            literal(self.name),
            literal(op),
            literal(threshold),
        )
        return Test(
            None, args, relation=self.relation, levels=self.levels,
            label=f"{self.name} {op} {threshold}",
        )

    def gt(self, threshold) -> "Test":
        """Fails when the quantity is greater than `threshold`."""
        return self._test(">", threshold)

    def ge(self, threshold) -> "Test":
        """Fails at `threshold` and above. premethyst's read cut is this one,
        which is why the boundary has to be sayable."""
        return self._test(">=", threshold)

    def lt(self, threshold) -> "Test":
        """Fails when the quantity is less than `threshold`."""
        return self._test("<", threshold)

    def le(self, threshold) -> "Test":
        """Fails at `threshold` and below."""
        return self._test("<=", threshold)

    def eq(self, threshold) -> "Test":
        """Fails on exactly this value."""
        return self._test("=", threshold)

    def ne(self, threshold) -> "Test":
        """Fails on anything but this value."""
        return self._test("<>", threshold)


def _refuse_unsatisfiable(quantity: "Quantity", op: str, threshold) -> None:
    """Refuse a depth threshold that no row can ever meet.

    A read, block or site has a row only because it has at least one hit, so
    its total is never below 1. `total.lt(1)` therefore binds, runs, and
    matches nothing -- which is the failure this library works hardest to
    prevent, because it looks exactly like a filter that found nothing to
    reject.

    What such a call means is "it saw none of this measure at all", and that is
    `missing()`, which reads `src` instead and so can see a read that the
    counting views never gave a row to.
    """
    if quantity.name != "total" or not isinstance(threshold, (int, float)):
        return
    reachable = {
        "<": threshold > 1,
        "<=": threshold >= 1,
        "=": threshold >= 1,
    }.get(op, True)
    if reachable:
        return
    raise ScriptError(
        f"total {op} {threshold} can never match: a row exists only because "
        "there was at least one hit, so the total is never below 1. A read that "
        "saw none of the measure has no row here at all -- "
        f"{quantity.relation.measure.name}.missing() is the filter for that."
    )


class Test:
    """A quantity and a fail condition, not yet attached to a grouping.

    `drop_reads`, `drop_blocks` and `drop_sites` each render the same test
    against their own level, which is why "blocks whose depth is below 1000"
    needs no name of its own.
    """

    def __init__(self, macro, args, *, relation=None, measure=None,
                 levels=("read", "block", "site"), label=None):
        self.macro = macro
        self.args = args
        self.relation = relation
        self.measure = measure or (relation.measure if relation else None)
        self.levels = levels
        self.label = label or (macro or "")

    def __repr__(self) -> str:
        return f"Test({self.label!r})"

    def at(self, level: str) -> "Fails":
        """Render against one grouping: `read_fails(...)` and its siblings."""
        if level not in self.levels:
            where = " or ".join(self.levels)
            raise ScriptError(
                f"{self.label} is not a {level}-level test; it applies to "
                f"{where}, because that quantity is not defined per {level}"
            )
        macro = f"{level}_fails" + (f"_{self.macro}" if self.macro else "")
        named = self.relation or self.measure
        name = named.name if named else ""
        return Fails(
            f"{macro}({', '.join(self.args)})",
            " ".join(p for p in (level, name, self.label) if p),
            level=level,
            measure=self.measure,
            relation=self.relation,
        )


class _Both:
    """The intersection of several tests, rendered at whichever level it is
    handed to. Holds no SQL of its own: each part renders itself and the parts
    are intersected, so `both` gains every quantity the parts have."""

    def __init__(self, tests):
        self.tests = tests
        self.levels = set.intersection(*(set(t.levels) for t in tests))
        self.label = " and ".join(t.label for t in tests)
        self.relation = next((t.relation for t in tests if t.relation), None)
        self.measure = next((t.measure for t in tests if t.measure), None)
        self.macro = None

    def __repr__(self) -> str:
        return f"both({self.label!r})"

    def at(self, level: str) -> "Fails":
        if level not in self.levels:
            where = " or ".join(sorted(self.levels)) or "no level in common"
            raise ScriptError(
                f"{self.label} is not a {level}-level test; its parts apply to "
                f"{where}"
            )
        column = {"read": "read", "block": "block", "site": "site"}[level]
        parts = [f"SELECT {column} FROM {t.at(level).sql}" for t in self.tests]
        return Fails(
            "(" + " INTERSECT ".join(parts) + ")",
            f"{level} {self.label}",
            level=level,
            measure=self.measure,
            relation=self.relation,
        )


def as_relation(item):
    """A `Relation` from whatever was passed, or None.

    A `Measure` with two categories is its only relation, which is what lets
    the common case never mention the word.
    """
    if isinstance(item, Relation):
        return item
    if isinstance(item, Measure):
        return item.relate()
    return None


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------


class Run:
    """A run under construction. Every method returns `self`, so calls chain.

    The object accumulates lines of a script and nothing else: it holds no
    connection and executes nothing until `run()`.

    Statements come out in dependency order rather than call order. The
    source, the measure table and `lib/filters.sql` are assembled at the top of
    the script when `sql()` renders it, so a measure first mentioned by the
    last format call still reaches the declarations block, where the filter
    library needs it.
    """

    def __init__(self, hits=None, *, read=READ, block=None, block_sep=".",
                 include=(SQL_DIR,), **values):
        #: Script lines, in order, exactly as they will be written.
        self._lines: list[str] = []
        #: Session variables the last format call set, so the next one can
        #: clear the ones it does not set. See `_assign` on why only formats.
        self._volatile: set[str] = set()
        self._include = tuple(Path(d) for d in include)
        self._hits = hits
        self._read = read
        self._block = block
        self._block_sep = block_sep
        self._where: list[str] = []
        self._measures: dict[str, Measure] = {}
        self._read_filters: list[Fails] = []
        self._block_filters: list[Fails] = []
        self._relations: dict[str, "Relation"] = {}
        self._site_filters: list[Fails] = []
        #: Lines that precede the source, for variables the source itself reads.
        self._head: list[str] = []
        #: Library files already read, so a requirement is met once.
        self._provided: set[str] = set()
        if values:
            self._head.append("\n".join(self._assign(values)))

    # -- the declarations ------------------------------------------------

    def source(self, hits, *, read=READ, block=None, block_sep=".") -> "Run":
        """The hit tables, what identifies one read, and what groups them.

        `read` is the expression a read-level filter judges -- the read, not
        the alignment record. `block` is the partition over reads that the
        outputs are split by: a cell barcode, a spatial spot, a library.
        `block=None` puts every read in one block, which is a bulk run, and
        every level below works unchanged.

        Both take several columns. A `read` becomes a struct, because nothing
        ever prints it -- it is only compared and counted. A `block` becomes
        the columns joined by `block_sep`, because a block *is* printed: it is
        the barcode column of a cellInfo file and the dataset name inside an
        Amethyst H5, and a struct there would be unreadable.

            block=["condition", "XB"]     -->  ctrl.AACGTT
        """
        self._refuse_if_resolved("source()")
        self._hits, self._read = hits, read
        self._block, self._block_sep = block, block_sep
        return self

    def where(self, *conditions: str) -> "Run":
        """The hit-level filter, applied before anything is counted.

        Conditions in one call are OR'd; separate calls are AND'd. So a list of
        independent cuts is a list of calls, and the one place you write a
        boolean yourself is where you actually meant one:

            run.where("qual >= 20")
            run.where("mapq >= 30")
            run.where("read_5p >= 10", "read_3p >= 2")   # either end is enough

        This is where per-position artefacts go: end repair sits a fixed
        distance from the 3' end and random priming the same distance from the
        5' end, both measured on the read as sequenced, including hard-clipped
        bases.
        """
        self._refuse_if_resolved("where()")
        if not conditions:
            return self
        joined = " OR ".join(conditions)
        self._where.append(f"({joined})")
        return self

    def measures(self, *items) -> "Run":
        """Declare measures or relations the chain does not otherwise mention.

        Rarely needed: anything used in a filter or handed to a format declares
        itself. This is for the run whose measures are found by name rather
        than passed, such as `cgmap`, which collects every relation named for a
        trinucleotide.
        """
        for item in items:
            self._register(item)
        return self

    def drop_reads(self, *filters) -> "Run":
        """Reads to throw away. Each argument is a reason; they compose by union."""
        for item in filters:
            self._read_filters.append(self._filter(item, "read"))
        return self

    def drop_blocks(self, *filters) -> "Run":
        """Blocks to throw away, for the same reason and in the same way."""
        for item in filters:
            self._block_filters.append(self._filter(item, "block"))
        return self

    def drop_sites(self, *filters) -> "Run":
        """Positions to throw away, one block at a time.

        A site is one position in one block, so this judges a position per
        block: "this cell did not see this CG often enough" rather than
        "nobody did". A filter meant to apply everywhere, such as `in_regions(...)`
        over a blacklist, rejects the position in every block.

        Not the same thing as a format's `min_depth`, which is what one file
        shows. These are what the run keeps, and they appear in `report()`
        beside every other filter.
        """
        if self._has("CREATE OR REPLACE VIEW site_relation"):
            raise ScriptError(
                "drop_sites() has to come before anything that reads sites, "
                "and this run has already decided which ones survive"
            )
        for item in filters:
            self._site_filters.append(self._filter(item, "site"))
        return self

    def resolve(self) -> "Run":
        """`resolved`: `src` less everything `drop_reads` and `drop_blocks` named.

        Called automatically by the first call that needs it, so a chain rarely
        says it.
        """
        if self._has("CREATE OR REPLACE VIEW resolved"):
            return self
        clauses = []
        if self._read_filters:
            clauses.append(_not_in("read", self._read_filters))
        if self._block_filters:
            clauses.append(_not_in("block", self._block_filters))
        body = "\n  AND ".join(clauses)
        where = f"\nWHERE {body}" if body else ""
        self._sql(f"CREATE OR REPLACE VIEW resolved AS\nSELECT * FROM src{where};")
        return self

    def report(self) -> "Run":
        """A QC query: reads seen, reads kept, and what each filter rejected.

        Rendered from the filters recorded above rather than written out a
        second time, which is the point. A hand-written script repeats each
        filter and its threshold here, and nothing makes the two copies agree.
        """
        self.resolve()
        parts = [
            "SELECT 'reads seen' AS what, count(DISTINCT read) AS n FROM src",
            "SELECT 'reads kept', count(DISTINCT read) FROM resolved",
        ]
        for f in self._read_filters + self._block_filters + self._site_filters:
            parts.append(f"SELECT {literal(f.label)}, count(*) FROM {f.sql}")
        self._sql("\nUNION ALL\n".join(parts) + ";")
        return self

    # -- reading files --------------------------------------------------

    def read(self, path: str, *args, **values) -> "Run":
        """`.read path`, with `values` set first and stale variables cleared.

        The general form. `run.sites()` is `run.read("lib/count_sites.sql")` and
        `run.bismark_cov(mCG, out=...)` is `run.read("formats/bismark_cov.sql",
        mCG, out=...)`, so this is the only method that does anything; the rest
        are spellings of it found by `__getattr__`.

        A `Measure` anywhere in the arguments is declared into the script and
        passed on as its name. One given positionally is the `measure`
        variable, which is what every format that takes one calls it.
        """
        self._provided.add(path)
        values = dict(values)
        for item in args:
            relation = as_relation(item)
            if relation is None:
                raise ScriptError(
                    f"{path}: a positional argument is the relation to write, "
                    f"not {type(item).__name__}; everything else is a keyword"
                )
            values["relation"] = relation
        for key, value in list(values.items()):
            if isinstance(value, (Measure, Relation)):
                self._register(value)
                values[key] = value.name
        self._satisfy(path)
        block = self._assign(values, sticky=not path.startswith("formats/"))
        block.append(f".read {path}")
        companion = _companion(path, self._include)
        if companion:
            block.append(f".python {companion}")
        self._lines.append("\n".join(block))
        self._restrict_sites()
        return self

    def python(self, path: str, **values) -> "Run":
        """`.python path`. Needed only for a format whose Python is not named
        after its SQL, since `read` follows that pairing on its own."""
        block = self._assign(values)
        block.append(f".python {path}")
        self._lines.append("\n".join(block))
        return self

    def sql_statement(self, text: str) -> "Run":
        """Any SQL, written out. The escape hatch: a chain never has to be
        abandoned because the library has no method for something."""
        self._sql(text if text.rstrip().endswith(";") else text + ";")
        return self

    def __getattr__(self, name: str):
        """`run.bismark_cov(...)` for `.read formats/bismark_cov.sql`.

        A name is looked up in `formats/` and then in `lib/`, across the
        include path, so a forked file shadows the shipped one here exactly as
        it does for a `.read` directive. A name that matches no file raises,
        naming where it looked, which is the closest thing to a list of methods
        this class can honestly give.
        """
        if name.startswith("_"):
            raise AttributeError(name)
        target = _find(name, self._include)
        if target is None:
            raise AttributeError(
                f"no method {name!r} and no {name}.sql in formats/ or lib/ on "
                f"the include path ({', '.join(str(d) for d in self._include)}). "
                f"`alnbase-agg lib --list` names the files."
            )

        def call(*args, **values) -> "Run":
            return self.read(target, *args, **values)

        return call

    # -- output ---------------------------------------------------------

    def sql(self) -> str:
        """The script, as text. An ordinary one: `alnbase-agg run` executes it."""
        self.resolve()
        return "\n\n".join(self._head + self._prelude() + self._lines) + "\n"

    def run(self, con=None, *, show=print, trace=None):
        """Write the script to a temporary file and execute it.

        Through `script.run`, so a chained run and a hand-written one take the
        same path through the same code, and an error is reported the same way
        -- quoted with its statement, against a file that exists on disk while
        the run lasts.
        """
        import tempfile

        from alnbase_agg import script

        if con is None:
            import duckdb

            con = duckdb.connect()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "run.sql"
            path.write_text(self.sql())
            script.run(con, path, include=self._include, trace=trace, show=show)
        return con

    # -- internals ------------------------------------------------------

    def _sql(self, text: str) -> None:
        self._lines.append(text)

    def _has(self, prefix: str) -> bool:
        return any(line.startswith(prefix) for line in self._lines)

    def _refuse_if_resolved(self, what: str) -> None:
        if self._has("CREATE OR REPLACE VIEW resolved"):
            raise ScriptError(
                f"{what} has to come before anything that reads `src`, and this "
                "run has already decided which observations survive"
            )

    def _filter(self, item, level: str) -> Fails:
        if isinstance(item, (Test, _Both)):
            item = item.at(level)
        if isinstance(item, str):
            item = fails(item, level=level)
        if not isinstance(item, Fails):
            raise ScriptError(
                f"a filter is a quantity compared against something -- "
                f"m.total.lt(3) -- or a string of SQL, not "
                f"{type(item).__name__}"
            )
        if item.level != level:
            verb = {"read": "drop_reads", "block": "drop_blocks",
                    "site": "drop_sites"}[level]
            raise ScriptError(
                f"{verb}() was given a {item.level}-level filter ({item.label}); "
                f"it judges {level}s, and the two never match"
            )
        if item.relation is not None:
            self._register(item.relation)
        elif item.measure is not None:
            self._register(item.measure)
        return item

    def _satisfy(self, path: str) -> None:
        """Read whatever `path` declares it needs, first and once.

        This is why a chain says nothing about counting sites or loading contig
        order: the format that needs them says so in its own header, so the
        line never has to appear in a run that could not do without it.
        """
        for need in _requires(path, self._include):
            if need == "resolved":
                self.resolve()
            elif need not in self._provided:
                self.read(need)

    def _restrict_sites(self) -> None:
        """`site_relation`: `site_all` less everything `drop_sites` named.

        Emitted once the file that defines `site_all` has been read and before
        anything that reads sites, which is the same shape as `resolved` one
        level down. Redefining the view rather than passing an argument is what
        lets every format follow without being told.
        """
        if not self._site_filters or self._has("CREATE OR REPLACE VIEW site_relation"):
            return
        if "lib/count_sites.sql" not in self._provided:
            return
        self._sql(
            "CREATE OR REPLACE VIEW site_relation AS\nSELECT * FROM site_all\n"
            "WHERE " + _not_in("site", self._site_filters, key="site_key(block, "
                               "contig, refr_pos, strand)") + ";"
        )

    def _register(self, item) -> None:
        """Declare a measure, or a relation and the measure behind it.

        A two-category measure also declares its only relation, which takes
        the measure's own name. That is what lets a chain hand a measure
        straight to a format and never mention the word relation.
        """
        if isinstance(item, Relation):
            self._register(item.measure)
            existing = self._relations.get(item.name)
            if existing is not None and (
                existing.main_categories,
                existing.other_categories,
            ) != (item.main_categories, item.other_categories):
                raise ScriptError(
                    f"two different relations are both called {item.name!r}; a "
                    "relation's name is what every filter and output file calls "
                    "it, so it has to mean one thing"
                )
            self._relations.setdefault(item.name, item)
            return

        existing = self._measures.get(item.name)
        if existing is not None and existing.categories != item.categories:
            raise ScriptError(
                f"two different measures are both called {item.name!r}; a "
                "measure's name is what the script and every output file call "
                "it, so it has to mean one thing"
            )
        if existing is None:
            self._measures[item.name] = item
        if len(item.categories) == 2:
            self._relations.setdefault(item.name, item.relate())

    def _prelude(self) -> list[str]:
        """`src`, the measure table and the filter library, as text.

        Built when the script is rendered rather than when `source()` is
        called, so that a `where` or a measure first mentioned halfway down the
        chain still reaches the top of the script, where the library needs it.
        """
        blocks = []
        if self._hits is None:
            raise ScriptError(
                "this run has no hits; pass them to Run(...) or call source()"
            )
        fields = {self._read: self._read} if isinstance(self._read, str) else self._read
        if not isinstance(fields, dict):
            fields = {k: k for k in fields}
        struct = ", ".join(f"{literal(k)}: {v}" for k, v in fields.items())
        source = self._hits.sql if isinstance(self._hits, Var) else literal(self._hits)
        block_key = _block_expression(self._block, self._block_sep)
        clause = "\nWHERE " + "\n  AND ".join(self._where) if self._where else ""
        aliases = "".join(
            f"       {column} AS {alias},\n" for alias, column in ALIASES.items()
        )
        blocks.append(
            f"CREATE OR REPLACE VIEW src AS\n"
            f"SELECT *,\n"
            f"{aliases}"
            f"       {{{struct}}} AS read,\n"
            f"       {block_key} AS block\n"
            f"FROM read_parquet({source}){clause};"
        )
        if self._measures:
            rows = [row for m in self._measures.values() for row in m.rows()]
            blocks.append(
                "CREATE OR REPLACE VIEW measure (measure, name, category) AS VALUES\n"
                "  " + ",\n  ".join(rows) + ";"
            )
            rows = [r.row() for r in self._relations.values()]
            blocks.append(
                "CREATE OR REPLACE VIEW relation (relation, measure, main, other) "
                "AS VALUES\n  " + ",\n  ".join(rows) + ";"
            )
            blocks.append(".read lib/filters.sql")
        return blocks

    def _assign(self, values: dict, *, sticky: bool = True) -> list[str]:
        """Set these variables, and clear what a previous format call set.

        Clearing is what makes a format call independent of the calls before
        it. Without it `min_depth`, passed to one output, silently applies to
        every later output in the script as well -- a hand-written script has
        no defence against this at all, because a session variable is global
        and stays set.

        Only what a *format* set is cleared, which is the part worth knowing.
        A format consumes its variables at once: it reads them, writes a file,
        and is done. A `lib/` file does not -- it defines a view whose body
        calls `getvariable`, and a view is resolved every time it is queried,
        so `merge`, set for `lib/merge_strands.sql`, is read again by every format
        that follows and must still be set when they run. Hence `sticky`:
        variables set for a library file persist, variables set for a format
        do not.
        """
        lines = []
        if not sticky:
            for name in sorted(self._volatile - set(values)):
                lines.append(f"RESET VARIABLE {name};")
            self._volatile = set(values)
        for name, value in values.items():
            lines.append(f"SET VARIABLE {name} = {literal(value)};")
        return lines


def _block_expression(block, sep: str) -> str:
    """The SQL for the block key.

    Several columns are joined rather than made into a struct, because a block
    is written out: it names a dataset inside an Amethyst H5 and fills the
    barcode column of a cellInfo file. `concat_ws` skips a NULL rather than
    poisoning the whole key, which means a row missing one part lands in the
    block named by the rest; a run that would rather see that fail writes the
    expression itself, which is always accepted.
    """
    if block is None:
        return "'bulk'"
    if isinstance(block, str):
        return block
    columns = list(block)
    if not columns:
        return "'bulk'"
    if len(columns) == 1:
        return columns[0]
    return f"concat_ws({literal(sep)}, {', '.join(columns)})"


def _not_in(column: str, filters: list[Fails], key: str | None = None) -> str:
    """`key NOT IN (everything the filters reject)`.

    `key` differs from `column` only for sites, where what identifies the row
    is a struct built from four columns rather than one column's value.
    """
    branches = f"\n        UNION ALL SELECT {column} FROM ".join(f.sql for f in filters)
    return f"{key or column} NOT IN (\n        SELECT {column} FROM {branches})"


def _requires(path: str, include) -> list[str]:
    """What a file says it needs, read off the `REQUIRES` line in its header.

    Declared rather than inferred. A format that reads `site_relation` needs
    `lib/count_sites.sql` to have run, and saying so in the file means a person
    reading it learns the same thing the builder does -- and that a forked file
    with different needs is believed about them.

    The one entry that is not a path is `resolved`, which no file defines: it
    is built from the run's own filters.

    Every `REQUIRES` line in the header counts, not only the first. A file that
    splits a long list over two lines would otherwise have half of it honoured
    and nothing said about the rest, which is the silent failure this whole
    library is built to avoid.
    """
    for directory in include:
        candidate = Path(directory) / path
        if not candidate.exists():
            continue
        needs: list[str] = []
        for line in candidate.read_text().splitlines():
            if not line.startswith("--"):
                break                       # past the header
            text = line[2:].strip()
            if text.startswith("REQUIRES "):
                needs += [p.strip()
                          for p in text[len("REQUIRES "):].split(",") if p.strip()]
        return list(dict.fromkeys(needs))
    return []


def _find(name: str, include) -> str | None:
    """The library-relative path for a method name, or None."""
    for directory in ("formats", "lib"):
        relative = f"{directory}/{name}.sql"
        if any((Path(d) / relative).exists() for d in include):
            return relative
    return None


def _companion(path: str, include) -> str | None:
    """The `.python` half of a format, if it has one.

    Paired by name: `formats/amethyst_h5.sql` is finished by
    `formats/amethyst_h5.py`. Naming the pair here means a format that grows a
    Python half does not have to be called differently.
    """
    if not path.endswith(".sql"):
        return None
    candidate = path[: -len(".sql")] + ".py"
    if any((Path(d) / candidate).exists() for d in include):
        return candidate
    return None
