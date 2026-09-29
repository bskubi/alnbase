"""The chained builder, against the script it is supposed to be a spelling of.

The load-bearing test is `the_two_front_ends_write_the_same_files`: a chain and
a hand-written script producing the same bytes is the whole claim, because it
says the builder added no behaviour of its own. Everything else here pins one
rendering decision each.
"""

from __future__ import annotations

import pytest

from alnbase_agg.pipeline import (Measure, Run, at_position, both, fails,
                                  in_regions, listed, literal, outside_regions,
                                  var)
from alnbase_agg.script import ScriptError

mCG = Measure("mCG", methylated="CG", unmethylated="TG")
mCH = Measure("mCH", methylated=["CHG", "CHH"], unmethylated=["THG", "THH"])


def build(hits, exclude, fai, sample, min_total=2):
    """The example run, as a chain. Kept beside the SQL below, which is it."""
    run = Run(hits, block="XB", fai=fai)
    run.where("qual >= 20")
    run.where("mapq >= 30")
    run.where("off_5p >= 10")
    run.where("off_3p >= 2")
    run.drop_reads(mCH.ratio.gt(0.4), mCH.consecutive.gt(3))
    run.drop_reads(listed(exclude))
    run.drop_blocks(mCG.total.lt(min_total))
    run.merge_strands(merge=mCG)
    run.bismark_cov(mCG, out=f"{sample}.CG.cov")
    run.bismark_cov(mCH, out=f"{sample}.CH.cov")
    return run


SCRIPT = """
CREATE OR REPLACE VIEW src AS
SELECT *, ref_name AS contig,
       {{'qname': qname, 'mate': is_first_in_template}} AS read, XB AS block
FROM read_parquet('{hits}')
WHERE (qual >= 20) AND (mapq >= 30) AND (off_5p >= 10) AND (off_3p >= 2);

CREATE OR REPLACE VIEW measure (measure, name, category) AS VALUES
  ('mCG', 'CG', 'methylated'), ('mCG', 'TG', 'unmethylated'),
  ('mCH', 'CHG', 'methylated'), ('mCH', 'CHH', 'methylated'),
  ('mCH', 'THG', 'unmethylated'), ('mCH', 'THH', 'unmethylated');

CREATE OR REPLACE VIEW relation (relation, measure, main, other) AS VALUES
  ('mCG', 'mCG', ['methylated'], ['unmethylated']),
  ('mCH', 'mCH', ['methylated'], ['unmethylated']);

.read lib/filters.sql

CREATE OR REPLACE VIEW resolved AS SELECT * FROM src
WHERE read NOT IN (SELECT read FROM read_fails('mCH', 'ratio', '>', 0.4)
                   UNION ALL SELECT read FROM read_fails_consecutive('mCH', '>', 3)
                   UNION ALL SELECT read FROM read_fails_listed('{exclude}'))
  AND block NOT IN (SELECT block FROM block_fails('mCG', 'total', '<', 2));

.read lib/count_sites.sql
SET VARIABLE merge = 'mCG';
.read lib/merge_strands.sql
SET VARIABLE fai = '{fai}';
.read lib/load_contigs.sql

SET VARIABLE relation = 'mCG';
SET VARIABLE out = '{sample}.CG.cov';
.read formats/bismark_cov.sql
SET VARIABLE relation = 'mCH';
SET VARIABLE out = '{sample}.CH.cov';
.read formats/bismark_cov.sql
"""


def test_the_two_front_ends_write_the_same_files(
    tmp_path, hits, exclude, fai, con, write
):
    """The claim the module exists to make: a chain is a way of writing the
    script, not a second implementation of it."""
    import duckdb

    from alnbase_agg import script

    hand = tmp_path / "hand"
    chain = tmp_path / "chain"
    hand.mkdir()
    chain.mkdir()

    path = write(
        SCRIPT.format(hits=hits, exclude=exclude, fai=fai, sample=hand / "s")
    )
    script.run(duckdb.connect(), path, show=None)

    build(hits, exclude, fai, str(chain / "s")).run(con, show=None)

    written = sorted(p.name for p in hand.iterdir())
    assert written == ["s.CG.cov", "s.CH.cov"]
    for name in written:
        assert (hand / name).read_bytes() == (chain / name).read_bytes(), name


# --- the hit-level filter --------------------------------------------------


def test_conditions_in_one_call_are_or_and_separate_calls_are_and(hits):
    """The reason `where` takes several: a list of independent cuts is a list
    of calls, and you write a boolean yourself only where you meant one."""
    run = Run(hits, block="XB")
    run.where("qual >= 20")
    run.where("off_5p >= 10", "off_3p >= 2")
    text = run.sql()
    assert "WHERE (qual >= 20)\n  AND (off_5p >= 10 OR off_3p >= 2);" in text


def test_a_run_with_no_where_has_no_where_clause(hits):
    assert "WHERE" not in Run(hits).sql().split("CREATE OR REPLACE VIEW resolved")[0]


def test_a_where_after_the_filters_are_fixed_is_refused(hits):
    """A measure may arrive at any point, because the declarations block is
    assembled last. A `where` may not, once `resolved` has said what survives:
    it would read as a later filter and act as an earlier one."""
    run = Run(hits, block="XB")
    run.drop_reads(mCG.total.lt(2))
    run.bismark_cov(mCG, out="/dev/null")
    with pytest.raises(ScriptError) as caught:
        run.where("qual >= 20")
    assert "before anything that reads `src`" in str(caught.value)


def test_a_measure_first_used_by_the_last_call_still_reaches_the_top(hits):
    """The reason the script is assembled at the end rather than as it goes."""
    run = Run(hits, block="XB")
    run.bismark_cov(mCG, out="/dev/null")
    text = run.sql()
    assert text.index("VIEW measure") < text.index("VIEW resolved")


def test_the_read_key_defaults_to_qname_and_mate(hits):
    """Not the alignment record: a chimeric read is several records."""
    assert "{'qname': qname, 'mate': is_first_in_template} AS read" in Run(hits).sql()


def test_a_run_with_no_block_is_bulk(hits):
    assert "'bulk' AS block" in Run(hits).sql()


# --- measures --------------------------------------------------------------


def test_a_measure_declares_itself_when_it_is_first_used(hits):
    """There is no registration step, so a filter on an undeclared measure is
    not a state the script can reach."""
    run = Run(hits, block="XB")
    run.drop_reads(mCH.ratio.gt(0.4))
    text = run.sql()
    assert "('mCH', 'CHG', 'methylated')" in text
    assert text.index("VIEW measure") < text.index("lib/filters.sql")
    assert text.index("lib/filters.sql") < text.index("VIEW resolved")


def test_a_two_category_measure_is_its_own_relation():
    """The common case never says the word. A measure with two categories has
    one possible relation, and it takes the measure's own name."""
    assert mCH.ratio.gt(0.4).at("read").sql == "read_fails('mCH', 'ratio', '>', 0.4)"
    assert mCH.relate().name == "mCH"
    assert mCH.relate().main_categories == ["methylated"]
    assert mCH.relate().other_categories == ["unmethylated"]


def test_a_filter_is_a_grouping_a_quantity_and_a_condition():
    """The same test rendered against three groupings. "Blocks whose total is
    below 1000" needs no macro of its own."""
    test = mCG.total.lt(1000)
    assert test.at("read").sql == "read_fails('mCG', 'total', '<', 1000)"
    assert test.at("block").sql == "block_fails('mCG', 'total', '<', 1000)"
    assert test.at("site").sql == "site_fails('mCG', 'total', '<', 1000)"


def test_every_comparison_is_spelled_out():
    """`> 0.4` and `>= 0.4` are different filters, and premethyst's read cut
    uses the second; with the direction in a macro name it was unaskable."""
    for method, op in [("gt", ">"), ("lt", "<"), ("ge", ">="),
                       ("le", "<="), ("eq", "="), ("ne", "<>")]:
        rendered = getattr(mCH.ratio, method)(0.4).at("read").sql
        assert rendered == f"read_fails('mCH', 'ratio', '{op}', 0.4)"


def test_a_quantity_that_is_not_a_plain_aggregate_keeps_its_own_macro():
    """A run along the read, and a count across blocks: neither is an
    aggregate over one group's own totals, so neither fits the uniform shape
    and each is refused at the levels where it means nothing."""
    assert mCH.consecutive.gt(3).at("read").sql == \
        "read_fails_consecutive('mCH', '>', 3)"
    with pytest.raises(ScriptError) as caught:
        mCH.consecutive.gt(3).at("block")
    assert "not a block-level test" in str(caught.value)

    assert mCG.blocks.lt(3).at("site").sql == "site_fails_blocks('mCG', '<', 3)"
    with pytest.raises(ScriptError):
        mCG.blocks.lt(3).at("read")


def test_a_measure_with_more_categories_has_to_say_which_relation():
    snp = Measure("snp", A="A>A", C="A>C", G="A>G", T="A>T")
    with pytest.raises(ScriptError) as caught:
        snp.ratio.gt(0.9)
    assert "no single relation to mean" in str(caught.value)
    assert snp.relate("G", "A").name == "snp.G_vs_A"


def test_a_relation_is_two_disjoint_sides():
    """`other` never repeats `main`, exactly as the categories do not overlap;
    the denominator is their sum."""
    snp = Measure("snp", A="A>A", C="A>C", G="A>G", T="A>T")
    assert snp.relate("G", "A").other_categories == ["A"]
    assert sorted(snp.relate("G").other_categories) == ["A", "C", "T"]
    assert snp.relate(["G", "A"], "C").main_categories == ["G", "A"]
    with pytest.raises(ScriptError) as caught:
        snp.relate("G", ["G", "A"])
    assert "both sides" in str(caught.value)


def test_a_category_a_measure_does_not_have_is_refused():
    snp = Measure("snp", A="A>A", C="A>C")
    with pytest.raises(ScriptError) as caught:
        snp.relate("G", "A")
    assert "no category 'G'" in str(caught.value)


def test_a_query_name_in_two_categories_of_one_measure_is_refused():
    """Within a measure the categories partition the hits, so a denominator can
    be a plain sum. Two independent readouts are two measures instead."""
    with pytest.raises(ScriptError) as caught:
        Measure("m", methylated=["CG"], unmethylated=["CG", "TG"])
    assert "at most one category" in str(caught.value)


def test_a_query_name_may_belong_to_several_measures(hits):
    """NOMe-seq: accessibility and methylation are separate readouts over
    separate queries, and a name shared between two measures is fine."""
    a = Measure("acc", open=["GCH"], closed=["GTH"])
    b = Measure("meth", methylated=["GCH"], unmethylated=["GTH"])
    run = Run(hits, block="XB")
    run.drop_reads(a.ratio.gt(0.9), b.ratio.gt(0.9))
    text = run.sql()
    assert "('acc', 'GCH', 'open')" in text
    assert "('meth', 'GCH', 'methylated')" in text


def test_a_measure_with_many_categories_reaches_a_format_through_a_relation(hits):
    snp = Measure("snp", A="A>A", C="A>C", G="A>G", T="A>T")
    run = Run(hits, block="XB")
    run.bismark_cov(snp.relate("G", "A"), out="/dev/null")
    text = run.sql()
    assert "('snp.G_vs_A', 'snp', ['G'], ['A'])" in text
    assert "SET VARIABLE relation = 'snp.G_vs_A';" in text


def test_one_name_cannot_mean_two_measures(hits):
    run = Run(hits, block="XB")
    run.drop_reads(mCG.total.lt(3))
    with pytest.raises(ScriptError) as caught:
        run.drop_reads(Measure("mCG", methylated="CA", unmethylated="TA").total.lt(3))
    assert "has to mean one thing" in str(caught.value)


def test_a_forked_macro_is_reachable_without_a_method():
    """Code intelligence is worth having, but not at the price of a library
    you cannot extend without editing Python. The level still comes from the
    `drop_*` call, as it does for everything else."""
    assert mCG.wobbly(7).at("read").sql == "read_fails_wobbly('mCG', 7)"
    assert mCG.wobbly(7).at("block").sql == "block_fails_wobbly('mCG', 7)"


# --- filters ---------------------------------------------------------------


def test_a_site_only_quantity_is_refused_by_drop_reads(hits):
    """A count across blocks means nothing about one read, and saying so here
    beats a filter that binds and quietly drops nothing."""
    run = Run(hits, block="XB")
    with pytest.raises(ScriptError) as caught:
        run.drop_reads(mCG.blocks.lt(3))
    assert "not a read-level test" in str(caught.value)


def test_a_filter_may_be_written_as_sql(hits):
    """The escape hatch. A ladder is a CASE, and must stay writable."""
    f = fails("(SELECT read FROM read_relation WHERE den < 3)", label="too few")
    run = Run(hits, block="XB")
    run.drop_reads(f)
    assert "WHERE den < 3" in run.sql()


def test_a_filter_of_the_wrong_kind_is_refused(hits):
    run = Run(hits, block="XB")
    with pytest.raises(ScriptError):
        run.drop_reads(42)


def test_the_report_names_every_filter_that_was_applied(hits, exclude, fai, con):
    """`report()` is rendered from the recorded filters, so a filter cannot be
    applied without appearing in the report."""
    run = build(hits, exclude, fai, "s")
    run.report()
    text = run.sql()
    for label in ("read mCH ratio > 0.4", "read mCH consecutive > 3",
                  "read listed", "block mCG total < 2"):
        assert f"'{label}'" in text


# --- reading files ---------------------------------------------------------


def test_a_format_argument_does_not_reach_the_next_format(hits):
    """`min_depth`, given to one output, must not filter the next one.

    A hand-written script has no defence here, because a session variable is
    global and stays set. This is the one behaviour the builder adds.
    """
    run = Run(hits, block="XB")
    run.bismark_cov(mCG, out="/dev/null", min_depth=99)
    run.bismark_cov(mCG, out="/dev/null")
    assert "RESET VARIABLE min_depth;" in run.sql()


def test_a_library_argument_does_reach_the_next_format(hits):
    """The other half of that rule. `merge` is read by a *view* body, so it is
    read again by every format that follows and must stay set."""
    run = Run(hits, block="XB")
    run.merge_strands(merge=mCG)
    run.bismark_cov(mCG, out="/dev/null")
    assert "RESET VARIABLE merge;" not in run.sql()


def test_a_format_carries_its_python_half(hits, fai):
    """`amethyst_h5` is two files paired by name, and calling it reads both."""
    run = Run(hits, block="XB")
    run.amethyst_h5(mCG, context="CG", out="/dev/null")
    text = run.sql()
    assert ".read formats/amethyst_h5.sql" in text
    assert ".python formats/amethyst_h5.py" in text


def test_a_method_with_no_file_behind_it_says_where_it_looked(hits):
    run = Run(hits)
    with pytest.raises(AttributeError) as caught:
        run.bismarck_cov
    assert "bismarck_cov" in str(caught.value)
    assert "formats/" in str(caught.value)


def test_a_run_with_no_filters_still_defines_resolved(hits):
    """`lib/count_sites.sql` reads `resolved`, so it must exist even when nothing is
    dropped; the view is emitted by whichever call first needs it."""
    run = Run(hits, block="XB")
    run.drop_reads(mCG.total.lt(2))
    run.bismark_cov(mCG, out="/dev/null")
    text = run.sql()
    assert text.index("CREATE OR REPLACE VIEW resolved") < text.index("lib/count_sites.sql")


def test_every_requires_line_in_a_header_is_read(tmp_path, hits):
    """Not only the first. A forked file may spell a long list over two lines,
    and honouring half of it would be silent: the missing include is a view
    that is simply not there when the format asks for it."""
    from alnbase_agg.pipeline import _requires

    (tmp_path / "formats").mkdir()
    (tmp_path / "formats" / "two.sql").write_text(
        "-- A format that declares its needs twice.\n"
        "-- REQUIRES lib/count_sites.sql, lib/load_contigs.sql\n"
        "-- REQUIRES lib/merge_strands.sql, lib/count_sites.sql\n"
        "SELECT 1;\n"
    )
    assert _requires("formats/two.sql", [tmp_path]) == [
        "lib/count_sites.sql", "lib/load_contigs.sql", "lib/merge_strands.sql"
    ]


def test_a_requires_line_after_the_header_is_not_one(tmp_path):
    """The header ends at the first line that is not a comment, so a `REQUIRES`
    inside the body is prose about the file rather than a declaration."""
    from alnbase_agg.pipeline import _requires

    (tmp_path / "formats").mkdir()
    (tmp_path / "formats" / "late.sql").write_text(
        "-- REQUIRES lib/count_sites.sql\n"
        "SELECT 1;\n"
        "-- REQUIRES lib/load_contigs.sql\n"
    )
    assert _requires("formats/late.sql", [tmp_path]) == ["lib/count_sites.sql"]


# --- rendering -------------------------------------------------------------


@pytest.mark.parametrize(
    "value,expected",
    [
        ("plain", "'plain'"),
        ("it's", "'it''s'"),
        (5, "5"),
        (0.4, "0.4"),
        (True, "true"),
        (None, "NULL"),
        (var("block"), "getvariable('block')"),
        (var("block") + ".cov", "(getvariable('block') || '.cov')"),
    ],
)
def test_a_value_renders_as_the_literal_for_it(value, expected):
    assert literal(value) == expected


def test_the_script_is_an_ordinary_one(hits, exclude, fai, tmp_path, con):
    """Generated, saved, and run by the ordinary runner with no Python left.

    This is what makes the builder removable rather than a second toolchain.
    """
    import duckdb

    from alnbase_agg import script

    out = tmp_path / "s"
    path = tmp_path / "generated.sql"
    path.write_text(build(hits, exclude, fai, str(out)).sql())
    script.run(duckdb.connect(), path, show=None)
    assert (tmp_path / "s.CG.cov").exists()


def test_a_many_category_measure_filters_against_real_rows(hits, con):
    """A relation drawn from the middle of a measure, executed rather than
    rendered: CHH against CHG, which a two-outcome model cannot express.

    Per read that is rA 1/1, rB 1/2, rC 3/4, rD 1/1; rE and rF have no such
    hits at all, so they have no row and no ratio.
    """
    ctx = Measure("ctx", hh=["CHH"], hg=["CHG"], converted=["THH", "THG"])
    run = Run(hits, block="XB")
    run.drop_reads(ctx.relate("hh", "hg").ratio.gt(0.6))
    run.run(con, show=None)
    failed = con.execute(
        "SELECT read.qname FROM read_fails('ctx.hh_vs_hg', 'ratio', '>', 0.6)"
    ).fetchall()
    assert sorted(q for (q,) in failed) == ["rA", "rC", "rD"]


def test_sites_can_be_dropped(hits, con):
    """The level the run had no filters for: a position in a block, judged on
    its own coverage. Redefining the view is what makes every format follow."""
    run = Run(hits, block="XB")
    run.drop_sites(mCG.total.lt(2))
    run.bismark_cov(mCG, out="/dev/null")
    text = run.sql()
    assert "CREATE OR REPLACE VIEW site_relation AS" in text
    assert "site_fails('mCG', 'total', '<', 2)" in text
    assert text.index("lib/count_sites.sql") < text.index("VIEW site_relation")
    run.run(con, show=None)
    kept = con.execute("SELECT count(*) FROM site_relation WHERE relation = 'mCG'")
    assert kept.fetchone()[0] < con.execute(
        "SELECT count(*) FROM site_all WHERE relation = 'mCG'").fetchone()[0]


def test_a_bed_file_drops_sites_everywhere(hits, con, tmp_path):
    """A blacklist is a site filter that judges the position and rejects it in
    every block."""
    bed = tmp_path / "black.bed"
    bed.write_text("chr2\t200\t205\n")
    run = Run(hits, block="XB")
    run.drop_sites(in_regions(str(bed)))
    run.bismark_cov(mCG, out="/dev/null")
    assert f"site_fails_in_regions('{bed}')" in run.sql()
    run.run(con, show=None)
    left = con.execute(
        "SELECT count(*) FROM site_relation WHERE contig = 'chr2' AND refr_pos = 200"
    ).fetchone()[0]
    assert left == 0


def test_a_total_threshold_no_row_can_meet_is_refused():
    """A read has a row only because it had a hit, so a total is never below 1.
    `total.lt(1)` would bind, run and match nothing -- the exact failure the
    library works hardest to prevent -- and what it means is `missing()`."""
    for test in (lambda: mCG.total.lt(1), lambda: mCG.total.le(0),
                 lambda: mCG.total.eq(0)):
        with pytest.raises(ScriptError) as caught:
            test()
        assert "missing()" in str(caught.value)
    mCG.total.lt(2)        # reachable, and must stay allowed
    mCG.total.ge(0)


def test_missing_catches_what_no_threshold_can(hits, con):
    """rC covers no CG position at all, so it has no row in read_relation and
    every threshold passes it by default."""
    run = Run(hits, block="XB")
    run.drop_reads(mCG.missing())
    run.run(con, show=None)
    absent = con.execute("SELECT read.qname FROM read_fails_missing('mCG')").fetchall()
    assert [q for (q,) in absent] == ["rC"]
    # rE and rF have a total of 1, so the threshold catches them -- but never rC,
    # whatever the number, because rC has no row for it to compare.
    thin = con.execute(
        "SELECT read.qname FROM read_fails('mCG', 'total', '<', 2)"
    ).fetchall()
    assert "rC" not in [q for (q,) in thin]


# --- the filters the target formats need -----------------------------------
#
# Each of these closes a gap found by going through what Bismark, MethylDackel,
# premethyst, ScaleMethyl, Amethyst, scbs and the hyper-editing pipeline
# actually filter on, and checking it against the syntax.


def test_a_block_can_be_judged_on_reads_and_on_positions(hits, con):
    """`unique_reads >= N` and a minimum of covered CG *positions* are the two
    commonest block cuts in this ecosystem, and neither is a hit count."""
    run = Run(hits, block="XB")
    run.drop_blocks(mCG.reads.lt(2))
    run.bismark_cov(mCG, out="/dev/null")
    run.run(con, show=None)
    extent = con.execute(
        "SELECT block, reads, sites FROM block_extent WHERE relation = 'mCG' "
        "ORDER BY 1"
    ).fetchall()
    # cellA is rA and rB; cellB is rC alone, because rD has no CG at all.
    assert extent == [("cellA", 2, 2), ("cellB", 1, 2), ("cellC", 2, 2)]
    assert mCG.sites.lt(3).at("block").sql == "block_fails_sites('mCG', '<', 3)"
    with pytest.raises(ScriptError):
        mCG.reads.lt(3).at("read")          # a read is one read


def test_a_bed_file_works_in_both_directions(hits, con, tmp_path):
    """A blacklist and a capture panel are different filters, and one is not
    reachable by inverting the other without complementing the BED by hand."""
    bed = tmp_path / "r.bed"
    bed.write_text("chr2\t200\t205\n")
    assert in_regions(str(bed)).at("site").sql == \
        f"site_fails_in_regions('{bed}')"
    assert outside_regions(str(bed)).at("site").sql == \
        f"site_fails_outside_regions('{bed}')"

    import duckdb

    def kept(test):
        c = duckdb.connect()
        run = Run(hits, block="XB")
        run.drop_sites(test)
        run.bismark_cov(mCG, out="/dev/null")
        run.run(c, show=None)
        return sorted({r[0] for r in c.execute(
            "SELECT refr_pos FROM site_relation WHERE relation = 'mCG'"
        ).fetchall()})

    inside, outside = kept(in_regions(str(bed))), kept(outside_regions(str(bed)))
    assert inside and outside and set(inside).isdisjoint(outside)


def test_a_position_is_vetoed_on_pooled_evidence(hits, con):
    """The SNP veto. The evidence is a different relation at the same position,
    and it is summed across blocks before it is compared -- testing block by
    block would veto a position because one cell saw one variant read."""
    other = Measure("other", hit="CHH", miss="THH")
    run = Run(hits, block="XB")
    run.drop_sites(at_position(other.ratio.gt(0.5)))
    run.bismark_cov(mCG, out="/dev/null")
    run.run(con, show=None)
    assert "site_fails_position('other', 'ratio', '>', 0.5)" in run.sql()


def test_both_is_the_other_connective(hits, con):
    """Filters compose by union, which is right. `both` is the conjunction, for
    a reason that is itself a conjunction -- scbs drops a site that is neither
    wholly methylated nor wholly unmethylated."""
    mixed = both(mCG.main.gt(0), mCG.other.gt(0))
    rendered = mixed.at("site").sql
    assert rendered.startswith("(SELECT site FROM site_fails('mCG', 'main', '>', 0)")
    assert " INTERSECT " in rendered
    assert mixed.at("read").sql.count("SELECT read FROM") == 2
    with pytest.raises(ScriptError):
        both(mCG.main.gt(0))                       # one test is not a conjunction
    with pytest.raises(ScriptError):
        both(mCG.main.gt(0), "not a test")


def test_both_makes_a_ladder_expressible(hits, con):
    """premethyst's read cut is a threshold that depends on how much evidence
    there is. It is two conjunctions unioned, and needs no mechanism of its own.

    What is lost with it: the recipe language this replaced required a ladder's
    rungs to tile the range exactly once. Here an overlap is harmless and a gap
    just means those reads are kept, and nothing says so.
    """
    run = Run(hits, block="XB")
    run.drop_reads(
        both(mCH.total.ge(3), mCH.total.le(5), mCH.main.gt(1)),
        both(mCH.total.gt(5), mCH.ratio.ge(0.4)),
    )
    run.bismark_cov(mCG, out="/dev/null")
    run.run(con, show=None)
    # rC has ten CH calls, four of them methylated: rung three, 0.4 >= 0.4.
    dropped = {r[0] for r in con.execute(
        "SELECT DISTINCT read.qname FROM src WHERE read NOT IN "
        "(SELECT read FROM resolved)"
    ).fetchall()}
    assert "rC" in dropped
