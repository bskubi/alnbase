"""Verifying the coordinates against the reference genome.

A coordinate convention is the one thing in this pipeline that can be wrong by
exactly one and produce output that looks perfect: every file writes, every
count is plausible, every plot has the right shape, and every cytosine is
reported one base to the right of where it is. Nothing downstream can catch it,
because nothing downstream knows where the cytosines were supposed to be.

The reference does. alnbase's hit rows carry `refr_base`, the reference base at
the column the row is anchored on, so the check needs no interpretation of what
a context means: take the position as an aggregation will report it, look the
base up in the FASTA, and see whether it is the base alnbase read there. An
off-by-one in a CG lands on the G, which is the wrong letter, so a few thousand
sampled positions settle the question.

This lives beside the aggregation rather than beside any one output format
because it is about hit rows and a FASTA, and touches neither the run script nor
HDF5. Every writer benefits from it for the same reason: a shifted coordinate is
wrong in a bedGraph exactly as it is in an Amethyst file.

Two checks, because two different things can be wrong and they want different
instruments:

*The reference is not the one the run used.* Checked exactly and for nothing:
the manifest records every `@SQ` contig with its length, so comparing that
against the FASTA index catches a wrong assembly, and catches `chr1` versus `1`
completely rather than statistically. No base is read.

*The positions are shifted.* Checked on a sample, and the sample is scored at
several candidate shifts rather than only at the one we intend to apply. A
report saying "0 of 2000 match as written, 2000 of 2000 match one base to the
left" names the bug; one saying only "0 of 2000 match" leaves the reader to
guess between a shift, the wrong assembly and a bug in this code.

One wrinkle deserves stating plainly, because getting it backwards would make
the check pass on shifted data. alnbase walks a read 5'->3' along its conversion
strand, so for a minus-strand read the walk runs backwards along the reference
*with both sides reverse complemented*: `refr_base` on such a row is the
complement of what the FASTA holds at `refr_pos`. The comparison therefore
complements the expected base for those rows, which it identifies from the hit
table's `strand` column -- the walk's own direction, not the FLAG's `is_reverse`,
which differs from it for exactly one mate of every pair. Without that column the
check accepts either orientation, which weakens it precisely where it matters
most: a CG's C and the G opposite it are complements, so a one-base shift within
a CG cannot be told from a correct call. The shift probe still has power, since
CH contexts are not palindromic, but the report says the check was weakened.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from .manifest import expand


def sql_string(text: str) -> str:
    """A SQL single-quoted literal, for the one place DDL forces one."""
    return "'" + str(text).replace("'", "''") + "'"

# Reported through the stdlib, because this package's dependencies are DuckDB and
# click and a logging framework is not worth a third. `alnbase-agg` attaches a
# handler; a library caller gets whatever its own configuration says.
logger = logging.getLogger("alnbase_agg.verify")

# Only the four unambiguous bases are compared. alnbase spells a column that is
# not a plain reference base with punctuation or a letter IUPAC leaves free (a
# pad past the read end, a gap, a skipped intron), and an ambiguity code in the
# reference is not evidence either way, so both are filtered out in SQL.
DNA = ("A", "C", "G", "T")
COMPLEMENT = str.maketrans("ACGT", "TGCA")

# Candidate shifts the sample is scored at, in the order they are reported.
# 0 is the shift this package intends; the rest are what a bug would look like.
SHIFTS = (0, -1, 1, -2, 2)

# Below this fraction of sampled bases matching, the run is refused.
DEFAULT_MIN_MATCH = 0.99


class VerifyError(Exception):
    """The reference and the hit tables do not agree about where the bases are."""


# ---------------------------------------------------------------------------
# the reference
# ---------------------------------------------------------------------------


class Reference:
    """A FASTA opened for single-base lookups, by 1-based position.

    pyfaidx rather than our own indexer: it is the mainstream library for this,
    it reads bgzipped FASTA, and a base lookup through it is one seek.
    """

    def __init__(self, path: str | Path):
        """Open the FASTA and read its index, without reading a base.

        pyfaidx is imported here rather than at the top of the module because it
        is the one dependency this package does not require. A user who never
        passes `--reference` never needs it, so the import failing has to be a
        message about that flag rather than a traceback on `import alnbase_agg`.

        The `Fasta` options are each load-bearing. `as_raw` returns a plain
        string, so a comparison is between two strings rather than between a
        string and a pyfaidx sequence object. `sequence_always_upper` keeps a
        soft-masked repeat from reading as a mismatch, which would otherwise show
        up as a scattering of failures in exactly the regions that are lowercase.
        `read_ahead` buffers, because the sample walks many positions on one
        contig.

        `lengths` is taken from the index, so constructing a `Reference` costs a
        `.fai` read and nothing else -- which is what lets the contig check run
        before any position is sampled.
        """
        try:
            from pyfaidx import Fasta
        except ImportError:
            raise VerifyError(
                "verifying against the reference needs pyfaidx, which is not "
                "installed. `pip install pyfaidx`, or omit --reference."
            ) from None

        self.path = Path(path)
        if not self.path.exists():
            raise VerifyError(f"no reference FASTA at {self.path}")
        try:
            # as_raw gives a plain string; sequence_always_upper keeps a
            # soft-masked (lowercase) repeat from reading as a mismatch.
            self.fasta = Fasta(
                str(self.path), as_raw=True, sequence_always_upper=True, read_ahead=4096
            )
        except Exception as exc:
            raise VerifyError(
                f"{self.path} could not be opened as a FASTA ({exc}). An index is "
                f"written beside it on first use, so the directory has to be "
                f"writable, or a .fai has to exist already."
            ) from None
        # From the index, not from the sequence: no base is read here.
        self.lengths: dict[str, int] = {
            name: len(self.fasta[name]) for name in self.fasta.keys()
        }

    def base(self, chrom: str, pos: int) -> str | None:
        """The reference base at a 1-based position, or None if there is none."""
        length = self.lengths.get(chrom)
        if length is None or pos < 1 or pos > length:
            return None
        return self.fasta[chrom][pos - 1 : pos] or None

    def close(self) -> None:
        """Release the FASTA's file handles. Failures are ignored.

        The check has already produced its verdict by the time anything closes
        this, so an error here could only replace a real answer with an unrelated
        one. Worth calling anyway: a long-running caller that verifies many runs
        would otherwise hold a handle per reference until the process ends.
        """
        try:
            self.fasta.close()
        except Exception:  # pragma: no cover - closing is best-effort
            pass


# ---------------------------------------------------------------------------
# check one: is this the reference the run used?
# ---------------------------------------------------------------------------


@dataclass
class ContigReport:
    """Whether the FASTA holds the contigs the run saw, at the lengths it saw."""

    reference: Path
    compared: int = 0
    absent: list[str] = field(default_factory=list)
    wrong_length: list[tuple[str, int, int]] = field(default_factory=list)
    run_reference: str | None = None

    @property
    def ok(self) -> bool:
        """Whether every contig the run saw is in the FASTA at the same length.

        Stricter than `fatal`, and the two are not opposites: a check can be not
        `ok` and not `fatal`, which is the ordinary case of verifying a run
        against a primary-only FASTA. That state is reported and allowed to
        proceed.
        """
        return not self.absent and not self.wrong_length

    @property
    def fatal(self) -> bool:
        """Whether this is bad enough that no position can be verified.

        A contig missing from the FASTA is not, on its own: a run aligned to a
        full assembly is often verified against a primary-only FASTA, and the
        positions on the chromosomes both have are checkable and worth checking.
        It becomes fatal when *nothing* matches, which is what a naming
        difference looks like, or when a shared contig has a different length,
        which means a different assembly and makes every coordinate suspect.
        """
        return self.compared == 0 or bool(self.wrong_length)

    def message(self) -> str:
        """One line stating what was compared and what disagreed.

        Both failures are reported, not just the first, because they mean
        different things together than apart -- `advice` reads the combination.
        Each list is cut to five examples with an ellipsis: a whole-assembly
        mismatch names thousands of contigs, and the count already says how many
        there are.
        """
        if self.ok:
            return (
                f"{self.compared} contig(s) of the run are in {self.reference.name} "
                f"at the same lengths"
            )
        parts = []
        if self.absent:
            shown = ", ".join(self.absent[:5])
            parts.append(
                f"{len(self.absent)} contig(s) the run aligned to are not in "
                f"{self.reference.name}: {shown}"
                f"{' ...' if len(self.absent) > 5 else ''}"
            )
        if self.wrong_length:
            shown = "; ".join(
                f"{name}: run {ran:,}, reference {ref:,}"
                for name, ran, ref in self.wrong_length[:5]
            )
            parts.append(
                f"{len(self.wrong_length)} contig(s) have a different length: {shown}"
                f"{' ...' if len(self.wrong_length) > 5 else ''}"
            )
        return "; ".join(parts)

    def advice(self) -> str:
        """What the pattern of failure usually means."""
        if self.ok:
            return ""
        lines = []
        if self.absent and not self.wrong_length:
            lines.append(
                "Contigs missing but none of the shared ones disagreeing usually "
                "means a naming difference (Ensembl's `1` against UCSC's `chr1`) "
                "or a FASTA of primary chromosomes only."
            )
        elif self.wrong_length:
            lines.append(
                "A contig of a different length means a different assembly, not a "
                "different naming scheme; no coordinate in this file refers to "
                "the same place in that FASTA."
            )
        if self.run_reference:
            lines.append(f"The run itself read {self.run_reference}.")
        return " ".join(lines)


def check_contigs(manifest, reference: Reference) -> ContigReport:
    """Compare the run's `@SQ` contigs with the FASTA index.

    Exact, and free: the manifest records every contig the BAM header declared,
    with its length, so this settles the "is this even the right genome" question
    without reading a base. It runs first for that reason -- a sampled base
    mismatch has several possible causes, and this rules the loudest one out.
    """
    lengths = manifest.contig_lengths
    report = ContigReport(
        reference=reference.path, run_reference=manifest.reference
    )
    for name, length in sorted(lengths.items()):
        if name not in reference.lengths:
            report.absent.append(name)
            continue
        report.compared += 1
        if length and reference.lengths[name] != length:
            report.wrong_length.append((name, length, reference.lengths[name]))
    return report


# ---------------------------------------------------------------------------
# check two: do the positions we write land on the right base?
# ---------------------------------------------------------------------------


@dataclass
class BaseReport:
    """How a sample of hit rows scored against the FASTA, at several shifts."""

    reference: Path
    offset: int
    sampled: int = 0
    checked: int = 0
    matched: dict[int, int] = field(default_factory=lambda: {s: 0 for s in SHIFTS})
    oriented: bool = True
    strand_source: str | None = None
    minus_rows: int = 0
    out_of_range: int = 0
    absent_contigs: Counter = field(default_factory=Counter)
    examples: list[tuple[str, int, str, str]] = field(default_factory=list)

    def rate(self, shift: int = 0) -> float:
        """The fraction of checked positions matching, at one candidate shift.

        The denominator is `checked`, not `sampled`: rows whose contig is absent
        from the FASTA, or whose position falls off the end of it, were never
        looked up and are counted separately. Dividing by `sampled` would let a
        FASTA missing half the assembly report a low match rate, which reads as a
        coordinate bug and is not one.

        Zero checked gives 0.0 rather than raising. A caller that cares tests
        `checked` -- `ok` does -- and the reporting path wants a number it can
        print beside an explanation.
        """
        return self.matched[shift] / self.checked if self.checked else 0.0

    @property
    def best_shift(self) -> int:
        """The candidate shift that matched most often.

        The diagnostic, not the verdict. When shift 0 fails and this returns -1
        or +1, the sample is saying the coordinates are systematically off by
        that much, which is a far more useful thing to report than a bare match
        rate. `SHIFTS` is ordered so that ties resolve towards 0.
        """
        return max(SHIFTS, key=lambda s: self.matched[s])

    def ok(self, min_match: float = DEFAULT_MIN_MATCH) -> bool:
        """Whether the sample passes at shift 0, the only shift that may pass.

        A run matching 100% at -1 is a failure, not a success at another offset:
        the shift this package intends is 0, and the others exist to name the bug
        rather than to admit it.

        `checked > 0` is part of passing. An empty sample has a match rate of 0.0
        by `rate`, but the requirement is stated here too, because "nothing was
        verified" must never read as "everything verified".
        """
        return self.checked > 0 and self.rate(0) >= min_match

    def summary(self) -> str:
        """One line, short enough to record in the file's /metadata."""
        return (
            f"{self.reference.name}: {self.matched[0]:,}/{self.checked:,} sampled "
            f"positions hold the expected base"
            + ("" if self.oriented else " (orientation not checked)")
        )

    def message(self) -> str:
        """The full account: the verdict's number, then everything that shapes it.

        The other shifts are always printed, including on success, because their
        rates are what show the check had power to fail. A sample where every
        shift matches 99% is not evidence of correct coordinates; it is evidence
        of a low-complexity or heavily sampled region, and only seeing the row
        reveals that.

        The remaining clauses appear only when non-zero, and each names a reason
        the headline number might understate the truth rather than a failure of
        its own.
        """
        parts = [
            f"{self.matched[0]:,} of {self.checked:,} sampled positions hold the "
            f"base alnbase read there ({self.rate(0):.1%})"
        ]
        others = ", ".join(
            f"{s:+d}: {self.rate(s):.1%}" for s in SHIFTS if s != 0
        )
        parts.append(f"at other shifts -- {others}")
        if self.minus_rows:
            parts.append(
                f"{self.minus_rows:,} of the sample sit on the minus strand, "
                f"where the expected base is complemented"
            )
        if self.out_of_range:
            parts.append(f"{self.out_of_range:,} fell outside their contig")
        if self.absent_contigs:
            worst = ", ".join(f"{c} ({n:,})" for c, n in self.absent_contigs.most_common(3))
            parts.append(f"contigs not in the reference: {worst}")
        return "; ".join(parts)

    def advice(self) -> str:
        """Name the bug when the numbers name it."""
        if self.checked == 0:
            return (
                "No position could be looked up, so nothing was verified. Either "
                "the contig names disagree or the sampled rows carry no plain "
                "reference base."
            )
        best = self.best_shift
        if best != 0 and self.rate(best) > self.rate(0) + 0.5:
            written = "1-based" if self.offset else "unshifted"
            return (
                f"The sample matches the reference {abs(best)} base(s) to the "
                f"{'left' if best < 0 else 'right'} of where these positions would "
                f"be written. Positions are being shifted by {self.offset:+d} to "
                f"make them {written}, taken from the manifest's coordinate_base; "
                f"a shift of {self.offset + best:+d} would put them on the "
                f"reference. Do not paper over this with a different offset until "
                f"it is clear which side is wrong: the manifest's coordinate_base "
                f"may not describe these files, this may not be the reference the "
                f"run read, or -- on CG data especially -- the rows' orientation "
                f"may have been read backwards, since a C's complement is the G "
                f"beside it and mistaking one for the other looks exactly like a "
                f"shift of one."
            )
        if not self.oriented:
            return (
                "No shift explains the mismatch. The hit table has no `strand` "
                "column, so the comparison accepted either orientation and cannot "
                "distinguish a one-base shift inside a palindromic context; write "
                "the hits with `-F strand` for the strict check. Otherwise this "
                "looks like the wrong reference."
            )
        return (
            "No shift explains the mismatch, so this is most likely not the "
            "reference the reads were aligned to, even though its contig names "
            "and lengths matched."
        )


def table_columns(con, path: Path) -> tuple[str, ...]:
    """The columns a hit table actually has.

    Read from the file rather than from the manifest's `fields`: that list names
    the record columns the run selected, while a query's captures become columns
    of their own, so trusting it would mean generating SQL for a column that is
    not there.
    """
    # `DESCRIBE` rather than `parquet_schema`, which reports a list column as a
    # group with `list`/`element` children and would hide `capture_refr` behind
    # the parquet encoding of a list.
    rows = con.execute(
        f"DESCRIBE SELECT * FROM read_parquet({sql_string(str(path))})"
    ).fetchall()
    return tuple(r[0] for r in rows)


def _strand_expression(columns) -> tuple[str, str | None]:
    """SQL that is TRUE where the walk ran backwards along the reference.

    `strand` is the walk's own direction, which is what `refr_base` was
    complemented against. `is_reverse` is FLAG 0x10 and is deliberately not used:
    for a pair it differs from the walk direction for exactly one mate, so it
    would mis-orient half the rows and turn a clean check into a coin flip.
    """
    if "strand" in columns:
        return "strand = '-'", "strand"
    if "conv_strand" in columns:
        return "conv_strand IN ('OB', 'CTOB')", "conv_strand"
    return "NULL", None


def _sample_sql(shard: Path, columns, rows: int, seed: int) -> str:
    """Rows of (contig, position, reference base, walked backwards?).

    Two layouts to read. A run where every query captures only its anchor writes
    the anchor's bases as scalars (`refr_base`); one where any query captures
    more writes every capture as lists, and then each capture carries its own
    position, so they are unnested and all of them are checked.
    """
    reversed_sql, _ = _strand_expression(columns)
    if "refr_base" in columns:
        source = f"""
        SELECT "ref_name" AS chr, "refr_pos" AS pos, "refr_base" AS base,
               ({reversed_sql}) AS reversed
        FROM read_parquet('{shard}')
        """
    elif "capture_refr" in columns:
        source = f"""
        SELECT "ref_name" AS chr,
               unnest("capture_refr_pos") AS pos,
               unnest("capture_refr") AS base,
               ({reversed_sql}) AS reversed
        FROM read_parquet('{shard}')
        """
    else:
        raise VerifyError(
            "the hit table carries no reference base column (`refr_base` or "
            "`capture_refr`), so there is nothing to compare the reference "
            "against. Hit tables written by alnbase always have one; this file "
            "may come from `alnbase extract`, which reads no reference."
        )

    bases = ", ".join(f"'{b}'" for b in DNA)
    return f"""
        SELECT chr, pos, base, reversed
        FROM ({source})
        WHERE chr IS NOT NULL AND pos IS NOT NULL AND base IN ({bases})
        USING SAMPLE reservoir({rows} ROWS) REPEATABLE ({seed})
    """


def sample_bases(
    con,
    shards,
    reference: Reference,
    offset: int,
    columns,
    sites: int = 2000,
    max_shards: int = 4,
    seed: int = 20260921,
) -> BaseReport:
    """Score a sample of hit rows against the FASTA, at each candidate shift.

    The sample is drawn from a bounded number of shards rather than from the
    whole run: a shift is systematic, so more shards buy accuracy nobody needs,
    while every extra shard is another pass over a column of a large parquet
    file. Which contigs exist is settled exactly by `check_contigs`, so breadth
    is not this check's job either.
    """
    shards = list(shards)
    if not shards or sites <= 0:
        return BaseReport(reference=reference.path, offset=offset)

    # Evenly spaced rather than the first few, so a run whose shards are split by
    # contig is not sampled entirely from chromosome 1.
    chosen = shards
    if len(shards) > max_shards:
        step = len(shards) / max_shards
        chosen = [shards[int(i * step)] for i in range(max_shards)]
    per_shard = max(1, sites // len(chosen))

    _, strand_source = _strand_expression(columns)
    report = BaseReport(
        reference=reference.path,
        offset=offset,
        oriented=strand_source is not None,
        strand_source=strand_source,
    )

    for shard in chosen:
        rows = con.execute(_sample_sql(Path(shard), columns, per_shard, seed)).fetchall()
        for chrom, pos, base, reversed_walk in rows:
            report.sampled += 1
            if chrom not in reference.lengths:
                report.absent_contigs[chrom] += 1
                continue

            expected = base.upper()
            complemented = expected.translate(COMPLEMENT)
            if report.oriented:
                wanted = {complemented} if reversed_walk else {expected}
                if reversed_walk:
                    report.minus_rows += 1
            else:
                # No strand column: either orientation is accepted, and the
                # report says so rather than pretending to the strict check.
                wanted = {expected, complemented}

            written = pos + offset
            hits = {}
            for shift in SHIFTS:
                found = reference.base(chrom, written + shift)
                hits[shift] = found
                if found is not None and found.upper() in wanted:
                    report.matched[shift] += 1

            if hits[0] is None:
                report.out_of_range += 1
                continue
            report.checked += 1
            if hits[0].upper() not in wanted and len(report.examples) < 5:
                report.examples.append(
                    (chrom, written, "/".join(sorted(wanted)), hits[0].upper())
                )

    return report


# ---------------------------------------------------------------------------
# the two checks together
# ---------------------------------------------------------------------------


def verify_reference(
    con,
    manifest,
    hits: str,
    reference_path: str | Path,
    sites: int = 2000,
    max_shards: int = 4,
    min_match: float = DEFAULT_MIN_MATCH,
    strict: bool = True,
) -> str:
    """Check the run against a reference FASTA, and describe the result.

    Returns a one-line summary for a writer to record beside the data, so that an
    output carries the evidence its coordinates were checked rather than leaving
    it in a terminal that has since been closed.

    Raises `VerifyError` when the data and the reference disagree, before any
    output is written -- an hour of aggregation that produces shifted files is
    worse than a refusal. `strict=False` downgrades that to a warning, for the
    case where the user knows the reference is a near relation of the run's.

    The positions checked are the ones an aggregation reports, so the shift comes
    from the same `Manifest.position_sql` the generated SQL uses: checking a
    number this module worked out for itself would verify this module against
    itself.
    """
    reference = Reference(reference_path)
    try:
        contigs = check_contigs(manifest, reference)
        if contigs.fatal:
            _fail(f"{contigs.message()}. {contigs.advice()}", strict)
            # A wrong assembly makes the sampled check meaningless noise, so with
            # the refusal downgraded the sample is skipped rather than reported.
            return f"{Path(reference_path).name}: contigs do not match the run"
        if contigs.ok:
            logger.info(f"Reference check: {contigs.message()}.")
        else:
            # Some contigs absent, every shared one the right length: carry on,
            # and the sample reports how many rows land on the missing ones.
            logger.warning(
                f"Reference check: {contigs.message()}. {contigs.advice()} The "
                f"positions on the {contigs.compared} shared contig(s) are still "
                f"checked."
            )

        if sites <= 0:
            # Asked for explicitly, so not a failure: the contig check is the
            # exact half of the verification, and it passed.
            logger.warning(
                "--verify-sites 0: the contigs were checked but no position was, "
                "so a shifted coordinate would not have been noticed."
            )
            return (
                f"{reference.path.name}: {contigs.compared} contig(s) match; "
                f"positions not sampled"
            )

        files = expand(con, hits)
        report = sample_bases(
            con, files, reference, manifest.shift(1),
            table_columns(con, files[0]), sites=sites, max_shards=max_shards,
        )
        if report.ok(min_match):
            logger.info(f"Reference check: {report.message()}.")
            if not report.oriented:
                logger.warning(
                    "The hit table has no `strand` column, so the reference check "
                    "accepted either orientation of each base. A one-base shift "
                    "inside a CG cannot be seen that way, because the C and the G "
                    "opposite it are complements. Write the hits with "
                    "`-F strand` for the strict check."
                )
            return report.summary()

        examples = "".join(
            f"\n  {chrom}:{pos} -- alnbase read {want}, the reference holds {got}"
            for chrom, pos, want, got in report.examples
        )
        _fail(f"{report.message()}.{examples}\n{report.advice()}", strict)
        return report.summary()
    finally:
        reference.close()


def _fail(message: str, strict: bool) -> None:
    """Refuse the run, or warn and go on, depending on `strict`.

    One function so that every way the check can fail offers the same choice, and
    so the non-strict path cannot quietly become silent: the message written is
    the message that would have been raised, and it says plainly that the run is
    continuing after a failed check.
    """
    if strict:
        raise VerifyError(message)
    logger.warning(f"Reference check failed, continuing anyway: {message}")
