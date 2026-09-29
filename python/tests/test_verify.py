"""Tests for the reference check.

The point of the check is to catch a shift of one base, so most of what is
tested here is a deliberately shifted or mis-oriented input that must *fail*. A
verifier that cannot be made to fail verifies nothing, so every case below either
takes a correct table and breaks one thing about it, or takes a broken one and
fixes it.

Expected bases are computed from the fixture sequence rather than written out, so
a change to the sequence cannot leave a stale letter behind that passes for the
wrong reason.

These rows are written with pyarrow rather than through `conftest.hits`, because
the cases need columns that fixture does not have -- `refr_base`, `strand`, and a
list-valued capture -- and need to control the manifest's contigs, which is what
a reference is compared against.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from alnbase_agg import Manifest
from alnbase_agg import verify

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")
pytest.importorskip("pyfaidx")

# 79 bases with CG, CH and CHH contexts in it. Lowercase at the end stands for a
# soft-masked repeat, which a naive comparison would read as a mismatch.
CHR1 = "ACGTTGCAGGCCATTACGATCGATCGGATCCAAGGTTCCAAGGATCTAGCTAGGCTTACGaacgttggccaacgttaac"
CHR2 = "CCGGTTAACCGGTTAACCGG"
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

# Both names the run scripts below use, so that a script naming `TG` is not
# refused for asking about a query the run never declared.
QUERY_FILE = '[query.CG]\npat = "cg"\n[query.TG]\npat = "cg"\n'


def write_fasta(tmp_path, chr1=CHR1, chr2=CHR2, name="ref.fa"):
    path = tmp_path / name
    path.write_text(f">chr1 a test contig\n{chr1}\n>chr2\n{chr2}\n")
    return path


def base_at(pos_0based, seq=CHR1):
    """The reference base at an alnbase (0-based) position, as the walk saw it."""
    return seq[pos_0based].upper()


def contigs(chr1_len=len(CHR1), chr2_len=len(CHR2)):
    return [
        {"name": "chr1", "length": chr1_len, "md5": None, "in_reference": None},
        {"name": "chr2", "length": chr2_len, "md5": None, "in_reference": None},
    ]


def manifest_doc(files, contig_list=None, reference=None, coordinate_base=0):
    return {
        "alnbase_version": "0.1.25",
        "format_version": "9",
        "coordinate_base": coordinate_base,
        "run_id": "t",
        "command": "query",
        "reference": reference,
        "output": {
            "files": list(files),
            "partition_by": ["XB"],
            "fields": ["qname", "record_id", "XB", "ref_name", "refr_pos", "name"],
        },
        "query_files": [{"path": "queries.toml", "text": QUERY_FILE}],
        "contigs": contig_list if contig_list is not None else contigs(),
    }


def write_hits(tmp_path, rows, stem="calls", extra=None, contig_list=None, **kw):
    """Write hit rows as alnbase writes them: 0-based `refr_pos`, one row per call.

    `rows` are (qname, record_id, barcode, chr, pos, query_name). The manifest
    goes in the parquet footer as well as beside it, because that is where
    `Manifest.read` looks: a hit table has to be able to describe itself even when
    the run that wrote it was interrupted.
    """
    table = pa.table(
        {
            "record_id": pa.array([r[1] for r in rows], pa.uint64()),
            "qname": pa.array([r[0] for r in rows]),
            "ref_name": pa.array([r[3] for r in rows]),
            "XB": pa.array([r[2] for r in rows]),
            "name": pa.array([r[5] for r in rows]),
            "refr_pos": pa.array([r[4] for r in rows], pa.int64()),
            **(extra or {}),
        }
    )
    name = f"{stem}_0_0.parquet"
    doc = manifest_doc([name], contig_list=contig_list, **kw)
    table = table.replace_schema_metadata({"alnbase_manifest": json.dumps(doc)})
    pq.write_table(table, tmp_path / name)
    (tmp_path / f"{stem}.manifest.json").write_text(json.dumps(doc))
    return str(tmp_path / f"{stem}_*_*.parquet")


def hit_rows(positions, strand="+", shift=0, seq=CHR1):
    """Hit rows whose `refr_base` is the reference base at their position.

    `shift` moves the recorded *position* away from the base, which is exactly
    the bug the check exists to find: everything is internally consistent and
    every row is one base off the genome.
    """
    rows, bases, strands = [], [], []
    for i, pos in enumerate(positions):
        rows.append((f"r{i}", i, "AAA", "chr1", pos + shift, "CG"))
        base = base_at(pos, seq)
        # The walk reverse complements both sides for a minus-strand read, so a
        # row on that strand records the complement of what the FASTA holds.
        bases.append(base.translate(COMPLEMENT) if strand == "-" else base)
        strands.append(strand)
    return rows, {"refr_base": pa.array(bases), "strand": pa.array(strands)}


class Table:
    """A written hit table and what the check needs to know about it."""

    def __init__(self, con, glob: str):
        self.glob = glob
        self.manifest = Manifest.read(con, glob)
        self.shards = [Path(p) for p in __import__("glob").glob(glob)]
        self.columns = verify.table_columns(con, self.shards[0])
        self.offset = self.manifest.shift(1)


def table(con, tmp_path, positions=(1, 5, 9, 17, 21), strand="+", shift=0, **kw):
    rows, extra = hit_rows(positions, strand=strand, shift=shift)
    return Table(con, write_hits(tmp_path, rows, extra=extra, **kw))


def written(con, tmp_path, rows, extra, **kw) -> Table:
    return Table(con, write_hits(tmp_path, rows, extra=extra, **kw))


def scored(con, t: Table, fasta, **kw):
    return verify.sample_bases(
        con, t.shards, verify.Reference(fasta), t.offset, t.columns, **kw
    )


# --- the reference itself -----------------------------------------------------

def test_contig_lengths_come_from_the_manifest(con, tmp_path):
    t = table(con, tmp_path)
    assert t.manifest.contig_lengths == {"chr1": len(CHR1), "chr2": len(CHR2)}


def test_a_matching_reference_passes_the_contig_check(con, tmp_path):
    t = table(con, tmp_path)
    report = verify.check_contigs(t.manifest, verify.Reference(write_fasta(tmp_path)))

    assert report.ok
    assert report.compared == 2


def test_a_contig_the_reference_lacks_is_named(con, tmp_path):
    """Ensembl's `1` against UCSC's `chr1` is this case, and it is total."""
    t = table(con, tmp_path)
    path = tmp_path / "ensembl.fa"
    path.write_text(f">1\n{CHR1}\n>2\n{CHR2}\n")

    report = verify.check_contigs(t.manifest, verify.Reference(path))

    assert not report.ok
    assert report.absent == ["chr1", "chr2"]
    assert "naming" in report.advice()


def test_a_contig_of_a_different_length_is_a_different_assembly(con, tmp_path):
    t = table(con, tmp_path)
    report = verify.check_contigs(
        t.manifest, verify.Reference(write_fasta(tmp_path, chr1=CHR1[:-10]))
    )

    assert not report.ok
    assert report.wrong_length == [("chr1", len(CHR1), len(CHR1) - 10)]
    assert "assembly" in report.advice()


def test_the_run_s_own_reference_is_reported_when_it_differs(con, tmp_path):
    t = table(con, tmp_path, reference="/refs/mm10.aref")
    report = verify.check_contigs(
        t.manifest, verify.Reference(write_fasta(tmp_path, chr1="ACGT"))
    )

    assert "/refs/mm10.aref" in report.advice()


def test_a_missing_reference_is_reported_not_raised_as_an_os_error(tmp_path):
    with pytest.raises(verify.VerifyError, match="no reference FASTA"):
        verify.Reference(tmp_path / "absent.fa")


def test_lookups_are_one_based_and_bounded(tmp_path):
    ref = verify.Reference(write_fasta(tmp_path))

    assert ref.base("chr1", 1) == CHR1[0]
    assert ref.base("chr1", len(CHR1)) == CHR1[-1].upper()
    assert ref.base("chr1", len(CHR1) + 1) is None
    assert ref.base("chr1", 0) is None, "position 0 must not wrap to the end"
    assert ref.base("chrZ", 1) is None


def test_a_soft_masked_base_is_not_a_mismatch(con, tmp_path):
    """Lowercase in a FASTA marks a repeat; it is the same base."""
    masked = CHR1.index("aacgtt")
    t = table(con, tmp_path, positions=(masked, masked + 2, masked + 3))
    report = scored(con, t, write_fasta(tmp_path))

    assert report.checked == 3
    assert report.matched[0] == 3


# --- the sampled positions ----------------------------------------------------

def test_correct_positions_match_as_written(con, tmp_path):
    report = scored(con, table(con, tmp_path), write_fasta(tmp_path))

    assert report.checked == 5
    assert report.matched[0] == 5
    assert report.ok()
    assert report.oriented, "the fixture has a strand column, so the check is strict"


def test_the_positions_checked_are_the_ones_the_query_reports(con, tmp_path):
    """The shift comes from the manifest, through the same method the SQL uses.

    A check that worked the offset out for itself would pass on a run whose
    generated SQL shifted differently, which is the one disagreement that matters.
    """
    t = table(con, tmp_path)
    assert t.manifest.position_sql("refr_pos") == "refr_pos + 1"
    assert t.offset == 1


def test_a_one_base_shift_is_caught_and_named(con, tmp_path):
    """The whole reason this check exists: output like this looks perfect."""
    report = scored(con, table(con, tmp_path, shift=1), write_fasta(tmp_path))

    assert not report.ok()
    assert report.matched[0] == 0
    # The position was recorded one base to the right of the base itself, so the
    # base is one to the left of where the position now points.
    assert report.best_shift == -1
    assert report.matched[-1] == 5
    assert "1 base(s) to the left" in report.advice()


def test_the_other_direction_is_caught_too(con, tmp_path):
    report = scored(con, table(con, tmp_path, shift=-1), write_fasta(tmp_path))

    assert report.best_shift == 1
    assert "1 base(s) to the right" in report.advice()


def test_a_minus_strand_row_expects_the_complemented_base(con, tmp_path):
    """The walk runs backwards along the reference for a minus-strand read, with
    both sides reverse complemented, so `refr_base` is the complement there."""
    report = scored(con, table(con, tmp_path, strand="-"), write_fasta(tmp_path))

    assert report.matched[0] == 5
    assert report.minus_rows == 5


def test_getting_the_orientation_backwards_is_caught(con, tmp_path):
    """The same complemented bases labelled `+`. If the check ignored the strand
    column, or complemented the wrong rows, this would pass.

    The positions are chosen so that no base within two of them is the
    complement of the base itself. On CG sites that is not true -- a C's
    complement is the G beside it -- so a mis-orientation there reads as a
    one-base shift instead, which `advice` says out loud.
    """
    rows, extra = hit_rows((31, 32, 33, 34, 39), strand="-")
    extra["strand"] = pa.array(["+"] * 5)
    report = scored(con, written(con, tmp_path, rows, extra), write_fasta(tmp_path))

    assert not report.ok()
    assert report.matched[0] == 0
    assert "not the reference" in report.advice()


def test_without_a_strand_column_either_orientation_is_accepted(con, tmp_path):
    """A weaker check, and it says so: a CG's C and the G opposite it are
    complements, so a shift inside a CG is invisible this way."""
    rows, extra = hit_rows((1, 5, 9, 17, 21), strand="-")
    del extra["strand"]
    report = scored(con, written(con, tmp_path, rows, extra), write_fasta(tmp_path))

    assert report.ok(), "complemented bases still match, unstranded"
    assert not report.oriented
    assert "orientation not checked" in report.summary()


def test_the_conversion_strand_column_serves_when_strand_is_absent(con, tmp_path):
    rows, extra = hit_rows((1, 5, 9), strand="-")
    del extra["strand"]
    extra["conv_strand"] = pa.array(["OB", "CTOB", "OB"])
    report = scored(con, written(con, tmp_path, rows, extra), write_fasta(tmp_path))

    assert report.oriented and report.strand_source == "conv_strand"
    assert report.matched[0] == 3


def test_a_pad_or_gap_column_is_not_compared(con, tmp_path):
    """alnbase spells a pad, a gap or a skipped intron with a character that is
    not a base; those rows carry no reference base to check."""
    rows, extra = hit_rows((1, 5, 9))
    extra["refr_base"] = pa.array(["_", ".", base_at(9)])
    report = scored(con, written(con, tmp_path, rows, extra), write_fasta(tmp_path))

    assert report.checked == 1 and report.matched[0] == 1


def test_a_position_past_the_end_of_its_contig_is_counted_not_compared(con, tmp_path):
    rows, extra = hit_rows((1, 5))
    rows.append(("rX", 9, "AAA", "chr1", len(CHR1) + 50, "CG"))
    extra["refr_base"] = pa.array(list(extra["refr_base"].to_pylist()) + ["C"])
    extra["strand"] = pa.array(["+", "+", "+"])
    report = scored(con, written(con, tmp_path, rows, extra), write_fasta(tmp_path))

    assert report.out_of_range == 1
    assert report.checked == 2 and report.ok()


def test_list_captures_are_unnested_and_all_checked(con, tmp_path):
    """A run where any query captures more than its anchor writes every capture
    as a list, with a position per capture."""
    rows = [("r0", 0, "AAA", "chr1", 1, "CG")]
    extra = {
        "capture_refr_pos": pa.array([[1, 2, 3]], pa.list_(pa.int64())),
        "capture_refr": pa.array([[base_at(1), base_at(2), base_at(3)]]),
        "strand": pa.array(["+"]),
    }
    report = scored(
        con, written(con, tmp_path, rows, extra), write_fasta(tmp_path), sites=100
    )

    assert report.checked == 3 and report.matched[0] == 3


def test_a_table_with_no_reference_base_at_all_is_refused(con, tmp_path):
    rows, _ = hit_rows((1, 5))
    t = written(con, tmp_path, rows, None)

    with pytest.raises(verify.VerifyError, match="no reference base column"):
        scored(con, t, write_fasta(tmp_path))


def test_sampling_nothing_checks_nothing(con, tmp_path):
    report = scored(con, table(con, tmp_path), write_fasta(tmp_path), sites=0)

    assert report.checked == 0 and not report.ok()
    assert "nothing was verified" in report.advice()


# --- the two checks together --------------------------------------------------

def verified(con, t: Table, fasta, **kw) -> str:
    return verify.verify_reference(con, t.manifest, t.glob, fasta, **kw)


def test_verify_reference_returns_a_summary_for_the_file(con, tmp_path):
    t = table(con, tmp_path)
    summary = verified(con, t, write_fasta(tmp_path))

    assert "5/5" in summary and "ref.fa" in summary


def test_verify_reference_raises_before_anything_is_written(con, tmp_path):
    t = table(con, tmp_path, shift=1)
    with pytest.raises(verify.VerifyError, match="to the left"):
        verified(con, t, write_fasta(tmp_path))


def test_verify_reference_can_be_downgraded_to_a_warning(con, tmp_path):
    t = table(con, tmp_path, shift=1)
    summary = verified(con, t, write_fasta(tmp_path), strict=False)

    assert "0/5" in summary


def test_a_wrong_assembly_skips_the_sample_rather_than_reporting_noise(con, tmp_path):
    t = table(con, tmp_path)
    summary = verified(
        con, t, write_fasta(tmp_path, chr1=CHR1[:-10]), strict=False
    )

    assert "contigs do not match" in summary


def test_verify_sites_zero_checks_the_contigs_only(con, tmp_path):
    """Documented as contigs-only, so it must not fail for want of a block."""
    summary = verified(con, table(con, tmp_path), write_fasta(tmp_path), sites=0)

    assert "contig(s) match" in summary and "not sampled" in summary


def test_the_suggested_shift_is_the_one_that_would_work(con, tmp_path):
    """`offset + best_shift`, not `offset - best_shift`. Naming the wrong remedy
    in a message about an off-by-one would be a poor joke.
    """
    report = scored(con, table(con, tmp_path, shift=1), write_fasta(tmp_path))

    # The fixture records each position one base too high, so with the +1 that
    # makes positions 1-based they land two bases out; +0 is what would work.
    assert report.best_shift == -1
    assert "a shift of +0 would put them on the reference" in report.advice()


def test_a_primary_only_reference_still_verifies_what_it_covers(con, tmp_path):
    """A run against a full assembly is often checked against a primary-only
    FASTA. The contigs it does have are worth checking, so a missing alt contig
    warns rather than refusing the run."""
    path = tmp_path / "primary.fa"
    path.write_text(f">chr1\n{CHR1}\n")
    t = table(con, tmp_path)

    report = verify.check_contigs(t.manifest, verify.Reference(path))
    assert not report.ok and not report.fatal
    assert report.absent == ["chr2"]

    assert "5/5" in verified(con, t, path)


def test_no_contig_in_common_is_fatal(con, tmp_path):
    t = table(con, tmp_path)
    path = tmp_path / "ensembl.fa"
    path.write_text(f">1\n{CHR1}\n")

    assert verify.check_contigs(t.manifest, verify.Reference(path)).fatal


def test_a_length_disagreement_is_fatal_even_beside_matching_contigs(con, tmp_path):
    t = table(con, tmp_path)
    report = verify.check_contigs(
        t.manifest, verify.Reference(write_fasta(tmp_path, chr1=CHR1[:-10]))
    )

    assert report.compared == 2 and report.fatal


def test_the_columns_are_read_from_the_file_not_the_manifest(con, tmp_path):
    """A query's captures are columns of the file that `output.fields` never names.

    Trusting the manifest's list would mean generating SQL for `refr_base` when
    the file holds `capture_refr`, which fails at the point where a clear message
    matters most.
    """
    rows = [("r0", 0, "AAA", "chr1", 1, "CG")]
    extra = {
        "capture_refr_pos": pa.array([[1]], pa.list_(pa.int64())),
        "capture_refr": pa.array([[base_at(1)]]),
    }
    t = written(con, tmp_path, rows, extra)

    assert "capture_refr" not in t.manifest.fields
    assert "capture_refr" in t.columns


# --- the command --------------------------------------------------------------
#
# The check is its own subcommand rather than a flag on `run`, because it is
# about a hit table and a FASTA and has nothing to do with any script: the
# question "are these coordinates where they say they are" is asked once per
# alnbase run, not once per output format.


def run_cli(*args):
    from click.testing import CliRunner

    from alnbase_agg.cli import main

    return CliRunner().invoke(main, [str(a) for a in args])


def check_case(tmp_path, shift=0):
    """A written hit table and the FASTA it should agree with."""
    rows, extra = hit_rows((1, 5, 9, 17, 21), shift=shift)
    glob = write_hits(tmp_path, rows, extra=extra)
    return glob, write_fasta(tmp_path)


def test_check_passes_on_an_unshifted_table(tmp_path):
    glob, fasta = check_case(tmp_path)

    result = run_cli("check", "--hits", glob, "--reference", fasta)

    assert result.exit_code == 0, result.output


def test_check_fails_on_a_shifted_table(tmp_path):
    """The whole point: one base of shift, found before anything is written."""
    glob, fasta = check_case(tmp_path, shift=1)

    result = run_cli("check", "--hits", glob, "--reference", fasta)

    assert result.exit_code != 0
    assert "alnbase-agg" in result.output


def test_check_can_be_downgraded_to_a_warning(tmp_path):
    glob, fasta = check_case(tmp_path, shift=1)

    result = run_cli("check", "--hits", glob, "--reference", fasta, "--no-strict")

    assert result.exit_code == 0, result.output


def test_check_refuses_a_glob_that_matches_nothing(tmp_path):
    """A mistyped glob otherwise reports a clean run over zero rows."""
    fasta = write_fasta(tmp_path)

    result = run_cli("check", "--hits", str(tmp_path / "nope_*.parquet"),
                     "--reference", fasta)

    assert result.exit_code != 0
    assert "no files match" in result.output
