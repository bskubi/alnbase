"""The toy fixture, and the tour that states what it should contain.

The tour prints `ok` or `DIFFERS` per claim, so this test is one assertion: no
claim disagrees with the data. That keeps the numbers in the tutorial honest as
the library changes underneath them, which is the failure a hand-written
example always has.
"""

from __future__ import annotations

import io
import pathlib
import runpy
import subprocess
import sys

TOY = pathlib.Path(__file__).resolve().parents[1] / "examples" / "toy"


def test_the_toy_is_rebuilt_exactly(tmp_path):
    """Regenerating must give the same bytes: the fixture has no randomness in
    it, so a diff here means someone changed the data without meaning to."""
    before = (TOY / "toy_hits.parquet").read_bytes()
    fasta = (TOY / "toy.fa").read_text()
    subprocess.run([sys.executable, str(TOY / "make_toy.py")], check=True,
                   capture_output=True)
    assert (TOY / "toy_hits.parquet").read_bytes() == before
    assert (TOY / "toy.fa").read_text() == fasta


def test_every_claim_the_tour_makes_holds():
    captured = io.StringIO()
    stdout, sys.stdout = sys.stdout, captured
    try:
        runpy.run_path(str(TOY / "tour.py"), run_name="__main__")
    finally:
        sys.stdout = stdout
    text = captured.getvalue()
    assert "DIFFERS" not in text, text
    assert text.count("ok ") == 11


def test_the_toy_is_small_enough_to_check_by_hand():
    """The whole point. If it grows past a screenful it has stopped teaching."""
    import duckdb

    hits = duckdb.connect().execute(
        f"SELECT count(*) FROM '{TOY / 'toy_hits.parquet'}'"
    ).fetchone()[0]
    bases = sum(
        len(line.strip())
        for line in (TOY / "toy.fa").read_text().splitlines()
        if not line.startswith(">")
    )
    assert hits <= 40, f"{hits} hits is too many to check by hand"
    assert bases <= 100, f"{bases} bases is too many to index by eye"
