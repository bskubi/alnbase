#!/usr/bin/env python3
"""Check alnbase's column coordinates against an independent derivation.

alnbase's walk decides, for every aligned column, its reference position, the
read and reference bases (complemented on the bottom strand), the base quality,
and the read offsets `off_5p` / `off_3p`. Unit tests inside alnbase reuse its own
arithmetic, so they cannot catch that arithmetic being wrong. This test derives
the same values with pysam (`get_aligned_pairs(with_cigar=True)`) and pyfaidx,
for random records covering every flag case, every CIGAR operation, hard and
soft clips, and missing qualities, and requires every value to agree.

    python walk_oracle.py --alnbase /path/to/alnbase [--records 3000] [--seed 1]

Needs samtools-free Python only: pysam, pyfaidx, duckdb. Exits 0 when everything
agrees, 1 with the first disagreements otherwise.

What is checked, per record, at several --end-context / --splice-context settings:
  - aligned bases (M, =, X): refr_pos, read_base, refr_base, qual, off_5p, off_3p
  - deleted bases (D): refr_pos, refr_base, off_5p, off_3p
  - inserted bases (I, with --insertions emit): read_base, qual, off_5p, off_3p
  - the flank, end_context columns past each end of the aligned span: refr_pos,
    refr_base (`_` off the contig); clips (`:`) nearest the aligned span, one per
    soft-clipped base up to end_context, with that base's off_5p, off_3p; pads
    (`_`) beyond, with null off_5p, off_3p
  - intron columns (N): the whole intron when it is at most twice the splice context,
    otherwise that many bases at each end and one `,@,` marker; refr_pos, refr_base,
    off_5p, off_3p
  - soft-clipped bases produce no rows
  - the record fields soft_clip_5p and soft_clip_3p

Conventions under test (docs/reference/03-walk-and-matching.md):
  - Positions are 0-based.
  - A record is walked on the bottom strand unless (read 2) == (reverse), and then
    both bases are reported complemented.
  - off_5p counts from the first base the sequencer read, over the whole read:
    SEQ plus any hard-clipped bases. A reverse record's 5' end is SEQ's last base.
    off_5p + off_3p == (SEQ length + hard clips) - 1.
  - A '*' QUAL gives a null quality; columns without a read base have null quality.
  - A deleted or intron column takes the offset of the aligned base before it in
    walk order: SEQ's preceding base on the top strand, its following base on the
    bottom strand (03 §4.3, §4.7).
  - A clip column stands for a soft-clipped base: the i-th column left of the
    aligned span is SEQ index soft_left - i, the i-th right of it
    len(SEQ) - soft_right - 1 + i. Pads are past the read and have no offsets
    (03 §5.2).
  - soft_clip_5p is the soft clip at the read's sequenced 5' end: the CIGAR's first
    non-hard operation on a forward record, its last on a reverse one.
  - The intron marker sits at the first elided position in walk order: intron start
    + context on the top strand, intron end - context - 1 on the bottom (03 §4.3).
"""

import argparse
import random
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import duckdb
import pysam
from pyfaidx import Fasta

COMPLEMENT = str.maketrans("ACGT", "TGCA")
HARD_CLIP = pysam.CHARD_CLIP
ALIGNED = (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF)
# (--end-context, --splice-context) settings to check: none, narrow, wide.
CONTEXTS = [(0, 0), (2, 3), (5, 1)]
QUERY = """
# Fires once on every column the read observed, so each row is one column.
[query.column]
read = "~"
refr = "~"
"""


def random_cigar(rng: random.Random) -> list[tuple[int, int]]:
    """A valid CIGAR: optional H, optional S, M blocks joined by I/D/N, optional S, H."""
    ops = []
    if rng.random() < 0.3:
        ops.append((pysam.CHARD_CLIP, rng.randint(1, 30)))
    if rng.random() < 0.3:
        ops.append((pysam.CSOFT_CLIP, rng.randint(1, 5)))
    for block in range(rng.randint(1, 4)):
        if block:
            op = rng.choice([pysam.CINS, pysam.CDEL, pysam.CREF_SKIP])
            ops.append((op, rng.randint(1, 20 if op == pysam.CREF_SKIP else 3)))
        ops.append((rng.choice([pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF]), rng.randint(1, 15)))
    if rng.random() < 0.3:
        ops.append((pysam.CSOFT_CLIP, rng.randint(1, 5)))
    if rng.random() < 0.3:
        ops.append((pysam.CHARD_CLIP, rng.randint(1, 30)))
    return ops


def write_inputs(workdir: Path, n_records: int, rng: random.Random) -> None:
    contigs = {f"c{i}": "".join(rng.choice("ACGT") for _ in range(rng.randint(150, 400))) for i in range(3)}
    with open(workdir / "ref.fa", "w") as fa:
        for name, seq in contigs.items():
            fa.write(f">{name}\n{seq}\n")

    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": n, "LN": len(s)} for n, s in contigs.items()]})
    with pysam.AlignmentFile(str(workdir / "reads.bam"), "wb", header=header) as out:
        written = 0
        while written < n_records:
            name = rng.choice(list(contigs))
            cigar = random_cigar(rng)
            ref_len = sum(n for op, n in cigar if op in (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF, pysam.CDEL, pysam.CREF_SKIP))
            seq_len = sum(n for op, n in cigar if op in (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF, pysam.CINS, pysam.CSOFT_CLIP))
            if ref_len > len(contigs[name]):
                continue
            rec = pysam.AlignedSegment(header)
            rec.query_name = f"r{written}"
            rec.reference_id = header.get_tid(name)
            rec.reference_start = rng.randint(0, len(contigs[name]) - ref_len)
            rec.cigartuples = cigar
            rec.flag = rng.choice([0, 16, 0x41, 0x51, 0x81, 0x91]) | rng.choice([0, 0, 0x100, 0x800])
            rec.query_sequence = "".join(rng.choice("ACGT") for _ in range(seq_len))
            if rng.random() < 0.9:
                rec.query_qualities = pysam.qualitystring_to_array("".join(chr(33 + rng.randint(2, 40)) for _ in range(seq_len)))
            rec.mapping_quality = 60
            out.write(rec)
            written += 1


def expected_rows(workdir: Path, end_context: int, splice_context: int) -> Counter:
    """The oracle: every row alnbase should write, derived from pysam and pyfaidx."""
    ref = Fasta(str(workdir / "ref.fa"), sequence_always_upper=True)
    rows = Counter()
    with pysam.AlignmentFile(str(workdir / "reads.bam"), "rb", check_sq=False) as bam:
        for rec in bam:
            contig = ref[rec.reference_name]
            seq = rec.query_sequence
            quals = rec.query_qualities  # None when QUAL is '*'
            reverse, read2 = rec.is_reverse, rec.is_read2
            bottom = read2 != reverse
            oriented = (lambda b: b.translate(COMPLEMENT)) if bottom else (lambda b: b)

            def refr_base(pos: int) -> str:
                return oriented(contig[pos].seq) if 0 <= pos < len(contig) else "_"

            cig = [c for c in rec.cigartuples if c[0] != HARD_CLIP]
            hard_first = rec.cigartuples[0][1] if rec.cigartuples[0][0] == HARD_CLIP else 0
            hard_last = rec.cigartuples[-1][1] if len(rec.cigartuples) > 1 and rec.cigartuples[-1][0] == HARD_CLIP else 0
            hard_5p, hard_3p = (hard_last, hard_first) if reverse else (hard_first, hard_last)
            read_len = len(seq) + hard_5p + hard_3p
            soft_left = cig[0][1] if cig[0][0] == pysam.CSOFT_CLIP else 0
            soft_right = cig[-1][1] if cig[-1][0] == pysam.CSOFT_CLIP else 0

            def offsets(q: int) -> tuple[int, int]:
                """Offsets of SEQ index q; q may lie outside SEQ, for a pad."""
                five_in_seq = len(seq) - 1 - q if reverse else q
                off_5p = five_in_seq + hard_5p
                return off_5p, read_len - 1 - off_5p

            pairs = rec.get_aligned_pairs(with_cigar=True)
            aligned_q = [q for q, r, op in pairs if op in ALIGNED]

            def walk_previous(before: int) -> tuple[int, int]:
                """Offsets of the aligned base before a gap in walk order; `before` is how
                many aligned bases precede the gap in SEQ order."""
                return offsets(aligned_q[before if bottom else before - 1])

            aligned_so_far = 0
            for q, r, op in pairs:
                qual = None if quals is None or q is None else quals[q]
                if op in ALIGNED:
                    rows[("aligned", rec.query_name, r, oriented(seq[q]), oriented(contig[r].seq), qual, *offsets(q))] += 1
                    aligned_so_far += 1
                elif op == pysam.CDEL:
                    rows[("deleted", rec.query_name, r, refr_base(r), *walk_previous(aligned_so_far))] += 1
                elif op == pysam.CINS:
                    rows[("inserted", rec.query_name, oriented(seq[q]), qual, *offsets(q))] += 1
                # CSOFT_CLIP positions must produce no rows; introns are handled below.

            # Introns, from the CIGAR: each N and the aligned bases before it.
            pos, aligned_so_far = rec.reference_start, 0
            for op, length in cig:
                if op == pysam.CREF_SKIP:
                    prev = walk_previous(aligned_so_far)
                    intron = range(pos, pos + length)
                    shown = intron if length <= 2 * splice_context else \
                        [*intron[:splice_context], *intron[length - splice_context:]]
                    for r in shown:
                        rows[("intron", rec.query_name, r, refr_base(r), *prev)] += 1
                    if length > 2 * splice_context:
                        marker = pos + splice_context if not bottom else pos + length - splice_context - 1
                        rows[("marker", rec.query_name, marker, *prev)] += 1
                if op in ALIGNED:
                    aligned_so_far += length
                if op in (*ALIGNED, pysam.CDEL, pysam.CREF_SKIP):
                    pos += length

            # Flank: end_context reference positions either side of the aligned span,
            # clips for the soft-clipped bases nearest it, pads beyond.
            for i in range(1, end_context + 1):
                for pos, clipped, q in ((rec.reference_start - i, i <= soft_left, soft_left - i),
                                        (rec.reference_end - 1 + i, i <= soft_right, len(seq) - soft_right - 1 + i)):
                    if clipped:
                        rows[("clip", rec.query_name, pos, refr_base(pos), *offsets(q))] += 1
                    else:
                        rows[("pad", rec.query_name, pos, refr_base(pos), None, None)] += 1

            rows[("clip fields", rec.query_name, *((soft_right, soft_left) if reverse else (soft_left, soft_right)))] += 1
    return rows


def alnbase_rows(workdir: Path, alnbase: str, end_context: int, splice_context: int) -> Counter:
    (workdir / "q.toml").write_text(QUERY)
    run = lambda *args: subprocess.run([alnbase, *args], check=True, capture_output=True, text=True)
    run("index", str(workdir / "ref.fa"), str(workdir / "ref.aref"))
    for old in workdir.glob("hits_*.parquet"):
        old.unlink()
    run("query", "--query-file", str(workdir / "q.toml"), "--insertions", "emit", "--max-discordance", "1",
        "--end-context", str(end_context), "--splice-context", str(splice_context),
        "--parquet", "--only-hits", "-F", "qname,soft_clip_5p,soft_clip_3p", str(workdir / "reads.bam"), str(workdir / "ref.aref"),
        str(workdir / "hits.parquet"))
    hits = f"read_parquet('{workdir}/hits_*.parquet')"
    table = duckdb.sql(f"select qname, refr_pos, read_base, refr_base, qual, off_5p, off_3p from {hits}").fetchall()
    clips = duckdb.sql(f"select distinct qname, soft_clip_5p, soft_clip_3p from {hits}").fetchall()
    rows = Counter(("clip fields", *c) for c in clips)
    for qname, refr_pos, read_base, refr_base, qual, off_5p, off_3p in table:
        if read_base == "_":
            rows[("pad", qname, refr_pos, refr_base, off_5p, off_3p)] += 1
        elif read_base == ":":
            rows[("clip", qname, refr_pos, refr_base, off_5p, off_3p)] += 1
        elif read_base == "," and refr_base == ",":
            rows[("marker", qname, refr_pos, off_5p, off_3p)] += 1
        elif read_base == ",":
            rows[("intron", qname, refr_pos, refr_base, off_5p, off_3p)] += 1
        elif read_base == ".":
            rows[("deleted", qname, refr_pos, refr_base, off_5p, off_3p)] += 1
        elif refr_base == ".":
            rows[("inserted", qname, read_base, qual, off_5p, off_3p)] += 1
        else:
            rows[("aligned", qname, refr_pos, read_base, refr_base, qual, off_5p, off_3p)] += 1
        if read_base not in "ACGT" and qual is not None:
            rows[("quality on a column without a read base", qname, refr_pos)] += 1
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alnbase", required=True)
    parser.add_argument("--records", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    failed = False
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        write_inputs(workdir, args.records, random.Random(args.seed))
        for end_context, splice_context in CONTEXTS:
            want = expected_rows(workdir, end_context, splice_context)
            got = alnbase_rows(workdir, args.alnbase, end_context, splice_context)
            missing, extra = want - got, got - want
            kinds = dict(sorted(Counter(k[0] for k in want.elements()).items()))
            print(f"--end-context {end_context} --splice-context {splice_context}: "
                  f"{args.records} records, expected rows {kinds}")
            if not missing and not extra:
                print("  alnbase agrees with the oracle on every column")
                continue
            failed = True
            print(f"  DISAGREEMENT: {sum(missing.values())} expected rows missing, {sum(extra.values())} unexpected rows")
            for label, rows in (("missing (oracle has, alnbase lacks)", missing), ("extra (alnbase has, oracle lacks)", extra)):
                print(f"    {label}:")
                for row in sorted(rows, key=str)[:10]:
                    print(f"      {row}")
    return 1 if failed else 0

if __name__ == "__main__":
    sys.exit(main())
