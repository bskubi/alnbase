"""A tour of the toy data: every concept once, with the answer beside it.

    python examples/toy/tour.py

Runs against `toy_hits.parquet`, which is 30 hits you can print in full. Each
step prints what it computed and what it should be, so a step that disagrees
with its own comment is visible rather than plausible.

Nothing here writes a file. The point is the numbers.
"""

from __future__ import annotations

import pathlib

import duckdb

from alnbase_agg.pipeline import Measure, Run

HERE = pathlib.Path(__file__).resolve().parent
HITS = str(HERE / "toy_hits.parquet")

# Six measures over one hit table. Only `edit` and `snp` share a query name --
# they read the same chr2 hits two different ways -- and none of them knows the
# others exist.
mCG = Measure("mCG", methylated="CG",  unmethylated="TG")
mCH = Measure("mCH", methylated="CHH", unmethylated="THH")
edit = Measure("edit", edited="A>G", unedited="A>A")
snp = Measure("snp", A="A>A", C="A>C", G="A>G", T="A>T")
acc = Measure("acc", open="GCH", closed="GTH")
meth = Measure("meth", methylated="HCG", unmethylated="HTG")

# A four-category measure has no single relation, so this one says which pair.
gva = snp.relate("G", "A")


def show(con, title, sql, expected):
    rows = con.execute(sql).fetchall()
    got = ", ".join(" ".join(str(v) for v in r) for r in rows) or "(nothing)"
    flag = "ok " if got == expected else "DIFFERS"
    print(f"  {flag} {title}")
    print(f"        got      {got}")
    if got != expected:
        print(f"        expected {expected}")


def main() -> None:
    con = duckdb.connect()
    run = Run(HITS, block="XB")
    run.measures(mCG, mCH, edit, snp, acc, meth, gva)
    run.count_sites()
    run.run(con, show=None)

    print(__doc__.strip().splitlines()[0])
    print()

    print("1. A measure is query names in categories. Nothing else knows the assay.")
    print("   (The views expose `num` and `den`, which is what a file format")
    print("    writes. `main`, `other`, `total` and `ratio` are the same numbers")
    print("    under the names a filter uses. The layers differ on purpose.)")
    show(con, "mCG per read: rA and rB saw two CGs, both methylated",
         "SELECT read.qname, num, den FROM read_relation "
         "WHERE relation = 'mCG' ORDER BY 1",
         "rA 2 2, rB 2 2, rC 1 3, rE 1 1, rF 1 1, rG 0 1")
    print()

    print("2. A name in no category of a measure simply does not join.")
    show(con, "edit per read: e3's A>C and A>T are in no category of `edit`",
         "SELECT read.qname, num, den FROM read_relation "
         "WHERE relation = 'edit' ORDER BY 1",
         "e1 1 3, e2 3 3, e3 1 1")
    show(con, "the same hits under `snp`, which has a category for each base",
         "SELECT read.qname, category, n FROM read_category "
         "WHERE measure = 'snp' AND read.qname = 'e3' ORDER BY 2",
         "e3 C 1, e3 G 1, e3 T 1")
    print()

    print("3. Absence is not zero. rD covers no CG, so it has no mCG row at all.")
    show(con, "a total threshold catches thin reads but never rD",
         "SELECT read.qname FROM read_fails('mCG', 'total', '<', 3) ORDER BY 1",
         "rA, rB, rE, rF, rG")
    show(con, "missing() does, because it reads src rather than the counts",
         "SELECT read.qname FROM read_fails_missing('mCG') "
         "WHERE read.qname LIKE 'r%' ORDER BY 1", "rD")
    print("     (without the LIKE it also returns every chr2 and chr3 read,")
    print("      which have no CG either. That is the price of putting four")
    print("      assays in one table, and a real single-assay run never sees it.)")
    print()

    print("4. The read is the read, not the record. rF is one read in two.")
    show(con, "one row, not two",
         "SELECT count(*) FROM read_relation "
         "WHERE relation = 'mCG' AND read.qname = 'rF'", "1")
    print()

    print("5. A CG is two cytosines. chr1:2 and chr1:3 are the same CG.")
    show(con, "unmerged: two rows, one base apart, opposite strands",
         "SELECT refr_pos, strand, num, den FROM site_relation "
         "WHERE relation = 'mCG' AND contig = 'chr1' AND refr_pos IN (2, 3) "
         "AND block = 'cellB' ORDER BY 1",
         "2 + 0 1, 3 - 1 1")
    print()

    print("6. Two independent readouts over one contig: that is NOMe-seq.")
    show(con, "accessibility at chr3:3, methylation at chr3:7, sharing no name",
         "SELECT relation, refr_pos, sum(num), sum(den) FROM site_relation "
         "WHERE relation IN ('acc', 'meth') AND contig = 'chr3' "
         "GROUP BY 1, 2 ORDER BY 1, 2", "acc 3 1 2, meth 7 1 2")
    print()

    print("7. A site can disagree with itself. cellA at chr1:2 is 2 of 3.")
    show(con, "the one mixed site, which scbs would drop",
         "SELECT block, refr_pos, num, den FROM site_relation "
         "WHERE relation = 'mCG' AND num > 0 AND den - num > 0 ORDER BY 1, 2",
         "cellA 2 2 3")
    print()

    print("8. A block is judged on reads and on positions, not on hits.")
    show(con, "cellA saw three reads over two CGs; cellB one read over three",
         "SELECT block, reads, sites FROM block_extent "
         "WHERE relation = 'mCG' ORDER BY 1",
         "cellA 3 2, cellB 1 3, cellC 2 1")
    print("     (`mCG.reads.lt(3)` drops cellB and cellC; `mCG.sites.lt(2)`")
    print("      drops cellC alone. Neither number is a count of hits.)")
    print()

    print("9. Evidence about a position is pooled across every block.")
    show(con, "G against A on chr2: only position 6 is mostly G",
         "SELECT refr_pos, sum(num), sum(den) FROM site_relation "
         "WHERE relation = 'snp.G_vs_A' GROUP BY 1 ORDER BY 1",
         "2 1 2, 4 1 2, 6 3 3")
    print("     (`at_position(gva.ratio.gt(0.5))` vetoes position 6 in every")
    print("      block. Judging each block alone would veto all three, because")
    print("      in the block that saw a G, every read there saw one.)")
    print()

    print("Read `make_toy.py` for the table itself; it is 30 lines of data.")


if __name__ == "__main__":
    main()
