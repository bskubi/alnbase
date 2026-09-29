"""Build the toy reference and the toy hit table the tutorials run on.

    python examples/toy/make_toy.py

Writes `toy.fa`, `toy.fa.fai` and `toy_hits.parquet` beside this file.

It is deliberately tiny: three contigs of 24, 12 and 12 bases, and 30 hits
written out one per line below. Nothing is generated at random. You can print
the whole table, count any number in any output by hand, and disagree with the
software when it is wrong -- which is the only way a fixture teaches anything.

WHY ONE FILE HOLDS SEVERAL ASSAYS
---------------------------------
A real run is one assay. This table is not, deliberately: it carries bisulfite
methylation, RNA editing, a genotype scan and NOMe-seq accessibility together,
because the thing a newcomer most has to believe is that a measure is only a
set of query names and that measures over different names do not interact.

    chr1  AACGTTAACATTAACGTTAACATT   bisulfite: 2 CG pairs, 2 CHH
    chr2  GGATAGATAGGT               editing and genotype, over the same names
    chr3  TTGCATACGTTA               NOMe: one GCH beside one HCG

THE SIX METHYLATION SITES ON chr1
---------------------------------
    pos  2  +  CG    pos  3  -  CG     <- one CG, both strands, one apart
    pos  8  +  CHH
    pos 14  +  CG    pos 15  -  CG     <- the second CG pair
    pos 20  +  CHH

A methylated call keeps its context as the query name (`CG`, `CHH`); an
unmethylated one is the same position read as a T (`TG`, `THH`). That pairing
is the whole of the biology, and a run declares it rather than the library
knowing it.

WHAT IS PLANTED, AND WHAT IT IS FOR
-----------------------------------
    rA   cellA   ordinary: CG methylated, CH not
    rB   cellA   every CH methylated -- bisulfite failed on this read
    rC   cellB   the minus strand, so merge_strands has something to fold
    rD   cellB   no CG at all: invisible to every threshold, seen by missing()
    rE   cellC   one hit, so cellC is the thin block
    rF   cellC   ONE read, TWO records -- the reason a read is never a record
    rG   cellA   disagrees at chr1:2, so cellA has one mixed site
    e1   cellA   1 of 3 adenines edited
    e2   cellB   3 of 3 edited: hyper-edited
    e3   cellC   one of each non-reference base: a four-category measure
    n1   cellA   open chromatin, methylated CG
    n2   cellB   closed chromatin, unmethylated CG

NUMBERS TO CHECK IT AGAINST
---------------------------
    mCG over chr1, per read:  rA 2/2   rB 2/2   rC 1/3   rE 1/1   rF 1/1
                              rG 0/1
    mCH over chr1, per read:  rA 0/2   rB 2/2   rD 0/2   rF 0/1
    edited over chr2:         e1 1/3   e2 3/3   e3 1/1  (e3's A>C and A>T are
                              in no category of the editing measure)
    merge_strands on mCG:     chr1:2 and chr1:3 become one row at 2
    cellA at chr1:2:          2 of 3 -- the one mixed site
"""

from __future__ import annotations

import pathlib

HERE = pathlib.Path(__file__).resolve().parent

#: The reference. Short enough to index by eye; the positions named in the
#: hit table below are all real bases of it.
GENOME = {
    "chr1": "AACGTTAACATTAACGTTAACATT",
    "chr2": "GGATAGATAGGT",
    "chr3": "TTGCATACGTTA",
}

#: One row per hit: read, block, query name, contig, position, strand.
#: Positions are 0-based, which is what alnbase writes.
HITS = [
    # -- chr1, bisulfite ---------------------------------------------------
    # rA: converted properly. Both CGs methylated, neither CHH.
    ("rA", "cellA", "CG",  "chr1",  2, "+"),
    ("rA", "cellA", "THH", "chr1",  8, "+"),
    ("rA", "cellA", "CG",  "chr1", 14, "+"),
    ("rA", "cellA", "THH", "chr1", 20, "+"),

    # rB: both CHH methylated too. Bisulfite failed, so its CGs are worthless.
    ("rB", "cellA", "CG",  "chr1",  2, "+"),
    ("rB", "cellA", "CHH", "chr1",  8, "+"),
    ("rB", "cellA", "CG",  "chr1", 14, "+"),
    ("rB", "cellA", "CHH", "chr1", 20, "+"),

    # rC: the minus-strand cytosines, plus one plus-strand CG unmethylated.
    # chr1:2 and chr1:3 are the two halves of one CG.
    ("rC", "cellB", "TG",  "chr1",  2, "+"),
    ("rC", "cellB", "CG",  "chr1",  3, "-"),
    ("rC", "cellB", "TG",  "chr1", 15, "-"),

    # rD: covers only CHH. No CG row at all, so no CG threshold can see it.
    ("rD", "cellB", "THH", "chr1",  8, "+"),
    ("rD", "cellB", "THH", "chr1", 20, "+"),

    # rE: one hit, and the only read of cellC on chr1 besides rF.
    ("rE", "cellC", "CG",  "chr1",  2, "+"),

    # rF: one read name, two alignment records, far apart. Judging the records
    # separately would filter half a molecule.
    ("rF", "cellC", "CG",  "chr1",  2, "+"),
    ("rF", "cellC", "THH", "chr1", 20, "+"),

    # rG: disagrees with rA and rB at chr1:2, making cellA's site there 2 of 3
    # -- the only site in the file that is neither wholly methylated nor
    # wholly unmethylated. scbs drops exactly these.
    ("rG", "cellA", "TG",  "chr1",  2, "+"),

    # -- chr2, editing and genotype ---------------------------------------
    # The adenines are at 2, 4 and 6. `A>G` is an edited adenine, `A>A` an
    # unedited one; `A>C` and `A>T` belong to the genotype measure only.
    ("e1", "cellA", "A>A", "chr2",  2, "+"),
    ("e1", "cellA", "A>A", "chr2",  4, "+"),
    ("e1", "cellA", "A>G", "chr2",  6, "+"),

    ("e2", "cellB", "A>G", "chr2",  2, "+"),
    ("e2", "cellB", "A>G", "chr2",  4, "+"),
    ("e2", "cellB", "A>G", "chr2",  6, "+"),

    ("e3", "cellC", "A>C", "chr2",  2, "+"),
    ("e3", "cellC", "A>T", "chr2",  4, "+"),
    ("e3", "cellC", "A>G", "chr2",  6, "+"),

    # -- chr3, NOMe --------------------------------------------------------
    # chr3:3 is a GCH (accessibility) and chr3:7 an HCG (methylation). The two
    # readouts share no query name, which is the point.
    ("n1", "cellA", "GCH", "chr3",  3, "+"),
    ("n1", "cellA", "HCG", "chr3",  7, "+"),

    ("n2", "cellB", "GTH", "chr3",  3, "+"),
    ("n2", "cellB", "HTG", "chr3",  7, "+"),
]

COLUMNS = (
    "qname VARCHAR, is_first_in_template BOOLEAN, XB VARCHAR, "
    "name VARCHAR, ref_name VARCHAR, refr_pos BIGINT, strand VARCHAR, "
    "qual UTINYINT, mapq UTINYINT, off_5p BIGINT, off_3p BIGINT"
)

#: What every row shares. Constant here so that a hit-level filter in a
#: tutorial keeps everything, and varying one is the tutorial's own business.
QUAL, MAPQ, OFF_5P, OFF_3P = 35, 40, 20, 20


def write_fasta(path: pathlib.Path) -> None:
    """The FASTA and its `.fai`, both written here so a tutorial needs no
    samtools. Each contig fits on one line, so the offsets are easy to check."""
    lines, index, offset = [], [], 0
    for name, seq in GENOME.items():
        header = f">{name}\n"
        lines.append(header)
        offset += len(header)
        index.append(f"{name}\t{len(seq)}\t{offset}\t{len(seq)}\t{len(seq) + 1}")
        lines.append(seq + "\n")
        offset += len(seq) + 1
    path.write_text("".join(lines))
    (path.parent / (path.name + ".fai")).write_text("\n".join(index) + "\n")


def main() -> None:
    import duckdb

    write_fasta(HERE / "toy.fa")
    con = duckdb.connect()
    con.execute(f"CREATE TABLE h ({COLUMNS})")
    con.executemany(
        "INSERT INTO h VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [(q, True, b, n, c, p, s, QUAL, MAPQ, OFF_5P, OFF_3P)
         for q, b, n, c, p, s in HITS],
    )
    con.execute("COPY h TO ? (FORMAT parquet)", [str(HERE / "toy_hits.parquet")])

    print(f"toy.fa            {len(GENOME)} contigs, "
          f"{sum(len(s) for s in GENOME.values())} bases")
    print(f"toy_hits.parquet  {len(HITS)} hits\n")
    for row in con.execute(
        "SELECT ref_name, count(*), count(DISTINCT qname), count(DISTINCT XB), "
        "       list_sort(list(DISTINCT name)) "
        "FROM h GROUP BY ref_name ORDER BY 1"
    ).fetchall():
        print("  %-5s %2d hits  %d reads  %d blocks  %s" % row)


if __name__ == "__main__":
    main()
