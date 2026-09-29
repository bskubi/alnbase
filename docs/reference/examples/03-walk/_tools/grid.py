"""Print every column the walk emitted, one block per record.

Usage: grid.py 'prefix_*.parquet'

The run must use a query named `col` that matches every column ("~" on both
sides, span 1) and must add `-F qname,flags,cigar`. Because `col` fires once
per column, its rows *are* the column list, in emission order
(file_row_number): read bases, deletions, inserted bases when emitted, and also
the pad columns past each read end, intron-context columns and the `,@,` marker
(a query fires on any window it matches; alnbase 0.1.3 to 0.1.6 fired only on
windows holding a read column, and this tool used a probe run then). A column that no pattern can match, not
even `~` (a read symbol that is the empty set), has no row and is not shown.
A '.' in the grid is a SQL NULL (qual is NULL wherever the column has no read base).
"""
import sys
import duckdb

rows = duckdb.sql(f"""
    select qname, flags, cigar, strand, off_5p, off_3p, refr_pos, qual, read_base, refr_base
    from read_parquet('{sys.argv[1]}', file_row_number = true)
    where name = 'col'
    order by record_id, file_row_number""").fetchall()

records = {}
for r in rows:
    records.setdefault(r[:4], []).append(r[4:])

labels = ["off_5p", "off_3p", "refr_pos", "qual", "read_base", "refr_base"]
for (qname, flags, cigar, strand), cols in records.items():
    print(f"{qname}  flags={flags}  cigar={cigar}  strand={strand}  columns={len(cols)}")
    grid = [["." if c[i] is None else str(c[i]) for c in cols] for i in range(len(labels))]
    width = [max(len(grid[i][j]) for i in range(len(labels))) for j in range(len(cols))]
    for i, label in enumerate(labels):
        print("  " + label.ljust(10) + " ".join(grid[i][j].rjust(width[j]) for j in range(len(cols))))
    print()
