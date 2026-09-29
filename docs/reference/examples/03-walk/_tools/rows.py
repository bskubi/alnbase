"""Print hit rows as an aligned table, in file order.

Usage: rows.py 'prefix_*.parquet' 'SQL select list' ['SQL where clause']
"""
import sys
import duckdb

glob, cols = sys.argv[1], sys.argv[2]
where = sys.argv[3] if len(sys.argv) > 3 else "true"
rel = duckdb.sql(f"""
    select {cols} from read_parquet('{glob}', file_row_number = true)
    where {where} order by record_id, file_row_number""")
names = [d[0] for d in rel.description]
data = [["NULL" if v is None else str(v) for v in r] for r in rel.fetchall()]
width = [max([len(n)] + [len(r[i]) for r in data]) for i, n in enumerate(names)]
print("  ".join(n.rjust(width[i]) for i, n in enumerate(names)))
for r in data:
    print("  ".join(v.rjust(width[i]) for i, v in enumerate(r)))
