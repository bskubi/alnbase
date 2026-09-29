import sys, glob, pyarrow.parquet as pq
for f in sys.argv[1:]:
    s = pq.read_schema(f)
    print(f"-- {f.split('/')[-1]}  metadata={ {k.decode(): v.decode() for k, v in (s.metadata or {}).items()} }")
    for fld in s:
        print(f"   {fld.name:18s} {str(fld.type):22s} nullable={fld.nullable}")
