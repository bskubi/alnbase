#!/usr/bin/env bash
# Partitioning, file naming, shard column, record_id order, determinism across thread counts.
set -uo pipefail
cd "$(dirname "$0")"
# The run's strand rule: how this aligner records the strand a read was converted on.
# alnbase ships one file per surveyed aligner; this is the repository's own copy.
strand=../../../../../queries/strand/directional.toml
exec > >(sed -E "s#/tmp/[^ /\t]+/#TMP/#g") 2>&1
A=${ALNBASE:-alnbase}
out=$(mktemp -d)
samtools view --no-PG -b -o "$out/reads.bam" reads.sam
"$A" index ref.fa "$out/ref.aref"
echo "## 1. --partition-by CB:Z -@ 4 --shards-per-worker 2 (fragment warning on stderr)"
"$A" query --query-file "$strand" --parquet -@ 4 --shards-per-worker 2 --partition-by CB:Z -f qname,mapq -F CB:Z \
  --query-file queries.toml "$out/reads.bam" "$out/ref.aref" "$out/cb.parquet"; echo "exit=$?"
ls "$out" | grep '^cb_'
python3 - "$out" <<'PY'
import sys, glob, os, duckdb, pyarrow.parquet as pq
d = sys.argv[1]
files = sorted(glob.glob(f"{d}/cb_*.parquet"))
con = duckdb.connect()
q = lambda s: con.sql(s.replace("FILES", f"read_parquet('{d}/cb_*.parquet', filename=true)")).fetchall()
print("files:", len(files), "schemas identical:", len({str(pq.read_schema(f)) for f in files}) == 1)
print("rows per file:", [(os.path.basename(f), pq.ParquetFile(f).metadata.num_rows) for f in files])
print("CB values in more than one file:", q("select CB, count(distinct filename) n from FILES group by CB having n > 1"))
print("CB -> file:", *sorted((str(cb), os.path.basename(fn)) for cb, fn in q("select distinct CB, filename from FILES")), sep="\n  ")
# shard column == worker * shards_per_worker + shard, from the file name
bad = [(fn, s) for fn, s in q("select distinct filename, shard from FILES")
       if int(os.path.basename(fn).split("_")[1]) * 2 + int(os.path.basename(fn).split("_")[2].split(".")[0]) != s]
print("shard column disagrees with file name:", bad)
ids = [r[0] for r in q("select distinct record_id from FILES order by 1")]
print("record_ids dense 0..1499:", ids == list(range(1500)))
names_ok = q("select count(*) from FILES where qname != 'r' || lpad(record_id::varchar, 4, '0')")
print("record_id == input ordinal (qname rNNNN):", names_ok[0][0] == 0)
mono = True
for f in files:
    r = pq.read_table(f, columns=["record_id"]).column(0).to_pylist()
    mono &= all(a <= b for a, b in zip(r, r[1:]))
print("record_id non-decreasing within every file:", mono)
print("unmapped rows (null name, null ref) count:", q("select count(*) from FILES where record_id % 13 = 0 and name is null")[0][0], "expected", len(range(0,1500,13)))
PY
echo "## 2. an empty partition: 1 record into 8 files still writes 8 readable files"
head -4 reads.sam | samtools view --no-PG -b -o "$out/one.bam" -
"$A" query --parquet -@ 2 --shards-per-worker 4 --query-file queries.toml --query-file "$strand" "$out/one.bam" "$out/ref.aref" "$out/one"
ls "$out" | grep '^one_'
python3 -c "
import glob, pyarrow.parquet as pq
print(sorted((f.split('/')[-1], pq.ParquetFile(f).metadata.num_rows) for f in glob.glob('$out/one_*')))"
echo "## 3. same rows whatever the thread / shard count (default qname key)"
for t in 1 3; do
  "$A" query --parquet -@ $t --shards-per-worker 2 -f qname --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/t$t.parquet" 2>/dev/null
done
python3 -c "
import duckdb
a = duckdb.sql(\"select * exclude (shard) from '$out/t1_*.parquet' order by all\").fetchall()
b = duckdb.sql(\"select * exclude (shard) from '$out/t3_*.parquet' order by all\").fetchall()
print('rows', len(a), len(b), 'identical except shard:', a == b)"
echo "## 4. tagged BAM: -@ 1 vs -@ 4"
"$A" query -@ 1 --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/t1.bam"
"$A" query -@ 4 --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/t4.bam"
cmp -s "$out/t1.bam" "$out/t4.bam" && echo "files byte-identical" || echo "files differ as bytes (at least the @PG CL, which records -@)"
diff <(samtools view --no-PG "$out/t1.bam") <(samtools view --no-PG "$out/t4.bam") >/dev/null && echo "records identical"
diff <(samtools view --no-PG "$out/reads.bam" | cut -f1-11) <(samtools view --no-PG "$out/t4.bam" | cut -f1-11) >/dev/null && echo "input order and fields 1-11 preserved"
diff <(samtools view --no-PG -H "$out/t1.bam") <(samtools view --no-PG -H "$out/t4.bam")
rm -rf "$out"
