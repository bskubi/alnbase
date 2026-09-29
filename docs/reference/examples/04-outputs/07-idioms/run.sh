#!/usr/bin/env bash
# Pipeline idioms: cell barcode in rows and as the partition key; MAPQ/flags for filtering;
# context by query name versus by capture.
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
echo "## 1. -F CB:Z,mapq,is_duplicate --partition-by CB:Z --only-hits"
"$A" query --query-file "$strand" --parquet -@ 2 --only-hits -F CB:Z,mapq,is_duplicate --partition-by CB:Z \
  --query-file by_name.toml "$out/reads.bam" "$out/ref.aref" "$out/calls.parquet"
python3 -c "
import duckdb
con = duckdb.connect()
f = \"read_parquet('$out/calls_*.parquet', filename=true)\"
print(con.sql(f'select CB, count(distinct filename) files from {f} group by CB order by CB').fetchall())
print(con.sql(f'''select CB, split_part(name, '_', 1) ctx, refr_pos,
       count(*) filter (where name in ('CG_meth', 'CHG_meth')) meth, count(*) calls
  from {f} where mapq >= 20 and not is_duplicate
  group by all order by all'''))"
echo "## 2. the same contexts from captures"
"$A" query --parquet -f qname --query-file by_capture.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/cap.parquet"
python3 -c "
import duckdb
print(duckdb.sql(\"select qname, refr_pos, read_base_list[1] read_c, capture_refr, case when capture_refr[2]='G' then 'CG' when capture_refr[3]='G' then 'CHG' else 'CHH' end ctx from (select *, capture_read as read_base_list from '$out/cap_0_0.parquet') where qname = 'a2' order by refr_pos\"))"
rm -rf "$out"
