#!/usr/bin/env bash
# Parquet schema, record fields, aux typing, hitless rows, footer metadata.
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
echo "## 1. non-UTF-8 aux text selected as a column, without --permissive"
"$A" query --parquet -F XZ:Z --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/bad.parquet"; echo "exit=$?"
echo "the failed run removed its output files:"
ls "$out" | grep bad || echo "no bad_* files"
echo "## 2. flat layout, default fields (core), --permissive"
"$A" query --parquet --permissive --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/flat.parquet"; echo "exit=$?"
python3 schema.py "$out"/flat_0_0.parquet
duckdb -c "select record_id, name, off_5p, off_3p, refr_pos, qual, read_base, refr_base, ref_name, strand from '$out/flat_0_0.parquet'" 2>/dev/null || python3 -c "
import duckdb; print(duckdb.sql(\"select record_id, name, off_5p, off_3p, refr_pos, qual, read_base, refr_base, ref_name, strand from '$out/flat_0_0.parquet'\"))"
echo "## 3. list layout"
"$A" query --parquet --permissive --query-file lists.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/list.parquet"
python3 schema.py "$out"/list_0_0.parquet
python3 -c "
import duckdb; print(duckdb.sql(\"select record_id, name, off_5p, capture_col, capture_read, capture_refr, capture_qual, capture_off_5p, capture_off_3p, capture_refr_pos from '$out/list_0_0.parquet'\"))"
echo "## 4. -f all plus aux columns (and type mismatches), --only-hits off"
"$A" query --parquet --permissive -f all -F NM:i,NM:Z,XA:A,XA:i,XC:i,XC:f,XI:i,XF:f,XF:i,XH:H,XB:B,XB:i,BF:B,XZ:Z,YY:i --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/all.parquet"; echo "exit=$?"
python3 schema.py "$out"/all_0_0.parquet
python3 -c "
import duckdb
con = duckdb.connect()
print(con.sql(\"select distinct record_id, qname, flags, tid, ref_name, pos, end_pos, unclipped_start, unclipped_end, mapq, cigar, read_len, strand, insert_size, mate_tid, mate_ref_name, mate_pos from '$out/all_0_0.parquet' order by record_id\"))
print(con.sql(\"select distinct record_id, seq_ascii, qual_phred, is_paired, is_unmapped, is_reverse, is_first_in_template, is_last_in_template from '$out/all_0_0.parquet' order by record_id\"))
print(con.sql(\"select distinct record_id, NM, XA, XC, XI, XF, XH, XB, BF, XZ, YY from '$out/all_0_0.parquet' order by record_id\"))
"
echo "## 4b. the first spelling of a tag wins; a later spelling of the same tag is dropped. Type mismatches read as null:"
"$A" query --parquet --permissive -f qname -F NM:Z,XA:i,XC:f,XI:f,XF:i,XB:i,XH:i,XZ:i --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/mis.parquet"
python3 -c "
import duckdb; print(duckdb.sql(\"select distinct record_id, NM, XA, XC, XI, XF, XB, XH, XZ from '$out/mis_0_0.parquet' order by record_id\"))"
echo "## 5. aliases: -f qual selects qual_phred (not the hit column qual); is_ot selects strand"
"$A" query --parquet --permissive -f qual,is_ot,rname,SEQ,is-qcfail --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/alias.parquet"
python3 schema.py "$out"/alias_0_0.parquet | head -8
echo "## 6. --only-hits drops hitless rows (nohit, unmapped)"
"$A" query --parquet --permissive --only-hits -f qname --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/hits.parquet"
python3 -c "
import duckdb; print(duckdb.sql(\"select record_id, qname, count(*) n from '$out/hits_0_0.parquet' group by all order by record_id\"))"
echo "## 7. IPC stream output carries the same schema; the extension is not checked"
"$A" query --parquet --permissive --output-format ipc -f qname --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/s.parquet"
python3 -c "
import pyarrow.ipc as ipc
r = ipc.open_stream('$out/s_0_0.parquet'); t = r.read_all()
import json
m = {k.decode(): v.decode() for k, v in t.schema.metadata.items()}
print(t.schema.names, t.num_rows, {k: m[k] for k in ('format_version', 'coordinate_base')})
print('manifest keys', list(json.loads(m['alnbase_manifest'])))"
echo "## 8. parquet footer: compression, row groups, key-value metadata"
"$A" query --parquet --permissive --row-group-rows 2 --output-profile small -f qname --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/rg.parquet"
python3 -c "
import pyarrow.parquet as pq
m = pq.ParquetFile('$out/rg_0_0.parquet').metadata
print('row_groups', m.num_row_groups, [m.row_group(i).num_rows for i in range(m.num_row_groups)])
print('codec', m.row_group(0).column(0).compression, 'kv keys', sorted(k.decode() for k in m.metadata))
print({c: m.row_group(0).column(i).encodings for i, c in enumerate(m.schema.names)})"
echo "## 9. row order within a record follows the column where each query completes"
printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:20\nord\t0\tchr1\t1\t60\t7M\t*\t0\t0\tTTTCGTT\tIIIIIII\n' | samtools view --no-PG -b -o "$out/ord.bam" -
"$A" query --parquet -f qname --query-file order.toml --query-file "$strand" "$out/ord.bam" "$out/ref.aref" "$out/ord.parquet" 2>/dev/null
python3 -c "
import pyarrow.parquet as pq
t = pq.read_table('$out/ord_0_0.parquet', columns=['name', 'off_5p', 'refr_pos'])
print(list(zip(*[t.column(c).to_pylist() for c in t.column_names])))"
rm -rf "$out"
