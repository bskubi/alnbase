#!/usr/bin/env bash
# extract: rows from a tag, compared with query --parquet.
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
"$A" query --query-file calls.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/tagged.bam"
samtools view "$out/tagged.bam" | cut -f1,2,6,10,12-
echo "## 1. extract without --tag: two bases tags are defined"
"$A" extract "$out/tagged.bam" "$out/x.parquet"; echo "exit=$?"
echo "## 2. extract --tag XM (definition from the header) vs query --parquet"
"$A" extract --tag XM "$out/tagged.bam" "$out/xm.parquet"; echo "exit=$?"
"$A" query --parquet --query-file calls.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/direct.parquet"
python3 - "$out" <<'PY'
import sys, duckdb
d = sys.argv[1]
cols = "record_id, name, off_5p, off_3p, refr_pos, qual, read_base, refr_base, ref_name, strand"
x = duckdb.sql(f"select {cols} from '{d}/xm_*.parquet' order by all").fetchall()
x = [r for r in x if r[1] is not None]
y = duckdb.sql(f"select {cols} from '{d}/direct_*.parquet' where name in ('CG','TG') order by all").fetchall()
print(duckdb.sql(f"select {cols} from '{d}/xm_*.parquet' order by record_id, off_5p"))
print("extract XM hit rows == query --parquet rows named CG/TG:", x == y)
print("hitless record_ids  extract:", [r[0] for r in duckdb.sql(f"select record_id from '{d}/xm_*.parquet' where name is null order by 1").fetchall()],
      " query --parquet:", [r[0] for r in duckdb.sql(f"select record_id from '{d}/direct_*.parquet' where name is null order by 1").fetchall()])
print(duckdb.sql(f"select name, count(*) from '{d}/direct_*.parquet' group by name order by name"))
PY
echo "## 3. extract --tag XX: refr_base per query (eqG -> read base, mm -> null) vs the walked value"
"$A" extract --tag XX "$out/tagged.bam" "$out/xx.parquet" 2>&1 | tail -1
python3 -c "
import duckdb
print(duckdb.sql(\"select x.record_id, x.name, x.off_5p, x.read_base, x.refr_base as extract_refr, d.refr_base as walked_refr from '$out/xx_*.parquet' x left join '$out/direct_*.parquet' d using (record_id, name, off_5p) where x.name is not null order by all\"))"
echo "## 4. extract a strand tag"
"$A" extract --tag XR "$out/tagged.bam" "$out/xr.parquet"; echo "exit=$?"
echo "## 5. standard input, and a foreign tag read with --query-file (header not consulted)"
samtools view --no-PG -b "$out/tagged.bam" | "$A" extract --query-file foreign.toml --query-file "$strand" - "$out/stdin.parquet"; echo "exit=$?"
echo "## 6. a character the definition lacks, and a tag of the wrong length"
printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:36\nodd\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTTTC\tIIII\tXM:Z:..u.\nshort\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTTTC\tIIII\tXM:Z:...\n' | samtools view --no-PG -b -o "$out/odd.bam" -
"$A" extract --query-file foreign.toml --query-file "$strand" "$out/odd.bam" "$out/odd.parquet"; echo "exit=$?"
ls "$out" | grep '^odd_' || echo "no odd_* files: a failed run removes its output"
samtools view --no-PG -b -o "$out/short.bam" <(printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:36\nshort\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTTTC\tIIII\tXM:Z:...\n')
"$A" extract --query-file foreign.toml --query-file "$strand" "$out/short.bam" "$out/short.parquet"; echo "exit=$?"
echo "## 7. extract needs no reference: a record on a contig absent from any index is decoded"
printf '@HD\tVN:1.6\n@SQ\tSN:chrZ\tLN:36\nz\t0\tchrZ\t1\t60\t4M\t*\t0\t0\tTTTC\tIIII\tXM:Z:.Z..\n' | samtools view --no-PG -b -o "$out/z.bam" -
"$A" extract --query-file foreign.toml --query-file "$strand" "$out/z.bam" "$out/z.parquet"
python3 -c "
import duckdb; print(duckdb.sql(\"select record_id, name, off_5p, refr_pos, read_base, refr_base, ref_name from '$out/z_*.parquet'\").fetchall())"
echo "## 8. a tag on an unmapped record is not decoded (hitless row)"
printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:36\nu\t4\t*\t0\t0\t*\t*\t0\t0\tTTTC\tIIII\tXM:Z:.Z..\n' | samtools view --no-PG -b -o "$out/u.bam" -
"$A" extract --query-file foreign.toml --query-file "$strand" "$out/u.bam" "$out/u.parquet"
python3 -c "
import duckdb; print(duckdb.sql(\"select record_id, name from '$out/u_*.parquet'\").fetchall())"
echo "## 9. extract never checks aux text: non-UTF-8 XZ is null, not an error and not counted (compare 03 section 1)"
samtools view --no-PG -b -o "$out/bad.bam" ../03-parquet-schema/reads.sam
"$A" extract -f qname -F XZ:Z --query-file foreign.toml --query-file "$strand" "$out/bad.bam" "$out/xbad.parquet"; echo "exit=$?"
python3 -c "
import duckdb; print(duckdb.sql(\"select record_id, qname, name, XZ from '$out/xbad_0_0.parquet'\").fetchall())"
rm -rf "$out"
