#!/usr/bin/env bash
# How query files are stored in the header, dumped, and chosen by extract.
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
"$A" query --query-file run1_calls.toml --query-file run1_end.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/run1.bam"
"$A" query --overwrite-tags --query-file run2_calls.toml --query-file "$strand" "$out/run1.bam" "$out/ref.aref" "$out/run2.bam"
echo "## 1. header of run2.bam (CL cut)"
samtools view --no-PG -H "$out/run2.bam" | cut -c1-110 | cat -A | sed 's/\$$//'
echo "## 2. records"
samtools view "$out/run2.bam" | cut -f1,12-
echo "## 3. dump-query"
"$A" dump-query --file 0 "$out/run2.bam" | tail -n +3 | cmp - run2_calls.toml && echo "latest run: byte-identical"
"$A" dump-query --pg alnbase "$out/run2.bam"; echo "exit=$?"
for k in 0 1; do
  f=$([ $k = 0 ] && echo run1_calls.toml || echo run1_end.toml)
  "$A" dump-query --pg alnbase --file $k "$out/run2.bam" | tail -n +3 | cmp - $f && echo "run 1 file $k: byte-identical"
done
"$A" dump-query --pg alnbase --file 3 "$out/run2.bam"; echo "exit=$?"
"$A" dump-query --pg nope "$out/run2.bam"; echo "exit=$?"
echo "## 4. extract picks each tag's most recent definition"
"$A" extract "$out/run2.bam" "$out/x.parquet"; echo "exit=$?"
"$A" extract --tag XM "$out/run2.bam" "$out/xm.parquet"; echo "exit=$?"
"$A" extract --tag XE "$out/run2.bam" "$out/xe.parquet"; echo "exit=$?"
python3 -c "
import duckdb
print(duckdb.sql(\"select 'XM' tag, name, off_5p, refr_pos, refr_base from '$out/xm_0_0.parquet' union all select 'XE', name, off_5p, refr_pos, refr_base from '$out/xe_0_0.parquet' order by all\").fetchall())"
echo "## 5. header comments survive samtools view / sort"
samtools sort --no-PG -o "$out/sorted.bam" "$out/run2.bam" 2>/dev/null
"$A" dump-query --pg alnbase --file 1 "$out/sorted.bam" | tail -n +3 | cmp - run1_end.toml && echo "after samtools sort: byte-identical"
echo "## 6. a lost @CO line is detected"
samtools view --no-PG -H "$out/run2.bam" | grep -v 'PG:alnbase.1	file:0	line:3	' > "$out/h.sam"
samtools reheader --no-PG "$out/h.sam" "$out/run2.bam" > "$out/lost.bam"
"$A" dump-query "$out/lost.bam"; echo "exit=$?"
"$A" extract --tag XE "$out/lost.bam" "$out/lost.parquet"; echo "exit=$?"
rm -rf "$out"
