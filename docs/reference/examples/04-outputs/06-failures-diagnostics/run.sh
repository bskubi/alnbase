#!/usr/bin/env bash
# Failure behaviour, leftovers, exit codes and odd option values.
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
echo "## 1. QUAL '*': qual and qual_phred are null in query --parquet, and qual is null in extract (0.1.3 gave 255)"
"$A" query --parquet -f qname,qual_phred --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/q.parquet"
"$A" extract -f qname --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/e.parquet"
python3 -c "
import duckdb
print(duckdb.sql(\"select qname, name, off_5p, qual, qual_phred[1:3] q3, unicode(qual_phred) cp from '$out/q_0_0.parquet' where qname='c'\").fetchall())
print(duckdb.sql(\"select qname, name, off_5p, qual from '$out/e_0_0.parquet' where qname='c'\").fetchall())"
echo "## 2. a failing parquet run reports the worker's own error and removes every output file"
python3 -c "
print('@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:20')
for i in range(80000):
    x = '..u.' if i == 10 else '.Z..'
    print(f'r{i}\t0\tchr1\t1\t60\t4M\t*\t0\t0\tTTTC\tIIII\tXM:Z:{x}')
" | samtools view --no-PG -b -o "$out/big.bam" -
samtools view --no-PG -h "$out/big.bam" | head -22 | samtools view --no-PG -b -o "$out/small.bam" -
echo "   20 records, record r10 bad, -@ 2:"
"$A" extract -@ 2 --query-file queries.toml --query-file "$strand" "$out/small.bam" "$out/small.parquet"; echo "exit=$?"
echo "   80000 records, record r10 bad, -@ 2 (the reader outlives the failed worker; the worker's error is still the one reported):"
"$A" extract -@ 2 --query-file queries.toml --query-file "$strand" "$out/big.bam" "$out/big.parquet"; echo "exit=$?"
ls "$out" | grep -E '^(small|big)_' || echo "   no small_* or big_* files left"
echo "## 3. a failing BAM run removes its output (stdout excepted)"
"$A" query --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/t.bam"; echo "exit=$?"; ls "$out" | grep '^t.bam' || echo "t.bam absent"
echo "## 4. degenerate option values (0 for --row-group-rows / --batch-rows is a usage error, exit 2)"
"$A" query --parquet --row-group-rows 0 --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/rg0.parquet" 2>&1 | grep -v '^note:' | sed -E "s/^thread '<unnamed>' \([0-9]+\)/thread '<unnamed>' (TID)/"; echo "exit=${PIPESTATUS[0]}"
"$A" query --parquet --shards-per-worker 0 --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/s0.parquet"; echo "exit=$?"
"$A" query --parquet -@ 0 --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/t0.parquet"; echo "exit=$?"; ls "$out" | grep '^t0_'
"$A" query --parquet --batch-rows 0 --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/b0.parquet"; echo "exit=$?"
(cd "$out" && "$A" query --parquet --query-file "$OLDPWD/queries.toml" --query-file "$OLDPWD/$strand" reads.bam ref.aref - ; echo "exit=$?"; ls | grep '^-')
"$A" query --parquet --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/nodir/x.parquet"; echo "exit=$?"
echo "## 5. dump-query on a BAM alnbase did not write, and from standard input"
"$A" dump-query "$out/reads.bam"; echo "exit=$?"
"$A" query --overwrite-tags --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" - 2>/dev/null | "$A" dump-query --file 0 - | tail -n +2
rm -rf "$out"
