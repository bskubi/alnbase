#!/usr/bin/env bash
# Edge cases of bases tags and of what a tagged-BAM run writes.
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
echo "## 1. off-reference record without --permissive"
"$A" query --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/t.bam"; echo "exit=$?"
test -e "$out/t.bam" && echo "t.bam exists" || echo "t.bam absent"
echo "## 2. with --permissive: unused-query warning, unplaced hits, every record written in order"
"$A" query --permissive --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/t.bam"; echo "exit=$?"
samtools view "$out/t.bam" | cut -f1,2,3,12-
echo "## 3. two codes marking one base"
"$A" query --permissive --query-file clash.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/c.bam"; echo "exit=$?"
echo "## 4. an anchor that can sit on a deletion"
"$A" query --query-file gap_anchor.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/g.bam"; echo "exit=$?"
echo "## 5. no tags declared"
"$A" query --query-file notags.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/n.bam"; echo "exit=$?"
echo "## 6. output name must end in .bam (or have no extension, or be -)"
"$A" query --permissive --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/t.sam"; echo "exit=$?"
"$A" query --permissive --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" - 2>/dev/null | samtools view -c -
echo "## 7. parquet options without --parquet, BAM options with it"
"$A" query --only-hits --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/x.bam"; echo "exit=$?"
"$A" query --parquet --writer-threads 2 --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/x.parquet"; echo "exit=$?"
echo "## 8. clap usage error exit code"
"$A" query --query-file "$strand" --no-such-flag 2>/dev/null; echo "exit=$?"
rm -rf "$out"
