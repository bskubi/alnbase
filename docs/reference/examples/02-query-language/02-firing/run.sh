#!/usr/bin/env bash
# Reproduce with: bash run.sh > expected.txt 2>&1
set -euo pipefail
cd "$(dirname "$0")"
# The run's strand rule: how this aligner records the strand a read was converted on.
# alnbase ships one file per surveyed aligner; this is the repository's own copy.
strand=../../../../../queries/strand/directional.toml
ALNBASE=${ALNBASE:-alnbase}
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
show() {  # print a parquet shard in file order
  python3 -c "import duckdb,sys; print(duckdb.sql(sys.argv[1]))" "$1"
}
samtools view -b -o "$work/reads.bam" reads.sam
"$ALNBASE" index ref.fa "$work/ref.aref"
"$ALNBASE" query --query-file queries.toml --query-file "$strand" --parquet -f qname --only-hits \
  "$work/reads.bam" "$work/ref.aref" "$work/hits.parquet"
show "select qname, name, off_5p, off_3p, refr_pos, read_base, refr_base from '$work/hits_0_0.parquet'"
echo '$ --trace-records 1 (the first record, grid shows the flank columns)'
"$ALNBASE" query --query-file queries.toml --query-file "$strand" --trace-records 1 "$work/reads.bam" "$work/ref.aref" "$work/unused.bam" | head -n 16
