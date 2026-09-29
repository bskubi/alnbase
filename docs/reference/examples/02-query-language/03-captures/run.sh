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
for f in flat list; do
  echo "== $f.toml"
  "$ALNBASE" query --query-file $f.toml --query-file "$strand" --parquet -f qname --only-hits \
    "$work/reads.bam" "$work/ref.aref" "$work/$f.parquet" 2>/dev/null
  show "select column_name, column_type from (describe select * from '$work/${f}_0_0.parquet') where column_name not in ('shard','record_id','qname')"
done
show "select qname, name, capture_col as col, capture_read as rd, capture_refr as rf, capture_off_5p as off5, capture_refr_pos as pos, capture_qual as q from '$work/list_0_0.parquet'"
