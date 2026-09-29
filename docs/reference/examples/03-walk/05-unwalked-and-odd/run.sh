#!/usr/bin/env bash
# Unwalked, malformed and unusual records.
# Needs samtools and a python with duckdb on PATH.
set -euo pipefail
ALNBASE=${ALNBASE:-alnbase}
here=$(cd "$(dirname "$0")" && pwd)
# The run's strand rule: how this aligner records the strand a read was converted on.
# alnbase ships one file per surveyed aligner; this is the repository's own copy.
strand=$here/../../../../../queries/strand/directional.toml
tools="$here/../_tools"
out=$(mktemp -d)
trap 'rm -rf "$out"' EXIT
samtools view -b -o "$out/reads.bam" "$here/reads.sam"
"$ALNBASE" index "$here/ref.fa" "$out/ref.aref" >/dev/null
echo "== without --permissive: the off-reference record stops the run"
"$ALNBASE" query --query-file "$here/queries.toml" --query-file "$strand" --parquet -F qname \
    "$out/reads.bam" "$out/ref.aref" "$out/strict.parquet" 2>&1 | sed 's/^/  /' || true
echo "== with --permissive"
"$ALNBASE" query --query-file "$here/queries.toml" --query-file "$strand" --parquet -F qname,flags,cigar --permissive \
    "$out/reads.bam" "$out/ref.aref" "$out/cols.parquet" 2>&1 | sed 's/^/  /'
echo "-- one row per record that got no column rows"
python "$tools/rows.py" "$out/cols_*.parquet" "record_id, qname, flags, cigar, name, off_5p, refr_pos, qual, read_base" "name is null"
echo "-- columns of the walked records"
python "$tools/grid.py" "$out/cols_*.parquet"
echo "-- eq_in_seq in a trace: read symbol 0 (empty set), no pattern progresses through it"
"$ALNBASE" query --query-file "$here/queries.toml" --query-file "$strand" --trace-records 5 \
    "$out/reads.bam" "$out/ref.aref" /dev/null 2>/dev/null | awk '/^eq_in_seq/{p=1} p&&/^  wide/{print; exit} p'
echo "== seq_star.sam: mapped record with SEQ '*'"
samtools view -b -o "$out/seq_star.bam" "$here/seq_star.sam"
if "$ALNBASE" query --query-file "$here/queries.toml" --query-file "$strand" --parquet \
    "$out/seq_star.bam" "$out/ref.aref" "$out/star.parquet" >"$out/star.log" 2>&1; then
    echo "  exit 0"
else
    echo "  exit non-zero"
fi
sed 's/^/  /' "$out/star.log"
python "$tools/rows.py" "$out/star_*.parquet" "record_id, name, off_5p, refr_pos, read_base"
