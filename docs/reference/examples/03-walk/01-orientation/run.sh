#!/usr/bin/env bash
# Walk orientation for every FLAG combination.
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
"$ALNBASE" query --query-file "$here/queries.toml" --query-file "$strand" --parquet -F qname,flags,cigar \
    "$out/reads.bam" "$out/ref.aref" "$out/cols.parquet" 2>/dev/null
echo "== every column of every record"
python "$tools/grid.py" "$out/cols_*.parquet"
echo "== per record: where off_5p = 0 sits, versus where sequencing started"
python "$tools/rows.py" "$out/cols_*.parquet" "qname, flags, strand, refr_pos as refr_pos_of_off5p_0,
    case when flags & 16 = 16 then 27 else 20 end as refr_pos_of_first_sequenced_base,
    (flags & 128 = 128) as is_read2" "name = 'col' and off_5p = 0"
echo "== off_5p is the sequencing-cycle offset for every FLAG, read 2 included"
python "$tools/rows.py" "$out/cols_*.parquet" "qname, refr_pos, off_5p, off_3p" \
    "name = 'col' and qual is not null and qname in ('f65','f81','f129','f145')"
