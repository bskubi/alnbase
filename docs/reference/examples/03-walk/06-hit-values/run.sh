#!/usr/bin/env bash
# How hit rows are derived from columns.
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
echo "== flat captures (every query records only its anchor)"
"$ALNBASE" query --query-file "$here/flat.toml" --query-file "$strand" --parquet -F qname,flags \
    "$out/reads.bam" "$out/ref.aref" "$out/flat.parquet" 2>/dev/null
python "$tools/rows.py" "$out/flat_*.parquet" "qname, flags, strand, name, off_5p, off_3p, refr_pos, qual, read_base, refr_base"
echo "== list captures"
"$ALNBASE" query --query-file "$here/lists.toml" --query-file "$strand" --parquet -F qname \
    "$out/reads.bam" "$out/ref.aref" "$out/lists.parquet" 2>/dev/null
python "$tools/rows.py" "$out/lists_*.parquet" "qname, off_5p, refr_pos, qual, capture_col, capture_read, capture_refr, capture_qual, capture_off_5p, capture_off_3p, capture_refr_pos" "refr_pos = 25 and qname in ('r0', 'ends_on_C', 'del_G')"
echo "== BAM tag (SEQ order) then extract"
"$ALNBASE" query --query-file "$here/tags.toml" --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/tagged.bam" 2>/dev/null
samtools view "$out/tagged.bam" | awk '{x=""; for(i=12;i<=NF;i++) if($i ~ /^XM:Z:/) x=$i; printf "  %-10s %4s %s %s\n", $1, $2, $10, x}'
"$ALNBASE" extract -F qname,flags "$out/tagged.bam" "$out/ext.parquet" >/dev/null 2>&1
"$ALNBASE" query --query-file "$here/tags.toml" --query-file "$strand" --parquet -F qname,flags "$out/reads.bam" "$out/ref.aref" "$out/direct.parquet" >/dev/null 2>&1
echo "-- extract"
python "$tools/rows.py" "$out/ext_*.parquet" "qname, flags, name, off_5p, off_3p, refr_pos, qual, read_base, refr_base" "qname in ('r0', 'r129')"
echo "-- direct parquet, same queries"
python "$tools/rows.py" "$out/direct_*.parquet" "qname, flags, name, off_5p, off_3p, refr_pos, qual, read_base, refr_base" "qname in ('r0', 'r129')"
