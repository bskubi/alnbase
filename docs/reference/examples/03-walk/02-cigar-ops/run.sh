#!/usr/bin/env bash
# Columns emitted for each CIGAR operation.
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
echo "== default: --insertions skip"
"$ALNBASE" query --query-file "$here/queries.toml" --query-file "$strand" --parquet -F qname,flags,cigar \
    "$out/reads.bam" "$out/ref.aref" "$out/skip.parquet" 2>/dev/null
python "$tools/grid.py" "$out/skip_*.parquet"
echo "== --insertions emit (insertion and soft-clip records; the sclip_ blocks are unchanged)"
"$ALNBASE" query --query-file "$here/queries.toml" --query-file "$strand" --parquet -F qname,flags,cigar \
    --insertions emit \
    "$out/reads.bam" "$out/ref.aref" "$out/emit.parquet" 2>/dev/null
python "$tools/grid.py" "$out/emit_*.parquet" | awk '/^(ins_|sclip_|lead_ins)/{p=1} /^$/{if(p)print; p=0} p'
echo "== soft-clipped bases have no option: --soft-clips was removed in 0.1.2"
"$ALNBASE" query --query-file "$here/queries.toml" --query-file "$strand" --parquet --soft-clips emit \
    "$out/reads.bam" "$out/ref.aref" "$out/clips.parquet" 2>&1 | head -n 1 || true
