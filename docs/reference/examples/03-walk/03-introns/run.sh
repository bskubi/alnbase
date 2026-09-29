#!/usr/bin/env bash
# Intron (CIGAR N) windows and markers.
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
for ctx in default 0 5; do
    echo "== --splice-context $ctx"
    opt=(); [ "$ctx" = default ] || opt=(--splice-context "$ctx")
    "$ALNBASE" query --query-file "$here/queries.toml" --query-file "$strand" --parquet -F qname,flags,cigar "${opt[@]}" \
        "$out/reads.bam" "$out/ref.aref" "$out/c$ctx.parquet" 2>/dev/null
    python "$tools/grid.py" "$out/c${ctx}_*.parquet"
done
