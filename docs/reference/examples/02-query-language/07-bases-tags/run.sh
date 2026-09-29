#!/usr/bin/env bash
# Reproduce with: bash run.sh > expected.txt 2>&1
set -euo pipefail
cd "$(dirname "$0")"
# The run's strand rule: how this aligner records the strand a read was converted on.
# alnbase ships one file per surveyed aligner; this is the repository's own copy.
strand=../../../../../queries/strand/directional.toml
ALNBASE=${ALNBASE:-alnbase}
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
samtools view -b -o "$work/reads.bam" reads.sam
"$ALNBASE" index ref.fa "$work/ref.aref"
echo '== docs/bismark-xm.toml (u/U: reference context not determined)'
"$ALNBASE" query --query-file ../../../../bismark-xm.toml --query-file "$strand" "$work/reads.bam" "$work/ref.aref" "$work/xm.bam"
samtools view "$work/xm.bam" | cut -f 1,10,12-
echo '== clash.toml'
"$ALNBASE" query --query-file clash.toml --query-file "$strand" "$work/reads.bam" "$work/ref.aref" "$work/clash.bam" 2>&1 | sed "s#$work#\$work#g" || echo "(exit status non-zero)"
echo '== unplaced.toml'
"$ALNBASE" query --query-file unplaced.toml --query-file "$strand" "$work/reads.bam" "$work/ref.aref" "$work/unplaced.bam"
samtools view "$work/unplaced.bam" | cut -f 1,10,12-
