#!/usr/bin/env bash
# `info` output, `info --seq` (1-based inclusive) versus the 0-based refr_pos
# of hit rows, and what lies past the ends of a contig (PAD, coordinates -1
# and LN emitted as-is).
set -uo pipefail
ALNBASE=${ALNBASE:-alnbase}
here=$(cd "$(dirname "$0")" && pwd)
# The run's strand rule: how this aligner records the strand a read was converted on.
# alnbase ships one file per surveyed aligner; this is the repository's own copy.
strand=$here/../../../../../queries/strand/directional.toml
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
cp "$here"/ref.fa "$here"/reads.sam "$here"/queries.toml "$work"; cd "$work"
run() { echo "\$ alnbase $*"; "$ALNBASE" "$@" 2>&1; echo "(exit $?)"; }

"$ALNBASE" index ref.fa ref.aref
run info ref.aref
echo "== --seq is 1-based inclusive, like samtools faidx"
run info --seq c1:1-5 ref.aref
samtools faidx ref.fa c1:1-5
run info --seq c1:18-40 ref.aref          # TO past the end is clipped
run info --seq c1:1,0-1,2 ref.aref        # commas are stripped from numbers
run info --seq c1:0-5 ref.aref
run info --seq c1:21-22 ref.aref
run info --seq c1:5 ref.aref
run info --seq chrX:1-2 ref.aref

echo "== abridged contig list: more than 11 contigs shows first 5 and last 5"
for i in $(seq 0 13); do printf '>s%d\nACGT\n' $i; done > many.fa
"$ALNBASE" index many.fa many.aref
run info many.aref

echo "== hit coordinates are 0-based; off-contig reference is PAD"
samtools view -b -o reads.bam reads.sam
"$ALNBASE" query --query-file queries.toml --query-file "$strand" --parquet -F qname reads.bam ref.aref out.parquet 2>&1
python -c "
import duckdb
duckdb.sql('''
  select qname, name, refr_pos as anchor_pos, off_5p,
         capture_col, capture_refr_pos, capture_read, capture_refr, capture_off_5p
  from 'out_*.parquet' where name is not null order by qname, name
''').show(max_width=200)
"
