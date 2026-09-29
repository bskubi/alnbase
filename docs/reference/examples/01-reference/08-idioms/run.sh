#!/usr/bin/env bash
# Idioms: selecting ambiguous / exactly-N reference symbols, and capturing
# reference context.
set -euo pipefail
ALNBASE=${ALNBASE:-alnbase}
here=$(cd "$(dirname "$0")" && pwd)
# The run's strand rule: how this aligner records the strand a read was converted on.
# alnbase ships one file per surveyed aligner; this is the repository's own copy.
strand=$here/../../../../../queries/strand/directional.toml
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
cp "$here"/ref.fa "$here"/reads.sam "$here"/queries.toml "$work"; cd "$work"

samtools view -b -o reads.bam reads.sam
"$ALNBASE" index ref.fa ref.aref
"$ALNBASE" query --query-file queries.toml --query-file "$strand" --parquet -F qname \
    reads.bam ref.aref out.parquet 2>&1
python -c "
import duckdb
duckdb.sql('''
  select qname, name, refr_pos, capture_refr_pos, capture_read, capture_refr
  from 'out_*.parquet' where name is not null
  order by name, qname, refr_pos
''').show(max_rows=100, max_width=200)
"
