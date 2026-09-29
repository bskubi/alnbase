#!/usr/bin/env bash
# Which single-column codes fire over which observed (read, reference) values.
# Shows: matching is subset-of-code; reference N/IUPAC match only codes that
# contain the whole ambiguity set; lowercase FASTA is folded; `=`/`/` need a
# plain base on both sides; a reverse read sees complemented reference codes.
set -euo pipefail
ALNBASE=${ALNBASE:-alnbase}
here=$(cd "$(dirname "$0")" && pwd)
# The run's strand rule: how this aligner records the strand a read was converted on.
# alnbase ships one file per surveyed aligner; this is the repository's own copy.
strand=$here/../../../../../queries/strand/directional.toml
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
cp "$here"/ref.fa "$here"/reads.sam "$here"/queries.toml "$work"; cd "$work"

samtools view -b -o reads.bam reads.sam
echo "== read odd after SAM->BAM conversion:"; samtools view reads.bam | awk '$1=="odd"{print $1, $10}'
"$ALNBASE" index ref.fa ref.aref
"$ALNBASE" query --query-file queries.toml --query-file "$strand" --parquet -F qname \
    reads.bam ref.aref out.parquet 2>&1

python -c "
import duckdb
duckdb.sql('''
  select qname, refr_pos, read_base, refr_base,
         string_agg(name, ' ' order by name) as fired
  from 'out_*.parquet'
  group by record_id, qname, refr_pos, read_base, refr_base
  order by record_id, refr_pos
''').show(max_rows=100, max_width=200)
"
