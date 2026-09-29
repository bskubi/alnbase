#!/usr/bin/env bash
# Which FASTA bytes the indexer accepts and how they are stored, printed by
# `info --seq`, and matched. Anything other than an IUPAC letter or U, in
# either case, stops `alnbase index` with its line, contig and position.
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
echo "== info (10 bases: whitespace inside the sequence line was dropped)"
"$ALNBASE" info ref.aref
echo "== info --seq: U is stored as T, lowercase is folded"
"$ALNBASE" info --seq c1:1-10 ref.aref
echo "== query"
"$ALNBASE" query --query-file queries.toml --query-file "$strand" --parquet \
    reads.bam ref.aref out.parquet 2>&1
python -c "
import duckdb
duckdb.sql('''
  with pos as (select unnest(range(10)) as refr_pos)
  select p.refr_pos, any_value(h.read_base) read_base, any_value(h.refr_base) refr_base,
         coalesce(string_agg(h.name, ' ' order by h.name), '(nothing fired)') as fired
  from pos p left join 'out_*.parquet' h using (refr_pos)
  group by p.refr_pos order by p.refr_pos
''').show(max_rows=100)
"
echo "== every other byte is refused, and no index is left behind"
for bad in 'X' '*' '7' '0' '.' '-' 'E' '%'; do
  printf '>c1\nACGT\n>c2 description\nAC%sT\n' "$bad" > bad.fa
  "$ALNBASE" index bad.fa bad.aref 2>&1 || true
  [ -e bad.aref ] && echo "  (an index was written)"
done
