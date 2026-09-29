#!/usr/bin/env bash
# Pads, the widest query, and how to keep a query's hits fixed.
#
# The walk adds pad columns past each read end: --end-context of them, by default
# the widest query's span. A query fires on any window it matches, including
# windows made only of pads. So a query that can match all-pad windows gains hits
# when a wider query is added to the run. Two ways to keep its hits fixed are shown:
# set --end-context explicitly, or exclude all-pad windows in the query.
set -euo pipefail
cd "$(dirname "$0")"
# The run's strand rule: how this aligner records the strand a read was converted on.
# alnbase ships one file per surveyed aligner; this is the repository's own copy.
strand=../../../../../queries/strand/directional.toml
ALNBASE="${ALNBASE:-alnbase}"
samtools view -b -o reads.bam reads.sam
"$ALNBASE" index ref.fa ref.aref > /dev/null

show() {  # $1 = label, rest = alnbase query options
    local label="$1"; shift
    rm -rf out && mkdir out
    "$ALNBASE" query --query-file "$strand" "$@" --parquet --only-hits reads.bam ref.aref out/hits.parquet 2> /dev/null
    echo "== $label"
    python -c "
import duckdb
print(duckdb.sql(\"\"\"select refr_pos as C_pos_0based, off_5p, read_base
                     from 'out/*.parquet' where name = 'CG_in_reference' order by refr_pos\"\"\"))"
}

show "1. alone: widest span 2, so 2 pads per end" \
    --query-file narrow.toml
show "2. beside a 6-column query: 6 pads per end, and the two all-pad CG windows fire" \
    --query-file narrow.toml --query-file wide.toml
show "3. beside the 6-column query, pads fixed with --end-context 1" \
    --query-file narrow.toml --query-file wide.toml --end-context 1
show "4. beside the 6-column query, all-pad windows excluded in the query" \
    --query-file not_all_pad.toml --query-file wide.toml
