#!/usr/bin/env bash
# Idioms: read ends, context past ends, indel anchors, reference-gap columns.
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
run() { "$ALNBASE" query --query-file "$here/$1" --query-file "$strand" --parquet -F qname,flags,cigar "${@:3}" \
          "$out/reads.bam" "$out/ref.aref" "$out/$2.parquet" >/dev/null 2>&1; }
cols="qname, flags, cigar, strand, name, off_5p, off_3p, refr_pos, read_base, refr_base"
echo "== 1+2. read ends and context past the end"
run ends.toml ends
python "$tools/rows.py" "$out/ends_*.parquet" "$cols" "qname in ('r1f','r1r','r2f','r2r','ends_on_C','clip_only')"
echo "== 3. beside a deletion"
run indel.toml indel
python "$tools/rows.py" "$out/indel_*.parquet" "$cols" "name is not null"
echo "== 4a. --insertions emit: every N@. column is an inserted base (clipped bases are clip columns, over the reference)"
run gapref.toml g1 --insertions emit
python "$tools/rows.py" "$out/g1_*.parquet" "$cols" "name = 'refgap'"
echo "== 4b. bounded and one-sided patterns over the same run: only ins_only fires"
python "$tools/rows.py" "$out/g1_*.parquet" "$cols" "name is not null and name <> 'refgap' and qname in ('ins_only','clip_only','clip_3p')"
echo "== 4c. the clipped records' first/last walked bases are aligned bases; offsets count the clip"
python "$tools/rows.py" "$out/ends_*.parquet" "$cols" "qname in ('clip_only','clip_3p') and name in ('walk_first','walk_last','unclipped_first','clipped_first')"
