#!/usr/bin/env bash
# How a BAM's @SQ contigs are matched to the index: by exact name, with a
# length check and, when the header carries @SQ M5, a sequence checksum check.
set -uo pipefail
ALNBASE=${ALNBASE:-alnbase}
here=$(cd "$(dirname "$0")" && pwd)
# The run's strand rule: how this aligner records the strand a read was converted on.
# alnbase ships one file per surveyed aligner; this is the repository's own copy.
strand=$here/../../../../../queries/strand/directional.toml
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
cp "$here"/*.fa "$here"/*.sam "$here"/queries.toml "$work"; cd "$work"
for s in *.sam; do samtools view -b -o "${s%.sam}.bam" "$s" 2>/dev/null; done
"$ALNBASE" index ref.fa ref.aref
"$ALNBASE" index other_assembly.fa other.aref
q() { # q BAM INDEX [extra args]: run a parquet query and print its hits
  local bam=$1 idx=$2; shift 2
  rm -f out_*.parquet
  echo "\$ alnbase query --parquet $* $bam $idx out.parquet"
  "$ALNBASE" query --query-file queries.toml --query-file "$strand" --parquet -F qname "$@" "$bam" "$idx" out.parquet 2>&1
  echo "(exit $?)"
  ls out_*.parquet >/dev/null 2>&1 && python -c "
import duckdb
duckdb.sql('''select qname, ref_name, refr_pos, read_base, refr_base, name
              from 'out_*.parquet' order by record_id, refr_pos''').show(max_rows=50)"
}

echo "== 1. order differs (BAM chr1,chr2,chr10; index chr1,chr10,chr2): resolved by name."
echo "      chr2's @SQ M5 agrees with the index; the unmapped record is skipped without error;"
echo "      chrRefOnly (only in the index) is never mentioned."
q reordered.bam ref.aref

echo "== 2. same names and lengths, DIFFERENT sequence (chr2 GGGGGTTTTT): the header's M5 catches it"
q reordered.bam other.aref

echo "   ... but without M5 in the header nothing can tell: it passes silently,"
echo "       and read G over reference T is reported as if it were real"
sed 's/\tM5:[0-9a-f]*//' reordered.sam | samtools view -b -o no_m5.bam - 2>/dev/null
q no_m5.bam other.aref

echo "== 3. a BAM-only contig with no records: warning only"
q extra_contig_no_records.bam ref.aref

echo "== 4. a record on a BAM-only contig: fatal by default ..."
q extra_contig_with_records.bam ref.aref
echo "   ... skipped and counted with --permissive (it still gets a null row)"
q extra_contig_with_records.bam ref.aref --permissive
echo "   ... and in BAM-tagging mode it is written through untagged"
cat > tag.toml <<'T'
[query.base]
read = "N"
refr = "~"
[tag.xb.bases]
fill = "."
b = "base"
T
"$ALNBASE" query --query-file tag.toml --query-file "$strand" --permissive extra_contig_with_records.bam ref.aref tagged.bam 2>&1
samtools view tagged.bam | cut -f1,3,4,12-
echo "   ... without --permissive, BAM mode removes its partial output (parquet mode, above, removed its files too)"
"$ALNBASE" query --query-file tag.toml --query-file "$strand" extra_contig_with_records.bam ref.aref failed.bam 2>&1 | grep -v 'cannot be scanned'
ls failed.bam 2>&1 | sed 's/^ls: //'

echo "== 5. a shared name with a different length: fatal"
q length_mismatch.bam ref.aref
echo "   header errors are raised before any output file exists:"
ls out_*.parquet 2>&1 | sed 's/^ls: //'

echo "== 6. no shared names at all (1 vs chr1): fatal, even with no mapped records"
q no_shared_names.bam ref.aref

echo "== 7. overlap --refr builds the same map (length mismatch fatal); --trace-records"
echo "      skips off-reference records SILENTLY, even without --permissive"
samtools sort -n -o length_mismatch.byname.bam length_mismatch.bam
"$ALNBASE" overlap --refr ref.aref --stale-tags recompute length_mismatch.byname.bam o.bam 2>&1; echo "(exit $?)"
"$ALNBASE" query --query-file queries.toml --query-file "$strand" --trace-records 5 --trace-grid false extra_contig_with_records.bam ref.aref - 2>&1 | grep -v '^ *$' | head -20; echo "(exit ${PIPESTATUS[0]})"
