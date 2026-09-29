#!/usr/bin/env bash
# Bismark-identical XM/XR/XG from docs/bismark-xm.toml plus two strand tables.
set -euo pipefail
cd "$(dirname "$0")"
# The run's strand rule: how this aligner records the strand a read was converted on.
# alnbase ships one file per surveyed aligner; this is the repository's own copy.
strand=../../../../../queries/strand/directional.toml
exec > >(sed -E "s#/tmp/[^ /\t]+/#TMP/#g") 2>&1
A=${ALNBASE:-alnbase}
out=$(mktemp -d)
samtools view --no-PG -b -o "$out/reads.bam" reads.sam
"$A" index ref.fa "$out/ref.aref"
echo "## query (tagged BAM), run summary on stderr"
"$A" query --query-file queries.toml --query-file "$strand" "$out/reads.bam" "$out/ref.aref" "$out/tagged.bam" 2>&1
echo "## tags written (qname, XM, XR, XG, and every aux field)"
samtools view "$out/tagged.bam" | cut -f1,2,12-
echo "## Bismark model (bismark_model.py) vs alnbase"
python3 bismark_model.py > "$out/model.txt"
samtools view "$out/tagged.bam" | awk -F'\t' '!and($2,4){print $1"\t"$12"\t"$13"\t"$14}' > "$out/alnbase.txt"
diff "$out/model.txt" "$out/alnbase.txt" && echo "IDENTICAL"
echo "## aux types (samtools view shows TAG:TYPE:VALUE)"
samtools view "$out/tagged.bam" | head -1 | tr '\t' '\n' | grep -E '^X[MRG]:'
echo "## header lines added"
samtools view --no-PG -H "$out/tagged.bam" | grep -E '^@PG|^@CO	alnbase' | head -5
samtools view --no-PG -H "$out/tagged.bam" | grep -c '^@CO	alnbase:v1:line'
echo "## dump-query round trip: --file picks one of the run's stored files"
"$A" dump-query --file 0 "$out/tagged.bam" | head -2 | cut -c1-60
"$A" dump-query --file 0 "$out/tagged.bam" | tail -n +3 | cmp - queries.toml && echo "BYTE-IDENTICAL"
echo "## re-tag the tagged BAM: existing tags stop the run (exit code shown)"
set +e
"$A" query --query-file queries.toml --query-file "$strand" "$out/tagged.bam" "$out/ref.aref" "$out/again.bam" 2>&1; echo "exit=$?"
ls "$out/again.bam" 2>&1 | sed "s#$out/##"
set -e
"$A" query --overwrite-tags --query-file queries.toml --query-file "$strand" "$out/tagged.bam" "$out/ref.aref" "$out/again.bam" 2>&1
samtools view --no-PG -H "$out/again.bam" | grep -E '^@PG|alnbase:v1:file' | cut -f1-5
"$A" dump-query --file 0 "$out/again.bam" | head -1 | cut -f1-2
"$A" dump-query --pg alnbase --file 0 "$out/again.bam" | tail -n +3 | cmp - queries.toml && echo "first run still dumpable"
samtools view "$out/again.bam" | cut -f1,12- | head -2
rm -rf "$out"
