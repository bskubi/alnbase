#!/usr/bin/env bash
# alnbase reference, section 05 (overlap), example 07-circles-and-two-copy. See README.txt.
set -uo pipefail
ALNBASE=${ALNBASE:-alnbase}
cd "$(dirname "$0")"
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
samtools view --no-PG -b -o "$T/in.bam" reads.sam
# run LABEL [overlap options...]: run overlap on $T/in.bam, print stderr and the output records.
run() {
  local label=$1; shift
  echo "### $label"
  local shown="$*"; echo "\$ alnbase overlap ${shown//$T\//} in.bam out.bam"
  rm -f "$T/out.bam"
  "$ALNBASE" overlap "$@" "$T/in.bam" "$T/out.bam" 2>"$T/err"; local rc=$?
  echo "exit status: $rc"
  sed 's/^/stderr| /' "$T/err"
  if [ -s "$T/out.bam" ]; then samtools view "$T/out.bam"; fi
  echo
}
run "defaults" --stats - --report-unresolved
echo "### samtools mpileup -A -B -Q 13 depth, input vs alnbase output"
"$ALNBASE" overlap "$T/in.bam" "$T/out.bam" 2>/dev/null
for spec in "self_circle chr1:36-55" "self_circle chr1:106-125" "sister_chromatid chr2:41-80"; do
  set -- $spec
  for f in in out; do
    samtools view -h "$T/$f.bam" | awk -v t=$1 '/^@/||$1==t' | samtools sort -o "$T/s.bam" - 2>/dev/null
    samtools index "$T/s.bam"
    for x in "" "-x"; do
      printf '%-18s %-15s %-3s %-3s ' "$1" "$2" "$f" "${x:-  }"
      samtools mpileup $x -A -B -Q 13 -a -r "$2" "$T/s.bam" 2>/dev/null | cut -f4 | tr '\n' ' '; echo
    done
  done
done
