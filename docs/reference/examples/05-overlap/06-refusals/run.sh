#!/usr/bin/env bash
# alnbase reference, section 05 (overlap), example 06-refusals. See README.txt.
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
run "defaults, --report-unresolved, statistics" --report-unresolved --stats -
run "--on-unresolved drop" --on-unresolved drop --stats -
run "relaxed: --min-support 2 --min-span-frac 0.5 --max-mismatch-frac 0.25" --min-support 2 --min-span-frac 0.5 --max-mismatch-frac 0.25 --report-unresolved
run "relaxed: --max-anchor-slack 20 (anchor_slack now fails the next test)" --max-anchor-slack 20 --report-unresolved --no-tag
run "relaxed: --length-tolerance 10" --length-tolerance 10 --report-unresolved --no-tag
