#!/usr/bin/env bash
# alnbase reference, section 05 (overlap), example 05-drops-promotion-readthrough. See README.txt.
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
run "defaults" --stats -
run "--keep r2" --keep r2
run "--max-past-end 10 (alone: now refused by the anchor test)" --max-past-end 10 --no-tag
run "--max-past-end 10 --max-anchor-slack 10" --max-past-end 10 --max-anchor-slack 10
run "--max-past-end 10 --max-anchor-slack 10 --keep-terminal-indels" --max-past-end 10 --max-anchor-slack 10 --keep-terminal-indels --no-tag
echo "### samtools flagstat of the default output (note primary counts)"
"$ALNBASE" overlap "$T/in.bam" "$T/out.bam" 2>/dev/null; samtools flagstat "$T/out.bam" | head -8
