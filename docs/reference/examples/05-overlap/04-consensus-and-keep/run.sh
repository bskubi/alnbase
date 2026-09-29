#!/usr/bin/env bash
# alnbase reference, section 05 (overlap), example 04-consensus-and-keep. See README.txt.
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
alnbase_index() { "$ALNBASE" index ref.fa "$T/ref.aref" >/dev/null 2>&1; }
run "defaults (subtract, sum-capped, cap 40, keep score, soft clip, strip)" --stats -
run "--mismatch-qual zero-both" --mismatch-qual zero-both
run "--mismatch-qual zero-loser" --mismatch-qual zero-loser
run "--mismatch-qual none" --mismatch-qual none
run "--mismatch-base set-n" --mismatch-base set-n
run "--match-qual max" --match-qual max
run "--match-qual none" --match-qual none
run "--qual-cap 60" --qual-cap 60
run "--keep r2" --keep r2
run "--clip-mode hard" --clip-mode hard
run "--keep r2 --clip-mode hard" --keep r2 --clip-mode hard
run "--score-qual-weight 0" --score-qual-weight 0
alnbase_index
run "--stale-tags recompute --refr ref.aref" --stale-tags recompute --refr "$T/ref.aref"
