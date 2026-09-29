#!/bin/bash
# Run GLORI's pileup step over reads that differ only in orientation, and count how many
# reads survive its read-end trim at each position on each strand.
#
#   GLORI=<GLORI-tools checkout> PY=<python with pysam and biopython> ./run.sh [workdir]
set -euo pipefail

GLORI=${GLORI:?set GLORI to a GLORI-tools checkout}
PY=${PY:-python}
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${1:-$(mktemp -d)}
mkdir -p "$WORK" && cd "$WORK"

"$PY" "$HERE/make_bam.py"

for trim in 1 2 3; do
  echo
  echo "=== --trim-head $trim --trim-tail $trim ==="
  "$PY" "$GLORI/pipelines/pileup_genome_multiprocessing.py" \
      -i reads.bam -f ref.fa -o "pileup.$trim.txt" \
      --trim-head "$trim" --trim-tail "$trim" > /dev/null 2>&1
  "$PY" "$HERE/summarize.py" "pileup.$trim.txt" "$trim"
done

echo
echo "=== which duplicate alignment's A-count survives the tie-break ==="
GLORI="$GLORI" "$PY" "$HERE/probe_tiebreak.py"
