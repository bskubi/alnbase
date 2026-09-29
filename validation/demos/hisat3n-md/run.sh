#!/bin/bash
# Align reads that carry a conversion immediately before a deletion, and compare the MD
# tag HISAT-3N writes against the one samtools recomputes from the same alignment.
#
#   HISAT3N=<dir with hisat-3n and hisat-3n-build> PY=<python with pysam> ./run.sh [workdir]
set -euo pipefail

HISAT3N=${HISAT3N:?set HISAT3N to a built hisat-3n checkout}
PY=${PY:-python}
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${1:-$(mktemp -d)}
mkdir -p "$WORK" && cd "$WORK"

"$PY" "$HERE/make_reads.py"

"$HISAT3N/hisat-3n-build" --base-change T,C ref.fa idx > hisat-3n-build.log 2>&1
"$HISAT3N/hisat-3n" -x idx -U r1.fq --base-change T,C -S out.sam 2> hisat-3n.log

echo
echo "=== what HISAT-3N wrote, and what samtools calmd says it should be ==="
samtools faidx ref.fa
samtools view -bS out.sam > out.bam 2>/dev/null
samtools calmd out.bam ref.fa 2> calmd.log > fixed.sam
paste <(grep -v '^@' out.sam   | sed -E 's/.*\t(MD:Z:[^\t]*).*/\1/') \
      <(grep -v '^@' fixed.sam | sed -E 's/.*\t(MD:Z:[^\t]*).*/\1/') \
      <(grep -v '^@' out.sam   | cut -f1,6) \
  | awk -F'\t' '{printf "  %-12s %-16s HISAT-3N %-16s samtools %-16s %s\n",
                        $3, $4, $1, $2, ($1==$2 ? "agree" : "DIFFER")}'
echo
echo "  calmd complained about:"
sed 's/^/    /' calmd.log | head -20

echo
echo "=== can pysam still read the reference bases back? ==="
"$PY" "$HERE/check_pysam.py" out.sam ref.fa

echo
echo "=== the repeat-index MD builder, fed to the same consumers ==="
"$PY" "$HERE/probe_repeat.py"
