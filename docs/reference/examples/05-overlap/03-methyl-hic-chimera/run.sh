#!/usr/bin/env bash
# alnbase reference, section 05 (overlap), example 03-methyl-hic-chimera. See README.txt.
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
echo "### samtools mpileup -A -B -Q 13 depth (column 4), per template, input vs alnbase output"
echo "### (mpileup's default overlap detection is on; -x turns it off)"
"$ALNBASE" overlap "$T/in.bam" "$T/out.bam" 2>/dev/null
for spec in "hic_sup_overlap chr2:101-125" "hic_junction_overlap chr1:86-100" "hic_junction_overlap chr2:101-125" \
            "hic_inverted chr2:116-140" "colinear_control chr1:156-185"; do
  set -- $spec
  for f in in out; do
    samtools view -h "$T/$f.bam" | awk -v t=$1 '/^@/||$1==t' | samtools sort -o "$T/s.bam" - 2>/dev/null
    samtools index "$T/s.bam"
    for x in "" "-x"; do
      printf '%-22s %-17s %-3s %-3s ' "$1" "$2" "$f" "${x:-  }"
      samtools mpileup $x -A -B -Q 13 -a -r "$2" "$T/s.bam" 2>/dev/null | cut -f4 | tr '\n' ' '; echo
    done
  done
done
