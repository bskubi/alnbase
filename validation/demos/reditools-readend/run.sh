#!/usr/bin/env bash
# REDItools3's read-end trim filters, -mbp and -Mbp, probed with hand-built BAMs.
#
# Three questions, each answered by a pair of files that differ in exactly one way:
#   1. Does a soft clip shift the position axis?      softclip_bams.py
#   2. Does the filter flip for reverse-strand reads? strand_bams.py
#   3. Does -Mbp trim as many bases as -mbp?          strand_bams.py, forward file
#
#   PY=<python with pysam and reditools installed> ./run.sh [workdir]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${1:-$PWD/reditools-readend-work}"
PY="${PY:-python3}"
RT=("$PY" -m reditools analyze)

mkdir -p "$WORK"
cd "$WORK"

# Report the substitution REDItools calls at a given reference position, or DROP.
call_at() { awk -F'\t' -v p="$2" '$2==p{print $8; f=1} END{if(!f) print "DROP"}' "$1"; }

echo "=== 1. Soft clip: does it shift the position axis? ==="
"$PY" "$HERE/softclip_bams.py" clip
echo
printf '%-10s %s\n' "" "-mbp 0   1   2   3   4"
for bam in endtoend clipped; do
    printf '%-10s ' "$bam"
    for mbp in 0 1 2 3 4; do
        "${RT[@]}" "clip/$bam.bam" -r clip/ref.fa -o c.tsv -mbp "$mbp" >/dev/null 2>&1
        printf '%-5s ' "$(call_at c.tsv 1003)"
    done
    echo
done
echo "The edit is read base 3 in both files. Identical thresholds => the clip does"
echo "not shift the axis; the filter counts from the read's own first base."

echo
echo "=== 2. Reverse strand: does the filter flip? ==="
"$PY" "$HERE/strand_bams.py" str
echo
printf '%-6s %-22s %s\n' "" "ref1003 (SEQ idx 2)" "ref1038 (SEQ idx 37)"
for bam in fwd rev; do
    for mbp in 0 3; do
        "${RT[@]}" "str/$bam.bam" -r str/ref.fa -o s.tsv -mbp "$mbp" >/dev/null 2>&1
        printf '%-6s -mbp %s -> %-15s %s\n' "$bam" "$mbp" \
            "$(call_at s.tsv 1003)" "$(call_at s.tsv 1038)"
    done
done
echo "For the FLAG-16 file, SEQ index 37 is read base 3 AS SEQUENCED and SEQ index 2"
echo "is read base 38. -mbp 3 drops index 2 and keeps index 37 in both files, so it"
echo "counts from SEQ position 0 and trims the wrong end of every reverse read."

echo
echo "=== 3. -mbp and -Mbp are not symmetric ==="
printf '%-8s %s\n' "" "3rd base from that end"
for m in 2 3; do
    "${RT[@]}" str/fwd.bam -r str/ref.fa -o a.tsv -mbp "$m" >/dev/null 2>&1
    printf -- '-mbp %-3s %s\n' "$m" "$(call_at a.tsv 1003)"
done
for M in 3 4; do
    "${RT[@]}" str/fwd.bam -r str/ref.fa -o b.tsv -Mbp "$M" >/dev/null 2>&1
    printf -- '-Mbp %-3s %s\n' "$M" "$(call_at b.tsv 1038)"
done
echo "-mbp N ignores the first N bases; -Mbp N ignores only the last N-1."
