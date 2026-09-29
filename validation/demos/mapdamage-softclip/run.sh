#!/usr/bin/env bash
# mapDamage measures damage from the first ALIGNED base, not the read's own 5' end.
#
# Case A shows the mechanism with hand-written CIGARs: identical reads, identical
# damage, one file aligned end to end and one soft-clipped by two bases.
# Case B shows what it costs on simulated ancient DNA aligned by a real aligner:
# bwa aln (global, never clips) versus bwa mem (local, clips terminal mismatches).
#
#   PY=<python with pysam>  MAPDAMAGE=<mapDamage checkout>  BWA_BIN=<dir with bwa+samtools> \
#     ./run.sh [workdir]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${1:-$PWD/mapdamage-softclip-work}"
PY="${PY:-python3}"
export PYTHONPATH="${MAPDAMAGE:?set MAPDAMAGE to a mapDamage checkout}:${PYTHONPATH:-}"
export PATH="${BWA_BIN:?set BWA_BIN to a directory holding bwa and samtools}:$PATH"

MD=("$PY" -m mapdamage --merge-libraries --no-stats --no-plot)

mkdir -p "$WORK"
cd "$WORK"

echo "=== Case A: same reads, two CIGARs ==="
"$PY" "$HERE/hardcoded_cigar.py" caseA
for bam in endtoend clipped; do
    "${MD[@]}" -i "caseA/$bam.bam" -r caseA/ref.fa -d "caseA/md_$bam" >/dev/null 2>&1
done
"$PY" "$HERE/summarise.py" caseA/md_endtoend caseA/md_clipped

echo
echo "=== Case B: simulated ancient DNA, bwa aln vs bwa mem ==="
for level in high:0.40:11 low:0.10:12; do
    IFS=: read -r name damage seed <<<"$level"
    "$PY" "$HERE/simulate_adna.py" "$name" "$damage" "$seed"
    bwa index "${name}_ref.fa" >/dev/null 2>&1
    samtools faidx "${name}_ref.fa"

    bwa aln -t 4 "${name}_ref.fa" "$name.fq" 2>/dev/null >"$name.sai"
    bwa samse "${name}_ref.fa" "$name.sai" "$name.fq" 2>/dev/null |
        samtools sort -o "${name}_aln.bam" -
    bwa mem -t 4 "${name}_ref.fa" "$name.fq" 2>/dev/null |
        samtools sort -o "${name}_mem.bam" -
    for kind in aln mem; do
        samtools index "${name}_${kind}.bam"
        "${MD[@]}" -i "${name}_${kind}.bam" -r "${name}_ref.fa" \
            -d "md_${name}_${kind}" >/dev/null 2>&1
    done

    clipped=$(samtools view -F 4 "${name}_mem.bam" | awk '$6 ~ /S/' | wc -l)
    total=$(samtools view -c -F 4 "${name}_mem.bam")
    echo
    echo "-- ${name} damage: bwa mem soft-clipped ${clipped} of ${total} mapped reads"
    "$PY" "$HERE/summarise.py" "md_${name}_aln" "md_${name}_mem" "${name}_truth.tsv"
done
