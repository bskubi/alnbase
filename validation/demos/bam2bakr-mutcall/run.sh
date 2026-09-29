#!/usr/bin/env bash
# bam2bakR's mutation caller, probed with hand-built BAMs.
#
# Two questions:
#   1. What does --minQual compare against? (config.yaml calls it
#      "Minimum base quality to call mutation" and ships the value 40)
#   2. Do the mutation numerator and the trials denominator apply the same
#      filters?
#
#   PY=<python with pysam> MUTCALL=<path to mut_call.py> ./run.sh [workdir]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${1:-$PWD/bam2bakr-mutcall-work}"
PY="${PY:-python3}"
MUTCALL="${MUTCALL:?set MUTCALL to bam2bakR workflow/scripts/mut_call.py}"

mkdir -p "$WORK"
"$PY" "$HERE/make_bams.py" "$WORK"
cd "$WORK"
: > snp.txt

# mut_call.py writes <input>_counts.csv beside the input.
field() { "$PY" -c "
import csv,sys
print(list(csv.DictReader(open(sys.argv[1])))[0][sys.argv[2]])" "$1" "$2"; }

call() {
    rm -f "${1%.bam}_counts.csv"
    "$PY" "$MUTCALL" -b "$1" --mutType TC --reads SE --minQual "$2" \
        --SNPs snp.txt --strandedness F >/dev/null 2>&1
}

echo "=== 1. What is --minQual measured in? ==="
echo "One T->C at Phred 10, one at Phred 5, both far from the read ends."
echo
printf '%-12s %-6s %-6s %s\n' "--minQual" "TC" "nT" "implied Phred cutoff"
for q in 30 33 34 40 44 76; do
    call qual.bam "$q"
    printf '%-12s %-6s %-6s %s\n' "$q" \
        "$(field qual_counts.csv TC)" "$(field qual_counts.csv nT)" \
        "q > $((q - 33))"
done
echo
echo "The comparison is 'phred + 33 > minQual', so the shipped value 40 is a"
echo "Phred cutoff of 8, and any value of 33 or less disables the filter."

echo
echo "=== 2. Numerator and denominator do not share filters ==="
echo "Every T in the read is mutated, so the true per-T rate is 1.000."
call ends.bam 40
nt=$(field ends_counts.csv nT); tc=$(field ends_counts.csv TC)
printf 'nT (trials) = %s   TC (mutations) = %s   reported rate = %s\n' \
    "$nt" "$tc" "$("$PY" -c "print('%.3f' % ($tc/$nt))")"
echo
echo "nT counts every T that passes the base-quality filter. TC additionally"
echo "requires the T to be more than --minDist bases from either end of the"
echo "aligned block, so the T's in the two end zones are trials that can never"
echo "become mutations. The rate is biased low by about 2*minDist/readlength."
