#!/usr/bin/env bash
# Reproduce the ALLCools bug candidates in README.md on hand-written reads.
#
#   ALLCOOLS_ENV=/path/to/allcools-env PY=/path/to/python-with-duckdb \
#   ALNBASE=/path/to/alnbase ./run.sh [workdir]
#
# ALLCOOLS_ENV is the environment from ../../envs/allcools.yaml (ALLCools c9f7be2). It must
# provide allcools, samtools, bgzip and tabix: bam-to-allc runs `samtools mpileup` from PATH.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${ALLCOOLS_ENV:?set ALLCOOLS_ENV}" "${PY:?set PY}" "${ALNBASE:?set ALNBASE}"
export PATH="$ALLCOOLS_ENV/bin:$PATH"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }
promise() { printf -- '-- promised: %s\n' "$*"; }
shows() { printf -- '-- %s:\n' "$*"; }

# Reference, 1-based positions (ALLC and SAM are both 1-based):
#   chrT  C at 6, G at 7 (a CG)
#   chrU  C at 6, R at 7, G at 8 (CRG: the G's context holds an IUPAC code)
#   chrV  C at 6, A at 7, G at 8 (CAG: the same layout without IUPAC)
#   chrW  CGs at 6-7 and 20-21
#   chrX  a CG at 6-7
"$PY" - <<'PYEOF'
def contig(name, sites, length=40):
    s = bytearray(b"A" * length)
    for pos, bases in sites:          # 1-based
        s[pos - 1:pos - 1 + len(bases)] = bases.encode()
    return f">{name}\n{s.decode()}\n"
with open("ref.fa", "w") as f:
    f.write(contig("chrT", [(6, "CG")]))
    f.write(contig("chrU", [(6, "CRG")]))
    f.write(contig("chrV", [(6, "CAG")]))
    f.write(contig("chrW", [(6, "CG"), (20, "CG")]))
    f.write(contig("chrX", [(6, "CG")]))
PYEOF
samtools faidx ref.fa
cut -f1,2 ref.fa.fai > chrom.sizes
"$ALNBASE" index ref.fa ref.aref >/dev/null 2>&1

# bam NAME: SAM records on stdin (space-separated) -> sorted, indexed NAME.bam
bam() {
  { printf '@HD\tVN:1.6\tSO:unsorted\n'
    awk '{printf "@SQ\tSN:%s\tLN:%s\n", $1, $2}' ref.fa.fai
    tr ' ' '\t'; } | samtools sort -o "$1.bam" - && samtools index "$1.bam"
}
# allc NAME: ALLCools bam-to-allc with its defaults (min MAPQ 10, min base quality 20)
allc() {
  allcools bam-to-allc --bam_path "$1.bam" --reference_fasta ref.fa \
    --output_path "$1.allc.tsv.gz" >/dev/null 2>&1
}
duck() { "$PY" -c 'import duckdb, sys; duckdb.sql(sys.argv[1]).show()' "$1"; }

say "AC-1: mpileup read-start marker counted as a base (bam-to-allc)"
promise "one read gives coverage 1. samtools mpileup writes '^' plus chr(MAPQ+33) before"
promise "the base of a read that starts at the column; ALLCools counts '.'/'T' (C column)"
promise "and ','/'a' (G column) without removing that MAPQ character."
printf 'samtools %s\n' "$(samtools --version | head -n 1 | cut -d' ' -f2)"
# name FLAG POS MAPQ SEQ; each read starts exactly at the C (POS 6) or the G (POS 7)
while read -r name flag pos mapq seq why; do
  bam "$name" <<SAM
$name $flag chrT $pos $mapq 10M * 0 0 $seq IIIIIIIIII
SAM
  allc "$name"
  printf '%-22s mpileup %-5s ALLC: ' "$name ($why)" "$(samtools mpileup -B -f ref.fa -r chrT:$((flag ? 7 : 6))-$((flag ? 7 : 6)) "$name.bam" 2>/dev/null | cut -f5)"
  zcat "$name.allc.tsv.gz" | cut -f1-6 | tr '\t' ' '
done <<'EOF2'
q12_C_unmeth 0 6 12 TGAAAAAAAA control:'-'
q13_C_unmeth 0 6 13 TGAAAAAAAA '.'=meth
q51_C_meth 0 6 51 CGAAAAAAAA 'T'=unmeth
q11_G_meth 16 7 11 GAAAAAAAAA ','=meth
q64_G_meth 16 7 64 GAAAAAAAAA 'a'=unmeth
EOF2
echo "expected: coverage 1 in every row; counts 0/1, 0/1, 1/1, 1/1, 1/1"
for name in q12_C_unmeth q13_C_unmeth q51_C_meth q11_G_meth q64_G_meth; do
  "$ALNBASE" query --query-file "$HERE/cg.toml" --parquet -F qname "$name.bam" ref.aref "ac1_$name.parquet" >/dev/null 2>&1
done
shows "alnbase, CG calls per read"
duck "select qname, count(*) filter (where name = 'CG_methylated') as m, count(*) as cov
      from 'ac1_*.parquet' group by qname order by qname"

say "AC-4: IUPAC code in a bottom-strand context (bam-to-allc)"
bam ac4 <<'SAM'
fwdU 0 chrU 1 60 10M * 0 0 AAAAATAGAA IIIIIIIIII
revU 16 chrU 1 60 10M * 0 0 AAAAACAGAA IIIIIIIIII
fwdV 0 chrV 1 60 10M * 0 0 AAAAATAGAA IIIIIIIIII
revV 16 chrV 1 60 10M * 0 0 AAAAACAGAA IIIIIIIIII
SAM
promise "a row per covered C on either strand, context from the reference: chrU and chrV"
promise "should each give a '+' row at 6 and a '-' row at 8."
allc ac4
shows "ALLCools ALLC"
zcat ac4.allc.tsv.gz | tr '\t' ' '

say "AC-10: extract-allc --strandness merge"
bam ac10 <<'SAM'
w6 0 chrW 6 60 10M * 0 0 CGAAAAAAAA IIIIIIIIII
w7 16 chrW 7 60 10M * 0 0 GAAAAAAAAA IIIIIIIIII
w21 16 chrW 21 60 10M * 0 0 GAAAAAAAAA IIIIIIIIII
x6 0 chrX 6 60 10M * 0 0 CGAAAAAAAA IIIIIIIIII
SAM
allc ac10
shows "input ALLC (bam-to-allc)"
zcat ac10.allc.tsv.gz | tr '\t' ' '
promise "help: 'merge: This will only merge the count on adjacent CpG in +/- strands'."
promise "Every CG keeps one row, on '+' at its C: chrW 6 (2/2), chrW 20 (1/1), chrX 6 (1/1)."
allcools extract-allc --allc_path ac10.allc.tsv.gz --output_prefix ac10 --mc_contexts CGN \
  --chrom_size_path chrom.sizes --strandness merge >/dev/null 2>&1
shows "ALLCools extract-allc --strandness merge"
zcat ac10.CGN-Merge.allc.tsv.gz | tr '\t' ' '
