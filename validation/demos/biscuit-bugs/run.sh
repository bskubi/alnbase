#!/usr/bin/env bash
# Reproduce the BISCUIT bug candidates in README.md on hand-written reads.
#
#   BISCUIT_ENV=/path/to/biscuit-env PY=/path/to/python-with-duckdb \
#   ALNBASE=/path/to/alnbase ./run.sh [workdir]
#
# BISCUIT_ENV is the environment from ../../envs/biscuit.yaml with BISCUIT built into it by
# ../../envs/build_biscuit.sh (commit 0a5ceae). It provides biscuit and samtools.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
: "${BISCUIT_ENV:?set BISCUIT_ENV}" "${PY:?set PY}" "${ALNBASE:?set ALNBASE}"
export PATH="$BISCUIT_ENV/bin:$PATH"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }
promise() { printf -- '-- promised: %s\n' "$*"; }
shows() { printf -- '-- %s:\n' "$*"; }

# Reference, 1-based positions (VCF and SAM are both 1-based):
#   chrH  30 bp of T with a CG at 11-12
#   chrE  30 bp of A with a C at 20 and a C at 30, the last base
#   chrG  10 A, 10 G, 10 C (21-30), 10 A
#   chrR  20 A, 10 C (21-30), 10 A: the same Cs with no G anywhere
#   chrF  60 bp of A with CGs at 5-6 and 55-56
#   chrP  200 bp of A with CGs at 150-151 and 160-161
#   chrI  200 bp of A with a CG at 120-121
#   chrJ  200 bp of A
#   chrS  40 bp of A with a CG at 21-22
#   chrQ  40 bp of A with a CG at 1-2, the contig's first base
"$PY" - <<'PYEOF'
contigs = {
    "chrH": "T" * 10 + "CG" + "T" * 18,
    "chrE": "A" * 19 + "C" + "A" * 9 + "C",
    "chrG": "A" * 10 + "G" * 10 + "C" * 10 + "A" * 10,
    "chrR": "A" * 20 + "C" * 10 + "A" * 10,
    "chrF": "A" * 4 + "CG" + "A" * 48 + "CG" + "A" * 4,
    "chrP": "A" * 149 + "CG" + "A" * 8 + "CG" + "A" * 39,
    "chrI": "A" * 119 + "CG" + "A" * 79,
    "chrJ": "A" * 200,
    "chrS": "A" * 20 + "CG" + "A" * 18,
    "chrQ": "CG" + "A" * 38,
}
with open("ref.fa", "w") as f:
    for name, seq in contigs.items():
        f.write(f">{name}\n{seq}\n")
PYEOF
samtools faidx ref.fa
"$ALNBASE" index ref.fa ref.aref >/dev/null 2>&1

# bam NAME [CONTIG...]: SAM records on stdin (space-separated) -> sorted, indexed NAME.bam.
# The header lists every contig in FASTA order, or only the given contigs in the given order.
bam() {
  local name="$1"; shift
  { printf '@HD\tVN:1.6\tSO:unsorted\n'
    if [ $# -eq 0 ]; then set -- $(cut -f1 ref.fa.fai); fi
    for contig in "$@"; do
      awk -v c="$contig" '$1 == c {printf "@SQ\tSN:%s\tLN:%s\n", $1, $2}' ref.fa.fai
    done
    tr ' ' '\t'; } | samtools sort -o "$name.bam" - 2>/dev/null && samtools index "$name.bam"
}
# refseq CONTIG START END: reference bases START-END (1-based, inclusive), for read SEQ
refseq() { samtools faidx ref.fa "$1:$2-$3" | tail -n +2 | tr -d '\n'; }
# quals N: N quality characters of Phred 40
quals() { printf 'I%.0s' $(seq "$1"); }
# chrS_reads PREFIX N FLAG YD SEQ: N copies of a 40M read over chrS 1-40, named PREFIX1..PREFIXN
chrS_reads() {
  local i; for i in $(seq "$2"); do echo "$1$i $3 chrS 1 60 40M * 0 0 $5 $(quals 40) YD:Z:$4"; done
}
# epirows FILE: one example epiBED row per distinct (strand, CG string, SNP string), as
# read name, strand, CG string, SNP string
epirows() { awk -F'\t' '!seen[$6 FS $7 FS $9]++ {print $4, $6, $7, $9}' "$1"; }
# pileup NAME [OPTIONS]: VCF body rows (CHROM POS REF ALT INFO SAMPLE). Filters are
# BISCUIT's defaults (MAPQ >= 40, base quality >= 20) except that the -5/-3 end
# filters are off, so every aligned base can be called.
pileup() {
  local name="$1"; shift
  biscuit pileup -5 0 -3 0 "$@" ref.fa "$name.bam" 2>/dev/null | grep -v '^#' | cut -f1,2,4,5,8,10 || true
}
# alnbase_hits NAME: C_retained / C_converted hits as NAME_*.parquet
alnbase_hits() {
  "$ALNBASE" query --query-file "$HERE/c.toml" --parquet -F qname "$1.bam" ref.aref "$1.parquet" >/dev/null 2>&1
}
duck() { "$PY" -c 'import duckdb, sys; duckdb.sql(sys.argv[1]).show()' "$1"; }
Q25=IIIIIIIIIIIIIIIIIIIIIIIII

say "BC-2: hard clip shifts the SEQ index (pileup)"
# Both reads align 25 bases at chrH 1-25 and retain the C at 11. hc drops its first 5
# bases as a hard clip (not in SEQ); sc keeps them in SEQ as a soft clip.
bam bc2_hc <<SAM
hc 0 chrH 1 60 5H25M * 0 0 TTTTTTTTTTCGTTTTTTTTTTTTT $Q25 YD:Z:f
SAM
bam bc2_sc <<SAM
sc 0 chrH 1 60 5S25M * 0 0 AAAAATTTTTTTTTTCGTTTTTTTTTTTTT IIIII$Q25 YD:Z:f
SAM
promise "SAM spec: H does not consume SEQ. Both reads retain the C at chrH 11, so CV=1 and"
promise "BT=1.000, and each read matches the reference everywhere else (no ALT rows)."
shows "BISCUIT, soft clip (control)"
pileup bc2_sc
shows "BISCUIT, hard clip"
pileup bc2_hc
alnbase_hits bc2_hc
shows "alnbase, hard clip (refr_pos is 0-based; off_5p counts the hard-clipped bases)"
duck "select qname, name, refr_pos, off_5p from 'bc2_hc_*.parquet' order by refr_pos"

say "BC-3: hard clip aborts epiread"
promise "one epiBED row for read hc."
shows "BISCUIT epiread on the BC-2 hard-clip BAM"
status=0
# The inner shell reports "Aborted (core dumped)" with a PID; its stderr is discarded.
bash -c 'biscuit epiread -5 0 -3 0 ref.fa bc2_hc.bam > bc3.epibed 2> bc3.err; exit $?' 2>/dev/null || status=$?
cat bc3.epibed; grep -v '^\[' bc3.err || true
echo "exit status $status"

say "BC-8: retention count uses the wrong strand (pileup -t)"
# Both reads are OT (YD:Z:f) over 10 reference Cs at 21-30. conv converts all ten (T);
# ret retains all ten (C). chrG has 10 Gs before the Cs; chrR has none.
bam bc8 <<SAM
conv 0 chrG 1 60 40M * 0 0 AAAAAAAAAAGGGGGGGGGGTTTTTTTTTTAAAAAAAAAA IIIIIIIIIIIIIII$Q25 YD:Z:f
ret 0 chrR 1 60 40M * 0 0 AAAAAAAAAAAAAAAAAAAACCCCCCCCCCAAAAAAAAAA IIIIIIIIIIIIIII$Q25 YD:Z:f
SAM
promise "help: '-t INT Maximum cytosine retention in a read'; docs: retained C's for OT reads."
promise "conv retains 0 Cs, so -t 5 keeps it; ret retains 10, so -t 5 drops it."
for t in 999999 5; do
  shows "BISCUIT -t $t, C calls per contig (chrG is read conv, chrR is read ret)"
  pileup bc8 -t "$t" | awk '$5 ~ /CX=/ {n[$1]++} END {printf "chrG %d  chrR %d\n", n["chrG"], n["chrR"]}'
done
alnbase_hits bc8
shows "alnbase, retained and converted reference Cs per read"
duck "select qname, count(*) filter (where name = 'C_retained') as retained,
             count(*) filter (where name = 'C_converted') as converted
      from 'bc8_*.parquet' group by qname order by qname"

say "BC-18: last base of a contig and region end (pileup)"
bam bc18 <<SAM
e1 0 chrE 1 60 30M * 0 0 AAAAAAAAAAAAAAAAAAACAAAAAAAAAC IIIII$Q25 YD:Z:f
SAM
promise "a row for every covered C: chrE 20 and chrE 30 (the contig's last base)."
promise "help: '-g STR Region'; a samtools-style region chrE:1-20 includes position 20."
shows "BISCUIT, whole BAM"
pileup bc18
shows "BISCUIT -g chrE:1-20"
pileup bc18 -g chrE:1-20
shows "BISCUIT -g chrE:1-21"
pileup bc18 -g chrE:1-21
alnbase_hits bc18
shows "alnbase (refr_pos is 0-based)"
duck "select name, refr_pos from 'bc18_*.parquet' order by refr_pos"

say "BC-10: -5/-3 measured along SEQ, not from the read's 5' end (pileup, epiread)"
# A single-end reverse-strand read (OB, YD:Z:r) over chrF 1-60 retaining both CG Gs, at 6
# and 56. Its 5' end is the right end: G 56 is 5 bases from it, G 6 is 55 bases from it.
bam bc10 <<SAM
rev 16 chrF 1 60 60M * 0 0 $(refseq chrF 1 60) $(quals 60) YD:Z:r
SAM
promise "help: '-5 INT Minimum distance to 5' end of a read', '-3 INT ... 3' end'."
promise "-5 10 drops G 56 (near the 5' end) and keeps G 6; -3 10 does the opposite."
for ends in "-5 10 -3 0" "-5 0 -3 10"; do
  shows "BISCUIT pileup $ends"
  biscuit pileup $ends ref.fa bc10.bam 2>/dev/null | grep -v '^#' | cut -f1,2,4,5,8,10 || true
done
shows "BISCUIT epiread -5 10 -3 0 (CG string from chrF 1: F = filtered, M = methylated)"
biscuit epiread -5 10 -3 0 ref.fa bc10.bam 2>/dev/null
alnbase_hits bc10
shows "alnbase, off_5p and off_3p from the read's 5' and 3' ends as sequenced"
duck "select qname, name, refr_pos, off_5p, off_3p from 'bc10_*.parquet' order by refr_pos"

say "BC-11: mate overlap without an MC tag (pileup)"
# Bismark-style pairs (XG:Z:CT, no MC tag). Pair p: R1 chrP 101-130 (30M), R2 111-190 (80M);
# the reads overlap at 111-130 and only R2 covers the C at 150. Pair q: R1 101-180 (80M),
# R2 151-180 (30M); both reads cover the C at 160.
bam bc11 <<SAM
p 99 chrP 101 60 30M = 111 90 $(refseq chrP 101 130) $(quals 30) XG:Z:CT
p 147 chrP 111 60 80M = 101 -90 $(refseq chrP 111 190) $(quals 80) XG:Z:CT
SAM
bam bc11_mc <<SAM
p 99 chrP 101 60 30M = 111 90 $(refseq chrP 101 130) $(quals 30) XG:Z:CT MC:Z:80M
p 147 chrP 111 60 80M = 101 -90 $(refseq chrP 111 190) $(quals 80) XG:Z:CT MC:Z:30M
SAM
bam bc11_rev <<SAM
q 99 chrP 101 60 80M = 151 80 $(refseq chrP 101 180) $(quals 80) XG:Z:CT
q 147 chrP 151 60 30M = 101 -80 $(refseq chrP 151 180) $(quals 30) XG:Z:CT
SAM
promise "help: '-d Double count cytosines in overlapping mate reads (avoided by default)'."
promise "Pair p: R2 alone covers 150, so CV=1 at 150 (and at 160). Pair q: one fragment, CV=1 at 160."
shows "BISCUIT, pair p without MC"
pileup bc11
shows "BISCUIT, pair p with MC (control)"
pileup bc11_mc
shows "BISCUIT, pair q without MC"
pileup bc11_rev

say "BC-13: overlap test ignores the mate's contig (pileup -p)"
# An improperly paired R2 (FLAG 177) on chrI 101-150 covering the C at 120. Its mate is at
# position 101 of chrJ; in the control it is at chrI 1001, not overlapping.
bam bc13 <<SAM
x 177 chrI 101 60 50M chrJ 101 0 $(refseq chrI 101 150) $(quals 50) XG:Z:CT MC:Z:50M
SAM
bam bc13_ctrl <<SAM
x 177 chrI 101 60 50M = 1001 0 $(refseq chrI 101 150) $(quals 50) XG:Z:CT MC:Z:50M
SAM
promise "-p keeps improper pairs; a mate on another contig cannot overlap, so CV=1 at chrI 120."
shows "BISCUIT -p, mate on chrJ"
pileup bc13 -p
shows "BISCUIT -p, mate on chrI 1001 (control)"
pileup bc13_ctrl -p

say "BP-5: multi-BAM pileup uses the first BAM's contig order for every BAM"
# Two BAMs over the same reference with the contigs listed in opposite orders. Only the
# second has a read: chrI 101-150, retaining the C at 120.
bam bp5_a chrI chrJ <<SAM
SAM
bam bp5_b chrJ chrI <<SAM
y 0 chrI 101 60 50M * 0 0 $(refseq chrI 101 150) $(quals 50) YD:Z:f
SAM
promise "usage: 'biscuit pileup [options] <ref.fa> <in1.bam> [in2.bam ...]', one sample column"
promise "per BAM; sample bp5_b has CV=1 at chrI 120, whichever order the BAMs are given in."
for order in "bp5_a.bam bp5_b.bam" "bp5_b.bam bp5_a.bam"; do
  shows "BISCUIT pileup ref.fa $order (columns CHROM POS REF ALT INFO, then one per BAM)"
  biscuit pileup -5 0 -3 0 ref.fa $order 2>/dev/null | grep -v '^#' | cut -f1,2,4,5,8,10- || true
done

# chrS read sequences: the reference (C at 21 retained), the C at 21 read as T (converted
# on a C-strand read, or a C>T SNP on a G-strand read).
CHRS_C=$(refseq chrS 1 40)
CHRS_T=$(refseq chrS 1 20)T$(refseq chrS 22 40)

say "BE-2: epiread -B ignores a multi-sample SNP BED"
# One BAM with a heterozygous C>T SNP at chrS 21: 10 C-strand reads retain the C, 10
# G-strand reads show T. It is piled up twice as a two-sample VCF.
{ chrS_reads m 10 0 f "$CHRS_C"; chrS_reads s 10 16 r "$CHRS_T"; } | bam be2
biscuit pileup -5 0 -3 0 -o be2.vcf ref.fa be2.bam be2.bam 2>/dev/null
biscuit vcf2bed -s ALL -t snp be2.vcf > be2_all.snp.bed 2>/dev/null
cut -f1-9 be2_all.snp.bed > be2_first9.snp.bed
promise "docs (methylextraction.md): for -t snp 'Columns 6-9 are repeated for each sample found in"
promise "the VCF file'; epiread help: -B takes a BISCUIT SNP BED (pileup -> vcf2bed -t snp), and"
promise "epibed_format.md: a C with a T allele at frequency >= 0.05 is not callable. So with either"
promise "BED, the SNP string shows the C and the CG string has no M at 21 for the m reads."
shows "vcf2bed -s ALL -t snp (13 columns)"
cat be2_all.snp.bed
for bed in none be2_all.snp.bed be2_first9.snp.bed; do
  shows "BISCUIT epiread with -B $bed (read, strand, CG string, SNP string from chrS 1)"
  if [ "$bed" = none ]; then biscuit epiread -5 0 -3 0 ref.fa be2.bam > be2.epibed 2>/dev/null
  else biscuit epiread -5 0 -3 0 -B "$bed" ref.fa be2.bam > be2.epibed 2>/dev/null; fi
  epirows be2.epibed
done

say "BE-3: SNP veto thresholds differ between pileup and epiread"
# 21 C-strand reads retain the C at chrS 21; 1 G-strand read shows T there. T/C = 1/21 =
# 0.048, T/(C+T) = 1/22 = 0.045.
{ chrS_reads m 21 0 f "$CHRS_C"; chrS_reads s 1 16 r "$CHRS_T"; } | bam be3
biscuit pileup -5 0 -3 0 -o be3.vcf ref.fa be3.bam 2>/dev/null
biscuit vcf2bed -t snp be3.vcf > be3.snp.bed 2>/dev/null
promise "epiread help: 'with an unfiltered BISCUIT SNP BED file ... should have the same results in"
promise "the epiBED as in pileup'. pileup calls the C at 21 (CV=21), so epiread gives M at 21."
shows "BISCUIT pileup (SP C21T1, AF1 = 1/22 printed to 2 decimals, then CV and BT)"
grep -v '^#' be3.vcf | awk '$2 == 21' | cut -f1,2,4,5,10
shows "vcf2bed -t snp"
cat be3.snp.bed
shows "BISCUIT epiread -B (read, strand, CG string, SNP string from chrS 1)"
biscuit epiread -5 0 -3 0 -B be3.snp.bed ref.fa be3.bam > be3.epibed 2>/dev/null
epirows be3.epibed

say "BP-4: vcf2bed -c and mergecg -c rebuild counts from a 3-decimal beta"
# 2001 C-strand reads over the C at chrS 21: 1 retains it, 2000 convert it.
{ chrS_reads m 1 0 f "$CHRS_C"; chrS_reads c 2000 0 f "$CHRS_T"; } | bam bp4
biscuit pileup -5 0 -3 0 -o bp4.vcf ref.fa bp4.bam 2>/dev/null
promise "help: vcf2bed '-c Output Beta-M-U', mergecg '-c Output Beta-M-U'; M=1 and U=2000."
shows "BISCUIT pileup (SP C1Y2000, then CV and BT)"
grep -v '^#' bp4.vcf | awk '$2 == 21' | cut -f1,2,10
shows "vcf2bed -c -t cg (chrom, start, end, percent, M, U)"
biscuit vcf2bed -c -t cg bp4.vcf 2>/dev/null
shows "vcf2bed -t cg | mergecg -c (chrom, start, end, percent, M, U, strand detail)"
biscuit vcf2bed -t cg bp4.vcf 2>/dev/null > bp4.bed
biscuit mergecg -c ref.fa bp4.bed 2>/dev/null

say "BP-8: mergecg exits on a row at a contig's first base"
echo "q1 0 chrQ 1 60 40M * 0 0 $(refseq chrQ 1 40) $(quals 40) YD:Z:f" | bam bp8
biscuit pileup -5 0 -3 0 -o bp8.vcf ref.fa bp8.bam 2>/dev/null
biscuit vcf2bed -t c bp8.vcf > bp8.bed 2>/dev/null
promise "mergecg help: in.bed is vcf2bed output; the C at chrQ 1 is written as one row."
shows "vcf2bed -t c"
cat bp8.bed
shows "BISCUIT mergecg"
status=0
biscuit mergecg ref.fa bp8.bed 2>&1 || status=$?
echo "exit status $status"
