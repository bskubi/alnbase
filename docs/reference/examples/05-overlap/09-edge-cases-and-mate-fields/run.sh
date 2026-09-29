#!/usr/bin/env bash
# alnbase reference, section 05 (overlap), example 09-edge-cases-and-mate-fields. See README.txt.
set -uo pipefail
ALNBASE=${ALNBASE:-alnbase}
cd "$(dirname "$0")"
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
samtools view --no-PG -b -o "$T/in.bam" reads.sam
# show BAM: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN, SEQ length, then MC/MQ/SA tags.
show() {
  samtools view "$1" | awk -F'\t' 'BEGIN{OFS="\t"} {t=""; for(i=12;i<=NF;i++) if($i ~ /^(MC|MQ|SA):/) t=t"\t"$i;
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,"SEQ "length($10)t}'
}
# key BAM: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN and the MC tag, one line per record.
key() {
  samtools view "$1" | awk -F'\t' 'BEGIN{OFS="\t"} {mc="-"; for(i=12;i<=NF;i++) if($i ~ /^MC:/) mc=$i;
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,mc}'
}
# fixmate_check BAM: run samtools fixmate -m on BAM and show every one of those fields it changed.
fixmate_check() {
  samtools fixmate -m "$1" "$T/fm.bam"
  if diff <(key "$1") <(key "$T/fm.bam") > "$T/diff"; then
    echo "  fixmate -m changed none of FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN MC"
  else
    sed 's/^/  /' "$T/diff"
  fi
}
run() {
  local label=$1; shift
  echo "### $label"
  echo "\$ alnbase overlap${*:+ $*} in.bam out.bam"
  rm -f "$T/out.bam"
  "$ALNBASE" overlap "$@" "$T/in.bam" "$T/out.bam" 2>"$T/err"; echo "exit status: $?"
  sed 's/^/stderr| /' "$T/err"
  show "$T/out.bam"
  echo "--- samtools fixmate -m on this output"
  fixmate_check "$T/out.bam"
  echo "--- samtools flagstat (primary lines, read 1, read 2, singletons)"
  samtools flagstat "$T/out.bam" | grep -E ' primary$| read1$| read2$|singletons'
  echo
}
echo "### input"
show "$T/in.bam"
echo
run "defaults"
run "--keep r2 (read_through: R1 is the read left unmapped)" --keep r2
run "--clip-mode hard (the promoted primary keeps no hard clip it did not make)" --clip-mode hard
echo "### unmapped read in full (defaults): SEQ reverse-complemented back to sequencing orientation"
"$ALNBASE" overlap "$T/in.bam" "$T/out.bam" 2>/dev/null
samtools view "$T/in.bam" | awk -F'\t' '$1=="read_through" && $2==147 {print "input : " $2, $6, $10}'
samtools view "$T/out.bam" | awk -F'\t' '$1=="read_through" && and($2,4) {print "output: " $2, $6, $10}'
