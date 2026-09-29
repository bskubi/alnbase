#!/usr/bin/env bash
# Build a subset of the genome and all its indexes, so a real-data run finishes
# today.
#
#   prep_chrom.sh WORKDIR NAME CHROMS [GTF]
#
#   NAME     what to call the subset; also the directory it lands in
#   CHROMS   comma-separated contig names, or one name
#
# `prep_real.sh` builds thirteen whole-genome indexes. bwa costs about half an
# hour each and hisat2-build costs hours each, so the whole set is most of a day
# and around 60 Gb. chr21 alone is four minutes and 1.3 Gb, and every arm of the
# comparison still sees exactly the same reference, which is what the comparison
# needs.
#
# Which chromosomes to take is a signal question, not a size one. Hyper-editing
# happens in inverted Alu pairs in expressed transcripts, and Alu density varies
# several-fold across the karyotype -- chr19 is the densest in the genome, chr21
# and chr18 among the sparsest. The subset used here is every chromosome shorter
# than 100 Mb (chr16 through chr22), which is a simple rule that happens to pick
# out the gene- and Alu-dense end: 15% of the genome for a good deal more than
# 15% of the editing.
#
# What the subset costs, stated plainly: a read whose true locus is outside it
# has nowhere correct to go, and some of those reads will find a
# spurious home. Two things keep that honest.
#
#   Stage 1 still runs against the whole genome. The untransformed bwa index in
#   WORKDIR is the real one, so the pool of unmapped reads that reaches the
#   three-letter stage is the pool the pipeline would really see -- not a pool
#   manufactured by hiding most of the genome from the first pass.
#
#   The pipeline's own controls measure what gets through. Of its twelve
#   transform combinations only a2g and t2c can carry an A-to-I event; a2c, a2t,
#   g2c and the rest are noise channels by construction. Their yield on this
#   reference is the false-positive floor, and it is reported alongside every
#   number rather than assumed away.
#
# So the subset supports comparisons between arms, which is what is being asked
# here, and not absolute editing rates for the library.
set -euo pipefail

WORK="${1:?usage: prep_chrom.sh WORKDIR NAME CHROMS [GTF]}"
NAME="${2:?usage: prep_chrom.sh WORKDIR NAME CHROMS [GTF]}"
CHROMS="${3:-$NAME}"
GTF="${4:-}"
HE="${HE:-$(dirname "$WORK")/Hyper-editing}"
ALNBASE="${ALNBASE:-alnbase}"
THREADS="${THREADS:-8}"
REF="$WORK/hg38.analysisSet.fa"
OUT="$WORK/$NAME"

mkdir -p "$OUT"
cd "$OUT"

[ -f "$REF.fai" ] || samtools faidx "$REF"

if [ ! -f "$NAME.fa" ]; then
  echo "[$(date +%T)] extracting $CHROMS"
  # shellcheck disable=SC2046
  samtools faidx "$REF" $(echo "$CHROMS" | tr ',' ' ') > "$NAME.fa"
  samtools faidx "$NAME.fa"
fi

# The untransformed reference, for the splice-aware first pass that proposes
# junctions. This is the pass that runs in four letters with GT..AG intact, so
# it is the one that can find a junction the annotation does not have.
[ -f "$NAME.bwt" ]      || { echo "[$(date +%T)] bwa index $NAME";      bwa index -p "$NAME" "$NAME.fa" 2>/dev/null; }
[ -f "$NAME.ht2.1.ht2" ] || { echo "[$(date +%T)] hisat2-build $NAME"; hisat2-build -q -p "$THREADS" "$NAME.fa" "$NAME.ht2"; }
if [ ! -f "$NAME.aref" ] && command -v "$ALNBASE" >/dev/null 2>&1; then
  echo "[$(date +%T)] alnbase index $NAME"
  "$ALNBASE" index "$NAME.fa" "$NAME.aref"
fi

for tt in a2g t2c a2c t2g g2c a2t; do
  [ -f "$NAME.$tt.fa" ] || "$HE/$tt.pl" "$NAME.fa" > "$NAME.$tt.fa"
  [ -f "$NAME.$tt.bwt" ]       || { echo "[$(date +%T)] bwa index $tt";      bwa index -p "$NAME.$tt" "$NAME.$tt.fa" 2>/dev/null; }
  [ -f "$NAME.$tt.ht2.1.ht2" ] || { echo "[$(date +%T)] hisat2-build $tt";   hisat2-build -q -p "$THREADS" "$NAME.$tt.fa" "$NAME.$tt.ht2"; }
done

# The annotation arm of the junction-list comparison. hisat2_extract_splice_sites.py
# walks every transcript in the GTF and emits every adjacent exon pair, so this
# is not one isoform per gene -- annotated skipping and alternative 5'/3' sites
# are already in here. That is the point of using the real annotation rather
# than the three-intron stand-in the simulation used.
if [ -n "$GTF" ] && [ ! -f "gencode.$NAME.txt" ]; then
  echo "[$(date +%T)] gencode splice sites for $CHROMS"
  hisat2_extract_splice_sites.py <(zcat "$GTF") \
    | awk -v cs="$CHROMS" 'BEGIN{n=split(cs,a,","); for(i=1;i<=n;i++) keep[a[i]]=1}
                           $1 in keep' \
    | sort -k1,1 -k2,2n -k3,3n -u > "gencode.$NAME.txt"
fi

echo "[$(date +%T)] done"
[ -f "gencode.$NAME.txt" ] && echo "  gencode junctions: $(wc -l < "gencode.$NAME.txt")"
du -sh "$OUT"
