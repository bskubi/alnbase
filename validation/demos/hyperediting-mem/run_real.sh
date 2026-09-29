#!/usr/bin/env bash
# The pipeline on a real library, over a subset of the genome.
#
#   run_real.sh WORKDIR NAME READS [CHROMS]
#
#   WORKDIR   holds hg38.analysisSet.fa, the whole-genome bwa index, and
#             NAME/ from prep_chrom.sh
#   NAME      the subset the three-letter stage runs against
#   READS     the FASTQ, Phred+33, single-end
#   CHROMS    the contigs in the subset, comma-separated; defaults to NAME
#
# Two references are in play and the split between them is the whole design.
#
#   Stage 1 runs against the WHOLE GENOME. It is the pass that decides which
#   reads ordinary alignment cannot place, and if it only saw one chromosome it
#   would declare most of the library unmappable and hand the three-letter stage
#   a pool made of nothing but other chromosomes' reads. So stage 1 uses the
#   real hg38 index, exactly as the published pipeline does.
#
#   Stage 2 runs against NAME. This is the subset, and the reason for it is
#   cost: hisat2-build on hg38 is hours per transform and there are six.
#
# The reads that survive stage 1 come from everywhere, so some of them will find
# a spurious home in it. The pipeline's own controls measure exactly that.
# Of its twelve transform combinations only a2g and t2c can carry an A-to-I
# event; the rest are noise channels by construction, and their yield here is
# the false-positive floor. Every number below is reported against it.
#
# Environment: PATH must carry bwa, hisat2, samtools, perl and xa2multi.pl;
# ALNBASE the binary; HE the patched Hyper-editing tree.
set -euo pipefail

WORK="${1:?usage: run_real.sh WORKDIR NAME READS [CHROMS]}"
NAME="${2:?usage: run_real.sh WORKDIR NAME READS [CHROMS]}"
READS="${3:?usage: run_real.sh WORKDIR NAME READS [CHROMS]}"
CHROMS="${4:-$NAME}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HE="${HE:?set HE to the patched Hyper-editing tree}"
ALNBASE="${ALNBASE:-alnbase}"
THREADS="${THREADS:-16}"
PY="${PY:-python3}"

WORK="$(cd "$WORK" && pwd)"
READS="$(cd "$(dirname "$READS")" && pwd)/$(basename "$READS")"
REF="$WORK/hg38.analysisSet.fa"
C="$WORK/$NAME/$NAME"
OUT="$WORK/run_$NAME"
mkdir -p "$OUT"/{stage1,lists}
cd "$OUT"

# --- stage 1: what ordinary alignment cannot place ---------------------------
# The published pre_unmapped.sh does bwa aln, then bwa mem -k 50 over what aln
# left behind, and keeps the reads neither could place. It is reproduced here
# rather than called because of two hardcoded values that only cost time on a
# real library. `-t 5` becomes $THREADS. And `samse -n 50`, which asks for up to
# fifty suboptimal hits in the XA tag, becomes plain samse: samse is
# single-threaded and XA is the expensive part of it, while the only thing read
# off this alignment is the unmapped flag. Nothing downstream sees XA.
if [ ! -s stage1/um.fastq ]; then
  echo "[$(date +%T)] stage 1: bwa aln against the whole genome"
  bwa aln -t "$THREADS" "$WORK/hg38" "$READS" -f stage1/aln.sai 2> stage1/aln.log
  bwa samse "$WORK/hg38" stage1/aln.sai "$READS" 2>> stage1/aln.log \
    | samtools view -b -o stage1/aln.bam -
  samtools fastq -f 4 stage1/aln.bam > stage1/aln.um.fastq 2>/dev/null
  echo "[$(date +%T)] stage 1: bwa mem -k 50 over what aln left behind"
  bwa mem -M -t "$THREADS" -k 50 "$WORK/hg38" stage1/aln.um.fastq 2> stage1/mem.log \
    | samtools view -b -o stage1/mem.bam -
  samtools fastq -f 4 stage1/mem.bam > stage1/um.fastq 2>/dev/null
  rm -f stage1/aln.sai stage1/aln.um.fastq
fi

# Reads that stage 1 placed inside the subset. These are the library's ordinary
# reads for it, and they are the input to junction discovery below.
#
# Both alignments have to be read, and the second one is the important one.
# `bwa aln` is gapless and aligns end to end, so a read crossing a junction
# cannot come out of it at all -- which is to say that the reads junction
# discovery exists to see are precisely the ones missing from aln.bam. They are
# in mem.bam, soft-clipped at the junction, because `bwa mem -k 50` places a
# read on the strength of one anchor. Taking only aln.bam here would have built
# the first-pass junction list out of the reads that have no junctions in them.
if [ ! -s "stage1/$NAME.fastq" ]; then
  echo "[$(date +%T)] stage 1: reads placed in $CHROMS"
  # By RNAME rather than by region, so neither BAM needs a sort or an index.
  for b in stage1/aln.bam stage1/mem.bam; do
    samtools view -h -F 4 "$b" \
      | awk -v cs="$CHROMS" 'BEGIN{n=split(cs,a,","); for(i=1;i<=n;i++) keep[a[i]]=1}
                             $1 ~ /^@/ || ($3 in keep)' \
      | samtools fastq - 2>/dev/null
  done > "stage1/$NAME.fastq"
fi

printf "stage 1: %s reads in, %s unplaced, %s on %s\n" \
  "$(( $(wc -l < "$READS") / 4 ))" \
  "$(( $(wc -l < stage1/um.fastq) / 4 ))" \
  "$(( $(wc -l < "stage1/$NAME.fastq") / 4 ))" "$CHROMS" | tee stage1/counts.txt

# --- the three junction lists ------------------------------------------------
# The simulation said the list should be the union of the annotation and the
# sample's own first pass, because the cost of a splice-site list is entirely in
# what it omits. This is that claim at real annotation density, where GENCODE
# already carries every transcript's junctions and the gap it leaves is much
# narrower than three simulated introns suggested.
if [ ! -s lists/firstpass.txt ]; then
  echo "[$(date +%T)] first pass: junctions from the library's own $NAME reads"
  # Four letters, untransformed reference, GT..AG intact -- the easy regime.
  # --novel-splicesite-outfile is hisat2's own de novo junction report.
  hisat2 -p "$THREADS" -x "$C.ht2" -U "stage1/$NAME.fastq" \
    --novel-splicesite-outfile lists/firstpass.raw \
    -S /dev/null 2> lists/firstpass.log
  sort -k1,1 -k2,2n -k3,3n -u lists/firstpass.raw > lists/firstpass.txt
fi
cp -f "$WORK/$NAME/gencode.$NAME.txt" lists/gencode.txt
sort -k1,1 -k2,2n -k3,3n -u lists/gencode.txt lists/firstpass.txt > lists/union.txt

printf "junction lists: gencode %s, first pass %s, union %s (overlap %s)\n" \
  "$(wc -l < lists/gencode.txt)" "$(wc -l < lists/firstpass.txt)" \
  "$(wc -l < lists/union.txt)" \
  "$(( $(wc -l < lists/gencode.txt) + $(wc -l < lists/firstpass.txt) - $(wc -l < lists/union.txt) ))" \
  | tee lists/counts.txt

echo "[$(date +%T)] stage 1 and the junction lists are done"

# --- stage 2: the arms -------------------------------------------------------
# Six of them. Two settle the aligner question on real data; four settle where
# the junction list should come from, sharing everything except the list.
#
# Each arm is a whole pipeline -- transform alignments, edit calling, detection
# -- and arms are independent, so they run concurrently. On a real library this
# is not a nicety. `Get_orig_read.pl` and `filter_sam.pl` are single-threaded
# perl over the entire unmapped pool, twelve times per arm, and they dominate
# everything else; run one arm at a time and six arms is most of a day. The
# limit is memory rather than cores, because each concurrent arm holds its own
# copy of the read pool in a perl hash, so $ARMS_AT_ONCE is small by default.
ARMS_AT_ONCE="${ARMS_AT_ONCE:-3}"

"$PY" "$HERE/bin/patch_analyse_mm.py" "$HE/analyse_mm.pl" "$OUT/analyse_mm_cigar.pl"
"$PY" "$HERE/bin/patch_detect_ue.py"  "$HE/detect_ue.pl"  "$OUT/detect_ue_ported.pl"
chmod +x "$OUT"/analyse_mm_cigar.pl "$OUT"/detect_ue_ported.pl
UE_ARGS="0.05 0.6 30 0.6 0.1 0.8 0.2"
STRAND_RULE="${STRAND_RULE:-$HERE/../../../queries/strand/unconverted.toml}"
mkdir -p stage2

# A marker file rather than the output directory, because the directory appears
# when the stage starts. Keying the skip on it would make an interrupted run
# look finished to the next one.
arm () {
  local arm=$1 aligner=$2 caller=$3 splice=${4:-}
  local trans="$OUT/stage2/$arm" out="$OUT/arm_$arm"
  [ -f "$out/.done" ] && { echo "[$(date +%T)] $arm already done"; return; }

  if [ ! -f "$trans/.done" ]; then
    rm -rf "$trans"
    echo "[$(date +%T)] $arm: twelve transform alignments"
    if [ "$aligner" = aln ]; then
      "$HE/TransRun.sh" "$C" "$OUT/stage1/um.fastq" "$trans" real \
        "$trans.stat" bwa samtools > "stage2/$arm.log" 2>&1
    else
      "$HERE/bin/TransRunAligner.sh" "$aligner" "$C" "$OUT/stage1/um.fastq" \
        "$trans" real "$trans.stat" "$HE" $splice > "stage2/$arm.log" 2>&1
    fi
    touch "$trans/.done"
  fi

  rm -rf "$out"; mkdir -p "$out/beds"; : > "$out/hits.txt"
  echo "[$(date +%T)] $arm: calling edits"
  for bam in "$trans"/bamFiles/*FilterReads*.bam; do
    case "$caller" in
      published) "$OUT/analyse_mm_cigar.pl" "$C.fa" "$bam" samtools ;;
      alnbase)   "$PY" "$HERE/bin/alnbase_mm.py" "$C.fa" "$bam" \
                   --aref "$C.aref" --strand-rule "$STRAND_RULE" \
                   --alnbase "$ALNBASE" --work "$out/alnbase_work" ;;
    esac >> "$out/hits.txt" 2>> "$out/call.log"
  done
  "$HE/sort_R_read.pl" "$out/hits.txt" "$out/sort.stat" > "$out/real.analyseMM"
  "$OUT/detect_ue_ported.pl" "$out/real.analyseMM" "$out/beds/real" \
    "$out/detect.stat" $UE_ARGS samtools 0 > "$out/detect.log" 2>&1
  touch "$out/.done"
  echo "[$(date +%T)] $arm: finished"
}

# aligner, caller, and for the hisat2 arms the junction list under test.
ARMS=(
  "aln-published    aln     published"
  "mem-alnbase      mem     alnbase"
  "hisat2-none      hisat2  alnbase"
  "hisat2-gencode   hisat2  alnbase  $OUT/lists/gencode.txt"
  "hisat2-firstpass hisat2  alnbase  $OUT/lists/firstpass.txt"
  "hisat2-union     hisat2  alnbase  $OUT/lists/union.txt"
)
# ARMS_ONLY narrows the run to named arms. The `mem` arm is the reason it
# exists: `bwa mem` soft-clips, so on a subset reference it places reads from
# outside the subset rather than rejecting them, and filter_sam.pl then grinds
# through the extra records for hours. That makes the aligner comparison
# subset-confounded in a way the junction-list comparison is not -- those four
# arms are all HISAT2 and differ only in the list -- so on a subset it is
# usually right to leave `mem` out and run it when the whole genome is indexed.
for spec in "${ARMS[@]}"; do
  if [ -n "${ARMS_ONLY:-}" ]; then
    case " $ARMS_ONLY " in *" ${spec%% *} "*) ;; *) continue ;; esac
  fi
  while [ "$(jobs -rp | wc -l)" -ge "$ARMS_AT_ONCE" ]; do wait -n || true; done
  # shellcheck disable=SC2086
  arm $spec &
done
wait || true

# --- report ------------------------------------------------------------------
# There is no ground truth here, so the comparison is against the pipeline's own
# controls. Of the six mismatch classes detect_ue.pl reports, only A2G can carry
# an A-to-I event. The other five are noise channels by construction, and their
# yield is the floor that every signal number has to be read against.
"$PY" "$HERE/score_real.py" "$OUT" --out "$OUT/real_$NAME.tsv"
