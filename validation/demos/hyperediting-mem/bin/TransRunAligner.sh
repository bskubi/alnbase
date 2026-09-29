#!/usr/bin/env bash
# The pipeline's transform-and-align stage, with the aligner as a parameter.
#
#   TransRunAligner.sh ALIGNER GENOME_PREFIX READS_FASTQ TRANSRUN_DIR OUT_PRE \
#                      STATFILE HE_DIR [SPLICESITES]
#
#   ALIGNER   mem     bwa mem
#             hisat2  HISAT2, splice-aware
#
# Same twelve read/genome transform pairs as the published TransRun.sh, same
# output file names, so the stages after it do not know the difference. What
# changes is the aligner and, with it, which reads can be placed at all.
#
# Why not `bwa aln`, which the pipeline uses: it is given `-o 0`, forbidding
# gaps, and it aligns end to end. That is what makes the published edit
# caller's flat window correct, and it is also why no junction-spanning,
# indel-bearing or clipped read survives to be examined.
#
# Notes that apply to both aligners here:
#
#   No `grep -v chrUn`. The original drops unplaced contigs with three greps
#   over the SAM text, which also deletes any read whose name or sequence
#   contains the string. On a real genome, filter the header instead.
#
#   Supplementary and secondary alignments keep their full SEQ. Get_orig_read.pl
#   substitutes the whole original read into every record, which is only valid
#   if the CIGAR consumes the whole read. `bwa mem -Y` soft-clips instead of
#   hard-clipping for exactly this reason; HISAT2 does not hard-clip.
#
# bwa mem only:
#
#   xa2multi.pl, from bwa's own distribution, expands the alternative positions
#   in XA into real records. `bwa aln` reports them in XA too and analyse_mm.pl
#   walks them, but mem's XA entries each carry their own CIGAR, which a flat
#   window cannot use. As records they are on the same footing as any other
#   alignment, and the >51-hit rejection in sort_R_read.pl keeps working: it
#   sums hits per read name and each record contributes one.
#
# HISAT2 only:
#
#   --known-splicesite-infile is not optional here. None of the six three-letter
#   transforms leaves GT..AG intact -- a2g turns AG into GG, t2c turns GT into
#   GC -- so a transformed genome has no canonical splice motifs. Junctions can
#   still be found without them, which is measured in probe_splice_motif.py: the
#   motif is worth about three bases of anchor, not the difference between
#   working and not. The reason the list is not optional is the other result
#   from the same probe. De novo splicing on a collapsed genome placed 499 gaps
#   at a wrong intron across 9,000 reads, and with the junctions supplied, zero
#   -- and a misplaced junction in this pipeline is not a visible gap, it is a
#   read aligned against the wrong sequence and therefore a burst of mismatches
#   that looks exactly like dense editing.
#
#   Where the list should come from is measured in probe_junction_list_source.py:
#   the union of the annotation and the sample's own first pass, because a list's
#   cost is entirely in what it omits. An omitted junction is worse than no list
#   at all, since the aligner snaps the read onto a listed junction sharing one
#   of its two sites and returns a confident alignment on the wrong exon.
set -euo pipefail

ALIGNER=$1
genome_ind=$2
source_file=$3
Trans_run_dir=$4
pre_out=$5
statistic=$6
HE=$7
SPLICE=${8:-}
# The published TransRun.sh hardcodes five threads. On the simulation that is
# the whole machine; on a real library it is the bottleneck, so it is settable.
THREADS="${THREADS:-5}"

Trans_dir="$Trans_run_dir/TransFiles"
bam_files="$Trans_run_dir/bamFiles"
reject_files="$Trans_run_dir/Reject"
filter_statistic="$reject_files/statistic"
mkdir -p "$Trans_run_dir" "$Trans_dir" "$bam_files" "$reject_files"

align_origReads="origReads.$pre_out"
align_FilterReads="FilterReads.$pre_out"

for pair in ag tc ac tg gc cg at ta; do
  from=${pair:0:1}; to=${pair:1:1}
  "$HE/fastq_transform.pl" "$from" "$to" "$source_file" \
    > "$Trans_dir/${from}2${to}.$pre_out.fastq"
done

echo -e "\nTrans Run ($ALIGNER)" >> "$statistic"

# combo name, genome transform, read transform
COMBOS="
a2g++ a2g a2g
a2g+- a2g t2c
a2g-+ t2c a2g
a2g-- t2c t2c
a2c++ a2c a2c
a2c+- a2c t2g
a2c-+ t2g a2c
a2c-- t2g t2g
g2c++ g2c g2c
g2c+- g2c c2g
a2t++ a2t a2t
a2t+- a2t t2a
"

echo "$COMBOS" | while read -r combo gtr rtr; do
  [ -n "${combo:-}" ] || continue
  echo "statistic of $combo" >> "$statistic"
  reads="$Trans_dir/$rtr.$pre_out.fastq"

  case "$ALIGNER" in
    mem)
      bwa mem -Y -t "$THREADS" "$genome_ind.$gtr" "$reads" \
        2>> "$Trans_run_dir/align.log" > "$bam_files/temp.sam"
      xa2multi.pl "$bam_files/temp.sam" > "$bam_files/temp2.sam"
      ;;
    hisat2)
      # No splice-site file means de novo splicing, which the probes measure and
      # advise against; it is allowed here so that "no list" can be an arm.
      hisat2 -p "$THREADS" -x "$genome_ind.$gtr.ht2" -U "$reads" \
        ${SPLICE:+--known-splicesite-infile "$SPLICE"} --no-unal \
        -S "$bam_files/temp2.sam" 2>> "$Trans_run_dir/align.log"
      ;;
    *) echo "unknown aligner: $ALIGNER" >&2; exit 1 ;;
  esac

  samtools view -F 4 -S -h -o "$bam_files/temp3.sam" "$bam_files/temp2.sam"

  "$HE/Get_orig_read.pl" "$bam_files/temp3.sam" "$source_file" \
    > "$bam_files/$combo.$align_origReads.sam"
  samtools view -bS "$bam_files/$combo.$align_origReads.sam" \
    -o "$bam_files/$combo.$align_origReads.bam"

  "$HE/filter_sam.pl" "$bam_files/$combo.$align_origReads.sam" "$reject_files" \
    "$filter_statistic" > "$bam_files/$combo.$align_FilterReads.sam"
  samtools view -bS "$bam_files/$combo.$align_FilterReads.sam" \
    -o "$bam_files/$combo.$align_FilterReads.bam"
  cat "$filter_statistic.$combo" >> "$statistic"
  echo "" >> "$statistic"

  # Every stage here is written as SAM and again as BAM, and only the BAMs are
  # read downstream. On the simulation the duplicate is invisible; on a real
  # library it is tens of gigabytes per arm, so the SAM goes as soon as its BAM
  # exists. The one exception is origReads.sam, which filter_sam.pl reads -- it
  # is dropped after that, not before.
  rm -f "$bam_files/$combo.$align_origReads.sam" \
        "$bam_files/$combo.$align_FilterReads.sam"
done

rm -f "$bam_files"/temp*.sam
