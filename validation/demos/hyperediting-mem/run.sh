#!/usr/bin/env bash
# Run the hyper-editing pipeline three ways over one simulation and score all
# three against ground truth.
#
#   run.sh WORKDIR SIMDIR
#
#   WORKDIR   holds the patched pipeline tree from setup.sh
#   SIMDIR    holds genome.fa, reads.fastq and the truth tables from simulate.py
#
# Environment: ALNBASE (path to the binary), PATH must carry bwa, samtools,
# perl, xa2multi.pl and a python with pysam and pyarrow.
#
# The three arms:
#
#   aln           the published pipeline, unmodified.
#
#   mem-flat      bwa mem, with the original edit caller still cutting a flat
#                 window out of the genome. This is the change on its own, and
#                 it is in the comparison because it is what a swap looks like
#                 when nobody notices what the caller assumed. Its output is
#                 not obviously broken; that is the finding.
#
#   mem-alnbase   bwa mem, with the caller replaced by one that gets the
#                 read-to-reference pairing from alnbase's CIGAR walk.
#
# Arms two and three share one alignment stage, so every difference between
# them comes from the caller and nothing else.
set -euo pipefail

WORK="${1:?usage: run.sh WORKDIR SIMDIR}"
SIM="${2:?usage: run.sh WORKDIR SIMDIR}"
HE="$WORK/Hyper-editing"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALNBASE="${ALNBASE:-alnbase}"
STRAND_RULE="${STRAND_RULE:-$HERE/../../../queries/strand/unconverted.toml}"
PY="${PY:-python3}"
UE_ARGS="0.05 0.6 30 0.6 0.1 0.8 0.2"

SIM="$(cd "$SIM" && pwd)"
cd "$SIM"
export PATH="$HE:$PATH"

# --- indexes -----------------------------------------------------------------
if [ ! -f chrS.bwt ]; then
  cp -f genome.fa chrS.fa
  bwa index -p chrS chrS.fa
  for tt in a2c a2g a2t g2c t2c t2g; do
    "$HE/$tt.pl" chrS.fa > "chrS.$tt.fa"
    bwa index -p "chrS.$tt" "chrS.$tt.fa"
  done
fi
[ -f chrS.aref ] || "$ALNBASE" index chrS.fa chrS.aref
for tt in a2c a2g a2t g2c t2c t2g; do
  [ -f "chrS.$tt.ht2.1.ht2" ] || hisat2-build -q -p 4 "chrS.$tt.fa" "chrS.$tt.ht2"
done

# --- the three callers and the one downstream patch --------------------------
# Every arm's caller now puts the alignment's CIGAR on its hit lines, and every
# arm runs the same patched detect_ue.pl, which uses it to convert an edit's
# offset along the read into a genomic coordinate. The conversion is the
# identity on a CIGAR that is a single run of M, which is every read the
# published pipeline can produce, so the aln arm is unchanged by it.
#
# This is the part of the port that does not fit in one box. The flat-window
# assumption is made twice: once by the caller when it cuts a genome window,
# and again by detect_ue.pl at line 538 when it adds a read offset to a start
# position. Fixing only the first would leave every recovered clipped or
# indel-bearing read placed at the wrong base, which is worse than not
# recovering it at all.
#
# The mem-flat arm needs one line more. analyse_mm.pl counts a read's alignment
# positions from X0 and X1, which only `bwa aln` writes; under mem the count is
# zero, and sort_R_read.pl -- which reads exactly that many hit lines after each
# record -- then falls out of step with the file.
"$PY" "$HERE/bin/patch_analyse_mm.py" "$HE/analyse_mm.pl" "$WORK/analyse_mm_cigar.pl"
"$PY" "$HERE/bin/patch_analyse_mm.py" "$HE/analyse_mm.pl" "$WORK/analyse_mm_memhits.pl" --mem-hits
"$PY" "$HERE/bin/patch_analyse_mm.py" "$HE/analyse_mm.pl" "$WORK/analyse_mm_memonly.pl" --mem-hits --no-cigar
"$PY" "$HERE/bin/patch_detect_ue.py" "$HE/detect_ue.pl" "$WORK/detect_ue_ported.pl"
"$PY" "$HERE/bin/patch_detect_ue.py" "$HE/detect_ue.pl" "$WORK/detect_ue_div0.pl" --div0-only
chmod +x "$WORK"/analyse_mm_*.pl "$WORK"/detect_ue_*.pl

# --- shared first stage ------------------------------------------------------
# Reads that ordinary alignment could not place. Identical for all three arms,
# so it runs once.
mkdir -p stage/unMap
if [ ! -s stage/unMap/sim.mem.um.fastq ]; then
  "$HE/pre_unmapped.sh" reads.fastq "$SIM/chrS" 33 0 "$SIM/stage/unMap" \
    "$SIM/stage/sim" "$SIM/stage/general" bwa bwa bam2fastx samtools \
    > stage/pre_unmapped.log 2>&1
fi
UM="$SIM/stage/unMap/sim.mem.um.fastq"

# --- alignment stages --------------------------------------------------------
if [ ! -d stage/TransRun_aln/bamFiles ]; then
  "$HE/TransRun.sh" "$SIM/chrS" "$UM" "$SIM/stage/TransRun_aln" sim \
    "$SIM/stage/stat_aln" bwa samtools > stage/transrun_aln.log 2>&1
fi
if [ ! -d stage/TransRun_mem/bamFiles ]; then
  "$HERE/bin/TransRunAligner.sh" mem "$SIM/chrS" "$UM" "$SIM/stage/TransRun_mem" sim \
    "$SIM/stage/stat_mem" "$HE" > stage/transrun_mem.log 2>&1
fi
if [ ! -d stage/TransRun_hisat2/bamFiles ]; then
  "$HERE/bin/TransRunAligner.sh" hisat2 "$SIM/chrS" "$UM" "$SIM/stage/TransRun_hisat2" sim \
    "$SIM/stage/stat_hisat2" "$HE" "$SIM/splicesites.txt" > stage/transrun_hisat2.log 2>&1
fi

# --- edit calling, then the untouched downstream stages ----------------------
call_and_detect () {
  local arm=$1 transrun=$2
  local out="$SIM/arm_$arm"
  rm -rf "$out"; mkdir -p "$out/beds"
  : > "$out/hits.txt"
  for bam in "$transrun"/bamFiles/*FilterReads*.bam; do
    case "${arm##*-}" in
      published) "$WORK/analyse_mm_cigar.pl" chrS.fa "$bam" samtools ;;
      naive)     "$WORK/analyse_mm_memonly.pl" chrS.fa "$bam" samtools ;;
      flat)      "$WORK/analyse_mm_memhits.pl" chrS.fa "$bam" samtools ;;
      alnbase)   "$PY" "$HERE/bin/alnbase_mm.py" chrS.fa "$bam" \
                   --aref chrS.aref --strand-rule "$STRAND_RULE" \
                   --alnbase "$ALNBASE" --work "$out/alnbase_work" ;;
      *) echo "unknown caller for arm $arm" >&2; exit 1 ;;
    esac >> "$out/hits.txt" 2>> "$out/call.log"
  done
  "$HE/sort_R_read.pl" "$out/hits.txt" "$out/sort.stat" > "$out/sim.analyseMM"
  local detect="$WORK/detect_ue_ported.pl"
  [ "${arm##*-}" = "naive" ] && detect="$WORK/detect_ue_div0.pl"
  "$detect" "$out/sim.analyseMM" "$out/beds/sim" "$out/detect.stat" \
    $UE_ARGS samtools 0 > "$out/detect.log" 2>&1
  printf "%-12s %6d hit lines  %6d called sites\n" "$arm" \
    "$(grep -vc $'\thits:' "$out/hits.txt" || true)" \
    "$(wc -l < "$out/beds/sim.ES.bed_files/A2G.bed")"
}

echo "running arms"
call_and_detect aln-published    "$SIM/stage/TransRun_aln"
call_and_detect mem-naive        "$SIM/stage/TransRun_mem"
call_and_detect mem-flat         "$SIM/stage/TransRun_mem"
call_and_detect mem-alnbase      "$SIM/stage/TransRun_mem"
call_and_detect hisat2-flat      "$SIM/stage/TransRun_hisat2"
call_and_detect hisat2-alnbase   "$SIM/stage/TransRun_hisat2"

# --- score -------------------------------------------------------------------
ARMS=(aln-published mem-naive mem-flat mem-alnbase hisat2-flat hisat2-alnbase)
SPECS=()
for a in "${ARMS[@]}"; do
  SPECS+=("$a=$SIM/arm_$a/beds/sim.ES.bed_files/A2G.bed")
done
"$PY" "$HERE/score.py" "$SIM" "${SPECS[@]}" --out "$SIM/score.tsv"
