#!/usr/bin/env bash
# Build everything the real-data run needs: the reference, its six three-letter
# transforms, and an index of each for every aligner in the comparison.
#
#   prep_real.sh WORKDIR [--no-hisat2]
#
# This is the expensive step -- thirteen whole-genome indexes -- so it is
# resumable: every stage checks for its own output and skips if it is there.
# Transformed FASTAs are deleted as soon as both indexes are built, because
# thirteen indexes plus six 3 Gb FASTAs does not fit on a laptop.
#
# Reference: the GRCh38 analysis set, which is the primary assembly with no alt
# contigs, under UCSC contig names. Alts are left out because bwa's backtrack
# mode has no alt-awareness and an alt-carrying reference would make every read
# over an alt region look like a multi-mapper -- which the pipeline rejects at
# 51 hits. UCSC names because the published TransRun.sh filters unplaced
# contigs by grepping for "chrUn", "_random" and "chrUextra" in the SAM text.
set -euo pipefail

WORK="${1:?usage: prep_real.sh WORKDIR [--no-hisat2]}"
DO_HISAT2=1
[ "${2:-}" = "--no-hisat2" ] && DO_HISAT2=0
HE="$(dirname "$WORK")/Hyper-editing"
[ -d "$HE" ] || HE="$WORK/../Hyper-editing"
ALNBASE="${ALNBASE:-alnbase}"
THREADS="${THREADS:-8}"
REF="$WORK/hg38.analysisSet.fa"

need_gb () {   # bail before starting a stage that cannot finish
  local want=$1 have
  have=$(df -BG --output=avail "$WORK" | tail -1 | tr -dc '0-9')
  if [ "$have" -lt "$want" ]; then
    echo "prep_real: ${have}G free, stage needs ${want}G -- stopping here" >&2
    exit 3
  fi
}

cd "$WORK"

if [ ! -f "$REF" ]; then
  need_gb 10
  echo "[$(date +%T)] unpacking reference"
  gunzip -c hg38.analysisSet.fa.gz > "$REF"
fi
[ -f "$REF.fai" ] || samtools faidx "$REF"

if [ ! -f hg38.bwt ]; then
  need_gb 12
  echo "[$(date +%T)] bwa index, untransformed (stage 1 needs this one)"
  bwa index -p hg38 "$REF"
fi

if [ ! -f hg38.aref ]; then
  need_gb 8
  echo "[$(date +%T)] alnbase index"
  "$ALNBASE" index "$REF" hg38.aref
fi

for tt in a2g t2c a2c t2g g2c a2t; do
  bwa_done=0; hisat_done=0
  [ -f "hg38.$tt.bwt" ] && bwa_done=1
  [ -f "hg38.$tt.ht2.1.ht2" ] && hisat_done=1
  [ "$DO_HISAT2" = 0 ] && hisat_done=1
  if [ "$bwa_done" = 1 ] && [ "$hisat_done" = 1 ]; then
    echo "[$(date +%T)] $tt already built"
    continue
  fi

  if [ ! -f "hg38.$tt.fa" ]; then
    need_gb 8
    echo "[$(date +%T)] transforming $tt"
    "$HE/$tt.pl" "$REF" > "hg38.$tt.fa"
  fi
  if [ "$bwa_done" = 0 ]; then
    need_gb 12
    echo "[$(date +%T)] bwa index $tt"
    bwa index -p "hg38.$tt" "hg38.$tt.fa"
  fi
  if [ "$hisat_done" = 0 ]; then
    need_gb 12
    echo "[$(date +%T)] hisat2-build $tt"
    hisat2-build -q -p "$THREADS" "hg38.$tt.fa" "hg38.$tt.ht2"
  fi
  rm -f "hg38.$tt.fa"
  echo "[$(date +%T)] $tt done; $(df -BG --output=avail "$WORK" | tail -1 | tr -d ' ') free"
done

echo "[$(date +%T)] all indexes built"
du -sh "$WORK"
