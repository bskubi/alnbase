#!/usr/bin/env bash
# The cross-species re-test, as one resumable job.
#
#   run_xspec.sh WORKDIR SPECIES [ONLY]
#
#   WORKDIR   scratch; everything is created under it and nothing outside it
#             is written. Needs ~170 GB for the opossum pair (see "Disk").
#   SPECIES   a species named in assemblies.tsv and reads.tsv, e.g. opossum
#   ONLY      optional: one assembly label, so a scheduler can give each
#             assembly its own job, or `score` to only re-score what exists
#
# What it runs is a 2x2: reference quality crossed with aligner sensitivity.
#
#                     | bwa aln -n 2 -o 0 -N  | HISAT2
#     ----------------+-----------------------+--------
#     old assembly    |  A  (= Porath 2014)   |  B
#     new assembly    |  C                    |  D
#
# A is not merely similar to the 2014 result, it is that result: the same
# assembly, the same pipeline at the same commit, the same reads. A->C isolates
# sixteen years of assembly improvement with the published pipeline untouched.
# A->B isolates the aligner. C->D asks the question that decides whether the
# aligner work still matters: if a modern assembly absorbs what the aligner swap
# buys, the two are substitutes; if the gains are additive, both confounders are
# real and both survive into the 2017 cross-species revision.
#
# THE CALLER IS HELD FIXED ACROSS THE ALIGNER AXIS. run_real.sh pairs its `aln`
# arm with the patched Perl caller and its HISAT2 arms with the alnbase caller,
# which is the right pairing when the caller is itself under test. Here it would
# confound the only contrast that matters, so both cells use the CIGAR-aware
# Perl caller. That also drops the alnbase binary from this script's
# dependencies, which is one less thing to build on a cluster.
#
# Resumability: every stage is guarded by a marker file, so a job killed by a
# wall-clock limit can be resubmitted unchanged and will pick up where it
# stopped. Markers are written only after the stage's output is complete.
#
# Disk. All fourteen indexes per assembly are built up front and kept, because
# on a cluster wall time binds and disk does not: that lets the six transforms
# index concurrently instead of serially, turning ~23 h of indexing into ~4 h.
# For a 3.6 Gb genome budget ~6.3 GB per bwa index and ~5.0 GB per hisat2 index,
# so ~80 GB per assembly and ~160 GB for a pair, plus reads and BAMs. On a
# disk-bound machine this script is the wrong shape -- serialise the transforms
# instead, at ~25 GB peak and roughly six times the wall clock.
#
# Memory. hisat2-build on a mammalian genome wants ~8 GB. INDEX_JOBS of them run
# at once, so the default of 3 asks for ~24 GB. Raise it only with the memory to
# match; it is the setting most likely to get a job killed.
#
# Environment: PATH must carry bwa, hisat2, hisat2-build,
# hisat2_extract_splice_sites.py, samtools, perl, python3 and curl.
set -euo pipefail

WORK="${1:?usage: run_xspec.sh WORKDIR SPECIES [ONLY]}"
SPECIES="${2:?usage: run_xspec.sh WORKDIR SPECIES [ONLY]}"
ONLY="${3:-}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEMO="$(cd "$HERE/.." && pwd)"

THREADS="${THREADS:-16}"
INDEX_JOBS="${INDEX_JOBS:-3}"
ARMS_AT_ONCE="${ARMS_AT_ONCE:-2}"
PY="${PY:-python3}"

# Where stage 2 does its churn. Stage 2 writes and rewrites the whole unmapped
# pool twelve times per arm as SAM and again as BAM, which is the heaviest I/O
# in the run and none of it is wanted afterwards -- only the beds are. On a
# cluster with node-local disk this belongs there rather than on the shared
# filesystem, so it is settable. It defaults to WORKDIR, which keeps the
# single-machine behaviour unchanged.
#
# The trade is resumability: node-local scratch is wiped between jobs, so an
# interrupted arm restarts its alignments from the beginning. That is a few
# hours against the twenty of index building, which stays on shared storage.
XSPEC_TMP="${XSPEC_TMP:-}"

mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
HE="$WORK/Hyper-editing"
LOG () { echo "[$(date +'%F %T')] $*"; }

# Shared preparation -- the git clone and the read download -- is done once per
# work directory, but several assemblies may be running as concurrent jobs
# against it. Without a lock two of them clone into the same directory and
# download to the same path at the same time.
#
# `mkdir` is the primitive rather than `flock` because this runs on shared
# storage: ARC's gscratch does not support POSIX locking, which rules flock out,
# while an atomic directory creation works on any filesystem worth using.
#
# The stale-lock timeout matters because a job killed on its wall clock leaves
# the directory behind. Ten minutes is comfortably longer than a clone plus a
# few-GB download and short enough not to waste a queue slot.
HELD=""
lock () {
  local d="$WORK/.lock.$1" waited=0
  while ! mkdir "$d" 2>/dev/null; do
    if [ -d "$d" ] && [ "$(( $(date +%s) - $(stat -c %Y "$d" 2>/dev/null || date +%s) ))" -gt 600 ]; then
      LOG "removing stale lock $1"; rm -rf "$d"; continue
    fi
    [ "$waited" -eq 0 ] && LOG "waiting for another job to finish $1"
    sleep 10; waited=$((waited + 10))
    [ "$waited" -gt 3600 ] && { echo "gave up waiting for $d" >&2; exit 1; }
  done
  HELD="$HELD $1"
  trap release EXIT
}
# Only ever removes locks this process actually holds. Removing one
# unconditionally would delete a concurrent job's lock in the window between it
# writing its completion marker and releasing.
release () {
  local n
  for n in $HELD; do rm -rf "$WORK/.lock.$n"; done
  HELD=""
}
unlock () {
  case " $HELD " in *" $1 "*) rm -rf "$WORK/.lock.$1"
    HELD="$(echo "$HELD" | tr ' ' '\n' | grep -vx "$1" | tr '\n' ' ')" ;;
  esac
}

# --- the pipeline itself, at its published commit ----------------------------
if [ ! -f "$WORK/hyper-editing-commit.txt" ]; then
  lock bootstrap
fi
if [ ! -f "$WORK/hyper-editing-commit.txt" ]; then
  LOG "bootstrap: fetching and patching the Hyper-editing tree"
  "$DEMO/setup.sh" "$WORK"
fi
unlock bootstrap
export PATH="$HE:$PATH"

# The two patched Perl stages. analyse_mm.pl reads the reference with a flat
# window from POS, which is correct only for a gapless end-to-end alignment;
# the patch makes it walk the CIGAR so a spliced or clipped record lands on the
# right reference bases. Without it the HISAT2 cells would not be wrong-ish,
# they would be nonsense -- every junction-spanning read compared against
# intron sequence.
if [ ! -x "$WORK/analyse_mm_cigar.pl" ]; then
  "$PY" "$DEMO/bin/patch_analyse_mm.py" "$HE/analyse_mm.pl" "$WORK/analyse_mm_cigar.pl"
  "$PY" "$DEMO/bin/patch_detect_ue.py"  "$HE/detect_ue.pl"  "$WORK/detect_ue_ported.pl"
  chmod +x "$WORK/analyse_mm_cigar.pl" "$WORK/detect_ue_ported.pl"
fi
UE_ARGS="0.05 0.6 30 0.6 0.1 0.8 0.2"
TRANSFORMS="a2g t2c a2c t2g g2c a2t"

# --- reads -------------------------------------------------------------------
# The FTP path comes from ENA rather than from the manifest: the directory
# convention depends on the accession's length and trailing digits, and a
# hand-written path that is wrong fails at download time, hours in.
RUN=$(awk -F'\t' -v s="$SPECIES" '$1==s && $1!~/^#/ {print $2; exit}' "$HERE/reads.tsv")
[ -n "$RUN" ] || { echo "no read run for '$SPECIES' in reads.tsv" >&2; exit 1; }
READS="$WORK/reads/$RUN.fastq.gz"
if [ ! -f "$WORK/reads/.done.$RUN" ]; then
  mkdir -p "$WORK/reads"
  lock "reads.$RUN"
fi
if [ ! -f "$WORK/reads/.done.$RUN" ]; then
  LOG "reads: resolving $RUN at ENA"
  # The column is located by header name, not by position: ENA prepends
  # run_accession to every filereport whether or not it was asked for, so
  # taking field 1 yields the accession back and the download fails later,
  # after the indexes have already been built.
  URL=$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN&result=read_run&fields=fastq_ftp" \
        | awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="fastq_ftp") c=i; next}
                      c {print $c; exit}' | cut -d';' -f1)
  [ -n "$URL" ] || { echo "ENA returned no fastq path for $RUN" >&2; exit 1; }
  LOG "reads: downloading $URL"
  # Downloaded beside the target and moved into place, so a reader never sees a
  # partial file under the final name even if this job dies mid-transfer.
  curl -fsSL -C - -o "$READS.part" "ftp://$URL"
  mv -f "$READS.part" "$READS"
  touch "$WORK/reads/.done.$RUN"
fi
unlock "reads.$RUN"

# --- one assembly, end to end ------------------------------------------------
# Called once per era. Everything it makes lives under $WORK/$LABEL.
prep_assembly () {
  local label=$1 fasta_url=$2 gtf_url=$3
  local d="$WORK/$label" ref="$WORK/$label/$label.fa" pre="$WORK/$label/$label"
  mkdir -p "$d"

  if [ ! -f "$d/.done.fasta" ]; then
    LOG "$label: fetching reference"
    curl -fsSL -C - -o "$d/src.gz" "$fasta_url"
    # goldenPath ships some older assemblies (mm9, dm3) as a tar of
    # per-chromosome files rather than one FASTA; both shapes end up as
    # $label.fa. `find` rather than a glob because the tars are not consistent
    # about whether the members sit at the top level or under a directory.
    if [[ "$fasta_url" == *.tar.gz ]]; then
      mkdir -p "$d/chroms" && tar -xzf "$d/src.gz" -C "$d/chroms"
      find "$d/chroms" -name '*.fa' -print0 | sort -z | xargs -0 cat > "$ref"
      rm -rf "$d/chroms"
    else
      gunzip -c "$d/src.gz" > "$ref"
    fi
    rm -f "$d/src.gz"
    touch "$d/.done.fasta"
  fi
  # Outside the guard: arc/fetch.sh may have staged the FASTA on a login node
  # that has no samtools, in which case the index is still missing here.
  [ -f "$ref.fai" ] || samtools faidx "$ref"

  if [ ! -f "$d/.done.gtf" ]; then
    LOG "$label: fetching annotation"
    curl -fsSL -C - -o "$d/genes.gtf.gz" "$gtf_url"
    touch "$d/.done.gtf"
  fi

  # A GTF whose sequence names do not intersect the FASTA's yields an empty
  # junction list and a HISAT2 arm quietly running without one. GenArk names
  # sequences by RefSeq accession and goldenPath uses chr1..chrX, so this is a
  # live risk whenever the two URLs come from different halves of UCSC.
  #
  # Checked outside the fetch guard, and on every run, because arc/fetch.sh
  # downloads these files on the login node and sets the same markers -- inside
  # the guard this would never execute in the job that actually uses them.
  if [ ! -f "$d/.done.names" ]; then
    local shared
    shared=$(comm -12 \
      <(cut -f1 "$ref.fai" | sort -u) \
      <(zcat "$d/genes.gtf.gz" | grep -v '^#' | cut -f1 | sort -u) | wc -l)
    if [ "$shared" -eq 0 ]; then
      echo "$label: GTF and FASTA share no sequence names -- wrong pairing" >&2
      exit 1
    fi
    LOG "$label: GTF and FASTA share $shared sequence names"
    touch "$d/.done.names"
  fi

  # The untransformed indexes: bwa for stage 1, hisat2 for junction discovery.
  if [ ! -f "$d/.done.index.base" ]; then
    LOG "$label: indexing the untransformed reference"
    bwa index -p "$pre" "$ref" 2> "$d/bwa.base.log"
    hisat2-build -q -p "$THREADS" "$ref" "$pre.ht2" 2> "$d/hisat2.base.log"
    touch "$d/.done.index.base"
  fi

  # The six three-letter references, and both aligners' indexes of each. These
  # are independent, so they are built concurrently; INDEX_JOBS bounds the
  # memory, not the cores.
  if [ ! -f "$d/.done.index.trans" ]; then
    LOG "$label: six transformed references, $INDEX_JOBS at a time"
    for tt in $TRANSFORMS; do
      while [ "$(jobs -rp | wc -l)" -ge "$INDEX_JOBS" ]; do wait -n || true; done
      (
        [ -f "$d/$label.$tt.fa" ] || "$HE/$tt.pl" "$ref" > "$d/$label.$tt.fa"
        [ -f "$pre.$tt.bwt" ]      || bwa index -p "$pre.$tt" "$d/$label.$tt.fa" 2> "$d/bwa.$tt.log"
        [ -f "$pre.$tt.ht2.1.ht2" ] || hisat2-build -q -p "$THREADS" "$d/$label.$tt.fa" "$pre.$tt.ht2" 2> "$d/hisat2.$tt.log"
        # The transformed FASTA is only an input to the two indexes. The
        # untransformed one stays: the edit caller reads reference bases from it.
        rm -f "$d/$label.$tt.fa"
        LOG "$label: $tt indexed"
      ) &
    done
    wait
    touch "$d/.done.index.trans"
  fi
}

# --- stage 1 and the junction lists, per assembly ----------------------------
# Stage 1 is the published pre_unmapped.sh reproduced rather than called, for
# the two reasons run_real.sh gives: `-t 5` is hardcoded there, and `samse -n 50`
# spends most of its time building an XA tag that nothing downstream reads.
stage1 () {
  local label=$1
  local d="$WORK/$label" pre="$WORK/$label/$label" out="$WORK/run_$label"
  mkdir -p "$out/stage1" "$out/lists"

  if [ ! -f "$out/stage1/.done" ]; then
    LOG "$label: stage 1, what ordinary alignment cannot place"
    bwa aln -t "$THREADS" "$pre" "$READS" -f "$out/stage1/aln.sai" 2> "$out/stage1/aln.log"
    bwa samse "$pre" "$out/stage1/aln.sai" "$READS" 2>> "$out/stage1/aln.log" \
      | samtools view -b -o "$out/stage1/aln.bam" -
    samtools fastq -f 4 "$out/stage1/aln.bam" > "$out/stage1/aln.um.fastq" 2>/dev/null
    bwa mem -M -t "$THREADS" -k 50 "$pre" "$out/stage1/aln.um.fastq" 2> "$out/stage1/mem.log" \
      | samtools view -b -o "$out/stage1/mem.bam" -
    samtools fastq -f 4 "$out/stage1/mem.bam" > "$out/stage1/um.fastq" 2>/dev/null
    rm -f "$out/stage1/aln.sai" "$out/stage1/aln.um.fastq"
    touch "$out/stage1/.done"
  fi

  # The denominator, written down explicitly. A better assembly places more
  # reads in stage 1, which raises "mapped reads" and can lower sites-per-read
  # even as the absolute number of sites rises. Reporting the ratio alone would
  # let the reference effect and the aligner effect cancel invisibly, so each
  # term is recorded separately and score_xspec.py prints all of them.
  if [ ! -s "$out/stage1/counts.tsv" ]; then
    local total mapped unplaced
    total=$(( $(zcat -f "$READS" | wc -l) / 4 ))
    unplaced=$(( $(wc -l < "$out/stage1/um.fastq") / 4 ))
    mapped=$(( total - unplaced ))
    printf "total_reads\t%s\nmapped_reads\t%s\nunmapped_pool\t%s\n" \
      "$total" "$mapped" "$unplaced" > "$out/stage1/counts.tsv"
  fi
  LOG "$label: $(tr '\n' ' ' < "$out/stage1/counts.tsv")"

  # Junction lists. The union of the annotation and the sample's own first pass,
  # because a list's cost is entirely in what it omits: an omitted junction is
  # worse than no list, since HISAT2 will snap the read onto a listed junction
  # sharing one of its two sites and return a confident alignment on the wrong
  # exon. The first pass is read-derived, so it partly covers for a thin
  # annotation -- which matters here, as the old and new assemblies cannot have
  # equal-quality gene sets and that inequality is part of the reference axis.
  if [ ! -s "$out/lists/union.txt" ]; then
    LOG "$label: junction lists"
    hisat2_extract_splice_sites.py <(zcat "$d/genes.gtf.gz") \
      | sort -k1,1 -k2,2n -k3,3n -u > "$out/lists/annotation.txt"
    hisat2 -p "$THREADS" -x "$pre.ht2" -U "$READS" \
      --novel-splicesite-outfile "$out/lists/firstpass.raw" \
      -S /dev/null 2> "$out/lists/firstpass.log"
    sort -k1,1 -k2,2n -k3,3n -u "$out/lists/firstpass.raw" > "$out/lists/firstpass.txt"
    sort -k1,1 -k2,2n -k3,3n -u "$out/lists/annotation.txt" "$out/lists/firstpass.txt" \
      > "$out/lists/union.txt"
  fi
  printf "junctions: annotation %s, first pass %s, union %s\n" \
    "$(wc -l < "$out/lists/annotation.txt")" \
    "$(wc -l < "$out/lists/firstpass.txt")" \
    "$(wc -l < "$out/lists/union.txt")" | tee "$out/lists/counts.txt"
}

# --- one cell of the 2x2 -----------------------------------------------------
arm () {
  local label=$1 arm=$2 aligner=$3 splice=${4:-}
  local pre="$WORK/$label/$label" out="$WORK/run_$label"
  local trans="${XSPEC_TMP:-$out/stage2}/$label.$arm" res="$out/arm_$arm"
  [ -f "$res/.done" ] && { LOG "$label/$arm: already done"; return; }
  mkdir -p "$(dirname "$trans")"

  if [ ! -f "$trans/.done" ]; then
    rm -rf "$trans"
    LOG "$label/$arm: twelve transform alignments"
    if [ "$aligner" = aln ]; then
      "$HE/TransRun.sh" "$pre" "$out/stage1/um.fastq" "$trans" real \
        "$trans.stat" bwa samtools > "$out/stage2/$arm.log" 2>&1
    else
      "$DEMO/bin/TransRunAligner.sh" "$aligner" "$pre" "$out/stage1/um.fastq" \
        "$trans" real "$trans.stat" "$HE" $splice > "$out/stage2/$arm.log" 2>&1
    fi
    touch "$trans/.done"
  fi

  rm -rf "$res"; mkdir -p "$res/beds"; : > "$res/hits.txt"
  LOG "$label/$arm: calling edits"
  for bam in "$trans"/bamFiles/*FilterReads*.bam; do
    "$WORK/analyse_mm_cigar.pl" "$WORK/$label/$label.fa" "$bam" samtools
  done >> "$res/hits.txt" 2>> "$res/call.log"
  "$HE/sort_R_read.pl" "$res/hits.txt" "$res/sort.stat" > "$res/real.analyseMM"
  "$WORK/detect_ue_ported.pl" "$res/real.analyseMM" "$res/beds/real" \
    "$res/detect.stat" $UE_ARGS samtools 0 > "$res/detect.log" 2>&1
  touch "$res/.done"
  LOG "$label/$arm: finished"

  # The transform FASTQs are the bulk of an arm's footprint and nothing reads
  # them after its alignments exist. Dropped only once the arm is complete, so
  # an interrupted arm still restarts cleanly. When stage 2 is on scratch that
  # will be torn down anyway, the whole tree goes rather than just the FASTQs:
  # leaving it would mean the next arm on the same node inherits a full disk.
  if [ -n "$XSPEC_TMP" ]; then rm -rf "$trans"; else rm -rf "$trans/TransFiles"; fi
}

# --- run it ------------------------------------------------------------------
LABELS=""
while IFS=$'\t' read -r sp label era fasta_url gtf_url; do
  case "$sp" in \#*|"") continue ;; esac
  [ "$sp" = "$SPECIES" ] || continue
  # LABELS is built from every row for this species, not only the ones this
  # invocation runs, so that the scoring step sees the whole 2x2 even when each
  # assembly was given its own job.
  LABELS="$LABELS $label"
  [ -z "$ONLY" ] || [ "$ONLY" = "$label" ] || continue
  LOG "=== $label ($era assembly) ==="
  prep_assembly "$label" "$fasta_url" "$gtf_url"
  stage1 "$label"
  mkdir -p "$WORK/run_$label/stage2"
  for spec in "aln aln" "hisat2-union hisat2 $WORK/run_$label/lists/union.txt"; do
    while [ "$(jobs -rp | wc -l)" -ge "$ARMS_AT_ONCE" ]; do wait -n || true; done
    # `< /dev/null` because this loop's stdin is the manifest: a background job
    # that read from it would eat assembly rows and silently skip an era.
    # shellcheck disable=SC2086
    arm "$label" $spec < /dev/null &
  done
  wait
done < "$HERE/assemblies.tsv"

[ -n "$LABELS" ] || { echo "no assemblies for '$SPECIES' in assemblies.tsv" >&2; exit 1; }

# A job that was handed a single assembly has only half the 2x2 and would print
# a table that invites the wrong reading, so it stops here. `score` re-scores
# whatever is finished, which is how the chained submission ends.
if [ -n "$ONLY" ] && [ "$ONLY" != score ]; then
  LOG "=== $ONLY done; score the pair with: run_xspec.sh $WORK $SPECIES score ==="
  exit 0
fi

LOG "=== scoring ==="
# shellcheck disable=SC2086
"$PY" "$HERE/score_xspec.py" "$WORK" --species "$SPECIES" \
  --out "$WORK/xspec_$SPECIES.tsv" $LABELS
