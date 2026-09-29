#!/usr/bin/env bash
# One ARC job: everything for a single assembly, or the final scoring step.
#
# Not submitted by hand -- `submit.sh` sets the sbatch options and the
# dependency chain. It is a plain script, so it also runs under `srun` or in an
# salloc session for debugging.
#
# Arguments arrive through the environment, because sbatch's own argument
# passing is awkward to combine with --wrap-free submission:
#
#   XSPEC_WORK     the gscratch directory prepared by fetch.sh
#   XSPEC_SPECIES  species name
#   XSPEC_ONLY     an assembly label, or `score`
#   XSPEC_ENV      conda/mamba environment prefix to put on PATH; when unset,
#                  read from $XSPEC_WORK/.xspec_env if that file exists
#   XSPEC_MODULES  optional, space-separated module names to load; when unset,
#                  read from $XSPEC_WORK/.xspec_modules if that file exists
set -euo pipefail

WORK="${XSPEC_WORK:?XSPEC_WORK is not set}"
SPECIES="${XSPEC_SPECIES:?XSPEC_SPECIES is not set}"
ONLY="${XSPEC_ONLY:?XSPEC_ONLY is not set}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
XSPEC="$(cd "$HERE/.." && pwd)"
LOG () { echo "[$(date +'%F %T')] $*"; }

LOG "job ${SLURM_JOB_ID:-local} on $(hostname): $SPECIES/$ONLY"

# --- software ----------------------------------------------------------------
# A mamba environment, not modules. ARC has no bwa, hisat2 or samtools module --
# `module avail` offers bowtie2, bamtools, bedtools2 and kallisto, but not these
# -- so setup_env.sh builds the environment once on gscratch and this puts its
# bin on PATH. Prepending rather than appending so its pinned versions win over
# anything the login profile happened to set up.
ENV_PREFIX="${XSPEC_ENV:-}"
[ -z "$ENV_PREFIX" ] && [ -f "$WORK/.xspec_env" ] && ENV_PREFIX=$(cat "$WORK/.xspec_env")
if [ -n "$ENV_PREFIX" ]; then
  [ -d "$ENV_PREFIX/bin" ] || { echo "no environment at $ENV_PREFIX/bin" >&2; exit 1; }
  export PATH="$ENV_PREFIX/bin:$PATH"
  LOG "environment: $ENV_PREFIX"
fi

# Modules stay available for anything the environment does not carry, but are
# optional and empty by default.
if [ -n "${XSPEC_MODULES:-}" ] || [ -f "$WORK/.xspec_modules" ]; then
  # shellcheck disable=SC1091
  [ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh || true
  mods="${XSPEC_MODULES:-}"
  [ -z "$mods" ] && [ -f "$WORK/.xspec_modules" ] && mods=$(cat "$WORK/.xspec_modules")
  for m in $mods; do
    LOG "module load $m"
    module load "$m" || { echo "module load $m failed" >&2; exit 1; }
  done
fi

# Every binary is checked by hand, whatever provided it. Failing in the first
# ten seconds beats failing after the queue wait and an hour of indexing.
missing=""
for b in bwa hisat2 hisat2-build hisat2_extract_splice_sites.py samtools perl python3; do
  command -v "$b" >/dev/null 2>&1 || missing="$missing $b"
done
if [ -n "$missing" ]; then
  cat >&2 <<EOF
missing from PATH:$missing

ARC has no module for bwa, hisat2 or samtools. Build the environment once:

  $HERE/setup_env.sh /home/exacloud/gscratch/<YourLab>/envs/hyperedit

then resubmit with --env pointing at it. Apptainer (1.4.1) is an alternative
if you would rather use a container.
EOF
  exit 1
fi
LOG "bwa $(bwa 2>&1 | awk '/^Version/{print $2}'), hisat2 $(hisat2 --version 2>&1 | head -1 | awk '{print $NF}'), samtools $(samtools --version 2>&1 | head -1 | awk '{print $2}')"

# --- node-local scratch for stage 2 ------------------------------------------
# Stage 2 writes the whole unmapped pool twelve times per arm, as SAM and again
# as BAM, and none of it is wanted afterwards -- only the beds are. That is the
# heaviest I/O in the run, so it goes on the node's own disk rather than across
# the network to gscratch.
#
# ACC's mkdir-scratch.sh makes /mnt/scratch/$SLURM_JOB_ID; the matching rmdir
# must run whatever happens, which is what the trap is for. The ARC storage
# guide asks that jobs clean up after themselves "regardless of job success",
# and data left there is deleted without warning anyway.
#
# Requires --gres disk:N on the submission, which submit.sh sets. Without the
# gres the job may land on a node with no /mnt/scratch at all, so falling back
# to gscratch rather than failing is deliberate: it is slower, not wrong.
cleanup () {
  if [ -n "${SCRATCH_MADE:-}" ]; then
    LOG "tearing down $SCRATCH_PATH"
    if command -v rmdir-scratch.sh >/dev/null 2>&1; then
      rmdir-scratch.sh || true
    else
      rm -rf "$SCRATCH_PATH" || true
    fi
  fi
}
trap cleanup EXIT INT TERM

if [ -n "${SLURM_JOB_ID:-}" ] && [ -d /mnt/scratch ]; then
  if command -v mkdir-scratch.sh >/dev/null 2>&1; then
    mkdir-scratch.sh && SCRATCH_MADE=1
  fi
  SCRATCH_PATH="/mnt/scratch/${SLURM_JOB_ID}"
  if [ -d "$SCRATCH_PATH" ]; then
    export XSPEC_TMP="$SCRATCH_PATH/stage2"
    mkdir -p "$XSPEC_TMP"
    LOG "stage 2 scratch: $XSPEC_TMP ($(df -BG --output=avail "$SCRATCH_PATH" | tail -1 | tr -dc '0-9') GB free)"
  else
    LOG "no /mnt/scratch/$SLURM_JOB_ID; stage 2 falls back to $WORK (slower)"
  fi
fi

# --- resources ---------------------------------------------------------------
# THREADS follows the allocation rather than a guess, so the same script is
# right whether it got 8 cores or 24.
export THREADS="${SLURM_CPUS_PER_TASK:-${THREADS:-8}}"

# INDEX_JOBS is a memory knob, not a CPU one, and it is the likeliest way to
# have the job killed for OOM. It is derived from the memory actually allocated
# rather than guessed, and capped at the six transforms.
#
# The 12 GB per concurrent build is deliberately above the ~8 GB usually quoted
# for hisat2-build on a mammalian genome: that figure is for the 3.1 Gb human
# assembly, and opossum's mMonDom1 is 3.59 Gb. 8 GB is held back for everything
# else in the job. Being wrong here costs a restart, so it errs low.
if [ -z "${INDEX_JOBS:-}" ]; then
  mem_mb="${SLURM_MEM_PER_NODE:-}"
  if [ -z "$mem_mb" ] && [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
    mem_mb=$(( SLURM_MEM_PER_CPU * THREADS ))
  fi
  if [ -n "$mem_mb" ]; then
    INDEX_JOBS=$(( (mem_mb / 1024 - 8) / 12 ))
    [ "$INDEX_JOBS" -lt 1 ] && INDEX_JOBS=1
    [ "$INDEX_JOBS" -gt 6 ] && INDEX_JOBS=6
  else
    INDEX_JOBS=2
  fi
  export INDEX_JOBS
fi
LOG "threads=$THREADS index_jobs=$INDEX_JOBS arms_at_once=${ARMS_AT_ONCE:-2}"

"$XSPEC/run_xspec.sh" "$WORK" "$SPECIES" "$ONLY"
LOG "job ${SLURM_JOB_ID:-local} finished cleanly"
