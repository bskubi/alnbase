#!/usr/bin/env bash
# Submit the 2x2 to ARC as a chain of jobs.
#
#   submit.sh WORKDIR SPECIES -A <account> [options]
#
#   -A, --account   slurm account (ARC project). Required: the ARC scheduler
#                   guide calls specifying it "a best practice" on every
#                   invocation, and without it the work is billed to whatever
#                   your default account happens to be.
#   -p, --partition default batch
#   -c, --cpus      default 16
#       --mem       default 64G
#       --time      default 36:00:00
#       --qos       default none; see "Why a chain" below
#       --disk      GB of node-local /mnt/scratch, default 400
#       --env       conda/mamba environment prefix built by setup_env.sh.
#                   Required unless the tools are already on PATH: ARC has no
#                   bwa, hisat2 or samtools module.
#       --modules   optional, quoted list passed to `module load`
#       --dry-run   print the sbatch commands and stop
#
# WHY A CHAIN, NOT ONE LONG JOB
#
# ARC's batch partition has a 36-hour limit, and the opossum pair needs roughly
# forty: about twenty hours of indexing per assembly plus alignment. One job
# cannot hold it without --qos long_jobs (10 days) or very_long_jobs (30 days,
# but capped at 24 CPUs).
#
# Splitting per assembly avoids needing a QOS at all. Each job is ~20 h, inside
# the default limit, and the two are independent up to the scoring step. It also
# schedules better: ARC's backfill scheduler starts shorter jobs sooner, and the
# guide notes that asking for more time than you need means waiting longer in
# queue. The chain is dependency-ordered only so the scoring job runs last --
# the assembly jobs themselves could run concurrently, and will if the cluster
# has room.
#
# `afterany` rather than `afterok` is deliberate. run_xspec.sh is resumable, so
# a job that dies on a wall-clock limit has still moved the work forward; the
# scoring job runs regardless and reports whichever cells are finished. Nothing
# is silently lost, and there is no chain to restart by hand.
#
# If a job does time out, resubmitting this same command is the fix: finished
# stages are skipped by their marker files.
set -euo pipefail

WORK="${1:?usage: submit.sh WORKDIR SPECIES -A <account> [options]}"; shift
SPECIES="${1:?usage: submit.sh WORKDIR SPECIES -A <account> [options]}"; shift
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
XSPEC="$(cd "$HERE/.." && pwd)"

# Site defaults, if arc/site.sh exists. Explicit flags below override them:
# site.sh uses `: "${VAR:=...}"`, so anything already in the environment wins,
# and the argument loop wins over both.
# shellcheck disable=SC1091
[ -f "$HERE/site.sh" ] && . "$HERE/site.sh"

ACCOUNT=""; PARTITION=batch; CPUS=16; MEM=64G; TIME=36:00:00; QOS=""
DISK=400; MODULES=""; ENVPREFIX=""; DRY=""
while [ $# -gt 0 ]; do
  case "$1" in
    -A|--account)   ACCOUNT=$2; shift 2 ;;
    -p|--partition) PARTITION=$2; shift 2 ;;
    -c|--cpus)      CPUS=$2; shift 2 ;;
    --mem)          MEM=$2; shift 2 ;;
    --time)         TIME=$2; shift 2 ;;
    --qos)          QOS=$2; shift 2 ;;
    --disk)         DISK=$2; shift 2 ;;
    --modules)      MODULES=$2; shift 2 ;;
    --env)          ENVPREFIX=$2; shift 2 ;;
    --dry-run)      DRY=1; shift ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

[ -z "$ACCOUNT" ] && ACCOUNT="${XSPEC_ACCOUNT:-}"
[ -z "$ENVPREFIX" ] && ENVPREFIX="${XSPEC_ENV:-}"
[ -n "$ENVPREFIX" ] && [ ! -d "$ENVPREFIX" ] && ENVPREFIX=""   # not built yet

if [ -z "$ACCOUNT" ]; then
  echo "-A/--account is required. Your accounts:" >&2
  sshare -U -u "$USER" 2>/dev/null >&2 || echo "  (run: sshare -U -u $USER)" >&2
  exit 1
fi

WORK="$(cd "$WORK" && pwd)"
# Compute nodes on ARC do have network access, so a job can fetch its own
# inputs. Doing it beforehand is still better -- downloading several GB while
# holding a 16-core allocation wastes the allocation, and a transfer that fails
# should not cost a queue wait -- but it is advice, not a requirement.
[ -f "$WORK/hyper-editing-commit.txt" ] || {
  echo "note: $WORK is not prepared; the first job will fetch its own inputs." >&2
  echo "      Cheaper to run first:  $HERE/fetch.sh $WORK $SPECIES" >&2
}
case "$WORK" in
  /home/users/*) echo "refusing to run under /home/users; use gscratch" >&2; exit 1 ;;
esac

LABELS=$(awk -F'\t' -v s="$SPECIES" '$1==s && $1!~/^#/ {print $2}' "$XSPEC/assemblies.tsv")
[ -n "$LABELS" ] || { echo "no assemblies for '$SPECIES' in assemblies.tsv" >&2; exit 1; }

mkdir -p "$WORK/logs"

# The module list goes in a file rather than through --export. Slurm's --export
# takes comma-separated KEY=VALUE pairs and its handling of values containing
# spaces varies between versions; a silently truncated list would show up as a
# missing-binary failure after the queue wait. The file also means a resubmit
# uses the same modules without repeating the flag.
if [ -n "$MODULES" ]; then printf '%s\n' "$MODULES" > "$WORK/.xspec_modules"; fi
if [ -n "$ENVPREFIX" ]; then
  ENVPREFIX="$(cd "$ENVPREFIX" && pwd)"
  [ -x "$ENVPREFIX/bin/hisat2" ] || {
    echo "$ENVPREFIX does not look like a built environment (no bin/hisat2)." >&2
    echo "Build it with: $HERE/setup_env.sh $ENVPREFIX" >&2
    exit 1
  }
  printf '%s\n' "$ENVPREFIX" > "$WORK/.xspec_env"
elif [ ! -f "$WORK/.xspec_env" ]; then
  echo "note: no --env given and no environment recorded for this workdir." >&2
  echo "      ARC has no bwa/hisat2/samtools module, so unless they are already" >&2
  echo "      on PATH the jobs will stop immediately. See arc/setup_env.sh." >&2
fi
common=( -A "$ACCOUNT" -p "$PARTITION" -c "$CPUS" --mem "$MEM" --time "$TIME"
         --gres "disk:$DISK" )
[ -n "$QOS" ] && common+=( --qos "$QOS" )

submit () {  # submit ONLY [dependency]
  local only=$1 dep=${2:-}
  local args=( "${common[@]}"
               --job-name "xspec-$SPECIES-$only"
               --output "$WORK/logs/%x-%j.out"
               --error  "$WORK/logs/%x-%j.err" )
  [ -n "$dep" ] && args+=( --dependency "afterany:$dep" )
  args+=( --export "ALL,XSPEC_WORK=$WORK,XSPEC_SPECIES=$SPECIES,XSPEC_ONLY=$only" )
  if [ -n "$DRY" ]; then
    { printf 'sbatch'; printf ' %q' "${args[@]}" "$HERE/job.sh"; printf '\n'; } >&2
    echo "<jobid-$only>"
  else
    sbatch --parsable "${args[@]}" "$HERE/job.sh"
  fi
}

ids=""
for label in $LABELS; do
  id=$(submit "$label")
  echo "submitted $label -> $id"
  ids="${ids:+$ids:}$id"
done

# The scoring job is tiny -- it reads bed files and prints a table -- so it is
# submitted with one core and a short clock rather than inheriting the indexing
# allocation. Short, small jobs are exactly what ARC's backfill scheduler slots
# in early, so this costs essentially no extra queue time.
score_args=( -A "$ACCOUNT" -p "$PARTITION" -c 1 --mem 4G --time 20:00
             --job-name "xspec-$SPECIES-score"
             --output "$WORK/logs/%x-%j.out" --error "$WORK/logs/%x-%j.err"
             --export "ALL,XSPEC_WORK=$WORK,XSPEC_SPECIES=$SPECIES,XSPEC_ONLY=score" )
[ -n "$QOS" ] && score_args+=( --qos "$QOS" )
[ -n "$ids" ] && score_args+=( --dependency "afterany:$ids" )
if [ -n "$DRY" ]; then
  printf 'sbatch'; printf ' %q' "${score_args[@]}" "$HERE/job.sh"; printf '\n'
else
  echo "submitted score -> $(sbatch --parsable "${score_args[@]}" "$HERE/job.sh")"
  echo
  echo "watch:   squeue -u $USER"
  echo "logs:    $WORK/logs/"
  echo "result:  $WORK/xspec_$SPECIES.tsv"
fi
