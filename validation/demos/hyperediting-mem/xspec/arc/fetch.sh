#!/usr/bin/env bash
# Pre-stage the inputs. OPTIONAL -- run_xspec.sh downloads whatever is missing,
# and ARC's compute nodes have network access, so a job can do this itself.
#
#   fetch.sh WORKDIR SPECIES
#
# Worth running anyway, for one reason: it is pure network and I/O with no
# computation in it, so doing it inside a batch job means a 16-core allocation
# sits idle through several GB of transfer. Running it on the login node first
# costs nothing and returns those core-hours. A second, smaller benefit is that
# a bad URL or a withdrawn ENA accession surfaces immediately instead of after
# a queue wait.
#
# It is small enough for the login node; if your session drops, run it under
# `screen` or `tmux`, or as an interactive job
# (`srun -p interactive --time=4:00:00 --pty bash -i`).
#
# WORKDIR must be on a filesystem the compute nodes can see and that survives
# between jobs, which on ARC means gscratch:
#
#   /home/exacloud/gscratch/<YourLab>/xspec
#
# Not /home/users -- the ARC storage guide says plainly that it is for shell
# profiles and scripts and "should not be used for any data processing in the
# cluster". Not /mnt/scratch either: that is node-local and wiped between jobs.
# /arc/scratch1 works and has no quota, but it is a technology preview that
# "will eventually be deleted", so it is a poor home for a week of indexing.
set -euo pipefail

WORK="${1:?usage: fetch.sh WORKDIR SPECIES}"
SPECIES="${2:?usage: fetch.sh WORKDIR SPECIES}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
XSPEC="$(cd "$HERE/.." && pwd)"
DEMO="$(cd "$XSPEC/.." && pwd)"

mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
LOG () { echo "[$(date +'%F %T')] $*"; }

case "$WORK" in
  /home/users/*) echo "refusing to stage data under /home/users; use gscratch" >&2; exit 1 ;;
  /mnt/scratch/*) echo "/mnt/scratch is node-local and wiped between jobs" >&2; exit 1 ;;
esac

# Space before rather than after. The indexes are the expensive thing here --
# about twenty hours per assembly -- and finding out they do not fit at hour
# nineteen is the worst available outcome.
NEED_GB="${NEED_GB:-200}"
avail=$(df -BG --output=avail "$WORK" | tail -1 | tr -dc '0-9')
LOG "space at $WORK: ${avail} GB available, ~${NEED_GB} GB wanted for a species pair"
if [ "${avail:-0}" -lt "$NEED_GB" ]; then
  echo "not enough room. gscratch quota is per-lab and purchased; check with" >&2
  echo "  df -h $WORK" >&2
  echo "and either free space or ask ACC (acc@ohsu.edu) to raise it." >&2
  exit 1
fi

# The pipeline itself, at its published commit.
if [ ! -f "$WORK/hyper-editing-commit.txt" ]; then
  LOG "cloning and patching the Hyper-editing tree"
  "$DEMO/setup.sh" "$WORK"
fi

# Reads. Resolved through ENA's API rather than a written-down path, because the
# directory convention depends on the accession's length and trailing digits.
RUN=$(awk -F'\t' -v s="$SPECIES" '$1==s && $1!~/^#/ {print $2; exit}' "$XSPEC/reads.tsv")
[ -n "$RUN" ] || { echo "no read run for '$SPECIES' in reads.tsv" >&2; exit 1; }
if [ ! -f "$WORK/reads/.done.$RUN" ]; then
  mkdir -p "$WORK/reads"
  LOG "resolving $RUN at ENA"
  URL=$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN&result=read_run&fields=fastq_ftp" \
        | awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="fastq_ftp") c=i; next}
                      c {print $c; exit}' | cut -d';' -f1)
  [ -n "$URL" ] || { echo "ENA returned no fastq path for $RUN" >&2; exit 1; }
  LOG "downloading $URL"
  curl -fsSL -C - -o "$WORK/reads/$RUN.fastq.gz" "ftp://$URL"
  touch "$WORK/reads/.done.$RUN"
fi
LOG "reads: $WORK/reads/$RUN.fastq.gz"

# References and annotations, one pair per assembly.
while IFS=$'\t' read -r sp label era fasta_url gtf_url; do
  case "$sp" in \#*|"") continue ;; esac
  [ "$sp" = "$SPECIES" ] || continue
  d="$WORK/$label"; ref="$d/$label.fa"
  mkdir -p "$d"

  if [ ! -f "$d/.done.fasta" ]; then
    LOG "$label ($era): fetching reference"
    curl -fsSL -C - -o "$d/src.gz" "$fasta_url"
    if [[ "$fasta_url" == *.tar.gz ]]; then
      mkdir -p "$d/chroms" && tar -xzf "$d/src.gz" -C "$d/chroms"
      find "$d/chroms" -name '*.fa' -print0 | sort -z | xargs -0 cat > "$ref"
      rm -rf "$d/chroms"
    else
      gunzip -c "$d/src.gz" > "$ref"
    fi
    rm -f "$d/src.gz"
    samtools faidx "$ref" 2>/dev/null || LOG "$label: no samtools here; .fai will be built in the job"
    touch "$d/.done.fasta"
  fi

  if [ ! -f "$d/.done.gtf" ]; then
    LOG "$label: fetching annotation"
    curl -fsSL -C - -o "$d/genes.gtf.gz" "$gtf_url"
    touch "$d/.done.gtf"
  fi
  LOG "$label: ready"
done < "$XSPEC/assemblies.tsv"

LOG "fetch complete. Submit with:"
echo "  $HERE/submit.sh $WORK $SPECIES -A <YourLabAccount>"
