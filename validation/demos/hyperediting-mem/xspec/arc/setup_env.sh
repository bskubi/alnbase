#!/usr/bin/env bash
# Build the tool environment once, on gscratch, where every compute node sees it.
#
#   setup_env.sh PREFIX
#
#   PREFIX   where the environment goes, e.g.
#            /home/exacloud/gscratch/<YourLab>/envs/hyperedit
#
# ARC has no bwa, hisat2 or samtools module -- `module avail` offers bowtie2,
# bamtools, bedtools2 and kallisto, but not these -- so the environment is
# built rather than loaded. It is ~1.5 GB and takes a few minutes.
#
# Uses whichever of mamba, micromamba or conda is on PATH, preferring mamba.
# Run it on the login node; it needs the network and no real compute.
#
# Put the result somewhere persistent and shared. A conda environment in a home
# directory is a common way to make a cluster job unreproducible: the ARC
# storage guide says /home/users is for shell profiles and scripts, not for
# data processing, and /mnt/scratch is wiped between jobs.
set -euo pipefail

PREFIX="${1:?usage: setup_env.sh PREFIX}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG () { echo "[$(date +'%F %T')] $*"; }

case "$PREFIX" in
  /mnt/scratch/*) echo "/mnt/scratch is node-local and wiped between jobs" >&2; exit 1 ;;
esac

MGR=""
for c in mamba micromamba conda; do
  command -v "$c" >/dev/null 2>&1 && { MGR=$c; break; }
done
[ -n "$MGR" ] || {
  echo "no mamba, micromamba or conda on PATH." >&2
  echo "If ARC provides one via a module, load it first; otherwise install" >&2
  echo "micromamba: https://mamba.readthedocs.io/en/latest/installation.html" >&2
  exit 1
}
LOG "using $MGR"

if [ -x "$PREFIX/bin/hisat2" ]; then
  LOG "environment already present at $PREFIX"
else
  LOG "creating $PREFIX from environment.yml"
  "$MGR" create -y -p "$PREFIX" -f "$HERE/environment.yml"
fi

# Verify by running the binaries, not by trusting the solve. A package can be
# present and still not execute -- a missing shared library shows up here and
# not in the install log.
LOG "verifying"
export PATH="$PREFIX/bin:$PATH"
fail=0
for b in bwa hisat2 hisat2-build hisat2_extract_splice_sites.py samtools; do
  if command -v "$b" >/dev/null 2>&1; then
    printf '  %-32s ok\n' "$b"
  else
    printf '  %-32s MISSING\n' "$b"; fail=1
  fi
done
bwa 2>&1 | awk '/^Version/{printf "  bwa      %s\n", $2}'
hisat2 --version 2>&1 | head -1 | awk '{printf "  hisat2   %s\n", $NF}'
samtools --version 2>&1 | head -1 | awk '{printf "  samtools %s\n", $2}'
[ "$fail" -eq 0 ] || { echo "environment is incomplete" >&2; exit 1; }

echo
echo "done. Submit with:"
echo "  $HERE/submit.sh <WORKDIR> <SPECIES> -A <YourLabAccount> --env $PREFIX"
