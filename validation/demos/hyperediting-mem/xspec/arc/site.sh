# Site settings for ARC. Sourced by submit.sh when present, and useful to
# source by hand so the long gscratch paths need not be retyped:
#
#   source arc/site.sh
#   arc/setup_env.sh "$XSPEC_ENV"
#   arc/fetch.sh     "$XSPEC_WORK" opossum
#   arc/submit.sh    "$XSPEC_WORK" opossum
#
# This is the one file in xspec/ that is specific to one person on one cluster.
# Edit it, or delete it and pass -A and --env explicitly.

# The slurm account work is billed to. ARC's scheduler guide calls naming it
# "a best practice" on every invocation, because without it jobs go to whatever
# your default account happens to be. Confirm yours with:  sshare -U -u $USER
: "${XSPEC_ACCOUNT:=YardimciLab}"

# Persistent, node-visible storage. gscratch rather than /home/users, which the
# ARC storage guide says is not for data processing, and rather than
# /mnt/scratch, which is node-local and wiped between jobs.
: "${XSPEC_WORK:=/home/exacloud/gscratch/YardimciLab/skubi/xspec}"

# The mamba environment built by setup_env.sh. Kept beside the work directory
# rather than inside it so that clearing a run's data does not destroy the
# toolchain, which takes a few minutes and a network round trip to rebuild.
: "${XSPEC_ENV:=/home/exacloud/gscratch/YardimciLab/skubi/envs/hyperedit}"

export XSPEC_ACCOUNT XSPEC_WORK XSPEC_ENV
