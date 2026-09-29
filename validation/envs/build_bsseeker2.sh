#!/usr/bin/env bash
# Put BS-Seeker2 into the bs-seeker2 environment at a pinned commit.
#   build_bsseeker2.sh ENV_PREFIX [WORKDIR]
# BS-Seeker2 is a set of Python 2 scripts with nothing to compile, so this is a
# checkout rather than a build; the scripts land in $PREFIX/share/BSseeker2.
set -euo pipefail
PREFIX="$1"; WORK="${2:-$(mktemp -d)}"
COMMIT=0976fee742a0d90a698a8985fb705c162c4aca7a   # v2.1.8, the tip of master
git clone -q https://github.com/BSSeeker/BSseeker2 "$WORK/BSseeker2"
cd "$WORK/BSseeker2" && git checkout -q "$COMMIT"
rm -rf "$PREFIX/share/BSseeker2"
mkdir -p "$PREFIX/share"
cp -r "$WORK/BSseeker2" "$PREFIX/share/BSseeker2"
echo "$COMMIT" > "$PREFIX/share/bsseeker2-commit.txt"
