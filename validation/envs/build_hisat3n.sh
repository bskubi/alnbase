#!/usr/bin/env bash
# Build HISAT-3N from source into the hisat-3n environment.
#   build_hisat3n.sh ENV_PREFIX [WORKDIR]
#
# HISAT-3N ships as the `hisat-3n` branch of the hisat2 repository rather than as its own
# project, and has no bioconda package, so the branch is the only source. The build produces
# hisat-3n and hisat-3n-build (perl wrappers) beside hisat2-align-s and hisat2-build-s, which
# are the binaries they invoke; all four are installed because the wrappers look for the
# others on PATH. Only the small-index binaries are built: the -l pair is for references
# above 4 Gbp, which no validation fixture here comes near.
set -euo pipefail
PREFIX="$1"; WORK="${2:-$(mktemp -d)}"
COMMIT=f5dda37bd1340f74ab91deace470aebc66e87a2d
git clone -q --branch hisat-3n https://github.com/DaehwanKimLab/hisat2 "$WORK/hisat-3n"
cd "$WORK/hisat-3n" && git checkout -q "$COMMIT"
make -j 4 hisat2-align-s hisat2-build-s
for f in hisat-3n hisat-3n-build hisat2-align-s hisat2-build-s; do
  install -m 755 "$f" "$PREFIX/bin/$f"
done
echo "$COMMIT" > "$PREFIX/share/hisat-3n-commit.txt"
