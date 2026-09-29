#!/usr/bin/env bash
# Build BSBolt from source into the premethyst environment.
#   build_bsbolt.sh ENV_PREFIX
# BSBolt's setup.py at this commit leaves bsbolt.GenotypeMatrix out of its package list,
# so a pip install cannot start `bsbolt`. A development install runs from the checkout,
# which has every module; the checkout is kept inside the environment for that reason.
set -euo pipefail
PREFIX="$1"
COMMIT=ea4870e975c546d1eb9cbb4d7021537e3cca05bd
SRC="$PREFIX/share/bsbolt-src"
git clone -q https://github.com/NuttyLogic/BSBolt "$SRC"
cd "$SRC" && git checkout -q "$COMMIT"
PATH="$PREFIX/bin:$PATH" "$PREFIX/bin/python" setup.py develop
echo "$COMMIT" > "$PREFIX/share/bsbolt-commit.txt"
