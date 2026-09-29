#!/usr/bin/env bash
# Build BISCUIT from source into the biscuit environment.
#   build_biscuit.sh ENV_PREFIX [WORKDIR]
set -euo pipefail
PREFIX="$1"; WORK="${2:-$(mktemp -d)}"
COMMIT=0a5ceaec623391a5df0a8d8f5d3735aa75f9bfc6
git clone -q https://github.com/huishenlab/biscuit "$WORK/biscuit"
cd "$WORK/biscuit" && git checkout -q "$COMMIT"
PATH="$PREFIX/bin:$PATH" cmake -S . -B build -DCMAKE_PREFIX_PATH="$PREFIX" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_C_COMPILER="$PREFIX/bin/cc"
# The source targets do not declare a dependency on the fetched libraries, so a parallel
# build can compile them first and fail; build the libraries before everything else.
PATH="$PREFIX/bin:$PATH" cmake --build build --target utils sgsl htslib
PATH="$PREFIX/bin:$PATH" cmake --build build -j 4
install -m 755 build/src/biscuit "$PREFIX/bin/biscuit"
echo "$COMMIT" > "$PREFIX/share/biscuit-commit.txt"
