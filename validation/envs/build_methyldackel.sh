#!/usr/bin/env bash
# Build MethylDackel from source into the methyldackel environment.
#   build_methyldackel.sh ENV_PREFIX [WORKDIR]
set -euo pipefail
PREFIX="$1"; WORK="${2:-$(mktemp -d)}"
COMMIT=3c77bda12141e99d80234d416e668a90ec70b3f7
git clone -q https://github.com/dpryan79/MethylDackel "$WORK/MethylDackel"
cd "$WORK/MethylDackel" && git checkout -q "$COMMIT"
make CC="$PREFIX/bin/cc" \
     CFLAGS="-Wall -g -O3 -pthread -I$PREFIX/include -I$PREFIX/include/libBigWig" \
     LIBS="-L$PREFIX/lib -Wl,-rpath,$PREFIX/lib" LIBBIGWIG="-lBigWig"
install -m 755 MethylDackel "$PREFIX/bin/MethylDackel"
echo "$COMMIT" > "$PREFIX/share/methyldackel-commit.txt"
