#!/usr/bin/env bash
# Cross-build alnbase for a cluster whose glibc is older than this machine's.
#
#   ./scripts/build-cluster.sh              # targets glibc 2.35
#   ./scripts/build-cluster.sh 2.31         # or whatever `ldd --version` says there
#
# Uses cargo-zigbuild, which borrows zig's bundled clang to link against an
# older glibc than the one installed here. Installs what is missing.

set -euo pipefail

GLIBC="${1:-2.35}"
TARGET="x86_64-unknown-linux-gnu.${GLIBC}"
OUT="target/x86_64-unknown-linux-gnu/release/alnbase"

say() { printf '\n\033[1m%s\033[0m\n' "$*"; }

# --- conda ----------------------------------------------------------------
# An active environment puts its own compiler, zlib and openssl ahead of the
# system's. The build then either fails inside htslib for no visible reason,
# or succeeds and produces a binary that needs the environment at runtime --
# which defeats the point of building something to copy to a cluster.
if [[ -n "${CONDA_DEFAULT_ENV:-}" && "${CONDA_DEFAULT_ENV}" != "base" ]]; then
    cat >&2 <<EOF
The conda environment '${CONDA_DEFAULT_ENV}' is active. Deactivate it and run
this again:

    conda deactivate
    $0 $*

(Deactivating matters even for installing ziglang below: cargo-zigbuild finds
zig through the python on PATH, so a zig installed into a conda environment
disappears the moment you deactivate.)
EOF
    exit 1
fi

# --- zig ------------------------------------------------------------------
if command -v zig >/dev/null 2>&1; then
    say "zig: $(command -v zig)"
elif python3 -c 'import ziglang' 2>/dev/null; then
    say "zig: the ziglang python package"
else
    say "Installing the ziglang python package (~50 MB)"
    # --user, not a virtualenv, so it stays available in a plain shell.
    python3 -m pip install --user ziglang
    python3 -c 'import ziglang' || {
        echo "ziglang installed but not importable; is ~/.local on PYTHONPATH?" >&2
        exit 1
    }
fi

# --- cargo-zigbuild -------------------------------------------------------
if command -v cargo-zigbuild >/dev/null 2>&1; then
    say "cargo-zigbuild: $(cargo-zigbuild --version 2>/dev/null || echo present)"
else
    say "Installing cargo-zigbuild (compiles from source, a few minutes)"
    cargo install --locked cargo-zigbuild
fi

# --- build ----------------------------------------------------------------
# target-cpu=x86-64-v3: hts-sys only compiles its AVX2 CRAM decoders when the
# Rust target features say AVX2 is available, and that is where a scan spends
# its time. v3 is Haswell (2013) and later. Not v4: hts-sys has no AVX-512
# path, so it would only narrow which nodes can run the binary.
#
# CFLAGS=-mevex512: zig's clang refuses libdeflate's AVX-512 CRC32 code under
# v3 alone. The flag is safe -- that code sits behind a runtime dispatch and
# only executes on CPUs that have AVX-512.
export RUSTFLAGS="${RUSTFLAGS:--C target-cpu=x86-64-v3}"
export CFLAGS="${CFLAGS:--mevex512}"

say "Building for glibc ${GLIBC}"
echo "  RUSTFLAGS=$RUSTFLAGS"
echo "  CFLAGS=$CFLAGS"
cargo zigbuild --release --target "$TARGET"

# --- check ----------------------------------------------------------------
say "Checking the result"
need=$(objdump -T "$OUT" | grep -o 'GLIBC_[0-9.]*' | sort -uV | tail -1)
echo "  highest glibc symbol required: ${need#GLIBC_}  (the cluster has ${GLIBC})"
if [[ "$(printf '%s\n%s\n' "${need#GLIBC_}" "$GLIBC" | sort -V | tail -1)" != "$GLIBC" ]]; then
    echo "  ERROR: that is newer than the target. The binary will not start there." >&2
    exit 1
fi

simd=$(nm -C "$OUT" 2>/dev/null | grep -c 'rans_uncompress_O1_32x16_avx2' || true)
echo "  AVX2 CRAM decoder: $([[ "$simd" -gt 0 ]] && echo present || echo 'MISSING -- check RUSTFLAGS')"
echo "  size: $(du -h "$OUT" | cut -f1)"

say "Done: $OUT"
cat <<EOF
Copy it over and check it there:

    scp $OUT CLUSTER:~/bin/alnbase
    ssh CLUSTER '~/bin/alnbase --version && ~/bin/alnbase --help | head -3'

To run the test suite on the cluster too:

    cargo zigbuild test --release --no-run --target $TARGET
    # then scp the test binary printed above and run it there
EOF
