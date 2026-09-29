#!/usr/bin/env bash
# Render every figures/*.bob (svgbob ASCII art) to a same-named .svg.
#
#   SVGBOB=/path/to/svgbob_cli ./render.sh        # cargo install svgbob_cli
#
# The .bob files are the source of truth and read as plain text; the .svg files are
# committed so the docs render on GitHub without a build step. Re-run after editing.
set -euo pipefail
cd "$(dirname "$0")"
SVGBOB="${SVGBOB:-svgbob_cli}"
for src in *.bob; do
    "$SVGBOB" "$src" -o "${src%.bob}.svg"
    echo "rendered ${src%.bob}.svg"
done
