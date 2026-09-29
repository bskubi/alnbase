#!/usr/bin/env bash
# Reproduce with: bash run.sh > expected.txt 2>&1
set -euo pipefail
cd "$(dirname "$0")"
ALNBASE=${ALNBASE:-alnbase}
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
echo "$ alnbase query --query-file queries.toml --trace 'TTC.GTT@TTCAGTT'"
"$ALNBASE" query --query-file queries.toml --trace 'TTC.GTT@TTCAGTT'
echo "$ ... --trace 'c__@CGA' --trace-grid false   (lowercase and pads accepted; no flank is added)"
"$ALNBASE" query --query-file queries.toml --trace 'c__@CGA' --trace-grid false
