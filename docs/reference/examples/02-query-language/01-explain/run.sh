#!/usr/bin/env bash
# Reproduce with: bash run.sh > expected.txt 2>&1
set -euo pipefail
cd "$(dirname "$0")"
ALNBASE=${ALNBASE:-alnbase}
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
echo '$ alnbase query --query-file queries.toml --explain'
"$ALNBASE" query --query-file queries.toml --explain
echo
echo '$ alnbase query --query-file queries.toml --list-codes   (tail)'
"$ALNBASE" query --query-file queries.toml --list-codes | tail -n 5
