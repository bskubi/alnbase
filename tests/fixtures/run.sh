#!/usr/bin/env bash
# Run every fixture in this directory and diff alnbase's calls against the
# expected table each one was written with.
#
#   PY=/path/to/python-with-duckdb ALNBASE=/path/to/alnbase ./run.sh [case ...]
#
# With no arguments, every case directory runs. Needs samtools on PATH.
# Exit status 0 means every case matched.
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
: "${PY:?set PY}" "${ALNBASE:?set ALNBASE}"
STRAND="$HERE/../../query/strand/directional.toml"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

sorted() {  # sorted OUT.bam IN.bam
  samtools sort -o "$1" "$2" 2>/dev/null && samtools index "$1"
}

# One alnbase query over one BAM, printing the columns the case asks for.
#   calls BAM EXPECTED.TSV PREFIX   ->  PREFIX.tsv, PREFIX.scan
#
# A case asks for the columns it is about by naming them in its expected table's
# header, so a case about offsets can show off_5p without every other case
# carrying a column it has nothing to say about. The sort is over every selected
# column, since two rows can share qname, position and query name -- which is
# exactly what an unresolved mate overlap looks like.
calls() {
  "$ALNBASE" query --query-file "$HERE/context.toml" --query-file "$STRAND" \
    --parquet --only-hits -F qname "$1" ref.aref "${3}_out" > /dev/null 2> "$3.scan"
  EXPECT_COLS="$(head -1 "$2" | tr '\t' ',')" PREFIX="${3}_out" "$PY" - <<'PYEOF' > "$3.tsv"
import duckdb, glob, os
cols = os.environ["EXPECT_COLS"]
rows = duckdb.sql(
    "select " + cols + " "
    "from read_parquet(" + repr(sorted(glob.glob(os.environ["PREFIX"] + "_*"))) + ") "
    "order by " + cols).fetchall()
print(cols.replace(",", "\t"))
for r in rows:
    print("\t".join(str(x) for x in r))
PYEOF
}

cases=("$@")
if [ ${#cases[@]} -eq 0 ]; then
  mapfile -t cases < <(cd "$HERE" && find . -mindepth 2 -name expected.tsv -printf '%h\n' | sed 's|^\./||' | sort)
fi

failed=0
for case in "${cases[@]}"; do
  src="$HERE/$case"
  dir="$WORK/$case"
  mkdir -p "$dir" && cp "$src/ref.fa" "$src/reads.sam" "$dir/" && cd "$dir"

  # samtools rejects a record whose CIGAR and SEQ lengths disagree, which is a
  # free check on every hand-typed record.
  samtools faidx ref.fa
  "$ALNBASE" index ref.fa ref.aref >/dev/null
  samtools view -b --no-PG -o grouped.bam reads.sam

  bad=""
  : > diff.txt

  # A case with an overlap.args file (its contents are extra options, possibly
  # none) is checked twice: expected_raw.tsv is what a caller sees while the
  # duplicated observations are still there, expected.tsv what it sees once
  # `alnbase overlap` has resolved them. Half of what such a case claims is the
  # double count itself, so the fixture has to show it.
  if [ -f "$src/overlap.args" ]; then
    sorted raw.bam grouped.bam
    calls raw.bam "$src/expected_raw.tsv" raw
    if ! diff -u "$src/expected_raw.tsv" raw.tsv >> diff.txt; then bad="$bad raw-calls"; fi
    # shellcheck disable=SC2046  # the args file is deliberately word-split
    "$ALNBASE" overlap $(cat "$src/overlap.args") grouped.bam resolved.bam 2>> diff.txt
    sorted reads.bam resolved.bam
  else
    sorted reads.bam grouped.bam
  fi

  calls reads.bam "$src/expected.tsv" got
  # The scan line is alnbase's own account of what it read, scanned and skipped.
  scan="$(grep -o 'read .*' got.scan | head -1)"

  if ! diff -u "$src/expected.tsv" got.tsv >> diff.txt; then bad="$bad calls"; fi
  if [ -f "$src/expected_scan.txt" ] && [ "$scan" != "$(cat "$src/expected_scan.txt")" ]; then
    bad="$bad scan"
    printf 'expected scan: %s\n     got scan: %s\n' "$(cat "$src/expected_scan.txt")" "$scan" >> diff.txt
  fi

  if [ -z "$bad" ]; then
    printf 'ok    %-26s %s\n' "$case" "$scan"
  else
    failed=1
    printf 'FAIL  %-26s (%s)\n' "$case" "${bad# }"
    sed 's/^/        /' diff.txt
  fi
done

exit "$failed"
