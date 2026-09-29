#!/usr/bin/env bash
# What `alnbase index` accepts in a FASTA and what it makes of it.
set -uo pipefail
ALNBASE=${ALNBASE:-alnbase}
here=$(cd "$(dirname "$0")" && pwd)
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT
cp "$here"/ref.fa "$work"; cd "$work"
run() { echo "\$ alnbase $*"; "$ALNBASE" "$@" 2>&1; echo "(exit $?)"; }

echo "== names, comments, blank lines"
run index ref.fa plain.aref
run info plain.aref
run info --seq chr1:1-14 plain.aref

echo "== duplicate contig names are refused (before 0.1.2 both were stored and lookup found the last)"
printf '>dup first record named dup\nAAAA\n>dup second record named dup\nCCCCCCCC\n' > dup.fa
run index dup.fa dup.aref
ls dup.aref dup.aref.partial 2>&1 | sed 's/^ls: //'

echo "== compression and line endings give byte-identical indexes"
gzip -c ref.fa > ref.fa.gz            # plain gzip, not BGZF
bgzip -c ref.fa > ref.bgz.fa.gz       # BGZF
sed 's/$/\r/' ref.fa > crlf.fa        # CRLF line endings
"$ALNBASE" index ref.fa.gz gz.aref
"$ALNBASE" index ref.bgz.fa.gz bgz.aref
"$ALNBASE" index crlf.fa crlf.aref
for f in gz bgz crlf; do cmp plain.aref $f.aref && echo "$f.aref identical to plain.aref"; done

echo "== old-Mac CR-only line endings: the file is ONE line; accepted silently"
echo "   (comment lines removed first; otherwise the whole file is one ';' comment)"
grep -v '^;' ref.fa | tr '\n' '\r' > cr.fa
run index cr.fa cr.aref
run info cr.aref

echo "== errors"
printf 'ACGT\n>c1\nACGT\n' > before.fa
run index before.fa before.aref
echo "-- index writes OUT.partial and renames it only on success, so nothing is left behind:"
ls before.aref before.aref.partial 2>&1 | sed 's/^ls: //'
echo "-- and a failed run over an existing index leaves that index intact:"
cp plain.aref existing.aref
run index before.fa existing.aref
cmp plain.aref existing.aref && echo "existing.aref unchanged"
printf '>\nACGT\n' > noname.fa
run index noname.fa noname.aref
printf '> c1\nACGT\n' > spacename.fa
run index spacename.fa spacename.aref
echo "-- a non-UTF-8 name is written, but the index then refuses to open"
printf '>chr\xe91\nACGT\n' > latin1.fa
run index latin1.fa latin1.aref
run info latin1.aref
run index missing.fa missing.aref
