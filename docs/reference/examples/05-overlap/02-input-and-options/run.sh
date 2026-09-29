#!/usr/bin/env bash
# alnbase reference, section 05 (overlap), example 02-input-and-options. See README.txt.
set -uo pipefail
ALNBASE=${ALNBASE:-alnbase}
cd "$(dirname "$0")"
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
samtools view --no-PG -b -o "$T/in.bam" reads.sam
# run LABEL [overlap options...]: run overlap on $T/in.bam, print stderr and the output records.
run() {
  local label=$1; shift
  echo "### $label"
  local shown="$*"; echo "\$ alnbase overlap ${shown//$T\//} in.bam out.bam"
  rm -f "$T/out.bam"
  "$ALNBASE" overlap "$@" "$T/in.bam" "$T/out.bam" 2>"$T/err"; local rc=$?
  echo "exit status: $rc"
  sed 's/^/stderr| /' "$T/err"
  if [ -s "$T/out.bam" ]; then samtools view "$T/out.bam"; fi
  echo
}
run "accepted: @HD SO:unsorted GO:query (as written)"
mk() { samtools view --no-PG -b -o "$T/in.bam" "$T/in.sam"; }
sed 's/^@HD.*/@HD\tVN:1.6\tSO:queryname/' reads.sam > "$T/in.sam"; mk
run "accepted: @HD SO:queryname"
sed 's/^@HD.*/@HD\tVN:1.6\tSO:coordinate/' reads.sam > "$T/in.sam"; mk
run "rejected: SO:coordinate"
sed 's/^@HD.*/@HD\tVN:1.6/' reads.sam > "$T/in.sam"; mk
run "rejected: @HD with neither SO nor GO"
grep -v '^@HD' reads.sam > "$T/in.sam"; mk
run "rejected: no @HD line"
sed 's/^@HD.*/@HD\tVN:1.6\tSO:coordinate\n@CO\tGO:query/' reads.sam > "$T/in.sam"; mk
run "rejected: GO:query on a @CO line does not count"
# Header claims grouping, but records are interleaved: pair_overlap R1, pair_b R1, pair_overlap R2, pair_b R2.
{ sed 's/^@HD.*/@HD\tVN:1.6\tSO:queryname/' reads.sam | grep '^@'; grep -v '^@' reads.sam | awk 'NR==1||NR==3'; grep -v '^@' reads.sam | awk 'NR==2||NR==4'; } > "$T/in.sam"; mk
run "NOT detected: header says SO:queryname but records are interleaved (only the header is checked)" --stats -
samtools view --no-PG -b -o "$T/in.bam" reads.sam
echo "=== option validation (each is refused before any record is read) ==="
for opts in "--stale-tags recompute" "--tag-prefix x" "--tag-prefix X" "--tag-prefix ov" "--tag-prefix 1" \
            "--length-tag XM" "--length-tag nm" "--length-tag X" "--length-tag 1L" "--length-tag oL" \
            "--qual-cap 94" "--min-support 0" "--min-overlap 0" "--min-span-frac 1.5" "--max-mismatch-frac=-0.1" "--keep r3" "--qual-cap 300"; do
  run "invalid: $opts" $opts
done
run "allowed: --tag-prefix x with --no-tag" --tag-prefix x --no-tag
run "allowed: --length-tag oL with --no-tag" --length-tag oL --no-tag
run "renamed tags: --tag-prefix q --length-tag ZF" --tag-prefix q --length-tag ZF
echo "### QUAL '*' on pair_overlap stays '*' (output piped through cat -v to show the bytes)"
awk -F'\t' 'BEGIN{OFS="\t"} /^@/{print;next} $1=="pair_overlap"{$11="*";print}' reads.sam | samtools view --no-PG -b -o "$T/q.bam" -
"$ALNBASE" overlap "$T/q.bam" - | samtools view - | cut -f1-6,11 | cat -v
echo
echo "### --compression-level sets the output compression (output sizes in bytes)"
for l in 0 3 9; do "$ALNBASE" overlap --compression-level $l "$T/in.bam" "$T/c$l.bam"; printf 'level %s: ' $l; wc -c < "$T/c$l.bam"; done
echo
echo "### --no-tag on an already-tagged file leaves the old tags in place"
"$ALNBASE" overlap "$T/in.bam" "$T/pass1.bam"
"$ALNBASE" overlap --no-tag "$T/pass1.bam" - | samtools view - | cut -f1,2,12-
echo
echo "### second pass over the first pass's output (default options)"
"$ALNBASE" overlap "$T/in.bam" "$T/pass1.bam"
"$ALNBASE" overlap --stats - "$T/pass1.bam" "$T/pass2.bam"
samtools view "$T/pass2.bam"
