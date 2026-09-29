#!/usr/bin/env bash
# Reproduce the resource-limit findings in README.md.
#
#   PY=/path/to/python-with-scalemethyl-deps \
#   MET_EXTRACT=/path/to/ScaleMethyl/bin/met_extract.py \
#   ALNBASE=/path/to/alnbase \
#   ./run.sh [workdir]
#
# PY needs ScaleMethyl's pinned stack (envs/scaleMethylTools.conda.yml:
# polars 0.20.18, pyarrow 16.1.0, python-duckdb 0.10.1, pysam 0.22.1,
# pyfaidx 0.8.1.4) plus psutil. Uses util-linux `prlimit` and `taskset`.
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
MONITOR="$HERE/../../tools/resource_monitor.py"
QUERIES="$HERE/../../../tests/methyldackel/query.toml"
WORK="$(realpath -m "${1:-$(mktemp -d)}")"
export WORK
: "${PY:?set PY}" "${MET_EXTRACT:?set MET_EXTRACT}" "${ALNBASE:?set ALNBASE}"
mkdir -p "$WORK" && cd "$WORK"

say() { printf '\n== %s\n' "$*"; }

# Leaves headroom of $1 threads above what this user is already running.
thread_limit_with_headroom() { echo $(( $(ps -L -u "$(id -un)" --no-headers | wc -l) + $1 )); }

say "input: 40 contigs x 2000 reads, 50 cells"
"$PY" "$HERE/make_scalemethyl_input.py" --outdir data

met_extract() {  # run met_extract.py in its own directory, as the pipeline does
    rm -rf "$1" && mkdir -p "$1" && cp data/reads.bam data/reads.bam.bai "$1"/
    (cd "$1" && shift && "$@" "$PY" "$MET_EXTRACT" reads.bam --sample s1 --threshold 0.5 \
        --subprocesses 4 --contexts CG,CH --aligner bwa-meth --ref "$WORK/data/ref.fa")
}

say "1. met_extract.py fan-out, 4 subprocesses, all cores visible"
met_extract me_open "$PY" "$MONITOR" --out samples.tsv --interval 0.05 --

say "2. the same inside a 5-core affinity mask (what a scheduler allocation looks like)"
met_extract me_taskset taskset -c 0-4 "$PY" "$MONITOR" --out samples.tsv --interval 0.05 --

say "3. which thread pools honour the affinity mask"
taskset -c 0-4 "$PY" - <<'EOF'
import os, polars, pyarrow, duckdb
print(f"cores visible to this process: {len(os.sched_getaffinity(0))} of {os.cpu_count()}")
print(f"polars thread pool:  {polars.thread_pool_size()}")
print(f"pyarrow cpu_count:   {pyarrow.cpu_count()}")
print(f"duckdb threads:      {duckdb.sql('select current_setting(%s)' % repr('threads')).fetchone()[0]}")
EOF

say "3b. three extractions at once, as a local executor runs per-sample tasks"
three_at_once() {
    for i in 1 2 3; do
        met_extract "me_concurrent/s$i" taskset -c 0-4 &
    done
    wait
}
export -f met_extract three_at_once
export PY MET_EXTRACT
rm -rf me_concurrent && mkdir -p me_concurrent
"$PY" "$MONITOR" --out me_concurrent/samples.tsv --interval 0.05 -- bash -c three_at_once

LIMIT=$(thread_limit_with_headroom 40)
say "4. met_extract.py under a per-user thread limit of $LIMIT (40 above current); expect a hang"
met_extract me_nproc timeout 120 prlimit --nproc="$LIMIT:$LIMIT" > me_nproc.log 2>&1
echo "exit status $? (124 = killed by the 120 s timeout)"
grep -m1 -o 'could not spawn threads.*' me_nproc.log
echo "partial parquet files left: $(ls me_nproc/*.parquet 2>/dev/null | wc -l)"

"$ALNBASE" index data/ref.fa data/ref.aref > /dev/null

say "5. alnbase parquet run into a directory that does not exist"
"$ALNBASE" query --query-file "$QUERIES" --parquet data/reads.bam data/ref.aref no_such_dir/calls.parquet
echo "exit status $?"

say "6. alnbase with an open-file limit of 32 and 100 output shards"
mkdir -p ab_nofile
prlimit --nofile=32:32 "$ALNBASE" query --query-file "$QUERIES" --parquet --shards-per-worker 100 \
    data/reads.bam data/ref.aref ab_nofile/calls.parquet
echo "exit status $?; partial files left: $(ls ab_nofile | wc -l)"

LIMIT=$(thread_limit_with_headroom 10)
say "7. alnbase -@ 16 under a per-user thread limit of $LIMIT (10 above current)"
mkdir -p ab_nproc
prlimit --nproc="$LIMIT:$LIMIT" "$ALNBASE" query --query-file "$QUERIES" -@ 16 \
    data/reads.bam data/ref.aref ab_nproc/out.bam
echo "exit status $?; left behind: $(ls ab_nproc)"
