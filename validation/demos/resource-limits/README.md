# Resource limits beyond CPU and memory

Operating systems cap more than CPU time and memory. On shared clusters the caps that
bite are usually per-user thread/process counts (`ulimit -u`, `RLIMIT_NPROC`, which
Linux counts **per thread**), per-process open files (`ulimit -n`, often 1024 or 4096 on
cluster nodes), and per-process memory maps (`vm.max_map_count`, default 65530). This
demo measures what ScaleMethyl's `met_extract.py` and alnbase hold of each, and shows
what happens when a cap is hit.

`run.sh` reproduces everything below (see its header for the environment it needs).
`../../tools/resource_monitor.py` samples a command's whole process tree: processes,
threads, open file descriptors, memory maps and RSS.

Measured on 2026-09-16 on a 20-core laptop, alnbase 0.1.1, ScaleMethyl's pinned Python
stack (polars 0.20.18, pyarrow 16.1.0, duckdb 0.10.1), synthetic input of 40 contigs,
80,000 reads, 50 cells.

## ScaleMethyl `met_extract.py`

The Nextflow module runs it with `--subprocesses task.cpus - 1` (4 for the default 5
CPUs). It starts a `multiprocessing.Pool` of that many processes, one task per
chromosome; each process runs polars and writes with pyarrow, and the parent then runs
DuckDB. Each library sizes its own thread pool.

| Run | Processes | Threads (tree) | Max threads in one process |
|---|---|---|---|
| all 20 cores visible | 5 | 108 | 40 |
| inside a 5-core affinity mask (`taskset -c 0-4`) | 5 | 48 | 24 |

Which pools respect the allocation, inside the 5-core mask:

| Library | Pool size |
|---|---|
| polars | 5 (honours the affinity mask) |
| pyarrow | 20 (all cores on the machine) |
| duckdb | 20 (all cores on the machine) |

**Consequence on a cluster.** SLURM and most schedulers confine a job with an affinity
mask or cpuset. pyarrow and DuckDB size their pools from the *machine's* core count,
so on a 128-core node a 5-CPU job gets roughly 128 pyarrow threads per process plus 128
DuckDB threads in the parent. They are mostly idle, but every one counts against the
user's `RLIMIT_NPROC`, and many concurrent per-sample jobs from one user add up.

**Concurrent extractions multiply it.** A Nextflow pipeline on a local executor runs
one extraction per sample at the same time on one node. Three concurrent runs sharing
the same 5-core mask peaked at **19 processes, 148 threads, 272 open file descriptors
and 1.8 GB RSS** here. The thread count scales with the machine's cores, not the
allocation, so on a 128-core node each run carries roughly four worker processes of
about 128 pyarrow threads plus the parent's DuckDB pool, and several samples at once
reach thousands of threads for one user.

**What happens at the cap.** Under a per-user thread limit 40 above the user's
existing threads (the run needs about 108), a polars worker panics with
`could not spawn threads: ... Resource temporarily unavailable`. The `multiprocessing`
pool does not recover: the job **hangs** until killed (here by a 120 s timeout), leaving
14 partial parquet files. In a pipeline this looks like a stuck task, not an error.

The fix is ordinary configuration (`POLARS_MAX_THREADS`, `pyarrow.set_cpu_count`,
DuckDB `SET threads`, sized to the allocation divided by the number of processes), and
a pool that fails the job when a worker dies. Neither is in the script.

## alnbase 0.1.1

| Run | Threads | Open fds | RSS |
|---|---|---|---|
| tagged BAM, `-@ 4` | 16 | 5 | 31 MB |
| parquet, `-@ 4 --shards-per-worker 2` (8 files) | 10 | 12 | 276 MB |

Open files grow with `threads × shards-per-worker` (one per output file) and memory with
`batch-rows × files`. Threads are bounded by the command line, not the machine.

Failures found, all fixed or to be fixed in the next patch release:

1. **The real error is hidden in parquet runs.** A missing output directory, or hitting
   the open-file limit, is reported only as
   `worker 0 stopped early (sending on a closed channel)`. The missing-directory case
   needs no resource limit at all.
2. **Partial parquet output is left behind** when a parquet run fails (97 of 128 files
   in the open-file test), and some of those files are valid parquet with a subset of
   the rows.
3. **Thread-spawn failure panics** (exit 101, a Rust panic message) instead of a clean
   error, and a tagged-BAM run that panics this way leaves its partial `out.bam`,
   although a failed BAM run otherwise removes it.

alnbase at least fails promptly rather than hanging.
