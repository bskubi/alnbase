#!/usr/bin/env python3
"""Run a command and record the OS resources its whole process tree holds.

Wall time and peak memory are the usual benchmark numbers. This also samples
the resources an operating system caps separately, which is where
multi-process, multi-threaded tools tend to fail on shared machines:

    processes     the command and every descendant
    threads       summed over the tree (RLIMIT_NPROC counts threads, not processes)
    open_fds      file descriptors summed over the tree (RLIMIT_NOFILE is per process,
                  so the per-process maximum is recorded too)
    mmaps         memory mappings summed over the tree (vm.max_map_count is per process)
    rss_mb        resident memory summed over the tree

Usage:
    resource_monitor.py --out samples.tsv [--interval 0.1] -- CMD [ARGS...]

Writes one row per sample to --out and a one-line summary of the peaks to stderr.
The command's own exit status is returned.
"""

import argparse
import subprocess
import sys
import time

import psutil


def tree(root: psutil.Process) -> list[psutil.Process]:
    try:
        return [root] + root.children(recursive=True)
    except psutil.NoSuchProcess:
        return []


def count_maps(pid: int) -> int:
    try:
        with open(f"/proc/{pid}/maps") as fh:
            return sum(1 for _ in fh)
    except OSError:
        return 0


def sample(root: psutil.Process) -> dict:
    row = dict(processes=0, threads=0, open_fds=0, max_fds_one_process=0,
               max_threads_one_process=0, mmaps=0, max_mmaps_one_process=0, rss_mb=0.0)
    for p in tree(root):
        try:
            with p.oneshot():
                threads = p.num_threads()
                fds = p.num_fds()
                rss = p.memory_info().rss
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
        maps = count_maps(p.pid)
        row["processes"] += 1
        row["threads"] += threads
        row["open_fds"] += fds
        row["mmaps"] += maps
        row["rss_mb"] += rss / 2**20
        row["max_fds_one_process"] = max(row["max_fds_one_process"], fds)
        row["max_threads_one_process"] = max(row["max_threads_one_process"], threads)
        row["max_mmaps_one_process"] = max(row["max_mmaps_one_process"], maps)
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", required=True, help="TSV of samples")
    parser.add_argument("--interval", type=float, default=0.1, help="seconds between samples")
    parser.add_argument("cmd", nargs=argparse.REMAINDER, help="command, after --")
    args = parser.parse_args()
    cmd = args.cmd[1:] if args.cmd[:1] == ["--"] else args.cmd
    if not cmd:
        parser.error("no command given")

    start = time.monotonic()
    child = subprocess.Popen(cmd)
    root = psutil.Process(child.pid)
    columns = ["seconds", "processes", "threads", "max_threads_one_process", "open_fds",
               "max_fds_one_process", "mmaps", "max_mmaps_one_process", "rss_mb"]
    peaks = {c: 0 for c in columns[1:]}

    with open(args.out, "w") as out:
        out.write("\t".join(columns) + "\n")
        while child.poll() is None:
            row = sample(root)
            row["seconds"] = round(time.monotonic() - start, 2)
            out.write("\t".join(str(round(row[c], 1)) for c in columns) + "\n")
            for c in peaks:
                peaks[c] = max(peaks[c], row[c])
            time.sleep(args.interval)

    elapsed = time.monotonic() - start
    summary = " ".join(f"{c}={round(v, 1)}" for c, v in peaks.items())
    print(f"exit={child.returncode} seconds={elapsed:.1f} peak: {summary}", file=sys.stderr)
    return child.returncode


if __name__ == "__main__":
    sys.exit(main())
