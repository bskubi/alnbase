#!/usr/bin/env python3
"""Compare alnbase's methylation calls against MethylDackel's, per position.

Both tools are given the *same* pre-filtered BAM, so a disagreement is a
disagreement about calling rather than about which reads to look at. That
matters more than it sounds: MethylDackel filters by flag and MAPQ by default
and alnbase filters nothing, so comparing them on a raw BAM mostly measures
the filters.

    ./compare_methyldackel.py sample.bam genome.fa genome.aref -o cmp

Writes a summary to stdout and, in the output directory:

    filtered.bam          what both tools actually read (and how big it is)
    alnbase_*.parquet     alnbase's per-read calls
    md_C*.bedGraph        MethylDackel's per-position calls
    differences.tsv       every position the two disagree about
    positions.tsv         the full joined table, with --keep-all

As well as the calls, it reports what each tool cost: wall clock, CPU, peak
and mean resident memory, and bytes written. Both read the same filtered BAM,
which is what makes those numbers comparable.

Requires: alnbase, MethylDackel, samtools on PATH, and pyarrow. Memory
sampling reads /proc, so the mean is Linux-only; peak works anywhere.
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil
import subprocess
import sys
import threading
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

# --- the queries ----------------------------------------------------------
#
# Six queries: for each of the three contexts, the cytosine read as C
# (protected) and as T (converted). The tag table is not used by
# `query --parquet`, but it is what makes the same file usable for tagging, so
# the comparison and the real run share one definition.

QUERY_FILE = """\
[query.TG]
read = "T~"
refr = "CG"

[query.CG]
read = "C~"
refr = "CG"

[query.THH]
read = "T~~"
refr = "CHH"

[query.CHH]
read = "C~~"
refr = "CHH"

[query.THG]
read = "T~~"
refr = "CHG"

[query.CHG]
read = "C~~"
refr = "CHG"

[tag.XM.bases]
fill = "."
z = "TG"
Z = "CG"
h = "THH"
H = "CHH"
x = "THG"
X = "CHG"
"""

# query name -> (context, was the cytosine protected)
CALL = {
    "CG": ("CpG", True),
    "TG": ("CpG", False),
    "CHG": ("CHG", True),
    "THG": ("CHG", False),
    "CHH": ("CHH", True),
    "THH": ("CHH", False),
}

# MethylDackel names its outputs after the context.
MD_SUFFIX = {"CpG": "CpG", "CHG": "CHG", "CHH": "CHH"}


# --- measuring what a run costs ------------------------------------------


@dataclass
class Cost:
    """What one command cost. Zeroes mean "not measured", not "free"."""

    wall: float = 0.0
    cpu: float = 0.0
    peak_rss: int = 0  # bytes
    mean_rss: float = 0.0  # bytes, time-weighted over the samples
    samples: int = 0
    outputs: int = 0  # bytes written, filled in by the caller

    def add(self, other: "Cost") -> "Cost":
        """Total two commands: time and bytes add, peak is the larger, and the
        mean is weighted by how long each ran."""
        wall = self.wall + other.wall
        mean = (
            (self.mean_rss * self.wall + other.mean_rss * other.wall) / wall if wall else 0.0
        )
        return Cost(
            wall=wall,
            cpu=self.cpu + other.cpu,
            peak_rss=max(self.peak_rss, other.peak_rss),
            mean_rss=mean,
            samples=self.samples + other.samples,
            outputs=self.outputs + other.outputs,
        )


def _rss_of_tree(pid: int) -> int:
    """Resident bytes for a process and its descendants, from /proc.

    Both tools are one multi-threaded process, so the tree walk is mostly
    insurance; samtools invoked through a shell would not be.
    """
    total, seen, stack = 0, set(), [pid]
    while stack:
        p = stack.pop()
        if p in seen:
            continue
        seen.add(p)
        try:
            with open(f"/proc/{p}/status") as fh:
                for line in fh:
                    if line.startswith("VmRSS:"):
                        total += int(line.split()[1]) * 1024
                        break
            for task in os.listdir(f"/proc/{p}/task"):
                with open(f"/proc/{p}/task/{task}/children") as fh:
                    stack.extend(int(c) for c in fh.read().split())
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue  # it exited between listing and reading
    return total


def run(cmd: list[str], measure: bool = False, interval: float = 0.05, **kw) -> Cost:
    """Run a command, optionally measuring what it cost.

    Peak RSS comes from the kernel via `wait4`, which is exact. The mean comes
    from sampling `/proc`, which is not: a spike between samples is missed.
    The two answer different questions -- "will this fit in memory" and "how
    much was it using most of the time" -- so both are reported.
    """
    argv = [str(c) for c in cmd]
    print("  $ " + " ".join(argv), file=sys.stderr)
    if not measure:
        subprocess.run(argv, check=True, **kw)
        return Cost()

    start = time.monotonic()
    proc = subprocess.Popen(argv, **kw)
    total, samples, stop = 0.0, 0, threading.Event()

    def sample() -> None:
        nonlocal total, samples
        while not stop.wait(interval):
            rss = _rss_of_tree(proc.pid)
            if rss:
                total += rss
                samples += 1

    watcher = threading.Thread(target=sample, daemon=True)
    watcher.start()
    try:
        _, status, usage = os.wait4(proc.pid, 0)
    finally:
        stop.set()
        watcher.join(timeout=1.0)
    proc.returncode = os.waitstatus_to_exitcode(status)
    wall = time.monotonic() - start
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, argv)

    return Cost(
        wall=wall,
        cpu=usage.ru_utime + usage.ru_stime,
        # ru_maxrss is kilobytes on Linux, bytes on macOS.
        peak_rss=usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024),
        mean_rss=total / samples if samples else 0.0,
        samples=samples,
    )


def dir_bytes(paths: list[Path]) -> int:
    return sum(p.stat().st_size for p in paths if p.exists())


def human(n: float) -> str:
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(n) < 1024 or unit == "TiB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{int(n)} B"
        n /= 1024
    return ""


def hms(sec: float) -> str:
    if sec < 60:
        return f"{sec:.1f}s"
    m, s = divmod(int(sec), 60)
    h, m = divmod(m, 60)
    return f"{h}h{m:02d}m{s:02d}s" if h else f"{m}m{s:02d}s"


def missing_from_index(bam_names: list[str], aref_names: list[str]) -> str | None:
    """Complain if the index does not cover the BAM.

    The failure this catches is the expensive one: the tools run for an hour
    and the comparison then covers half the genome. Only the index is checked
    -- what the FASTA does or does not contain is MethylDackel's business, and
    a FASTA that is deliberately a subset is a legitimate thing to hand it.
    """
    missing = set(bam_names) - set(aref_names)
    if not missing:
        return None
    shown = ", ".join(sorted(missing)[:5]) + (", ..." if len(missing) > 5 else "")
    return (
        f"the alnbase index does not cover the BAM:\n"
        f"  {len(missing)} of {len(set(bam_names))} BAM contigs are missing: {shown}\n"
        "  Build the index from the FASTA the BAM was aligned against:\n"
        "    alnbase index REF.fa REF.aref\n"
        "  Or pass --permissive-aln to let alnbase skip those records and count them."
    )


def ensure_fresh_fai(fasta: Path) -> None:
    """Rebuild the `.fai` if it is missing or older than the FASTA.

    A stale index is worse than none: `samtools faidx` keeps working and
    returns sequence from wrong offsets, so every position looks like a
    disagreement and the reference looks wrong. This costs seconds and
    removes a whole class of wrong answers.
    """
    fai = Path(str(fasta) + ".fai")
    if not fai.exists():
        print(f"  {fai.name} is missing; building it", file=sys.stderr)
    elif fai.stat().st_mtime < fasta.stat().st_mtime:
        print(
            f"  {fai.name} is older than {fasta.name}, so it may describe a different\n"
            "  file; rebuilding it",
            file=sys.stderr,
        )
        fai.unlink()
    else:
        # Cheap sanity check: the line width the .fai claims against the
        # file's own second line. A mismatch means they are different files.
        claimed = None
        with open(fai) as fh:
            first = fh.readline().split("\t")
            if len(first) >= 5:
                claimed = int(first[3])
        actual = None
        with open(fasta, "rb") as fh:
            fh.readline()  # the header
            line = fh.readline().rstrip(b"\r\n")
            actual = len(line)
        if claimed is not None and actual and claimed != actual:
            print(
                f"  {fai.name} says the lines are {claimed} bases but {fasta.name}'s are "
                f"{actual}:\n  they are different files. Rebuilding it",
                file=sys.stderr,
            )
            fai.unlink()
        else:
            return
    run(["samtools", "faidx", fasta])


def sort_order(bam: Path) -> str:
    """The `@HD SO:` of a BAM, or "unknown" if it does not say."""
    out = subprocess.run(
        ["samtools", "view", "-H", str(bam)], capture_output=True, text=True, check=True
    )
    for line in out.stdout.splitlines():
        if line.startswith("@HD"):
            for f in line.split("\t"):
                if f.startswith("SO:"):
                    return f[3:]
    return "unknown"


def ensure_sorted_and_indexed(bam: Path, need_index: bool) -> None:
    """Both tools need a coordinate-sorted BAM, and `--region` needs an index.

    Sorting is the caller's decision -- it rewrites their file and can take a
    while -- so that is an error with the command in it. An index is cheap,
    derived, and throwaway, so it is just made.
    """
    order = sort_order(bam)
    if order != "coordinate":
        sys.exit(
            f"{bam} is {'not marked as sorted' if order == 'unknown' else f'sorted by {order}'}.\n"
            "  MethylDackel piles up by position, so it needs coordinate order:\n"
            f"    samtools sort -o sorted.bam {bam} && samtools index sorted.bam\n"
            "  (If it really is coordinate-sorted and the header just does not say,\n"
            "  `samtools reheader` can fix the @HD line.)"
        )
    if not need_index:
        return
    if any(Path(str(bam) + ext).exists() for ext in (".bai", ".csi")):
        return
    print(f"  note: {bam} is not indexed; --region needs one", file=sys.stderr)
    run(["samtools", "index", bam])


def contig_names(bam: Path) -> list[str]:
    out = subprocess.run(["samtools", "view", "-H", str(bam)], capture_output=True, text=True, check=True)
    names = []
    for line in out.stdout.splitlines():
        if line.startswith("@SQ"):
            for f in line.split("\t"):
                if f.startswith("SN:"):
                    names.append(f[3:])
    return names


def aref_names(alnbase: str, aref: Path) -> list[str]:
    out = subprocess.run([alnbase, "info", "--all", str(aref)], capture_output=True, text=True, check=True)
    names, in_list = [], False
    for line in out.stdout.splitlines():
        if line.strip() == "contigs":
            in_list = True
            continue
        if in_list and line.startswith("  "):
            names.append(line.split()[0])
    return names


def need(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        sys.exit(f"{tool} is not on PATH")
    return path


# --- step 1: one input for both tools -------------------------------------


def filter_bam(bam: Path, out: Path, exclude: str, min_mapq: int, region: str | None) -> Path:
    """Apply the read-level filters once, so both tools see the same reads.

    MethylDackel applies its own by default; alnbase applies none. Rather than
    try to make the two agree flag for flag, the filtering happens here and
    both tools are then told not to filter.
    """
    cmd = ["samtools", "view", "-b", "-F", exclude, "-q", str(min_mapq), "-o", out, bam]
    if region:
        cmd.append(region)
    run(cmd)
    run(["samtools", "index", out])
    n = subprocess.run(
        ["samtools", "view", "-c", str(out)], capture_output=True, text=True, check=True
    ).stdout.strip()
    if n == "0":
        sys.exit(
            f"no records survived the filter{f' in {region}' if region else ''}.\n"
            "  Check the region name against the BAM's contigs, and that "
            "--min-mapq and --exclude are not excluding everything."
        )
    print(f"  {n} records")
    return out


# --- step 2: alnbase ------------------------------------------------------


def run_alnbase(
    alnbase: str,
    bam: Path,
    aref: Path,
    qfile: Path,
    outdir: Path,
    threads: int,
    via_tags: bool,
    permissive: bool,
) -> tuple[list[Path], Cost]:
    base = [alnbase] + (["--permissive"] if permissive else [])
    prefix = outdir / "alnbase.parquet"
    for stale in glob.glob(str(outdir / "alnbase_*.parquet")):
        os.remove(stale)

    outputs: list[Path] = []
    if via_tags:
        # The longer path: tag the BAM, then read the tag back. Exercises the
        # tag round trip as well as the walk, which is what a real run does,
        # and costs a whole extra BAM on disk -- which is the point of
        # measuring it separately.
        tagged = outdir / "tagged.bam"
        cost = run(base + ["query", "--query-file", qfile, "-@", threads, bam, aref, tagged],
                   measure=True)
        cost = cost.add(
            run(base + ["extract", "-@", threads, "-f", "ref_name,strand,mapq",
                 "--only-hits", tagged, prefix], measure=True)
        )
        outputs.append(tagged)
    else:
        cost = run([alnbase, "query", "--parquet", "--query-file", qfile, "-@", threads,
                    "-f", "ref_name,strand,mapq", "--only-hits", bam, aref, prefix],
                   measure=True)

    files = [Path(f) for f in sorted(glob.glob(str(outdir / "alnbase_*.parquet")))]
    if not files:
        sys.exit("alnbase wrote no parquet files")
    cost.outputs = dir_bytes(outputs + files)
    return files, cost


def load_alnbase(files: list[Path], min_qual: int) -> tuple[dict, dict]:
    """(chrom, pos, context) -> [meth, unmeth], from per-read calls."""
    import pyarrow.parquet as pq

    counts: dict = defaultdict(lambda: [0, 0])
    strands: dict = defaultdict(set)
    by_strand: dict = defaultdict(int)
    unknown: set = set()
    rows = 0
    for f in files:
        have = pq.read_schema(f).names
        cols = ["ref_name", "refr_pos", "name", "qual"] + (["strand"] if "strand" in have else [])
        d = pq.read_table(f, columns=cols).to_pydict()
        d.setdefault("strand", [None] * len(d["name"]))
        for chrom, pos, name, qual, strand in zip(
            d["ref_name"], d["refr_pos"], d["name"], d["qual"], d["strand"]
        ):
            # A record that matched nothing still has a row under some
            # configurations; it has no query name and no position.
            if name is None or pos is None or chrom is None:
                continue
            if qual is not None and qual < min_qual:
                continue
            call = CALL.get(name)
            if call is None:
                unknown.add(name)
                continue
            context, methylated = call
            counts[(chrom, int(pos), context)][0 if methylated else 1] += 1
            if strand:
                # Which strand's cytosine this is. A position is a cytosine on
                # exactly one strand, so this is what tells a reference
                # disagreement apart from a strand disagreement.
                strands[(chrom, int(pos))].add(strand)
                by_strand[strand] += 1
            rows += 1
    if unknown:
        print(f"note: ignored query names not in the six contexts: {sorted(unknown)}")
    print(f"  alnbase: {rows} calls over {len(counts)} positions")
    if by_strand:
        print("    by strand: " + ", ".join(f"{v} on {k}" for k, v in sorted(by_strand.items())))
    return counts, strands


# --- step 3: MethylDackel -------------------------------------------------


def run_methyldackel(
    md: str,
    bam: Path,
    fasta: Path,
    outdir: Path,
    threads: int,
    min_qual: int,
    extra: list[str],
) -> tuple[dict, Cost]:
    prefix = outdir / "md"
    cost = run([md, "extract", "--CHH", "--CHG", "-q", "0", "-p", str(min_qual),
                "-@", threads, *extra, "-o", prefix, fasta, bam], measure=True)

    counts: dict = defaultdict(lambda: [0, 0])
    for context, suffix in MD_SUFFIX.items():
        path = Path(f"{prefix}_{suffix}.bedGraph")
        if not path.exists():
            sys.exit(f"MethylDackel wrote no {path.name}")
        with open(path) as fh:
            for line in fh:
                if line.startswith(("track", "#")) or not line.strip():
                    continue
                f = line.split()
                # chrom  start  end  percent  nMeth  nUnmeth
                counts[(f[0], int(f[1]), context)] = [int(f[4]), int(f[5])]
    cost.outputs = dir_bytes([Path(f"{prefix}_{s}.bedGraph") for s in MD_SUFFIX.values()])
    total = sum(m + u for m, u in counts.values())
    print(f"  MethylDackel: {total} calls over {len(counts)} positions")
    return counts, cost


# --- step 4: compare ------------------------------------------------------


def compare(
    a: dict,
    m: dict,
    strands: dict,
    outdir: Path,
    keep_all: bool,
    examples: int,
    fasta: Path | None = None,
) -> int:
    """Join the two call sets on (chrom, position) and classify each position.

    Keyed on position, not on (position, context): the two tools can agree
    that a cytosine is methylated and disagree about whether its context is
    CHG or CHH, and that is one disagreement worth seeing, not two phantom
    positions with nothing in common.
    """
    by_pos: dict = defaultdict(lambda: {"aln": None, "md": None})
    for (chrom, pos, context), (meth, unmeth) in a.items():
        by_pos[(chrom, pos)]["aln"] = (context, meth, unmeth)
    for (chrom, pos, context), (meth, unmeth) in m.items():
        by_pos[(chrom, pos)]["md"] = (context, meth, unmeth)

    per_context: dict = defaultdict(lambda: defaultdict(int))
    diffs = []
    for (chrom, pos), sides in by_pos.items():
        aln, md = sides["aln"], sides["md"]
        # Attribute a position to whichever context is known; when they
        # disagree, to alnbase's, so the row appears once.
        context = (aln or md)[0]
        ac, am, au = aln if aln else ("-", 0, 0)
        mc, mm, mu = md if md else ("-", 0, 0)

        if aln is None:
            kind = "only_methyldackel"
        elif md is None:
            kind = "only_alnbase"
        elif ac != mc:
            kind = "context_differs"
        elif (am, au) == (mm, mu):
            kind = "identical"
        else:
            # Same position and context, different counts. Whether the
            # *direction* of the call agrees is usually the more interesting
            # question than whether the depths match exactly.
            same_call = (am > au) == (mm > mu) and (am == au) == (mm == mu)
            kind = "counts_differ" if same_call else "call_differs"

        per_context[context][kind] += 1
        if kind != "identical":
            # "none" rather than "-": a dash reads as the minus strand, which
            # made every only_methyldackel row look like a minus-strand
            # finding when it only meant alnbase had nothing there.
            seen = ",".join(sorted(strands.get((chrom, pos), []))) or "none"
            diffs.append((chrom, pos, ac, am, au, mc, mm, mu, kind, seen))

    header = (
        "chrom\tpos\taln_context\taln_meth\taln_unmeth\tmd_context\tmd_meth\tmd_unmeth"
        "\tkind\taln_strand\n"
    )
    diffs.sort(key=lambda d: (d[0], d[1]))
    with open(outdir / "differences.tsv", "w") as fh:
        fh.write(header)
        for d in diffs:
            fh.write("\t".join(str(x) for x in d) + "\n")

    if keep_all:
        with open(outdir / "positions.tsv", "w") as fh:
            fh.write(header.replace("\tkind\n", "\n"))
            for (chrom, pos), sides in sorted(by_pos.items()):
                ac, am, au = sides["aln"] if sides["aln"] else ("-", 0, 0)
                mc, mm, mu = sides["md"] if sides["md"] else ("-", 0, 0)
                fh.write(f"{chrom}\t{pos}\t{ac}\t{am}\t{au}\t{mc}\t{mm}\t{mu}\n")

    kinds = [
        "identical",
        "counts_differ",
        "call_differs",
        "context_differs",
        "only_alnbase",
        "only_methyldackel",
    ]
    contexts = sorted(per_context)
    width = max(len(k) for k in kinds)
    # Which strand alnbase saw a disagreeing position on. A position is a
    # cytosine on exactly one strand, so if the disagreements pile up on one
    # strand the two tools disagree about strand assignment; if they are even,
    # look at the reference instead.
    shapes: dict = defaultdict(int)
    for d in diffs:
        shapes[(d[8], d[9])] += 1
    if any(strand != "-" for _, strand in shapes):
        print("\ndisagreements by kind, and the strand alnbase saw")
        for (kind, strand), n in sorted(shapes.items(), key=lambda kv: -kv[1])[:12]:
            print(f"  {n:>10}  {kind:<18} strand {strand}")

    print("\nper position, by context")
    print(f"  {'':<{width}}  " + "".join(f"{c:>12}" for c in contexts) + f"{'total':>12}")
    for kind in kinds:
        total = sum(per_context[c][kind] for c in contexts)
        row = "".join(f"{per_context[c][kind]:>12}" for c in contexts)
        print(f"  {kind:<{width}}  {row}{total:>12}")

    total = len(by_pos)
    agree = sum(per_context[c]["identical"] for c in contexts)
    print(f"\n  {agree}/{total} positions identical", end="")
    print(f" ({100 * agree / total:.4f}%)" if total else "")
    shared = total - sum(
        per_context[c][k] for c in contexts for k in ("only_alnbase", "only_methyldackel")
    )
    if shared:
        print(f"  of the {shared} both tools call, {100 * agree / shared:.4f}% agree exactly")

    if diffs and examples and fasta is not None:
        reference_check(fasta, diffs, min(8, examples))

    if diffs and examples:
        print(f"\nfirst {min(examples, len(diffs))} differences (also in differences.tsv)")
        print("  " + header.strip().replace("\t", "  "))
        for d in diffs[:examples]:
            print("  " + "  ".join(str(x) for x in d))

    return len(diffs)


def context_from_reference(window: str, at: int) -> str:
    """The context of the cytosine at `at` in `window`, from the reference.

    A position is a cytosine on one strand only: C means the plus strand and
    the context runs rightwards; G means the minus strand, and the context is
    the complement of the two bases to the *left*, in that order. Whichever
    tool disagrees with this is the one that is wrong -- assuming the FASTA is
    the reference both tools were given, which is worth checking separately.
    """
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    base = window[at] if at < len(window) else "?"
    if base == "C":
        nxt = window[at + 1 : at + 3]
        if len(nxt) < 2:
            return "?"
        if nxt[0] == "G":
            return "CpG"
        return "CHG" if nxt[1] == "G" else "CHH"
    if base == "G":
        prev = window[max(0, at - 2) : at]
        if len(prev) < 2:
            return "?"
        b1, b2 = comp.get(prev[1], "?"), comp.get(prev[0], "?")
        if b1 == "G":
            return "CpG"
        return "CHG" if b2 == "G" else "CHH"
    return "not a cytosine"


def reference_check(fasta: Path, diffs: list, limit: int) -> None:
    """For a few disagreeing positions, print the reference itself.

    A position is a cytosine on exactly one strand, so its context follows
    from the reference alone. If the two tools disagree about the context,
    either they read different bases there or one of them has the strand
    wrong -- and the FASTA settles which. This asks samtools, so it is
    MethylDackel's view of the reference, not alnbase's.
    """
    interesting = [d for d in diffs if d[8] == "context_differs"][:limit]
    if not interesting:
        interesting = diffs[:limit]
    if not interesting:
        return
    # bedGraph and alnbase positions are 0-based; `samtools faidx` is 1-based
    # and inclusive. So 0-based p is 1-based p+1, and a window of 0-based
    # [p-2, p+2] is 1-based [p-1, p+3].
    regions = [f"{d[0]}:{max(1, d[1] - 1)}-{d[1] + 3}" for d in interesting]
    out = subprocess.run(
        ["samtools", "faidx", str(fasta), *regions], capture_output=True, text=True
    )
    if out.returncode != 0:
        print("\ncould not read the FASTA to check:", out.stderr.strip()[:200])
        return

    seqs, name = {}, None
    for line in out.stdout.splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            seqs[name] = ""
        elif name:
            seqs[name] += line.strip()

    print("\nthe reference at the first few disagreements (from the FASTA)")
    print("  the base at pos decides the strand: C is a plus-strand cytosine,")
    print("  G is a minus-strand one, anything else is neither")
    print(f"  {'position':<22} {'pos-2..pos+2':<14} {'base':<5} {'strand it can be':<18} "
          f"{'aln':<6} {'md':<6} {'from the FASTA'}")
    for d, region in zip(interesting, regions):
        seq = seqs.get(region, "")
        # The window starts at 0-based pos-2, so pos itself is two in --
        # unless it was clipped at the start of the contig.
        offset = 2 if d[1] >= 2 else d[1]
        base = seq[offset].upper() if len(seq) > offset else "?"
        can = {"C": "+", "G": "-"}.get(base, "neither")
        truth = context_from_reference(seq.upper(), offset)
        print(f"  {d[0] + ':' + str(d[1]):<22} {seq.upper():<14} {base:<5} {can:<18} "
              f"{d[2]:<6} {d[5]:<6} {truth}")
    print(
        "  If `base` is neither C nor G, the FASTA and the BAM's reference differ.\n"
        "  If it is C or G but a tool called the other strand's context, that tool\n"
        "  has the strand wrong."
    )


def report_cost(aln: Cost, md: Cost, threads: int, via_tags: bool, shared: Path) -> None:
    """What each tool cost. Comparable only because both read the same file.

    Peak is what has to fit in memory; mean is what it used most of the time,
    and the gap between them is usually a load-the-index spike. Output size
    counts only what each tool wrote itself -- the filtered BAM they share is
    reported separately, since it is not attributable to either.
    """
    rows = [
        ("alnbase" + (" (tag, then extract)" if via_tags else " (--parquet)"), aln),
        ("MethylDackel", md),
    ]
    name_w = max(len(n) for n, _ in rows)
    print(f"\nresources, at -@ {threads}")
    print(f"  {'':<{name_w}}  {'wall':>9}  {'cpu':>9}  {'peak rss':>10}  {'mean rss':>10}  {'output':>10}")
    for name, c in rows:
        mean = human(c.mean_rss) if c.samples else "n/a"
        print(
            f"  {name:<{name_w}}  {hms(c.wall):>9}  {hms(c.cpu):>9}  "
            f"{human(c.peak_rss):>10}  {mean:>10}  {human(c.outputs):>10}"
        )
    if shared.exists():
        print(f"  (both read {shared.name}, {human(shared.stat().st_size)}, not counted above)")
    if aln.samples < 3 or md.samples < 3:
        print("  note: too few memory samples to trust the mean; use a larger input")


def main() -> int:
    p = argparse.ArgumentParser(
        description="Compare alnbase and MethylDackel per position.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("bam", type=Path, nargs="?", help="coordinate-sorted, indexed BAM")
    p.add_argument(
        "fasta", type=Path, nargs="?",
        help="reference FASTA for MethylDackel; faidx'd if it is not already. "
             "Its contents are not checked against anything",
    )
    p.add_argument("aref", type=Path, nargs="?", help="the same reference as an alnbase index")
    p.add_argument("-o", "--outdir", type=Path, default=Path("md_compare"))
    p.add_argument("-@", "--threads", type=int, default=4)
    p.add_argument("--alnbase", default="alnbase", help="path to the alnbase binary")
    p.add_argument("--methyldackel", default="MethylDackel")
    p.add_argument(
        "--exclude", default="0xF04",
        help="samtools -F for the shared pre-filter: unmapped, secondary, "
             "QC-fail, duplicate, supplementary (default 0xF04)",
    )
    p.add_argument("--min-mapq", type=int, default=10, help="applied once, to both (default 10)")
    p.add_argument(
        "--region", metavar="CHR[:FROM-TO]",
        help="restrict to one region. Start here: a whole mammalian genome is "
             "hours per run, one chromosome is minutes, and every difference "
             "worth understanding shows up on one chromosome too",
    )
    p.add_argument(
        "--min-qual", type=int, default=5,
        help="minimum base quality, applied to both sides (default 5)",
    )
    p.add_argument(
        "--via-tags", action="store_true",
        help="tag the BAM and read the tag back, instead of writing rows "
             "directly; exercises the tag round trip too",
    )
    p.add_argument(
        "--permissive-aln", action="store_true",
        help="pass --permissive to alnbase, so a record on a contig the index "
             "lacks is skipped and counted rather than stopping the run",
    )
    p.add_argument(
        "--md-pair-defaults", action="store_true",
        help="let MethylDackel keep its default pair filtering, which drops "
             "singleton and discordant pairs. Off by default because a Hi-C "
             "ligation product is discordant by construction, and dropping "
             "those means the two tools read different reads",
    )
    p.add_argument(
        "--md-extra", nargs="*", default=[], metavar="ARG",
        help="extra arguments for MethylDackel extract, e.g. --md-extra --OT 5,5,5,5",
    )
    p.add_argument("--keep-all", action="store_true", help="also write positions.tsv")
    p.add_argument("--examples", type=int, default=20)
    p.add_argument(
        "--skip-run", action="store_true",
        help="reuse the outputs already in --outdir instead of running the tools",
    )
    args = p.parse_args()

    # --skip-run only re-reads what is already in --outdir, so it needs none
    # of the tools and none of the inputs.
    if not args.skip_run:
        for tool in (args.alnbase, args.methyldackel, "samtools"):
            need(tool)
        ensure_fresh_fai(args.fasta)

    if not args.skip_run and not (args.bam and args.fasta and args.aref):
        sys.exit("bam, fasta and aref are required unless --skip-run is given")

    args.outdir.mkdir(parents=True, exist_ok=True)
    qfile = args.outdir / "queries.toml"
    qfile.write_text(QUERY_FILE)

    filtered = args.outdir / "filtered.bam"
    if not args.skip_run:
        print("checking the input:")
        ensure_sorted_and_indexed(args.bam, need_index=args.region is not None)
        problem = missing_from_index(contig_names(args.bam), aref_names(args.alnbase, args.aref))
        if problem:
            sys.exit(problem)
        print("  ok: coordinate-sorted, and the index covers it")

        print("\nfiltering once, so both tools read the same records:")
        filter_bam(args.bam, filtered, args.exclude, args.min_mapq, args.region)
        print("\nalnbase:")
        _, aln_cost = run_alnbase(args.alnbase, filtered, args.aref, qfile, args.outdir,
                                  args.threads, args.via_tags, args.permissive_aln)
        print("\nMethylDackel:")
        md_extra = list(args.md_extra)
        if not args.md_pair_defaults:
            # MethylDackel drops singleton and discordant pairs by default.
            # For a Hi-C library that is most of the data -- a ligation
            # product is *supposed* to look discordant -- so the two tools
            # would be reading different halves of the BAM.
            md_extra = ["--keepDiscordant", "--keepSingleton"] + md_extra
        _, md_cost = run_methyldackel(args.methyldackel, filtered, args.fasta, args.outdir,
                                      args.threads, args.min_qual, md_extra)
        report_cost(aln_cost, md_cost, args.threads, args.via_tags, filtered)

    print("\nloading:")
    aln, strands = load_alnbase(
        [Path(f) for f in sorted(glob.glob(str(args.outdir / "alnbase_*.parquet")))], args.min_qual
    )
    md: dict = defaultdict(lambda: [0, 0])
    for context, suffix in MD_SUFFIX.items():
        path = args.outdir / f"md_{suffix}.bedGraph"
        if not path.exists():
            continue
        with open(path) as fh:
            for line in fh:
                if line.startswith(("track", "#")) or not line.strip():
                    continue
                f = line.split()
                md[(f[0], int(f[1]), context)] = [int(f[4]), int(f[5])]
    print(f"  MethylDackel: {sum(m + u for m, u in md.values())} calls over {len(md)} positions")

    n = compare(aln, md, strands, args.outdir, args.keep_all, args.examples, args.fasta)

    print(
        "\nExpected differences, before treating any of this as a bug:\n"
        "  - Overlapping mates are counted twice by both tools. Run\n"
        "    `alnbase overlap` first if you want them resolved, but then run\n"
        "    MethylDackel on the same resolved BAM.\n"
        "  - MethylDackel trims nothing by default; if you normally pass\n"
        "    --OT/--OB, pass the same trimming here or expect end-of-read\n"
        "    differences.\n"
        "  - A cytosine whose reference context contains an N is called by\n"
        "    neither of these six queries. Bismark would call it u/U.\n"
        "  - alnbase reads reference context past the end of a read, so a\n"
        "    cytosine at the last base still has its context.\n"
        "  Investigate a specific position with:\n"
        "    samtools view filtered.bam CHROM:POS-POS | head\n"
        "    alnbase query --query-file queries.toml --trace-records 1 \\\n"
        "        filtered.bam REF.aref out.bam"
    )
    return 1 if n else 0


if __name__ == "__main__":
    sys.exit(main())
