"""Score a hyper-editing run against the simulator's ground truth.

    score.py TRUTHDIR ARM_NAME=ES_BED [ARM_NAME=ES_BED ...] [--out report.tsv]

Each ES bed is the per-site output the pipeline writes to
`UEdetect.*/PREFIX.ES.bed_files/A2G.bed`: one line per called editing event,
carrying the genomic coordinate and the read it came from.

Three things are measured, and they fail independently:

  recall      of the editing events the simulation created, how many came back
  placement   of the events that came back, how many are on the right base
  bias        whether recall depends on something it should not -- the read's
              class, how heavily it was edited, or where along the read the
              event sits

A site is a (read_id, 1-based genomic position) pair. A called site whose read
is in the truth table but whose position is not one of that read's true sites
is *misplaced*: the pipeline found a real read and put its edit on the wrong
base. That is the failure the flat-window assumption produces, and it is
invisible without ground truth, because the call looks exactly like a good one.
"""
import argparse
import collections
import os
import re

BED_NAME = re.compile(r"^\([-+]\)(.*);(\d+)$")


def read_truth(truthdir):
    sites = collections.defaultdict(set)      # read_id -> {1-based pos}
    offsets = {}                              # (read_id, pos) -> read offset
    with open(os.path.join(truthdir, "truth_sites.tsv")) as f:
        next(f)
        for line in f:
            rid, pos, off, _rb, _gb = line.rstrip("\n").split("\t")
            sites[rid].add(int(pos))
            offsets[(rid, int(pos))] = int(off)

    reads = {}
    with open(os.path.join(truthdir, "truth_reads.tsv")) as f:
        hdr = next(f).rstrip("\n").split("\t")
        for line in f:
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            reads[row["read_id"]] = row
    return sites, offsets, reads


def read_calls(bed):
    """(read_id, 1-based pos) pairs from an ES bed."""
    out = set()
    if not os.path.exists(bed):
        return out
    with open(bed) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            m = BED_NAME.match(p[3])
            rid = m.group(1) if m else p[3]
            out.add((rid, int(p[1]) + 1))
    return out


def edit_bucket(n):
    n = int(n)
    if n <= 8:
        return "05-08"
    if n <= 12:
        return "09-12"
    if n <= 16:
        return "13-16"
    return "17+"


def score(truth_sites, offsets, reads, calls):
    truth_pairs = {(r, p) for r, ps in truth_sites.items() for p in ps}
    tp = calls & truth_pairs
    fn = truth_pairs - calls
    extra = calls - truth_pairs
    misplaced = {(r, p) for (r, p) in extra if r in truth_sites}
    spurious = extra - misplaced

    by_class = collections.Counter()
    by_class_tot = collections.Counter()
    by_level = collections.Counter()
    by_level_tot = collections.Counter()
    by_off = collections.Counter()
    by_off_tot = collections.Counter()

    for pair in truth_pairs:
        rid = pair[0]
        cls = reads[rid]["class"]
        lvl = edit_bucket(reads[rid]["n_edits"])
        ob = offsets[pair] // 10 * 10
        by_class_tot[cls] += 1
        by_level_tot[lvl] += 1
        by_off_tot[ob] += 1
        if pair in tp:
            by_class[cls] += 1
            by_level[lvl] += 1
            by_off[ob] += 1

    called_reads = {r for r, _ in calls}
    truth_reads_edited = {r for r in truth_sites if truth_sites[r]}

    return {
        "sites_true": len(truth_pairs),
        "sites_called": len(calls),
        "tp": len(tp),
        "fn": len(fn),
        "misplaced": len(misplaced),
        "spurious": len(spurious),
        "recall": len(tp) / len(truth_pairs) if truth_pairs else 0.0,
        "placement": len(tp) / len(calls) if calls else 0.0,
        "reads_true": len(truth_reads_edited),
        "reads_called": len(called_reads),
        "read_recall": len(called_reads & truth_reads_edited) / len(truth_reads_edited),
        "by_class": {c: (by_class[c], by_class_tot[c]) for c in by_class_tot},
        "by_level": {c: (by_level[c], by_level_tot[c]) for c in by_level_tot},
        "by_offset": {c: (by_off[c], by_off_tot[c]) for c in by_off_tot},
    }


def spread(d):
    """How uneven recall is across strata: max minus min, in percentage points.
    A pipeline with no stratum-dependent loss scores 0."""
    rates = [got / tot for got, tot in d.values() if tot >= 20]
    return (max(rates) - min(rates)) * 100 if rates else 0.0


def fmt_strata(d, order=None):
    keys = order or sorted(d)
    return "  ".join(
        f"{k}:{(d[k][0] / d[k][1] * 100 if d[k][1] else 0):5.1f}%" for k in keys if k in d
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("truthdir")
    ap.add_argument("arms", nargs="+", help="NAME=path/to/A2G.bed")
    ap.add_argument("--out")
    args = ap.parse_args()

    truth_sites, offsets, reads = read_truth(args.truthdir)
    results = {}
    for spec in args.arms:
        name, _, path = spec.partition("=")
        results[name] = score(truth_sites, offsets, reads, read_calls(path))

    classes = ["flat_edited", "junction_edited", "del_edited", "ins_edited",
               "clip_edited"]
    levels = ["05-08", "09-12", "13-16", "17+"]

    w = max(len(n) for n in results) + 2
    print(f"\n{'arm':<{w}}{'sites':>8}{'found':>8}{'recall':>9}{'misplc':>8}{'spur':>7}{'placed':>9}{'reads':>8}")
    for name, r in results.items():
        print(f"{name:<{w}}{r['sites_true']:>8}{r['tp']:>8}{r['recall'] * 100:>8.1f}%"
              f"{r['misplaced']:>8}{r['spurious']:>7}{r['placement'] * 100:>8.1f}%"
              f"{r['read_recall'] * 100:>7.1f}%")

    print("\nrecall by read class")
    for name, r in results.items():
        print(f"  {name:<{w}}{fmt_strata(r['by_class'], classes)}")
    print(f"  {'spread':<{w}}" + "  ".join(f"{spread(r['by_class']):.1f}pp ({n})"
                                           for n, r in results.items()))

    print("\nrecall by edits per read")
    for name, r in results.items():
        print(f"  {name:<{w}}{fmt_strata(r['by_level'], levels)}")
    print(f"  {'spread':<{w}}" + "  ".join(f"{spread(r['by_level']):.1f}pp ({n})"
                                           for n, r in results.items()))

    print("\nrecall by position along the read (10 bp bins)")
    offs = sorted(next(iter(results.values()))["by_offset"])
    for name, r in results.items():
        print(f"  {name:<{w}}" + " ".join(
            f"{(r['by_offset'][o][0] / r['by_offset'][o][1] * 100 if r['by_offset'][o][1] else 0):4.0f}"
            for o in offs))
    print(f"  {'bins':<{w}}" + " ".join(f"{o:>4}" for o in offs))
    print(f"  {'spread':<{w}}" + "  ".join(f"{spread(r['by_offset']):.1f}pp ({n})"
                                           for n, r in results.items()))

    if args.out:
        with open(args.out, "w") as f:
            f.write("arm\tmetric\tstratum\tfound\ttotal\trate\n")
            for name, r in results.items():
                for key in ("by_class", "by_level", "by_offset"):
                    for k, (got, tot) in sorted(r[key].items(), key=lambda x: str(x[0])):
                        f.write(f"{name}\t{key}\t{k}\t{got}\t{tot}\t{got / tot if tot else 0:.6f}\n")
                for k in ("sites_true", "tp", "misplaced", "spurious",
                          "reads_true", "reads_called"):
                    f.write(f"{name}\toverall\t{k}\t{r[k]}\t\t\n")
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
