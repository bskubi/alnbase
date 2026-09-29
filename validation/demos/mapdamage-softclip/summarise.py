"""Print mapDamage's 5' C>T table, pooled over strand and library, beside truth.

Usage: summarise.py <mapDamage-result-dir> [truth.tsv]
"""

import collections
import csv
import sys


def misincorporation(folder):
    agg = collections.defaultdict(lambda: [0, 0, 0])
    with open(folder + "/misincorporation.txt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["End"] != "5p":
                continue
            pos = int(row["Pos"])
            agg[pos][0] += int(row["C"])
            agg[pos][1] += int(row["C>T"])
            agg[pos][2] += int(row["S"])
    return agg


def main(folders, truth_path=None):
    truth = {}
    if truth_path:
        with open(truth_path, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                truth[int(row["Pos"])] = float(row["Freq"])

    tables = [(name, misincorporation(name)) for name in folders]
    header = ["Pos"] + (["truth"] if truth else []) + [n for n, _ in tables]
    print("\t".join("%-10s" % h for h in header))
    for pos in range(1, 6):
        cells = ["%-10d" % pos]
        if truth:
            cells.append("%-10.4f" % truth[pos])
        for _, agg in tables:
            c, t, _ = agg[pos]
            cells.append("%-10.4f" % (t / c if c else 0.0))
        print("\t".join(cells))


if __name__ == "__main__":
    args = sys.argv[1:]
    truth = args.pop() if args[-1].endswith(".tsv") else None
    main(args, truth)
