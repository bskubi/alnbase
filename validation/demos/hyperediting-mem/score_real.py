"""Score the real-data arms against the pipeline's own controls.

    score_real.py RUNDIR [--out TSV]

There is no ground truth for a real library, so nothing here is recall. What
takes its place is the pipeline's own design: `detect_ue.pl` reports six
mismatch classes, and only A2G can carry an A-to-I event. A2C, A2T, C2A, G2A
and G2C are noise channels by construction -- no biology produces clusters of
them -- so their yield is the floor that A2G has to clear. An arm that finds
more A2G and the same control yield has found more editing; an arm that raises
both has found more noise.

Three further checks, all computable from what the pipeline already writes, so
none of them needs an external truth set:

  the ADAR signature   ADAR is depleted for G immediately 5' of the edited
                       adenosine and enriched for G immediately 3'. The ES bed
                       carries the trinucleotide context in column five, so the
                       signature can be read straight off. Real editing shows
                       it; misalignment does not, which makes it a specificity
                       measure independent of the control channels.

  position along read  Where in the read the edits sit. A caller that gets the
                       read-to-reference pairing wrong on clipped or spliced
                       reads produces edits piled at the ends.

  junction proximity   Distance from each called site to the nearest junction
                       in the union list. Comparable between arms on one subset,
                       not between subsets: a denser junction list makes "within
                       10 bp of a junction" a bigger target.

One column is not about editing at all. `spliced` counts the alignments carrying
an N in their CIGAR that reached the edit caller. It is there because the
junction-list arms are identical in every other respect, so it is the only place
their difference is visible directly -- and putting it beside the called-site
count is the point: the list can move the alignments a great deal and move the
output not at all.
"""
import argparse
import bisect
import collections
import os
import re
import subprocess

COMBOS = ["A2G", "A2C", "A2T", "C2A", "G2A", "G2C"]
CONTROLS = [c for c in COMBOS if c != "A2G"]


def bed_lines(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [l.rstrip("\n").split("\t") for l in f if l.strip()]


def spliced_alignments(transrun):
    """Alignments with an N in the CIGAR across an arm's twelve filtered BAMs."""
    d = os.path.join(transrun, "bamFiles")
    if not os.path.isdir(d):
        return None
    n = 0
    for b in sorted(os.listdir(d)):
        if not b.endswith(".bam") or "FilterReads" not in b:
            continue
        out = subprocess.run(["samtools", "view", os.path.join(d, b)],
                             capture_output=True, text=True)
        n += sum(1 for l in out.stdout.split("\n")
                 if l and "N" in l.split("\t")[5])
    return n


def junction_sites(path):
    """Every donor and acceptor coordinate in a HISAT2 splice-site file, sorted.
    The file gives the last base of the left exon and the first base of the
    right exon, both 0-based; the two boundaries are what a site can be near."""
    out = []
    for line in open(path):
        f = line.split()
        if len(f) >= 3:
            out += [int(f[1]), int(f[2])]
    return sorted(set(out))


def nearest(sorted_xs, x):
    i = bisect.bisect_left(sorted_xs, x)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(sorted_xs):
            d = abs(sorted_xs[j] - x)
            best = d if best is None else min(best, d)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    arms = sorted(d[4:] for d in os.listdir(args.rundir) if d.startswith("arm_"))
    junc_path = os.path.join(args.rundir, "lists", "union.txt")
    juncs = junction_sites(junc_path) if os.path.exists(junc_path) else []

    rows = []
    for arm in arms:
        beds = os.path.join(args.rundir, f"arm_{arm}", "beds")
        es = {c: bed_lines(os.path.join(beds, f"real.ES.bed_files/{c}.bed")) for c in COMBOS}
        ue = {c: bed_lines(os.path.join(beds, f"real.UE.bed_files/{c}.bed")) for c in COMBOS}

        a2g = es["A2G"]
        ctrl_sites = sum(len(es[c]) for c in CONTROLS)
        ctrl_reads = sum(len(ue[c]) for c in CONTROLS)

        # The ADAR signature, off the trinucleotide in column five.
        up = collections.Counter(r[4][0] for r in a2g if len(r) > 4 and len(r[4]) == 3)
        dn = collections.Counter(r[4][2] for r in a2g if len(r) > 4 and len(r[4]) == 3)
        n_ctx = sum(up.values()) or 1

        # Where along the read. Column four is "(strand)name;offset".
        offs = [int(m.group(1)) for r in a2g
                if (m := re.search(r";(\d+)$", r[3]))]

        spliced = spliced_alignments(os.path.join(args.rundir, "stage2", arm))

        near = [nearest(juncs, int(r[1])) for r in a2g] if juncs else []
        near10 = sum(1 for d in near if d is not None and d <= 10)

        rows.append(dict(
            arm=arm,
            ue_a2g=len(ue["A2G"]), ue_ctrl=ctrl_reads,
            es_a2g=len(a2g), es_ctrl=ctrl_sites, spliced=spliced,
            per_read=len(a2g) / max(len(ue["A2G"]), 1),
            g_up=up["G"] / n_ctx, g_dn=dn["G"] / n_ctx,
            mean_off=sum(offs) / len(offs) if offs else 0.0,
            near_junction=near10 / max(len(a2g), 1),
        ))

    w = max(len(r["arm"]) for r in rows) + 2
    print(f"{'arm':<{w}}{'A2G reads':>11}{'ctrl reads':>12}{'A2G sites':>11}"
          f"{'ctrl sites':>12}{'sites/read':>12}{'signal:noise':>14}{'spliced':>10}")
    for r in rows:
        sn = r["es_a2g"] / max(r["es_ctrl"], 1)
        sp = f"{r['spliced']:,}" if r["spliced"] is not None else "-"
        print(f"{r['arm']:<{w}}{r['ue_a2g']:>11,}{r['ue_ctrl']:>12,}"
              f"{r['es_a2g']:>11,}{r['es_ctrl']:>12,}{r['per_read']:>12.1f}"
              f"{sn:>14.1f}{sp:>10}")

    print(f"\n{'arm':<{w}}{'G at -1':>10}{'G at +1':>10}{'mean offset':>13}"
          f"{'within 10bp of a junction':>27}")
    print(f"{'(ADAR: low at -1, high at +1)':<{w}}")
    for r in rows:
        print(f"{r['arm']:<{w}}{r['g_up'] * 100:>9.1f}%{r['g_dn'] * 100:>9.1f}%"
              f"{r['mean_off']:>13.1f}{r['near_junction'] * 100:>26.2f}%")

    if args.out:
        keys = list(rows[0])
        with open(args.out, "w") as f:
            f.write("\t".join(keys) + "\n")
            for r in rows:
                f.write("\t".join(str(r[k]) for k in keys) + "\n")
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
