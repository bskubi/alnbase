"""The 2x2, with every term of the rate shown separately.

    score_xspec.py WORKDIR --species NAME [--out TSV] LABEL [LABEL ...]

The quantity Porath 2014 compares across species is hyper-edited sites per
mapped read, and it is the wrong thing to report alone here, because both axes
of this experiment move its denominator and its numerator in different ways.

  A better assembly places more reads in stage 1. That raises mapped reads --
  the denominator -- so a modern reference can lower sites-per-read while
  finding strictly more sites.

  A more sensitive aligner in the transformed stage raises sites without
  touching stage 1 at all, because stage 1 is held at `bwa aln`/`bwa mem` in
  every cell. It moves the numerator only.

So the two effects can cancel in the ratio while both are real. Every cell
below prints mapped reads, the unmapped pool, hyper-edited reads and sites as
separate columns, and the ratio last, as a derived quantity rather than the
result.

Specificity travels with each number, because a cell that finds more sites by
aligning more badly is not a better cell. Two checks, neither needing a truth
set: the control channels (of the six mismatch classes detect_ue.pl reports
only A2G can carry an A-to-I event, so the other five are a noise floor by
construction) and the ADAR motif (G depleted immediately 5' of the edited
adenosine, enriched immediately 3', against a genome background near 21%).
"""
import argparse
import collections
import os

COMBOS = ["A2G", "A2C", "A2T", "C2A", "G2A", "G2C"]
CONTROLS = [c for c in COMBOS if c != "A2G"]
ARMS = [("aln", "bwa aln"), ("hisat2-union", "HISAT2")]


def bed_lines(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [l.rstrip("\n").split("\t") for l in f if l.strip()]


def counts(path):
    """stage1/counts.tsv, written by run_xspec.sh."""
    out = {}
    if os.path.exists(path):
        for line in open(path):
            k, _, v = line.partition("\t")
            out[k.strip()] = int(v)
    return out


def cell(work, label, arm):
    beds = os.path.join(work, f"run_{label}", f"arm_{arm}", "beds")
    es = {c: bed_lines(os.path.join(beds, f"real.ES.bed_files/{c}.bed")) for c in COMBOS}
    ue = {c: bed_lines(os.path.join(beds, f"real.UE.bed_files/{c}.bed")) for c in COMBOS}
    a2g = es["A2G"]

    up = collections.Counter(r[4][0] for r in a2g if len(r) > 4 and len(r[4]) == 3)
    dn = collections.Counter(r[4][2] for r in a2g if len(r) > 4 and len(r[4]) == 3)
    n_ctx = sum(up.values()) or 1
    all_sites = sum(len(es[c]) for c in COMBOS) or 1

    st = counts(os.path.join(work, f"run_{label}", "stage1", "counts.tsv"))
    return dict(
        assembly=label, aligner=arm,
        mapped_reads=st.get("mapped_reads", 0),
        unmapped_pool=st.get("unmapped_pool", 0),
        ue_a2g=len(ue["A2G"]),
        es_a2g=len(a2g),
        es_ctrl=sum(len(es[c]) for c in CONTROLS),
        a2g_share=len(a2g) / all_sites,
        g_up=up["G"] / n_ctx, g_dn=dn["G"] / n_ctx,
        sites_per_mread=len(a2g) / max(st.get("mapped_reads", 0), 1) * 1e3,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    ap.add_argument("labels", nargs="+")
    ap.add_argument("--species", default="")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    rows = [cell(args.work, lab, arm)
            for lab in args.labels for arm, _ in ARMS
            if os.path.isdir(os.path.join(args.work, f"run_{lab}", f"arm_{arm}"))]
    if not rows:
        raise SystemExit("no finished cells found under " + args.work)

    w = max(len(f"{r['assembly']}/{r['aligner']}") for r in rows) + 2
    print(f"\n{args.species}: hyper-editing across assembly and aligner\n")
    print(f"{'cell':<{w}}{'mapped reads':>14}{'unmapped':>11}{'A2G reads':>11}"
          f"{'A2G sites':>11}{'ctrl sites':>12}{'A2G share':>11}{'sites/Mread':>13}")
    for r in rows:
        print(f"{r['assembly'] + '/' + r['aligner']:<{w}}{r['mapped_reads']:>14,}"
              f"{r['unmapped_pool']:>11,}{r['ue_a2g']:>11,}{r['es_a2g']:>11,}"
              f"{r['es_ctrl']:>12,}{r['a2g_share'] * 100:>10.2f}%"
              f"{r['sites_per_mread']:>13.3f}")

    print(f"\n{'cell':<{w}}{'G at -1':>10}{'G at +1':>10}   (ADAR: low at -1, high at +1; ~21% background)")
    for r in rows:
        print(f"{r['assembly'] + '/' + r['aligner']:<{w}}"
              f"{r['g_up'] * 100:>9.1f}%{r['g_dn'] * 100:>9.1f}%")

    # The contrasts, as ratios of sites rather than of the rate, because the
    # rate's denominator is itself one of the things being manipulated.
    by = {(r["assembly"], r["aligner"]): r for r in rows}
    old, new = args.labels[0], args.labels[-1]

    def ratio(a, b, key="es_a2g"):
        ra, rb = by.get(a), by.get(b)
        if not ra or not rb or not rb[key]:
            return None
        return ra[key] / rb[key]

    print("\ncontrasts (A2G sites, ratio):")
    for name, a, b in [
        ("reference, bwa aln    (C/A)", (new, "aln"), (old, "aln")),
        ("aligner, old assembly (B/A)", (old, "hisat2-union"), (old, "aln")),
        ("aligner, new assembly (D/C)", (new, "hisat2-union"), (new, "aln")),
        ("both                  (D/A)", (new, "hisat2-union"), (old, "aln")),
    ]:
        r = ratio(a, b)
        print(f"  {name}  {'-' if r is None else format(r, '.2f') + 'x'}")

    # Whether the two fixes are independent is the question cell D exists to
    # answer, so it is stated rather than left to the reader to multiply.
    if len(args.labels) > 1 and all(by.get((l, a)) for l in (old, new) for a, _ in ARMS):
        ca = ratio((new, "aln"), (old, "aln"))
        ba = ratio((old, "hisat2-union"), (old, "aln"))
        da = ratio((new, "hisat2-union"), (old, "aln"))
        if None not in (ca, ba, da):
            expected = ca * ba
            gap = da / expected if expected else 0.0
            if gap < 0.9:
                verdict = ("short of it, so the better assembly already recovers "
                           "part of what the aligner recovers -- the two overlap")
            elif gap > 1.1:
                verdict = ("above it, so the two compound: reads the aligner can "
                           "only place once the assembly is contiguous")
            else:
                verdict = ("as expected, so the two fixes are independent and "
                           "address different reads")
            print(f"\n  independent-effects check: {ca:.2f}x (reference) x {ba:.2f}x "
                  f"(aligner) = {expected:.2f}x expected,\n  observed D/A is {da:.2f}x -- "
                  f"{verdict}.")

    if args.out:
        keys = list(rows[0])
        with open(args.out, "w") as f:
            f.write("\t".join(keys) + "\n")
            for r in rows:
                f.write("\t".join(str(r[k]) for k in keys) + "\n")
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
