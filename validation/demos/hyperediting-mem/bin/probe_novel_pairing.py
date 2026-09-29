"""Does protecting the motif at annotated sites recover novel pairings of them?

    probe_novel_pairing.py SIMDIR OUTDIR [--he DIR]

A splice-site file is a list of donor-acceptor *pairs*. Exon skipping and
alternative site usage produce junctions that join two sites both of which are
annotated, in a combination that is not. Those junctions are absent from the
list, so the aligner is back to finding them de novo -- on a three-letter genome
where the motif has been collapsed away.

Protecting the motif marks the sites individually rather than as pairs, so in
principle it should help exactly here, and only here. The four bases it protects
per intron (GT at the start, AG at the end) are all intronic, so a spliced read
never contains them and the protection is close to free: about 0.006% of hg38
against the 6.1% that made global protection collapse.

This measures it. Reads cross a skipping junction that joins the donor of one
annotated intron to the acceptor of another, which is a pairing no annotation
contains, and are scored against four setups: with and without the motif
protected, with and without the canonical pair list. A fifth adds the skipping
pair to the list, as the upper bound.
"""
import argparse
import os
import random
import re
import subprocess

READ_LEN = 100
ANCHORS = [6, 8, 10, 12, 14, 16, 20, 30]
PER_ANCHOR = 400


def read_fasta(path):
    name, chunks = None, []
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
        else:
            chunks.append(line.strip())
    return name, "".join(chunks)


def write_fasta(path, name, seq):
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("simdir")
    ap.add_argument("outdir")
    ap.add_argument("--he", default=None)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--edit-rate", type=float, default=0.35)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    chrom, genome = read_fasta(os.path.join(args.simdir, "genome.fa"))

    annotated = []
    for line in open(os.path.join(args.simdir, "splicesites.txt")):
        f = line.split()
        annotated.append((int(f[1]) + 1, int(f[2])))

    # Skipping junctions: the donor of one annotated intron joined to the
    # acceptor of a later one. Both sites are annotated; the pairing is not.
    novel = [(annotated[0][0], annotated[1][1]), (annotated[1][0], annotated[2][1])]

    # --- references ----------------------------------------------------------
    a2g = genome.replace("A", "G")
    protected = list(a2g)
    for s, e in annotated:
        protected[s:s + 2] = list("GT")
        protected[e - 2:e] = list("AG")
    protected = "".join(protected)
    kept = sum(1 for i in range(len(a2g)) if protected[i] != a2g[i])
    print(f"protected {kept} bases at {len(annotated)} annotated introns "
          f"({kept / len(genome) * 100:.4f}% of the genome), all of them intronic\n")

    refs = {"collapsed": a2g, "motif at annotated sites": protected}
    for tag, seq in refs.items():
        slug = tag.split()[0]
        fa = os.path.join(args.outdir, f"{slug}.fa")
        write_fasta(fa, chrom, seq)
        if not os.path.exists(os.path.join(args.outdir, slug) + ".1.ht2"):
            subprocess.run(["hisat2-build", "-q", "-p", "8", fa,
                            os.path.join(args.outdir, slug)], check=True)

    # --- splice-site files ---------------------------------------------------
    def write_sites(path, pairs):
        with open(path, "w") as f:
            for s, e in pairs:
                f.write(f"{chrom}\t{s - 1}\t{e}\t+\n")

    canon = os.path.join(args.outdir, "canonical.txt")
    withnovel = os.path.join(args.outdir, "with_novel.txt")
    write_sites(canon, annotated)
    write_sites(withnovel, annotated + novel)

    # --- reads across the skipping junctions ---------------------------------
    fq = os.path.join(args.outdir, "reads.fq")
    truth = {}
    with open(fq, "w") as f:
        for anchor in ANCHORS:
            for i in range(PER_ANCHOR):
                s, e = rng.choice(novel)
                if rng.random() < 0.5:
                    left, right = anchor, READ_LEN - anchor
                else:
                    left, right = READ_LEN - anchor, anchor
                seq = genome[s - left:s] + genome[e:e + right]
                a_pos = [k for k, b in enumerate(seq) if b == "A"]
                for k in rng.sample(a_pos, int(len(a_pos) * args.edit_rate)):
                    seq = seq[:k] + "G" + seq[k + 1:]
                name = f"skip{anchor}_{i:04d}"
                truth[name] = (anchor, e - s)
                f.write(f"@{name}\n{seq.replace('A', 'G')}\n+\n{'I' * len(seq)}\n")

    # --- five setups ---------------------------------------------------------
    setups = [
        ("collapsed, no list", "collapsed", None),
        ("collapsed, canonical list", "collapsed", canon),
        ("motif protected, no list", "motif", None),
        ("motif protected, canonical list", "motif", canon),
        ("skipping pair in the list", "collapsed", withnovel),
    ]

    results = {}
    for label, slug, sites in setups:
        cmd = ["hisat2", "-p", "8", "--no-unal", "-x", os.path.join(args.outdir, slug),
               "-U", fq, "-S", "/dev/stdout"]
        if sites:
            cmd += ["--known-splicesite-infile", sites]
        sam = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
        found = {a: 0 for a in ANCHORS}
        wrong = 0
        for line in sam.split("\n"):
            if not line or line.startswith("@"):
                continue
            fld = line.split("\t")
            if fld[0] not in truth:
                continue
            anchor, ilen = truth[fld[0]]
            ns = [int(x[:-1]) for x in re.findall(r"\d+N", fld[5])]
            if ns and ns[0] == ilen:
                found[anchor] += 1
            elif ns:
                wrong += 1
        results[label] = (found, wrong)

    total = {a: sum(1 for _, (aa, _) in truth.items() if aa == a) for a in ANCHORS}
    print(f"{'anchor':>7}" + "".join(f"{a:>7}" for a in ANCHORS) + "   wrong-intron gaps")
    for label, (found, wrong) in results.items():
        row = "".join(f"{found[a] / total[a] * 100:>6.0f}%" for a in ANCHORS)
        print(f"{label:>7}".replace(f"{label:>7}", f"{label:<34}") + row + f"{wrong:>10}")

    out = os.path.join(args.outdir, "novel_pairing.tsv")
    with open(out, "w") as f:
        f.write("setup\tanchor\tfound\ttotal\trate\n")
        for label, (found, _) in results.items():
            for a in ANCHORS:
                f.write(f"{label}\t{a}\t{found[a]}\t{total[a]}\t{found[a] / total[a]:.4f}\n")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
