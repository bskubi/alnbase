"""How much does losing GT..AG cost a splice-aware aligner on a three-letter genome?

    probe_splice_motif.py SIMDIR OUTDIR

None of the six three-letter transforms leaves a canonical splice motif intact:
a2g turns the acceptor AG into GG, t2c turns the donor GT into GC, and so on
for the rest. A splice-aware aligner scores a candidate junction partly on
those two dinucleotides, so the question is how much alignment it loses without
them.

The confound is that a transformed genome is also less unique, which hurts
alignment for reasons that have nothing to do with splicing. This probe removes
it with a counterfactual reference: take the transformed genome and write GT..AG
back in at the true intron boundaries, four bases per intron. Alphabet, length,
repeat content and everything the read sees are identical; the only difference
is whether the aligner finds motif support at the gap it is considering. The gap
between the two arms is the motif effect and nothing else.

The variable is anchor length -- how many bases of the read fall on the short
side of the junction. A read that overlaps a junction by forty bases is easy
whatever the motif says. A read that overlaps by six is the case where the
aligner has to decide between a short spurious extension and a real intron, and
that is where motif support should matter.
"""
import argparse
import os
import random
import subprocess
import sys

READ_LEN = 100
ANCHORS = [4, 6, 8, 10, 11, 12, 13, 14, 15, 16, 18, 20, 25, 30, 40]
PER_ANCHOR = 600
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    name, chunks = None, []
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
        else:
            chunks.append(line.strip())
    return name, "".join(chunks)


def read_introns(path):
    """splicesites.txt: chrom, last base of the left exon, first base of the
    right exon, strand -- both 0-based. The intron is the half-open span
    between them."""
    out = []
    for line in open(path):
        f = line.split()
        out.append((int(f[1]) + 1, int(f[2])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("simdir")
    ap.add_argument("outdir")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--edit-rate", type=float, default=0.35)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    chrom, genome = read_fasta(os.path.join(args.simdir, "genome.fa"))
    introns = read_introns(os.path.join(args.simdir, "splicesites.txt"))

    # --- the two references ---------------------------------------------------
    a2g = genome.replace("A", "G")
    restored = list(a2g)
    for s, e in introns:
        restored[s:s + 2] = list("GT")
        restored[e - 2:e] = list("AG")
    restored = "".join(restored)

    refs = {"collapsed": a2g, "motif-restored": restored}
    for tag, seq in refs.items():
        fa = os.path.join(args.outdir, f"{tag}.fa")
        with open(fa, "w") as f:
            f.write(f">{chrom}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
        idx = os.path.join(args.outdir, tag)
        if not os.path.exists(idx + ".1.ht2"):
            subprocess.run(["hisat2-build", "-q", "-p", "8", fa, idx], check=True)

    # --- reads: a fixed anchor on the short side of a junction ----------------
    # Written out already three-letter transformed, which is what the pipeline
    # hands the aligner. Editing is applied first, in four-letter space, so the
    # reads are the same shape as the ones that actually reach this stage.
    fq = os.path.join(args.outdir, "reads.fq")
    truth = {}
    with open(fq, "w") as f:
        for anchor in ANCHORS:
            for i in range(PER_ANCHOR):
                s, e = rng.choice(introns)
                left_first = rng.random() < 0.5
                if left_first:
                    left, right = anchor, READ_LEN - anchor
                else:
                    left, right = READ_LEN - anchor, anchor
                start = s - left
                if start < 0 or e + right > len(genome):
                    continue
                seq = genome[start:s] + genome[e:e + right]
                a_pos = [k for k, b in enumerate(seq) if b == "A"]
                for k in rng.sample(a_pos, int(len(a_pos) * args.edit_rate)):
                    seq = seq[:k] + "G" + seq[k + 1:]
                # Sense only. An antisense read carries the same event as T->C
                # and is collapsed by the t2c transform against the t2c genome,
                # a different one of the pipeline's twelve combinations; mixing
                # both here would just halve every number.
                name = f"anch{anchor}_{i:04d}"
                truth[name] = (anchor, e - s)
                f.write(f"@{name}\n{seq.replace('A', 'G')}\n+\n{'I' * len(seq)}\n")

    # --- align, count a junction as found only if it is the right one ---------
    intron_lens = {e - s for s, e in introns}
    rows = []
    for tag in refs:
        for annot in (False, True):
            cmd = ["hisat2", "-p", "8", "--no-unal", "-x",
                   os.path.join(args.outdir, tag), "-U", fq, "-S", "/dev/stdout"]
            if annot:
                cmd += ["--known-splicesite-infile",
                        os.path.join(args.simdir, "splicesites.txt")]
            sam = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout

            found = {a: 0 for a in ANCHORS}
            aligned = {a: 0 for a in ANCHORS}
            wrong = 0
            for line in sam.split("\n"):
                if not line or line.startswith("@"):
                    continue
                fld = line.split("\t")
                name, cigar = fld[0], fld[5]
                if name not in truth:
                    continue
                anchor, ilen = truth[name]
                aligned[anchor] += 1
                ns = [int(x[:-1]) for x in
                      __import__("re").findall(r"\d+N", cigar)]
                if ns and ns[0] == ilen:
                    found[anchor] += 1
                elif ns:
                    wrong += 1
            total = {a: sum(1 for n, (aa, _) in truth.items() if aa == a)
                     for a in ANCHORS}
            for a in ANCHORS:
                rows.append((tag, "annotated" if annot else "de novo", a,
                             found[a], total[a], aligned[a]))
            sys.stderr.write(f"{tag:16s} {'annotated' if annot else 'de novo':10s} "
                             f"junctions at a wrong intron: {wrong}\n")

    out = os.path.join(args.outdir, "splice_motif.tsv")
    with open(out, "w") as f:
        f.write("reference\tmode\tanchor\tfound\ttotal\taligned\trate\n")
        for tag, mode, a, found, total, aligned in rows:
            f.write(f"{tag}\t{mode}\t{a}\t{found}\t{total}\t{aligned}\t{found / total:.4f}\n")

    print(f"\n{'anchor':>7}" + "".join(f"{t + '/' + m:>22}" for t in refs
                                       for m in ("de novo", "annotated")))
    for a in ANCHORS:
        line = f"{a:>7}"
        for t in refs:
            for m in ("de novo", "annotated"):
                r = next(x for x in rows if x[0] == t and x[1] == m and x[2] == a)
                line += f"{r[3] / r[4] * 100:>21.1f}%"
        print(line)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
