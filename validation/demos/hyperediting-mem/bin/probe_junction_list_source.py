"""Where should the splice-site list come from: the annotation, or the sample?

    probe_junction_list_source.py SIMDIR OUTDIR [--he DIR]

probe_novel_pairing.py established that a junction missing from the list you
supply is worse off than if you had supplied no list at all: the aligner snaps
the read onto a listed junction that shares one of its two sites and returns a
confident alignment on the wrong exon. Omission is not neutral. It is harmful.

That makes "where does the list come from" a question about completeness, and
the two candidate sources fail in opposite directions.

  annotation     GENCODE and friends list every junction of every annotated
                 transcript, so alternative splicing is largely covered, down
                 to isoforms this sample never expresses. What they cannot have
                 is a junction nobody has catalogued yet.

  the first pass The pipeline already maps the whole library against the
                 untransformed genome before it does anything else. Junctions
                 the sample actually uses appear there, in four letters, with
                 the motif intact -- including novel ones. A junction only needs
                 one ordinary read to be discovered, and the read does not have
                 to be the hyper-edited one that failed to map. What this misses
                 is junctions too lowly expressed to be seen at all.

Since the cost is in omission, the two sources should be complementary rather
than rival, and the union should beat either. This measures that, with a gene
built to contain one instance of each failure:

  * two skipping junctions the sample uses and no annotation contains
  * one annotated junction the sample expresses too weakly for the first pass

and a sixth arm adding motif protection at annotated sites on top of the union,
which is the only thing left that can help a junction neither source listed.
"""
import argparse
import os
import random
import re
import subprocess
from collections import Counter

READ_LEN = 100
ANCHORS = [8, 10, 12, 14, 16, 20, 30]
PER_ANCHOR = 300

# How many ordinary, unedited reads the first-pass library carries across each
# junction. Intron 3 is the lowly expressed one and its depth is the swept
# variable (--low-depth), because how deep the first pass has to be before it
# stops needing the annotation is exactly the question the union is there to
# answer.
FIRST_PASS_DEPTH = {"intron1": 400, "intron2": 400, "intron3": None,
                    "skip1": 120, "skip2": 120}


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


def write_sites(path, chrom, pairs):
    """HISAT2 --known-splicesite-infile: last base of the left exon and first
    base of the right exon, both 0-based."""
    with open(path, "w") as f:
        for s, e in sorted(pairs):
            f.write(f"{chrom}\t{s - 1}\t{e}\t+\n")


def junctions_in_cigar(pos, cigar):
    """Reference spans of every N gap, as (first intron base, first base after),
    0-based. pos is the 1-based SAM POS."""
    r = pos - 1
    out = []
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op == "N":
            out.append((r, r + n))
            r += n
        elif op in "MD=X":
            r += n
    return out


def spliced_read(genome, s, e, anchor, left_first):
    """A read crossing the junction (s, e) with `anchor` bases on the short
    side."""
    left, right = (anchor, READ_LEN - anchor) if left_first \
        else (READ_LEN - anchor, anchor)
    return genome[s - left:s] + genome[e:e + right]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("simdir")
    ap.add_argument("outdir")
    ap.add_argument("--seed", type=int, default=23)
    ap.add_argument("--edit-rate", type=float, default=0.35)
    ap.add_argument("--low-depth", type=int, default=1,
                    help="first-pass reads across the lowly expressed intron")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    chrom, genome = read_fasta(os.path.join(args.simdir, "genome.fa"))

    annotated = []
    for line in open(os.path.join(args.simdir, "splicesites.txt")):
        f = line.split()
        annotated.append((int(f[1]) + 1, int(f[2])))

    # The sample's real junction set: the three annotated introns plus two
    # exon-skipping events that reuse annotated sites in unannotated pairings.
    skips = [(annotated[0][0], annotated[1][1]),
             (annotated[1][0], annotated[2][1])]
    names = {annotated[0]: "intron1", annotated[1]: "intron2",
             annotated[2]: "intron3", skips[0]: "skip1", skips[1]: "skip2"}
    truth_set = list(names)

    # --- the first pass -------------------------------------------------------
    # The pipeline's own first stage: the whole library against the untransformed
    # genome, four letters, motif intact, no annotation. Unedited reads only --
    # the edited ones are by definition the ones that fail here.
    fp_fq = os.path.join(args.outdir, "firstpass.fq")
    with open(fp_fq, "w") as f:
        depths = dict(FIRST_PASS_DEPTH, intron3=args.low_depth)
        for junction, depth in depths.items():
            s, e = next(k for k, v in names.items() if v == junction)
            for i in range(depth):
                anchor = rng.randint(8, READ_LEN - 8)
                seq = spliced_read(genome, s, e, anchor, rng.random() < 0.5)
                f.write(f"@fp_{junction}_{i:04d}\n{seq}\n+\n{'I' * len(seq)}\n")

    plain_idx = os.path.join(args.simdir, "chrS.ht2")
    if not os.path.exists(plain_idx + ".1.ht2"):
        plain_idx = os.path.join(args.outdir, "plain")
        if not os.path.exists(plain_idx + ".1.ht2"):
            write_fasta(os.path.join(args.outdir, "plain.fa"), chrom, genome)
            subprocess.run(["hisat2-build", "-q", "-p", "8",
                            os.path.join(args.outdir, "plain.fa"), plain_idx],
                           check=True)

    sam = subprocess.run(["hisat2", "-p", "8", "--no-unal", "-x", plain_idx,
                          "-U", fp_fq, "-S", "/dev/stdout"],
                         check=True, capture_output=True, text=True).stdout
    observed = Counter()
    for line in sam.split("\n"):
        if not line or line.startswith("@"):
            continue
        fld = line.split("\t")
        for j in junctions_in_cigar(int(fld[3]), fld[5]):
            observed[j] += 1

    # A junction seen once could be one misplaced read. Two is the usual floor.
    MIN_READS = 2
    discovered = {j for j, n in observed.items() if n >= MIN_READS}
    print(f"the first pass, run against the untransformed genome "
          f"({args.low_depth} read(s) across intron3):")
    for j in sorted(observed, key=lambda k: -observed[k]):
        label = names.get(j, "NOT A REAL JUNCTION")
        mark = "kept" if observed[j] >= MIN_READS else "dropped"
        print(f"  {label:<20} {observed[j]:>5} reads   {mark}")
    missed = [names[j] for j in truth_set if j not in discovered]
    print(f"  real junctions it missed: {', '.join(missed) or 'none'}")
    spurious = [j for j in discovered if j not in names]
    print(f"  junctions it invented:    {len(spurious)}\n")

    # --- the lists ------------------------------------------------------------
    lists = {
        "annotation only": set(annotated),
        "first pass only": set(discovered),
        "union": set(annotated) | set(discovered),
        "oracle": set(truth_set),
    }
    paths = {}
    for tag, pairs in lists.items():
        p = os.path.join(args.outdir, tag.replace(" ", "_") + ".txt")
        write_sites(p, chrom, pairs)
        paths[tag] = p
        have = sorted(names[j] for j in pairs if j in names)
        print(f"{tag:<16} {len(pairs)} junctions: {', '.join(have)}")
    print()

    # --- references -----------------------------------------------------------
    a2g = genome.replace("A", "G")
    protected = list(a2g)
    for s, e in annotated:
        protected[s:s + 2] = list("GT")
        protected[e - 2:e] = list("AG")
    protected = "".join(protected)

    for tag, seq in {"collapsed": a2g, "protected": protected}.items():
        fa = os.path.join(args.outdir, f"{tag}.fa")
        write_fasta(fa, chrom, seq)
        if not os.path.exists(os.path.join(args.outdir, tag) + ".1.ht2"):
            subprocess.run(["hisat2-build", "-q", "-p", "8", fa,
                            os.path.join(args.outdir, tag)], check=True)

    # --- the reads under test -------------------------------------------------
    # Hyper-edited, three-letter transformed: what actually reaches this stage.
    test_fq = os.path.join(args.outdir, "test.fq")
    truth = {}
    with open(test_fq, "w") as f:
        for anchor in ANCHORS:
            for i in range(PER_ANCHOR):
                s, e = rng.choice(truth_set)
                seq = spliced_read(genome, s, e, anchor, rng.random() < 0.5)
                a_pos = [k for k, b in enumerate(seq) if b == "A"]
                for k in rng.sample(a_pos, int(len(a_pos) * args.edit_rate)):
                    seq = seq[:k] + "G" + seq[k + 1:]
                name = f"t{anchor}_{i:04d}"
                truth[name] = (anchor, names[(s, e)], e - s)
                f.write(f"@{name}\n{seq.replace('A', 'G')}\n+\n{'I' * len(seq)}\n")

    # --- the arms -------------------------------------------------------------
    arms = [
        ("no list", "collapsed", None),
        ("annotation only", "collapsed", paths["annotation only"]),
        ("first pass only", "collapsed", paths["first pass only"]),
        ("union", "collapsed", paths["union"]),
        ("union + motif protection", "protected", paths["union"]),
        # The two arms that ask whether protection is a real belt: one where
        # the list is missing a junction, one where no list is given at all.
        ("first pass + protection", "protected", paths["first pass only"]),
        ("protection, no list", "protected", None),
        ("oracle", "collapsed", paths["oracle"]),
    ]

    by_anchor, by_junction, wrongs = {}, {}, {}
    for label, ref, sites in arms:
        cmd = ["hisat2", "-p", "8", "--no-unal", "-x",
               os.path.join(args.outdir, ref), "-U", test_fq, "-S", "/dev/stdout"]
        if sites:
            cmd += ["--known-splicesite-infile", sites]
        sam = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
        found_a, found_j, wrong = Counter(), Counter(), 0
        for line in sam.split("\n"):
            if not line or line.startswith("@"):
                continue
            fld = line.split("\t")
            if fld[0] not in truth:
                continue
            anchor, junction, ilen = truth[fld[0]]
            ns = [int(x[:-1]) for x in re.findall(r"\d+N", fld[5])]
            if ns and ns[0] == ilen:
                found_a[anchor] += 1
                found_j[junction] += 1
            elif ns:
                wrong += 1
        by_anchor[label], by_junction[label], wrongs[label] = found_a, found_j, wrong

    tot_a = Counter(a for a, _, _ in truth.values())
    tot_j = Counter(j for _, j, _ in truth.values())

    print(f"{'junction found, by anchor':<26}" +
          "".join(f"{a:>6}" for a in ANCHORS) + "   misplaced")
    for label, _, _ in arms:
        row = "".join(f"{by_anchor[label][a] / tot_a[a] * 100:>5.0f}%" for a in ANCHORS)
        print(f"{label:<26}{row}{wrongs[label]:>12}")

    order = ["intron1", "intron2", "intron3", "skip1", "skip2"]
    print(f"\n{'junction found, by junction':<26}" + "".join(f"{j:>9}" for j in order))
    for label, _, _ in arms:
        row = "".join(f"{by_junction[label][j] / tot_j[j] * 100:>8.0f}%" for j in order)
        print(f"{label:<26}{row}")

    out = os.path.join(args.outdir, "list_source.tsv")
    with open(out, "w") as f:
        f.write("arm\tstratum\tkey\tfound\ttotal\trate\n")
        for label, _, _ in arms:
            for a in ANCHORS:
                f.write(f"{label}\tanchor\t{a}\t{by_anchor[label][a]}\t{tot_a[a]}"
                        f"\t{by_anchor[label][a] / tot_a[a]:.4f}\n")
            for j in order:
                f.write(f"{label}\tjunction\t{j}\t{by_junction[label][j]}\t{tot_j[j]}"
                        f"\t{by_junction[label][j] / tot_j[j]:.4f}\n")
            f.write(f"{label}\tmisplaced\t-\t{wrongs[label]}\t{len(truth)}"
                    f"\t{wrongs[label] / len(truth):.4f}\n")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
