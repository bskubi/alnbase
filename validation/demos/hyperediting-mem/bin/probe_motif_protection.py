"""Can you keep the splice motifs by exempting them from the three-letter collapse?

    probe_motif_protection.py SIMDIR OUTDIR

The idea is reasonable: the transform is what destroys GT..AG, so transform
everything except the donor and acceptor dinucleotides and the motifs survive.
It comes in two forms, and they behave completely differently.

  targeted   exempt the motif only at annotated intron boundaries. Costs one
             base per intron and does restore motif support -- this is the
             `motif-restored` reference in probe_splice_motif.py. But it needs
             the annotation, and anyone who has the annotation can hand the
             aligner the junction list instead, which works better.

  global     exempt every occurrence of the motif anywhere in the genome. This
             needs no annotation, which is the whole appeal. This script
             measures it.

Global protection breaks the method, and the reason is worth stating as a rule.
The three-letter trick works because the transform is *context-free*: each base
maps to its image regardless of its neighbours, which is exactly what guarantees
that an edited A and an unedited A become the same letter. A rule like "leave an
A alone when the next base is G" is context-dependent, and the context can
itself be edited. A genomic AG whose adenosine was deaminated reads GG. The
genome keeps its A, the read has a G, and the position becomes a hard mismatch
-- at an edited base, which is precisely the base the pipeline exists to see.

Under a2g only the acceptor needs protecting, because the donor GT contains no
adenosine and already survives. That is the most favourable case for the idea,
and it is the one measured here.
"""
import argparse
import os
import subprocess

COMP = str.maketrans("ACGT", "TGCA")


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
    ap.add_argument("--he", default=None, help="the patched Hyper-editing tree")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    chrom, genome = read_fasta(os.path.join(args.simdir, "genome.fa"))
    he = args.he or os.path.join(args.simdir, "..", "Hyper-editing")

    plain = genome.replace("A", "G")

    # Global protection: leave an adenosine unconverted whenever the next base
    # is G, so every genomic AG keeps its acceptor.
    prot = []
    for i, b in enumerate(genome):
        if b == "A":
            prot.append("A" if (i + 1 < len(genome) and genome[i + 1] == "G") else "G")
        else:
            prot.append(b)
    prot = "".join(prot)

    kept = sum(1 for i in range(len(genome)) if genome[i] == "A" and prot[i] == "A")
    print(f"genome {len(genome):,} bp, {genome.count('A'):,} adenosines; "
          f"{kept:,} left unconverted by global protection "
          f"({kept / len(genome) * 100:.2f}% of all positions)\n")

    refs = {"plain": plain, "protected": prot}
    for tag, seq in refs.items():
        fa = os.path.join(args.outdir, f"{tag}.fa")
        write_fasta(fa, chrom, seq)
        idx = os.path.join(args.outdir, tag)
        if not os.path.exists(idx + ".1.ht2"):
            subprocess.run(["hisat2-build", "-q", "-p", "8", fa, idx], check=True)

    fq = os.path.join(args.outdir, "a2g.fq")
    with open(fq, "w") as f:
        subprocess.run([os.path.join(he, "fastq_transform.pl"), "a", "g",
                        os.path.join(args.simdir, "reads.fastq")],
                       check=True, stdout=f)

    print(f"{'reference':<12}{'alignment rate':>16}{'edited reads placed':>22}{'junction gaps':>16}")
    for tag in refs:
        sam = os.path.join(args.outdir, f"{tag}.sam")
        with open(sam, "w") as out:
            r = subprocess.run(["hisat2", "-p", "8", "-x", os.path.join(args.outdir, tag),
                                "-U", fq], check=True, stdout=out, stderr=subprocess.PIPE,
                               text=True)
        rate = next(l for l in r.stderr.split("\n") if "overall alignment rate" in l).split()[0]
        edited = junc = 0
        for line in open(sam):
            if line.startswith("@"):
                continue
            fld = line.split("\t")
            if int(fld[1]) & 4:
                continue
            if "_edited" in fld[0]:
                edited += 1
            if fld[0].startswith("junction") and "N" in fld[5]:
                junc += 1
        print(f"{tag:<12}{rate:>16}{edited:>22,}{junc:>16,}")


if __name__ == "__main__":
    main()
