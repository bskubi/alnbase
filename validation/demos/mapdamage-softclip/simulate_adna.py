"""Simulate an ancient-DNA library with a known damage profile, and write FASTQ.

Damage model: the standard one mapDamage itself fits. Cytosine deamination gives
C>T near the fragment's 5' end and G>A near its 3' end, with the per-position
probability decaying geometrically from the terminus. Ground truth is recorded
against the position in the ORIGINAL fragment, which is what an off_5p-style
query is supposed to report.

Usage: simulate_adna.py <prefix> <damage-at-position-1> [seed]
Writes <prefix>_ref.fa and <prefix>.fq, and prints the true C>T frequency per
5' position as a TSV to <prefix>_truth.tsv.
"""

import pathlib
import random
import sys

REF_LEN = 400_000
NREADS = 60_000
DECAY = 0.55
SEQ_ERR = 0.001
MIN_LEN, MAX_LEN = 35, 90
NPOS = 10
COMP = str.maketrans("ACGT", "TGCA")


def main(prefix, damage, seed):
    random.seed(seed)
    out = pathlib.Path(prefix)

    ref = "".join(random.choice("ACGT") for _ in range(REF_LEN))
    out.with_name(out.name + "_ref.fa").write_text(
        ">chr1\n" + "\n".join(ref[i : i + 60] for i in range(0, REF_LEN, 60)) + "\n"
    )

    truth = {i: [0, 0] for i in range(1, NPOS + 1)}
    with open(out.with_name(out.name + ".fq"), "w") as fq:
        for n in range(NREADS):
            length = random.randint(MIN_LEN, MAX_LEN)
            start = random.randrange(0, REF_LEN - length)
            frag = list(ref[start : start + length])

            # 5' end: C>T, recorded as ground truth
            for i in range(min(NPOS, length)):
                if frag[i] == "C":
                    truth[i + 1][0] += 1
                    if random.random() < damage * DECAY**i:
                        frag[i] = "T"
                        truth[i + 1][1] += 1
            # 3' end: G>A, the mirror image
            for i in range(min(NPOS, length)):
                j = length - 1 - i
                if frag[j] == "G" and random.random() < damage * DECAY**i:
                    frag[j] = "A"

            read = "".join(frag)
            if random.random() < 0.5:  # sequenced from the fragment's other end
                read = read.translate(COMP)[::-1]
            read = "".join(
                random.choice("ACGT") if random.random() < SEQ_ERR else b for b in read
            )
            fq.write("@f%d\n%s\n+\n%s\n" % (n, read, "I" * len(read)))

    with open(out.with_name(out.name + "_truth.tsv"), "w") as handle:
        handle.write("Pos\tC\tC_to_T\tFreq\n")
        for i in range(1, NPOS + 1):
            c, t = truth[i]
            handle.write("%d\t%d\t%d\t%.6f\n" % (i, c, t, t / c))


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 11)
