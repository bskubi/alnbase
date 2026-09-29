"""Build a reference and reads that align with a conversion immediately before a deletion.

Each read is 100 bp taken from the reference with two changes: the base at read offset 49
is changed T->C (a legal conversion under `--base-change T,C`), and the two reference bases
after it are deleted. Sites are skipped when a C appears in the four reference bases after
the conversion, because HISAT-3N would otherwise shift the deletion to absorb the C and
align the read with no mismatch at all.

The intended alignment of every read is POS = s+1, CIGAR 50M2D50M, MD 49T0^<deleted>50.
"""
import random

random.seed(11)
ref = "".join(random.choice("ACGT") for _ in range(20000))

with open("ref.fa", "w") as f:
    f.write(">chr1\n" + "\n".join(ref[i:i + 60] for i in range(0, len(ref), 60)) + "\n")

n = 0
with open("r1.fq", "w") as f:
    for s in range(200, 19000, 137):
        if ref[s + 49] != "T" or "C" in ref[s + 50:s + 54]:
            continue
        read = ref[s:s + 49] + "C" + ref[s + 52:s + 102]
        f.write(f"@s{s}_del{ref[s + 50:s + 52]}\n{read}\n+\n{'I' * len(read)}\n")
        n += 1

print(f"{n} reads written; each should align 50M2D50M with MD 49T0^<del>50")
