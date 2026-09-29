#!/usr/bin/env python
"""SlamDunk's SNP mask is keyed by chromosome name concatenated with position.

SNPtools.py:33 builds the key as `snp[0] + snp[1]` with no separator, so
"chr11" + "234" and "chr1" + "1234" are the same string. A T->C SNP on one
chromosome therefore masks genuine T->C conversions on the other.

Run with SLAMDUNK=<path to a slamdunk checkout>, so the real class is exercised.
"""
import os
import sys

sys.path.insert(0, os.environ.get("SLAMDUNK", "."))
from slamdunk.utils.SNPtools import SNPDictionary  # noqa: E402

# GRCh38 primary assembly.
LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}


def show_collision():
    """One declared SNP, queried on two different chromosomes."""
    snps = SNPDictionary(None)
    # A VCF row is CHROM POS ID REF ALT; SNPtools reads snp[3]/snp[4] as REF/ALT.
    snps._addSNP(["chr11", "234", ".", "T", "C"])
    print("declared: one T->C SNP at chr11:234 (1-based)\n")
    for chrom, pos in [("chr11", 233), ("chr1", 1233), ("chr12", 233), ("chr2", 1233)]:
        print("  isTCSnp(%-6s, %6d)  1-based %7d  -> %s"
              % (chrom, pos, pos + 1, snps.isTCSnp(chrom, pos)))
    print("\n  key is snp[0] + snp[1]:  'chr11' + '234' = %s"
          % ("chr11" + "234"))
    print("                           'chr1'  + '1234' = %s"
          % ("chr1" + "1234"))


def colliding_positions():
    """Every position whose key is shared with a position on another chromosome.

    Two names collide only when one is a prefix of the other AND the concatenated
    digits are a valid str(pos) -- no leading zero, so chr1/chr10 is safe.
    """
    pairs = []
    for a in LENGTHS:
        for b in LENGTHS:
            if a == b or not b.startswith(a):
                continue
            suffix = b[len(a):]
            if suffix[0] == "0":
                continue
            n = 0
            for digits in range(1, 10):
                lo = 1 if digits == 1 else 10 ** (digits - 1)
                hi = min(10 ** digits - 1, LENGTHS[b])
                if lo > hi:
                    continue
                first = int(suffix + str(lo))
                last = min(int(suffix + str(hi)), LENGTHS[a])
                if first <= last:
                    n += last - first + 1
            if n:
                pairs.append((a, b, n))
    return sorted(pairs, key=lambda p: -p[2])


def main():
    print("=== 1. One SNP, two chromosomes ===")
    show_collision()

    print("\n=== 2. Footprint on GRCh38 ===")
    pairs = colliding_positions()
    for a, b, n in pairs:
        print("  %-6s <-> %-6s  %11s  (%5.1f%% of %s)"
              % (a, b, format(n, ","), 100.0 * n / LENGTHS[a], a))

    per_chrom = {}
    for a, _, n in pairs:
        per_chrom[a] = per_chrom.get(a, 0) + n
    total = sum(per_chrom.values())
    print()
    for a in sorted(per_chrom, key=lambda c: -per_chrom[c]):
        print("  %-6s %11s of %11s positions affected  (%.1f%%)"
              % (a, format(per_chrom[a], ","), format(LENGTHS[a], ","),
                 100.0 * per_chrom[a] / LENGTHS[a]))
    print("\n  total %s bp, %.1f%% of the primary assembly"
          % (format(total, ","), 100.0 * total / sum(LENGTHS.values())))


if __name__ == "__main__":
    main()
