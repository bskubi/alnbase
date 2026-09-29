"""Bismark's XM/XR/XG, derived from its documented semantics, independent of alnbase.

XM: one character per SEQ base, in SEQ order (reference-forward, as stored in the BAM).
  Top-strand reads (OT, CTOT): a call at every reference C; context from the next two
  reference bases (CG -> z/Z, CHG -> x/X, CHH -> h/H, anything involving N -> u/U);
  read C = methylated (uppercase), read T = unmethylated (lowercase), else '.'.
  Bottom-strand reads (OB, CTOB): a call at every reference G; context from the two
  reference bases before it, complemented; read G = methylated, read A = unmethylated.
XR: read conversion  (OT CT, OB CT, CTOT GA, CTOB GA): CT on the two converted
  originals, which are T-rich as sequenced, GA on their complements.
XG: genome conversion (OT CT, OB GA, CTOT CT, CTOB GA).
Only CIGAR nM reads are modelled (no indels, no clips)."""
import sys
ref = open("ref.fa").read().split("\n")[1]
comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
for line in open("reads.sam"):
    if line.startswith("@"):
        continue
    f = line.rstrip("\n").split("\t")
    name, flag, pos, seq = f[0], int(f[1]), int(f[3]) - 1, f[9]
    if flag & 4:
        continue
    r2, rev = bool(flag & 0x80), bool(flag & 0x10)
    strand = {(False, False): "OT", (False, True): "OB", (True, True): "CTOT", (True, False): "CTOB"}[(r2, rev)]
    top = strand in ("OT", "CTOT")
    xm = []
    for i, b in enumerate(seq):
        p = pos + i
        if top:
            if ref[p] != "C":
                xm.append("."); continue
            g1, g2 = ref[p + 1], ref[p + 2]
            meth, unmeth = "C", "T"
        else:
            if ref[p] != "G":
                xm.append("."); continue
            g1, g2 = comp[ref[p - 1]], comp[ref[p - 2]]
            meth, unmeth = "G", "A"
        if g1 == "N" or (g1 != "G" and g2 == "N"):
            ctx = "u"
        elif g1 == "G":
            ctx = "z"
        elif g2 == "G":
            ctx = "x"
        else:
            ctx = "h"
        xm.append(ctx.upper() if b == meth else ctx if b == unmeth else ".")
    xr = {"OT": "CT", "OB": "CT", "CTOT": "GA", "CTOB": "GA"}[strand]
    xg = {"OT": "CT", "CTOT": "CT", "OB": "GA", "CTOB": "GA"}[strand]
    print(f"{name}\tXM:Z:{''.join(xm)}\tXR:Z:{xr}\tXG:Z:{xg}")
