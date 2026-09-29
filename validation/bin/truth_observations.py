#!/usr/bin/env python3
"""List every cytosine observation BSReadSim simulated, with its true state.

    truth_observations.py GOLDEN.bam VARIANTS.vcf.gz REF.fa OUT.parquet [--reads KEPT.bam]

GOLDEN.bam is BSReadSim's origin-annotated BAM: each read placed where it was
simulated, with a `zt` tag giving the truth for every base of SEQ (one character per
base from A-Z a-z 0-9 - _, values 0-63; bits 0-1 context 0 none / 1 CG / 2 CHG / 3 CHH,
bit 2 methylated, bit 3 converted, bit 4 changed by a simulated variant, bit 5
sequencing error). VARIANTS.vcf.gz is the phased VCF BSReadSim writes with
`--save-truth`; the `zr` tag's first value, masked with 3, is the read's haplotype, an
index into the VCF genotype (checked: SNVs inside reads carry the variant bit exactly
when that haplotype's allele is the alternate).

`zt` gives a context for cytosines on both strands, but a read measures only the strand
it was converted from: YS OT or CTOT is the top strand, where the reference has C; OB or
CTOB is the bottom strand, where it has G. One row is written per such base with context
CG, CHG or CHH:

    qname, mate (1 or 2), contig, pos (0-based; the G's position for a bottom-strand C),
    context, methylated, converted, variant, error,
    variant_nearby, deletion_nearby,  a variant, a deletion or an insertion in the simulated
    insertion_nearby                  read, within 2 reference bases: BSReadSim takes context
                                      from the simulated molecule, a reference-based caller
                                      from the reference, so these can disagree.
                                      `zt` marks variants only on bases of SEQ, so a
                                      variant past the read end, on the read's haplotype,
                                      is taken from the VCF
    sim_pos, sim_reverse, sim_cigar   where the read was simulated, to recognise reads an
                                      aligner placed elsewhere

A reference cytosine that a simulated variant replaced gets a row with context,
methylated and converted null, so that calls a reference-based caller makes there can be
attributed. A base where a variant created a cytosine (reference A or T) is left out.
`--reads` keeps only reads present in another BAM (for example after removing
overlapping pairs).
"""

import argparse
import string

import pyarrow as pa
import pyarrow.parquet as pq
import pyfaidx
import pysam

STATE64 = string.ascii_uppercase + string.ascii_lowercase + string.digits + "-_"
CONTEXT = {1: "CG", 2: "CHG", 3: "CHH"}
CYTOSINE_ON_TOP = {"OT": "C", "CTOT": "C", "OB": "G", "CTOB": "G"}  # reference base of a measured cytosine
NEARBY = 2


def indel_positions(r: pysam.AlignedSegment) -> tuple[set[int], set[int]]:
    """Reference positions deleted from the read, and the positions either side of an insertion."""
    deleted, inserted, pos = set(), set(), r.reference_start
    for op, length in r.cigartuples:
        if op == pysam.CDEL:
            deleted.update(range(pos, pos + length))
        elif op == pysam.CINS:
            inserted.update((pos - 1, pos))
        if op in (pysam.CMATCH, pysam.CDEL, pysam.CREF_SKIP, pysam.CEQUAL, pysam.CDIFF):
            pos += length
    return deleted, inserted


def haplotype_variant_positions(vcf: str) -> dict[tuple[str, int], set[int]]:
    """Reference positions each simulated variant changes, by (contig, haplotype).

    An SNV changes its base, a deletion the bases after its anchor, and an insertion is
    marked on the bases either side of it, as indel_positions does for the read."""
    positions = {}
    with pysam.VariantFile(vcf) as variants:
        for v in variants:
            ref, alt = len(v.ref), len(v.alts[0])
            if ref == alt:
                changed = range(v.start, v.stop)
            elif ref > alt:
                changed = range(v.start + alt, v.stop)
            else:
                changed = (v.start, v.start + 1)
            for haplotype, allele in enumerate(v.samples[0]["GT"]):
                if allele:
                    positions.setdefault((v.contig, haplotype), set()).update(changed)
    return positions


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("golden")
    p.add_argument("variants")
    p.add_argument("ref")
    p.add_argument("out")
    p.add_argument("--reads", help="keep only reads present in this BAM")
    a = p.parse_args()

    keep = None
    if a.reads:
        with pysam.AlignmentFile(a.reads) as bam:
            keep = {(r.query_name, 2 if r.is_read2 else 1) for r in bam}
    fasta = pyfaidx.Fasta(a.ref)
    haplotype_variants = haplotype_variant_positions(a.variants)
    sequences = {name: str(fasta[name][:]).upper() for name in fasta.keys()}

    cols = {k: [] for k in ("qname", "mate", "contig", "pos", "context", "methylated", "converted", "variant",
                            "error", "variant_nearby", "deletion_nearby", "insertion_nearby",
                            "sim_pos", "sim_reverse", "sim_cigar")}
    with pysam.AlignmentFile(a.golden) as bam:
        for r in bam:
            mate = 2 if r.is_read2 else 1
            if keep is not None and (r.query_name, mate) not in keep:
                continue
            measured = CYTOSINE_ON_TOP[r.get_tag("YS")]
            contig = sequences[r.reference_name]
            zt = r.get_tag("zt")
            states = {pos: STATE64.index(zt[i]) for i, pos in r.get_aligned_pairs(matches_only=True)}
            outside_read = {pos for pos in haplotype_variants.get((r.reference_name, r.get_tag("zr")[0] & 3), ())
                            if not r.reference_start <= pos < r.reference_end}
            variants = {pos for pos, v in states.items() if v & 16} | outside_read
            deleted, inserted = indel_positions(r)
            for pos, v in states.items():
                if contig[pos] != measured or not (v & 3 or v & 16):
                    continue
                window = set(range(pos - NEARBY, pos + NEARBY + 1))
                context = CONTEXT[v & 3] if v & 3 else None
                row = (r.query_name, mate, r.reference_name, pos, context, bool(v & 4) if context else None,
                       bool(v & 8) if context else None, bool(v & 16), bool(v & 32),
                       bool(window & variants - {pos}), bool(window & deleted), bool(window & inserted),
                       r.reference_start, r.is_reverse, r.cigarstring)
                for col, value in zip(cols.values(), row):
                    col.append(value)
    pq.write_table(pa.table(cols), a.out)
    print(f"{len(cols['qname'])} observations")


if __name__ == "__main__":
    main()
