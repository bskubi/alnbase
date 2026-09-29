#!/usr/bin/env python3
"""A small bwa-meth-style BAM for exercising ScaleMethyl's met_extract.py.

This is load, not ground truth: the point is to give met_extract.py many
chromosomes (it starts one task per chromosome with reads) and read names in
the ScaleBio form it parses barcodes from (`name:BC1+BC2+BC3`).

    make_scalemethyl_input.py --contigs 40 --reads-per-contig 2000 --outdir data

Writes data/ref.fa (+ .fai) and data/reads.bam (coordinate-sorted, indexed).
Reads are forward (original top strand), 100 bp, unmethylated outside CG:
every non-CG C in the reference reads as T, and CG cytosines read as C.
"""

import argparse
import random
from pathlib import Path

import pysam

READ_LEN = 100


def bisulfite_top_strand(ref: str) -> str:
    """Convert every C that is not followed by G, as a directional OT read would show."""
    out = []
    for i, base in enumerate(ref):
        next_base = ref[i + 1] if i + 1 < len(ref) else "N"
        out.append("T" if base == "C" and next_base != "G" else base)
    return "".join(out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--contigs", type=int, default=40)
    parser.add_argument("--contig-length", type=int, default=20_000)
    parser.add_argument("--reads-per-contig", type=int, default=2_000)
    parser.add_argument("--cells", type=int, default=50)
    parser.add_argument("--outdir", type=Path, default=Path("data"))
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    args.outdir.mkdir(parents=True, exist_ok=True)
    fasta = args.outdir / "ref.fa"
    unsorted = args.outdir / "unsorted.bam"
    bam = args.outdir / "reads.bam"

    contigs = {f"chr{i + 1}": "".join(rng.choice("ACGT") for _ in range(args.contig_length))
               for i in range(args.contigs)}
    with open(fasta, "w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n{seq}\n")
    pysam.faidx(str(fasta))

    barcodes = ["+".join("".join(rng.choice("ACGT") for _ in range(8)) for _ in range(3))
                for _ in range(args.cells)]
    header = {"HD": {"VN": "1.6", "SO": "unsorted"},
              "SQ": [{"SN": name, "LN": len(seq)} for name, seq in contigs.items()]}

    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for tid, (name, seq) in enumerate(contigs.items()):
            for n in range(args.reads_per_contig):
                start = rng.randrange(2, len(seq) - READ_LEN - 4)
                read = pysam.AlignedSegment(out.header)
                read.query_name = f"{name}_r{n}:{rng.choice(barcodes)}"
                read.reference_id = tid
                read.reference_start = start
                read.mapping_quality = 60
                read.cigarstring = f"{READ_LEN}M"
                read.query_sequence = bisulfite_top_strand(seq[start:start + READ_LEN + 1])[:READ_LEN]
                read.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
                read.set_tag("YD", "f")  # bwa-meth: read from the original top strand
                out.write(read)

    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    unsorted.unlink()
    print(f"wrote {fasta} and {bam}")


if __name__ == "__main__":
    main()
