#!/usr/bin/env python3
"""Compare two per-base methylation call strings in one BAM, and both with the truth.

    compare_xm_tags.py TAGGED.bam GOLDEN.bam OUT_PREFIX [--tool-tag XM] [--alnbase-tag xm]

TAGGED.bam is an aligner's BAM (here Bismark's, carrying `XM`) that alnbase has
tagged with its own call string (`xm`). GOLDEN.bam is BSReadSim's origin-annotated
BAM for the same reads, whose `zt` tag holds the truth for every base of SEQ.

Writes OUT_PREFIX.calls.parquet, one row per base where either call string makes a
call, and prints a summary. Both call strings use Bismark's letters: Z/z CG, X/x CHG,
H/h CHH, U/u unknown context; upper case methylated, lower case unmethylated.

Truth (BSReadSim docs/outputs/index.md, `zt`): one character per SEQ base from the
alphabet A-Z a-z 0-9 - _ (values 0-63). Bits 0-1 context (0 none, 1 CG, 2 CHG,
3 CHH), bit 2 methylated, bit 3 conversion succeeded, bit 4 affected by a variant,
bit 5 sequencing error. Truth is used only where the aligner placed the read where
it was simulated (same position, strand and CIGAR), so SEQ indices correspond.
"""

import argparse
import string
import sys

import pyarrow as pa
import pyarrow.parquet as pq
import pysam

STATE64 = string.ascii_uppercase + string.ascii_lowercase + string.digits + "-_"
CONTEXT = {0: None, 1: "CG", 2: "CHG", 3: "CHH"}
CALL_CONTEXT = {"Z": "CG", "X": "CHG", "H": "CHH", "U": "unknown"}


def decode(ch: str) -> dict:
    v = STATE64.index(ch)
    return {
        "truth_context": CONTEXT[v & 3],
        "truth_methylated": bool(v & 4),
        "truth_converted": bool(v & 8),
        "truth_variant": bool(v & 16),
        "truth_error": bool(v & 32),
    }


def golden_records(path: str) -> dict:
    """(read name, mate) -> (position, reverse, CIGAR, zt) for every simulated read."""
    out = {}
    with pysam.AlignmentFile(path) as bam:
        for r in bam:
            mate = 2 if r.is_read2 else 1
            out[(r.query_name, mate)] = (r.reference_start, r.is_reverse, r.cigarstring, r.get_tag("zt"))
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("tagged")
    p.add_argument("golden")
    p.add_argument("out_prefix")
    p.add_argument("--tool-tag", default="XM")
    p.add_argument("--alnbase-tag", default="xm")
    args = p.parse_args()

    truth = golden_records(args.golden)
    cols = {k: [] for k in ("qname", "mate", "contig", "refr_pos", "seq_index", "tool_call", "alnbase_call",
                            "placed_as_simulated", "truth_context", "truth_methylated", "truth_converted",
                            "truth_variant", "truth_error")}
    with pysam.AlignmentFile(args.tagged) as bam:
        for r in bam:
            if r.is_unmapped:
                continue
            tool, ours = r.get_tag(args.tool_tag), r.get_tag(args.alnbase_tag)
            if len(tool) != len(ours):
                sys.exit(f"{r.query_name}: call strings differ in length")
            # Bismark names reads "<name>_1" / "<name>_2" only for some inputs; strip that.
            name = r.query_name[:-2] if r.query_name.endswith(("_1", "_2")) else r.query_name
            mate = 2 if r.is_read2 else 1
            sim = truth.get((name, mate))
            placed = sim is not None and sim[:3] == (r.reference_start, r.is_reverse, r.cigarstring)
            positions = r.get_reference_positions(full_length=True)
            for i, (a, b) in enumerate(zip(tool, ours)):
                if a == "." and b == ".":
                    continue
                t = decode(sim[3][i]) if placed else dict.fromkeys(
                    ("truth_context", "truth_methylated", "truth_converted", "truth_variant", "truth_error"))
                row = dict(qname=name, mate=mate, contig=r.reference_name, refr_pos=positions[i], seq_index=i,
                           tool_call=a, alnbase_call=b, placed_as_simulated=placed, **t)
                for k, v in row.items():
                    cols[k].append(v)

    table = pa.table(cols)
    pq.write_table(table, f"{args.out_prefix}.calls.parquet")

    n = table.num_rows
    agree = sum(1 for a, b in zip(cols["tool_call"], cols["alnbase_call"]) if a == b)
    print(f"{n} called bases, {agree} identical, {n - agree} different")


if __name__ == "__main__":
    main()
