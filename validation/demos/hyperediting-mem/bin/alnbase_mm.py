"""Drop-in replacement for the hyper-editing pipeline's edit caller.

    alnbase_mm.py GENOME_FASTA BAM --aref REF.aref [--work DIR] > hit_lines

Emits analyse_mm.pl's format, so sort_R_read.pl runs unchanged:

    <the SAM record>\thits:N
    <qname>\t<rname>\t<pos>\t<ref_seq>\t<mm_count>\t<edit_type>\t<edit_count>\t<edit_locs>\t<mm_locs>\t<cigar>

with one field added. detect_ue.pl turns an edit's offset along the read into a
genomic coordinate by adding it to the alignment's start position, which is the
same flat assumption this file exists to remove, made a second time downstream.
It needs the CIGAR to do that conversion properly, and it cannot take it from
the SAM line it already holds: sort_R_read.pl groups a read's hits from several
transform combinations under one record, whose CIGAR need not be the selected
hit's. See bin/patch_detect_ue.py.

One thing is different, and it is the whole point. analyse_mm.pl decides which
reference base a read base sits on by cutting a window out of the genome at the
alignment's start position and pairing position i with position i:

    substr($ref_hash{$Fields[2]}, $Fields[3]-1, length($Fields[9]))

which is only true when the CIGAR is a single run of M. Here the pairing comes
from alnbase, which walks the CIGAR, so clips, insertions, deletions and splice
junctions are placed correctly and a read carrying any of them can be used at
all.

What is deliberately NOT changed: the classification arithmetic. Which of the
mismatches counts as the expected edit, how the two directions are compared,
how mm_locs is assembled, and the N-against-N quirk at analyse_mm.pl:60 are all
reproduced as written, including their bugs. The comparison this demo runs is
between arms that differ only in how read bases are paired with reference
bases; folding other fixes in at the same time would make the difference
unattributable. `--fix-nn` turns the one behavioural quirk off for anyone who
wants the corrected caller rather than the comparable one.

ref_seq is now a projection rather than a window: entry i is the reference base
alnbase paired with read base i, flanked by the genomic bases immediately
outside the alignment, so it keeps the +/-1 padding detect_ue.pl expects. Read
bases with no reference base -- soft clips, insertions -- are filled with the
read's own base, so they contribute neither a mismatch nor an edit. For a read
with a deletion the projection skips the deleted reference bases, so the
three-base context detect_ue.pl cuts around an edit sitting next to a deletion
is the read's neighbours rather than the genome's; every other context is exact.
"""
import argparse
import glob
import os
import subprocess
import sys
import tempfile

import pyarrow.parquet as pq
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
QUERY = os.path.join(HERE, os.pardir, "queries", "hyperedit.toml")

REVERSE_MM_TYPES = {"A2G": "T2C", "A2C": "T2G"}

EXTRACT_FIELDS = ("core,cigar,read_len,hard_clip_5p,hard_clip_3p,"
                  "is_reverse,pos,end_pos")


def expected_edit(comb):
    """analyse_mm.pl:39-48. Combos are like A2G++ or A2G-+; a '-' in the
    fourth character means the genome was transformed the other way round, and
    the edit is then seen as its complement."""
    if len(comb) >= 5 and comb[3] == "-":
        return REVERSE_MM_TYPES.get(comb[:3], "")
    return comb[:3]


def classify(comb, ref_seq, read, fix_nn=False):
    """analyse_mm.pl's analyse_mm(), same arithmetic, same output string."""
    exp = expected_edit(comb)
    mm_count = 0
    mm_loc = []
    edit_loc1, edit_loc2 = [], []

    for i, (r, d) in enumerate(zip(read.upper(), ref_seq.upper())):
        if r != d:
            mm_count += 1
            mm_loc.append(i)
        elif r == "N" and d == "N":
            # analyse_mm.pl:60-66 counts a position where both sides are N as a
            # mismatch and moves on. Reproduced so the arms stay comparable.
            if not fix_nn:
                mm_count += 1
                mm_loc.append(i)
            continue

        if exp and d == exp[0]:
            if r == exp[2]:
                edit_loc1.append(i)
        elif exp and d == exp[2]:
            if r == exp[0]:
                edit_loc2.append(i)

    if len(edit_loc2) > len(edit_loc1):
        edit_count, edit_type, edit_loc = len(edit_loc2), exp[::-1], edit_loc2
    else:
        edit_count, edit_type, edit_loc = len(edit_loc1), exp, edit_loc1

    chosen = set(edit_loc)
    mm_locs = ";".join(str(m) for m in mm_loc if m not in chosen) or ";"
    return (f"{mm_count}\t{edit_type}\t{edit_count}\t"
            f"{';'.join(str(e) for e in edit_loc)}\t{mm_locs}")


def run_alnbase(alnbase, bam, aref, strand_rule, workdir, threads):
    """Tag and extract, returning the parquet rows as column lists."""
    tagged = os.path.join(workdir, "tagged.bam")
    prefix = os.path.join(workdir, "hits.parquet")
    for stale in glob.glob(os.path.join(workdir, "hits_*.parquet")):
        os.remove(stale)
    subprocess.run(
        [alnbase, "query", "--query-file", QUERY, "--query-file", strand_rule,
         "-@", str(threads), bam, aref, tagged],
        check=True, stdout=sys.stderr)
    subprocess.run(
        [alnbase, "extract", "--tag", "XE", "-@", str(threads),
         "-f", EXTRACT_FIELDS, tagged, prefix],
        check=True, stdout=sys.stderr)

    parts = sorted(glob.glob(os.path.join(workdir, "hits_*.parquet")))
    if not parts:
        return {}
    table = pq.read_table(parts)
    cols = {n: table.column(n).to_pylist() for n in
            ("record_id", "off_5p", "read_base", "refr_base", "read_len",
             "hard_clip_5p", "is_reverse")}

    # record_id -> {index into SEQ: reference base}
    per_record = {}
    for k in range(table.num_rows):
        if cols["refr_base"][k] is None:
            continue                      # a record that produced no hit
        rid = cols["record_id"][k]
        off5 = cols["off_5p"][k] - cols["hard_clip_5p"][k]
        seq_len = cols["read_len"][k]
        idx = seq_len - 1 - off5 if cols["is_reverse"][k] else off5
        per_record.setdefault(rid, {})[idx] = cols["refr_base"][k]
    return per_record


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genome")
    ap.add_argument("bam")
    ap.add_argument("--aref", required=True)
    ap.add_argument("--strand-rule", required=True)
    ap.add_argument("--alnbase", default="alnbase")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--fix-nn", action="store_true")
    ap.add_argument("--work")
    args = ap.parse_args()

    comb = os.path.basename(args.bam)[:5].upper()

    contigs = {}
    name = None
    chunks = []
    with open(args.genome) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    contigs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name:
        contigs[name] = "".join(chunks)

    work = args.work or tempfile.mkdtemp(prefix="alnbase_mm.")
    os.makedirs(work, exist_ok=True)
    per_record = run_alnbase(args.alnbase, args.bam, args.aref,
                             args.strand_rule, work, args.threads)

    out = sys.stdout
    with pysam.AlignmentFile(args.bam, "rb", check_sq=False) as bam:
        for rid, rec in enumerate(bam):
            if rec.is_unmapped:
                # record_id is the record's ordinal in the file, and it is the
                # only key tying a parquet row back to a record, so the two
                # walks have to agree on what they are counting. The pipeline
                # hands this stage a `samtools view -F 4` file, so an unmapped
                # record here means an assumption has broken rather than a case
                # to handle quietly.
                sys.exit(f"{args.bam}: unmapped record {rec.query_name}; "
                         "this stage expects a -F 4 filtered file")
            seq = rec.query_sequence.upper()
            chrom = rec.reference_name
            pos = rec.reference_start + 1          # SAM POS, 1-based
            if rec.cigartuples and any(op == 5 for op, _ in rec.cigartuples):
                sys.exit(f"{args.bam}: hard-clipped record {rec.query_name}; "
                         "align with `bwa mem -Y` so SEQ stays whole")

            mismatches = per_record.get(rid, {})
            proj = list(seq)
            for idx, refb in mismatches.items():
                if 0 <= idx < len(proj):
                    proj[idx] = refb
            ref_seq = "".join(proj)

            contig = contigs[chrom]
            before = contig[pos - 2].upper() if pos - 2 >= 0 else "N"
            after_i = rec.reference_end
            after = contig[after_i].upper() if after_i < len(contig) else "N"
            padded = before + ref_seq + after

            mm_str = classify(comb, ref_seq, seq, args.fix_nn)
            out.write(rec.to_string() + "\thits:1\n")
            # The CIGAR rides along so that detect_ue.pl can convert a read
            # offset into a genomic coordinate for the right alignment: it
            # groups a read's hits from several transform combinations under
            # one SAM line, whose CIGAR need not be this hit's.
            out.write(f"{rec.query_name}\t{chrom}\t{pos}\t{padded}\t{mm_str}"
                      f"\t{rec.cigarstring}\n")


if __name__ == "__main__":
    main()
