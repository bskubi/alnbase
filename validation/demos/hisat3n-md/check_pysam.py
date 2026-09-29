"""Read each alignment's reference bases back through pysam and check them against the FASTA.

`get_aligned_pairs(with_seq=True)` is how bam2bakR's `mut_call.py` and SlamDunk recover the
reference base at each aligned position, so this is the question that decides whether a
malformed MD tag reaches a user or is absorbed by the parser.
"""
import sys
import pysam

sam, fasta = sys.argv[1], sys.argv[2]
ref = {}
for line in open(fasta):
    if line.startswith(">"):
        name = line[1:].split()[0]
        ref[name] = []
    else:
        ref[name].append(line.strip())
ref = {k: "".join(v) for k, v in ref.items()}

ok = bad = 0
for aln in pysam.AlignmentFile(sam):
    try:
        pairs = aln.get_aligned_pairs(matches_only=True, with_seq=True)
        got = "".join(p[2] for p in pairs).upper()
        want = "".join(ref[aln.reference_name][p[1]] for p in pairs)
        status = "ok" if got == want else "MISREAD"
    except Exception as exc:
        status = type(exc).__name__
    ok += status == "ok"
    bad += status != "ok"
    print(f"  {aln.query_name:<16} {aln.cigarstring:<10} "
          f"MD {aln.get_tag('MD'):<14} NM {aln.get_tag('NM')}  pysam {status}")
print(f"\n  pysam read {ok} correctly, {bad} not")
