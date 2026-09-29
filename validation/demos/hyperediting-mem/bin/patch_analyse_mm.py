"""Put the alignment's CIGAR on analyse_mm.pl's hit lines, and optionally make
it survive an aligner that writes no X0/X1.

    patch_analyse_mm.py IN.pl OUT.pl [--mem-hits]

The CIGAR field is not a behaviour change: it adds a column that the original
consumers ignore and the patched detect_ue.pl uses to convert a read offset
into a genomic coordinate. Every arm's caller emits it, so all three feed the
same downstream code and the comparison is between callers alone.

--mem-hits is needed only for the arm that runs the original caller over bwa
mem output. analyse_mm.pl counts a read's alignment positions from the X0 and
X1 tags, which `bwa aln` writes and mem does not. The count comes out zero, and
sort_R_read.pl reads exactly that many hit lines after each record, so it falls
out of step with the file and mis-parses everything after the first read.
"""
import argparse
import sys

PRIMARY_OLD = '    my $hits2print = "$Fields[0]\\t$Fields[2]\\t$Fields[3]\\t$ref_seq\\t$mm_str\\n";'
PRIMARY_NEW = '    my $hits2print = "$Fields[0]\\t$Fields[2]\\t$Fields[3]\\t$ref_seq\\t$mm_str\\t$Fields[5]\\n";'

XA_OLD = '\t\t$hits2print = "$hits2print$Fields[0]\\t$sub_optimal[0]\\t$sub_optimal[1]\\t$ref_seq\\t$mm_str\\n";'
XA_NEW = '\t\t$hits2print = "$hits2print$Fields[0]\\t$sub_optimal[0]\\t$sub_optimal[1]\\t$ref_seq\\t$mm_str\\t$Fields[5]\\n";'

HITS_ANCHOR = '    if ($hits > ($hits_amount+1))'
HITS_NEW = ('    $hits = 1 if $hits == 0; # bwa mem writes no X0/X1\n'
            '    if ($hits > ($hits_amount+1))')


def replace_once(text, old, new, label):
    if text.count(old) != 1:
        sys.exit(f"patch_analyse_mm: expected one {label}, found {text.count(old)}")
    return text.replace(old, new)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--mem-hits", action="store_true")
    ap.add_argument("--no-cigar", action="store_true",
                    help="leave the hit line in its published shape; for the arm "
                         "that models swapping the aligner and changing nothing else")
    args = ap.parse_args()

    s = open(args.src).read()
    if not args.no_cigar:
        s = replace_once(s, PRIMARY_OLD, PRIMARY_NEW, "primary hit line")
        s = replace_once(s, XA_OLD, XA_NEW, "alternative-position hit line")
    if args.mem_hits:
        s = replace_once(s, HITS_ANCHOR, HITS_NEW, "hit-count test")
    open(args.dst, "w").write(s)
    print(f"patched {args.src} -> {args.dst}")


if __name__ == "__main__":
    main()
