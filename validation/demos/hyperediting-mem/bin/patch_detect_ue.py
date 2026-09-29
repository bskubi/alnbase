"""Make detect_ue.pl convert read offsets to genomic coordinates through the CIGAR.

    patch_detect_ue.py IN.pl OUT.pl

The flat-window assumption lives in two places, not one. The edit caller makes
it when it cuts a genome window at the alignment's start position, and
detect_ue.pl makes it again when it turns an edit's offset along the read into
a genomic coordinate:

    my $es_bed_ref_base = $ref_base-1+$ed_site;

Fixing only the caller leaves every recovered clipped, indel-bearing or
junction-spanning read placed at the wrong base, which is worse than not
recovering it. So the port touches this file too.

Three changes, and each one is a no-op on a read whose CIGAR is a single run of
M -- which is every read the published pipeline can produce, so the patched
script reproduces the original's output exactly on the original's input.

  1. A `ref_offset_of` subroutine: walk the CIGAR, return how far along the
     reference a given offset into SEQ sits, or undef for a read base that has
     no reference base under it (an insertion, or a soft clip).

  2. The hit line gains a CIGAR field. It has to: sort_R_read.pl groups a
     read's alignments from different transform combinations under one SAM
     line, so the CIGAR in @Fields is not necessarily the CIGAR of the hit that
     was selected. Carrying it on the hit line is the only way the conversion
     can use the right one. The field goes before the strand mark that
     sort_R_read.pl appends, so that file needs no change at all.

  3. The three coordinate expressions go through the conversion.
"""
import sys

REF_OFFSET_SUB = r'''
# --- added by the alnbase port -----------------------------------------------
# How far along the reference the read base at index $i into SEQ sits, counting
# from the alignment's start. undef when that read base has no reference base
# under it: an inserted base, or one inside a soft clip. On a CIGAR that is a
# single run of M this returns $i, which is what the unported code assumed for
# every read.
sub ref_offset_of
{
    my ($cigar, $i) = @_;
    return $i if (!defined $cigar or $cigar eq "*");
    my ($q, $r) = (0, 0);
    while ($cigar =~ /(\d+)([MIDNSHP=X])/g)
    {
        my ($n, $op) = ($1, $2);
        if ($op eq "M" or $op eq "=" or $op eq "X")
        {
            return $r + ($i - $q) if ($i < $q + $n);
            $q += $n; $r += $n;
        }
        elsif ($op eq "I" or $op eq "S")
        {
            return undef if ($i < $q + $n);
            $q += $n;
        }
        elsif ($op eq "D" or $op eq "N")
        {
            $r += $n;
        }
    }
    return undef;
}
# -----------------------------------------------------------------------------
'''

SPLIT_OLD = ('($read_id,$ref_chr,$ref_base,$ref_seq,$mm_count,$edit_type,'
             '$edit_count,$edit_loc,$mm_loc,$read_sign) = split ("\\t",$hit);')
SPLIT_NEW = ('($read_id,$ref_chr,$ref_base,$ref_seq,$mm_count,$edit_type,'
             '$edit_count,$edit_loc,$mm_loc,$hit_cigar,$read_sign) = split ("\\t",$hit);')

DECL_OLD = ('    my ($read_id,$ref_chr,$ref_base,$ref_seq,$mm_count,$edit_type,'
            '$edit_count,$edit_loc,$mm_loc,$read_sign);')
DECL_NEW = ('    my ($read_id,$ref_chr,$ref_base,$ref_seq,$mm_count,$edit_type,'
            '$edit_count,$edit_loc,$mm_loc,$hit_cigar,$read_sign);')

# the hit line detect_ue.pl rebuilds for its own use, kept in the new shape
REBUILD_OLD = ('    $hit = "$read_id\\t$ref_chr\\t$ref_base\\t$ref_seq\\t$mm_count'
               '\\t$edit_type\\t$edit_count\\t$edit_loc\\t$mm_loc";')
REBUILD_NEW = ('    $hit = "$read_id\\t$ref_chr\\t$ref_base\\t$ref_seq\\t$mm_count'
               '\\t$edit_type\\t$edit_count\\t$edit_loc\\t$mm_loc\\t$hit_cigar";')

CLUSTER_OLD = '''    my $bed_ref_base = $ref_base-1+$edit_sites[0];
    my $end_pos=$ref_base+$edit_sites[scalar(@edit_sites)-1];'''
CLUSTER_NEW = '''    my @edit_ref_offs = grep { defined } map { ref_offset_of($hit_cigar,$_) } @edit_sites;
    if (!@edit_ref_offs) { next; } # every edit sits on a base with no reference base
    my @sorted_ref_offs = sort {$a <=> $b} @edit_ref_offs;
    my $bed_ref_base = $ref_base-1+$sorted_ref_offs[0];
    my $end_pos=$ref_base+$sorted_ref_offs[scalar(@sorted_ref_offs)-1];'''

SITE_OLD = '''	my $es_bed_ref_base = $ref_base-1+$ed_site;
	my $es_bed_end_base = $es_bed_ref_base+1;'''
SITE_NEW = '''	my $ed_ref_off = ref_offset_of($hit_cigar,$ed_site);
	next if (!defined $ed_ref_off); # inserted or clipped: no base to name
	my $es_bed_ref_base = $ref_base-1+$ed_ref_off;
	my $es_bed_end_base = $es_bed_ref_base+1;'''


def replace_once(text, old, new, label, count=1):
    n = text.count(old)
    if n != count:
        sys.exit(f"patch_detect_ue: expected {count} of {label}, found {n}")
    return text.replace(old, new)


def main():
    src, dst = sys.argv[1], sys.argv[2]
    div0_only = "--div0-only" in sys.argv[3:]
    s = open(src).read()

    if div0_only:
        # For the arm that swaps the aligner and changes nothing else. Only the
        # crash guard goes in, because without it that arm cannot finish; the
        # coordinate arithmetic is left exactly as published, which is the point
        # of running it.
        s = replace_once(s,
                         "\t$hits_ed_frac{$hit}=$edit_count/$mm_count;",
                         "\t$hits_ed_frac{$hit} = $mm_count ? $edit_count/$mm_count : 0;",
                         "edit fraction division")
        open(dst, "w").write(s)
        print(f"patched {src} -> {dst} (crash guard only)")
        return

    s = replace_once(s, SPLIT_OLD, SPLIT_NEW, "hit-line split", 2)
    s = replace_once(s, DECL_OLD, DECL_NEW, "hit-line declaration")
    s = replace_once(s, REBUILD_OLD, REBUILD_NEW, "hit-line rebuild")
    s = replace_once(s, CLUSTER_OLD, CLUSTER_NEW, "cluster coordinates")
    s = replace_once(s, SITE_OLD, SITE_NEW, "per-site coordinate")

    # HE-10: mm_count is zero once a read can be aligned without every base of
    # it being compared. See the note in run.sh.
    s = replace_once(s,
                     "\t$hits_ed_frac{$hit}=$edit_count/$mm_count;",
                     "\t$hits_ed_frac{$hit} = $mm_count ? $edit_count/$mm_count : 0;",
                     "edit fraction division")

    # the subroutine goes after the last `use` line, before any code runs
    marker = "use strict;\n"
    s = replace_once(s, marker, marker + REF_OFFSET_SUB, "use strict")

    open(dst, "w").write(s)
    print(f"patched {src} -> {dst}")


if __name__ == "__main__":
    main()
