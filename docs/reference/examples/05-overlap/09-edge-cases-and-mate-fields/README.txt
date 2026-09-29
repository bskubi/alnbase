09-edge-cases-and-mate-fields: the four cases drawn in the reference section "Edge cases and
how they are resolved", with the mate fields checked against samtools fixmate.

(a) pair_overlap            ordinary overlapping mates. R1 kept; R2 becomes 20S20M at chr1:81,
                            and R1's PNEXT and MC follow it.
(b) hic_sup_overlap         R1 split into primary + supplementary; R2 overlaps only the
                            supplementary, which is clipped to 60H5M15S. The primary's SA is
                            rebuilt to list the clipped supplementary.
(c) primary_inside_overlap  R1's primary has nothing left. Its supplementary's alignment moves
                            onto the primary record (20M60S at chr2:141, whole 80 bp read), and
                            no supplementary or unmapped R1 line remains.
(d) read_through            dovetailing mates: the molecule is shorter than the reads, so the
                            losing read gives up every aligned base and is written unmapped at
                            its mate's position; the mate gets 0x8 and MC:Z:*.

Each run is followed by samtools fixmate -m on alnbase's output; it changes no FLAG, RNAME,
POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN or MC, and flagstat counts one primary line per read.
TLEN is measured 5' end to 5' end, as fixmate measures it.
