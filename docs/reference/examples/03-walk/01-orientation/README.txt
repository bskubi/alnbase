Which way a record is walked, and what off_5p / off_3p / qual / read_base /
refr_base mean on each strand, for the six FLAG orientations plus secondary
and supplementary records. See reads.sam @CO lines for each record's purpose.

Key points shown in expected.txt:
- flags 0, 65, 145, 256, 2048: strand '+', refr_pos ascending.
- flags 16, 81, 129: strand '-', refr_pos descending, qual reversed with the
  bases, read_base AND refr_base complemented.
- off_5p / off_3p are measured from the ends of the read as sequenced
  (alnbase 0.1.2), whatever the walk strand: off_5p = 0 is the first base
  the sequencer read for every FLAG. For flag 129 (R2 forward) that is the
  leftmost base although the walk runs right to left; for flag 145 (R2
  reverse) the rightmost, although the walk runs left to right. So for read 2
  the walk emits columns from its sequenced 3' end, and off_5p counts down
  along the walk.
- The last table shows off_5p running 0..7 from the first sequenced base in
  every case: it is the sequencing cycle, with no FLAG-dependent correction.
