One record per CIGAR operation (M = X, D, I, S, H, P), top and bottom strand
for the ones whose coordinates depend on direction. Introns (N) are in
03-introns. The run is done twice: default policy, then with
--insertions emit. Soft-clipped bases have no policy: they are never aligned
columns (the --soft-clips option was removed in 0.1.2, and the run.sh shows
it is refused).

Things to read off expected.txt:
- D: one column per deleted reference base, read_base '.', qual NULL,
  off_5p = the offset of the read base emitted just before it.
- I skipped: no column; off_5p jumps (3 -> 6) so offsets stay read-based.
- I emitted: refr_base '.', qual real, refr_pos = the coordinate of the
  column emitted just before (the last aligned base in walk order).
- S: no aligned column, under every option (sclip_0 / sclip_16 are
  identical in both runs). The flank columns nearest the aligned part stand
  for the clipped bases: read_base ':', the real
  reference base, qual NULL, and the clipped base's own offset (sclip_0,
  2S6M3S: 0 1 before the first aligned base, 8 9 10 after the last; its first
  aligned base is off_5p 2). Past the clip the flank is pads, with no offsets
  ('.').
- lead_ins skipped: the first aligned base is off_5p 2.
- H: no column, but (since 0.1.4) offsets count hard-clipped bases, which
  were sequenced: hclip (3H6M2H) reports off_5p 3..8 for its aligned bases.
- P: no effect at all.
