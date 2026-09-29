Pads, the widest query, and keeping a query's hits fixed.

One forward read, 0-based 12..31. Reference CGs at 7-8 and 35-36 lie wholly outside
the read; the CG at 11-12 straddles its first base.

The walk adds pad columns past each read end (the read side is a pad, the reference
side is the real reference). There are --end-context of them, by default the widest
query span in the run minus 1. A query fires on every window it matches, including
windows made only of pads.

CG_in_reference (read "~~", refr "CG") accepts pads on the read side, so:
  1. alone (1 pad per end) it fires only at the straddling CG, 11;
  2. beside an unrelated 6-column query (5 pads per end) it also fires at 7 and 35,
     whose windows are all pads;
  3. with --end-context 1 fixed, the run shows 1 pad per end whatever else it holds,
     and the query fires at 11 only;
  4. with "cg_ref and not all_pad" the query itself refuses all-pad windows, and fires
     at 11 only whatever the pad count.

Queries whose read row needs a read base somewhere (for example "C~" over "CG") can
never match an all-pad window, so they are unaffected either way.

alnbase 0.1.3 to 0.1.6 instead refused to fire on windows with no observed column;
0.1.7 restored firing on any match and added --end-context.
