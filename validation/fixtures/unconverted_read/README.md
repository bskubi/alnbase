# A read that retains every cytosine, and its converted twin

`chrU` is `ACG` six times over. `unconv` matches the reference, so all six cytosines read
as C; `conv` reads T at every one. Both produce six calls, one per CG, differing only in
whether the call is `CG_met` or `CG_unmet`.

A fully unconverted read is the signature of a molecule the bisulfite treatment missed,
and premethyst and ScaleMethyl both discard reads by their retention rate. alnbase makes
the calls and leaves that decision to a `group by` over these rows, which is why a read
that is methylated at 6 of 6 sites appears here rather than disappearing.
