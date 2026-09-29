# A CG whose partner base is outside the alignment

`chrR` is 40 bp of A with a CG at 20-21 (1-based). `end_on_c` stops on the C, so the G it
needs is one base past its last aligned base; `start_on_g` begins on the G, and on the
bottom strand that G is the cytosine, whose context lies one base before its first aligned
base. `spans` and `spans_rev` cover both bases and are the controls.

The context has to come from the reference past the end of the alignment, which alnbase
supplies as pad columns. A caller that reads context out of the read sequence alone has
nothing to read; Bismark's rule reaches two reference bases past the alignment for exactly
this case, and MethylDackel works from the reference throughout. Both calls land, and
both land on the same CG from opposite strands at different coordinates -- merging that
pair into one site is the aggregation step's job, not this one.
