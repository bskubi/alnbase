"""The same run as `single_cell_methylation.sql`, built by chained calls.

    python examples/single_cell_methylation.py \
        --hits 'calls_*_*.parquet' --fai genome.fa.fai --sample pbmc

Run it with `--show` to print the script it would execute and stop. That script
is an ordinary one: `alnbase-agg run` takes it, so this file is a way of
writing the run rather than a second way of performing it.

Read it top to bottom. The chain is the pipeline in order: what the query names
mean, where the hits come from, what is thrown away, what is counted, and what
is written.
"""

import argparse

from alnbase_agg.pipeline import Measure, Run, both, listed

# WHAT THE QUERY NAMES MEAN
#
# The only place biology appears. A measure sorts alnbase query names into
# named categories; you choose the names. Everything downstream is arithmetic
# over these, so an assay this toolchain has never heard of is a different set
# of categories here and no change anywhere else.
#
# A measure may have any number of categories -- a genotype has one per base.
# What a filter tests and a file writes is an interpretation of a *relation*
# between two sets of them: `main` is the side in focus, `other` what it is set
# against, `total` is both together and `ratio` is `main` over `total`. A
# measure with exactly two categories has one possible relation, so the common
# case below never has to name it.
mCG = Measure("mCG", methylated="CG", unmethylated="TG")
mCH = Measure("mCH", methylated=["CHG", "CHH"], unmethylated=["THG", "THH"])


def build(hits, fai, sample, exclude=None, min_block_depth=1000):
    # 1. THE HITS, AND THE HIT-LEVEL FILTER
    #
    # `read` defaults to qname plus mate, which is the read and not the
    # alignment record: a chimeric Methyl-HiC read is several records, and
    # judging them apart would filter half a molecule. `block` is whatever the
    # outputs are split by -- a cell barcode here, a spot in a spatial run.
    run = Run(hits, block="XB", fai=fai)

    # Each `where` is AND'd with the others, and conditions inside one are
    # OR'd. They run before anything is counted, which is where per-position
    # artefacts go: end repair sits a fixed distance from the 3' end and random
    # priming the same distance from the 5' end, both measured on the read as
    # sequenced, including hard-clipped bases.
    run.where("qual >= 20")
    run.where("mapq >= 30")
    run.where("off_5p >= 10")
    run.where("off_3p >= 2")

    # 2. WHICH OBSERVATIONS SURVIVE
    #
    # A list of reasons to throw a read away.
    run.drop_reads(
        # Incomplete conversion: a read with too much CH methylation did not
        # convert, so its CG calls cannot be trusted either.
        mCH.ratio.gt(0.4),
        # A patch of adjacent CH calls, which a rate cannot see: four in a row
        # and four scattered give the same rate, and only the first looks like
        # a conversion failure.
        mCH.consecutive.gt(3),
        # premethyst's ladder, where the threshold depends on how much
        # evidence there is. Each rung is a conjunction, and the rungs are
        # unioned like any other two reasons. Nothing checks that they tile the
        # range: an overlap is harmless, and a gap quietly keeps the reads it
        # should have judged.
        both(mCH.total.ge(3), mCH.total.le(5), mCH.main.gt(1)),
        both(mCH.total.gt(5), mCH.ratio.ge(0.4)),
    )
    if exclude:
        # Reads named in a file. Regenerating the file and re-running changes
        # the answer without touching this script.
        run.drop_reads(listed(exclude))

    # And one level up: blocks with too little data to be worth keeping.
    run.drop_blocks(mCG.total.lt(min_block_depth))

    # 3. WHAT THE RUN THREW AWAY
    #
    # Rendered from the filters above, so the report and the filtering cannot
    # drift apart. `alnbase-agg run` prints the result of any statement that
    # returns rows.
    run.report()

    # 4. FOLDING THE TWO STRANDS OF EACH CG TOGETHER
    #
    # A CG is two cytosines, one on each strand, one base apart. Adding them
    # together is a choice, not a fact -- a run studying hemimethylation must
    # not do it -- so it is written down. `span` is how far apart the pair
    # sits, and 1 is a CG.
    run.merge_strands(merge=mCG)

    # 5. THE FILES
    #
    # One call per output, each naming everything it depends on. The measure is
    # the object declared at the top, so a typo is a Python error here rather
    # than an empty column in a file.
    #
    # Counting the hits into sites, and reading the contig order out of the
    # FASTA index, are not written here: each format names what it needs on a
    # REQUIRES line in its own header, and the run reads those first and once.
    run.bismark_cov(mCG, out=f"{sample}.CG.cov")
    run.bismark_cov(mCH, out=f"{sample}.CH.cov")
    run.cellinfo(cg_relation=mCG, ch_relation=mCH, out=f"{sample}.cellInfo.txt")

    # The per-block HDF5. The only output here that needs Python, and only
    # because the file is a dataset per block rather than one table. Called once
    # per context, because an Amethyst file holds a group for each; the
    # `.python` half is added automatically, being named after the `.sql`.
    run.amethyst_h5(mCG, context="CG", out=f"{sample}.h5")
    run.amethyst_h5(mCH, context="CH", out=f"{sample}.h5")

    return run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--exclude")
    parser.add_argument("--min-block-depth", type=int, default=1000)
    parser.add_argument("--show", action="store_true", help="print the script, run nothing")
    args = parser.parse_args()

    run = build(args.hits, args.fai, args.sample, args.exclude, args.min_block_depth)
    if args.show:
        print(run.sql())
    else:
        run.run()


if __name__ == "__main__":
    main()
