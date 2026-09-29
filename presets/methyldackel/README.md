# MethylDackel

`MethylDackel extract`, commit `3c77bda` (GitHub master; the bioconda release has
bugs this one fixes, so validation builds from source). Sources cited below are
that tree. Background: `alnbase-validation/docs/research/methylation-tools.md` §2.

Run it with `--insertions skip`, which is alnbase's default: MethylDackel walks a
pileup and never sees an inserted base, and skipping keeps the reference bases on
either side of an insertion adjacent, which is what its context rule reads.

## What counts as a call — `queries.toml`

MethylDackel pileups each position and, at a reference C (a reference G for reads
from the bottom strand), counts a read C as methylated and a read T as
unmethylated. The context comes from the **reference alone**, never from the read:

- next reference base `G` → CpG;
- otherwise the base after that decides: `G` → CHG, anything else → CHH;
- an `N` counts as "not G", and so does running off the end of the contig, so a C
  in the last base of a contig is a CHH;
- a deletion or an insertion in the read does not change the context.

alnbase's walk orients every read along its conversion strand, so a single pattern
pair (`read = "C~~"` / `refr = "C~~"`, and the same with `read = "T~~"`) covers
both strands, and each cytosine is reported at its own coordinate.

## What it writes — `run.sql`

One 0-based bedGraph per context, CpG only unless `--CHG` or `--CHH` is given, six
columns (`extract.c:38-52`):

```
chrom   start   end   percent   nMethylated   nUnmethylated
```

Two details of that line are reproduced exactly, because both are invisible when
wrong:

- **The percentage truncates.** It is a C cast, `(int)(100.0 * m / (m + u))`
  (`extract.c:46-51`), so two thirds methylated is `66` and 99.6% is `99`. Bismark,
  on the same counts, writes `66.6666666666667`.
- **The coordinates are 0-based half-open**, where Bismark's coverage file is
  1-based with `start == end`. The same six fields, one apart.

MethylDackel also writes a `track type="bedGraph" description=...` line first
(`extract.c:562-569`). That is a UCSC display directive rather than data, so
`methyldackel-bedgraph` does not write it; a diff against MethylDackel's own output
skips its first line.

`-d`/`--minDepth` defaults to 1, so positions with no coverage are absent. The
recipe gets that for free: a writer here leaves out a position where the context
was not observed, because an uncovered cytosine is not an unmethylated one.

## Filters, which are options and not definitions

`[filter.hit]` in the recipe carries MethylDackel's defaults
(`common.c:415-430`, `common.c:127`): `-q 10`, `-F 0xF00`, singletons and
discordant pairs dropped, and `-p 5` on the cytosine's own base quality. They are
in the recipe rather than in `queries.toml` because they are flags — a query file
says what a methylation call *is*, and these say which alignments a run chose to
believe. Comparing against a MethylDackel run with other thresholds means editing
this table to match, or turning them off on both sides.

`--ignoreNH` is not applied: Bismark writes no NH tag.

## Mate overlap

MethylDackel counts both mates independently in the overlapping region, so a
position covered by both contributes twice to its own pileup. The recipe's
`[unit] by = "fragment"` collapses that to one call. The two rules disagree wherever
mates overlap, which is why the comparison runs on fragments without overlap; see
`presets/README.md`.

## Known bugs

Confirmed by execution, with reproductions in
`alnbase-validation/demos/methyldackel-bugs/`. That directory's `run.sh` reads this
preset's `queries.toml` as the mirror it compares against.
