# Bismark

Bismark v0.25.1. The call is made at alignment time, in `bismark`
(`extract_corresponding_genomic_sequence_single_end`, `methylation_call`), and
written into the read's `XM` tag; `bismark_methylation_extractor` then reads that
tag rather than the reference. Background:
`docs/design/research/methylation-tools.md` §1.2.

Run it with `--insertions emit`. Bismark's context is read-projected, so it sees
inserted bases, and an inserted base is what makes a context unknown (rule 3
below).

## What counts as a call — `queries.toml`

Generated, not written by hand:

```sh
python generate.py > queries.toml
```

The file is 300-odd patterns because Bismark's context rule has to be enumerated:
each of the two context positions may be preceded by up to `MAX_DELETION`
deletion columns, and every combination is its own pattern. Generating it keeps the
rule in one readable place and the file exact; `generate.py`'s header states the
rule, and `queries.toml`'s header repeats it so the file can be read on its own.

The rule, in short:

1. The context is the next two characters of the **read-projected** reference — the
   reference bases the read aligns to, in read order, plus two reference bases past
   the end of the alignment.
2. A deleted reference base is skipped, not used as context.
3. An inserted or soft-clipped read base contributes an unknown character, and the
   call becomes `U`/`u`.
4. Next character `G` → CG (`Z`/`z`); `N` or unknown → `U`/`u`; otherwise the
   character after it decides: `G` → CHG (`X`/`x`), `N` or unknown → `U`/`u`,
   anything else → CHH (`H`/`h`).

Uppercase is methylated, lowercase unconverted, and the query names are those
letters, so `[tag.xm.bases]` reassembles Bismark's XM string from the hits — which
is the direct check that the preset is Bismark's rule and not an approximation.

One deliberate difference: an IUPAC ambiguity code in the reference is treated as
unknown here, where Bismark would use it as an H. Deletions longer than
`MAX_DELETION` are not covered.

## What it writes — `run.sql`

`.bismark.cov.gz`, via `bismark2bedGraph` and `coverage2cytosine`: one line per
cytosine, **1-based**, `start == end`, and the **methylated count first**:

```
chrom   start   end   percent   count methylated   count unmethylated
```

Two things about that line are reproduced exactly:

- **The column order.** Methylated before unmethylated. Swapping the last two
  produces a file every downstream tool reads happily and wrongly.
- **The percentage.** It is Perl's default stringification of `100 * m / (m + u)`,
  which is fifteen significant digits: `66.6666666666667`, not the seventeen a
  double prints by default, and not MethylDackel's truncated `66`.

Both strands appear as their own rows. `bismark2bedGraph` does not merge the two
cytosines of a CG; `coverage2cytosine --merge_CpG` is the separate step that does,
and strand merging is a downstream operation here too.

CG only, unless `bismark_methylation_extractor` is given `--CX`; the recipe has the
other two contexts commented out, which is that flag.

`U`/`u` calls — context unknown — reach no `.cov` file: Bismark writes them into
the XM string and its extractor discards them. They are declared in `queries.toml`
so the tag round-trips and are deliberately not a context in the recipe.

## Filters

None. `bismark_methylation_extractor` applies no MAPQ or base-quality threshold of
its own; filtering happens before it, in the aligner and in
`deduplicate_bismark`. So the recipe has no `[filter.hit]` — adding one would not
be reproducing Bismark.

`--ignore` and `--ignore_r2` (both default 0) trim positions from a read's ends.
Non-zero values are a `[filter.hit]` over `off_5p`, which is why alnbase defines
those offsets against the read as sequenced, including hard-clipped bases.

## Mate overlap

`bismark_methylation_extractor` defaults to `--no_overlap` for paired-end data:
where two mates overlap, the read-2 portion is ignored entirely. The recipe's
`[unit] by = "fragment"` instead collapses a doubly-covered position to one call and
drops it if the mates disagree. The two rules disagree wherever mates overlap, so
the comparison runs on fragments without overlap; see `presets/README.md`.
