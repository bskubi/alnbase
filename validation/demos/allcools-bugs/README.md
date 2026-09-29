# ALLCools bugs, reproduced

Three ALLCools behaviours that give wrong counts or drop rows, each shown on a few
hand-written reads. The candidates come from reading the source
(`docs/design/research/methylation-callers-hardcoded.md`, IDs AC-1, AC-4, AC-10); this
directory confirms each one by running the tool.

Tool: ALLCools 1.1.1, GitHub master `c9f7be2` (`../../envs/allcools.yaml`), with samtools
1.24. Checked on 2026-09-17 with alnbase 0.1.10.

```
ALLCOOLS_ENV=... PY=python-with-duckdb ALNBASE=... ./run.sh [workdir]
```

`expected.txt` is the output of that run. For AC-1, alnbase's calls for the same reads
follow, from the query file `cg.toml` in this directory.

| ID | Command | Expected | Observed | Cause |
|---|---|---|---|---|
| AC-1 | `bam-to-allc` | one read gives coverage 1 at a site | a read that **starts** at a CG C or G gets coverage 2 when its MAPQ is 13 or 51 (at the C) or 11 or 64 (at the G): MAPQ 13 adds a methylated count, MAPQ 51 an unmethylated one. MAPQ 12 is correct | samtools mpileup writes `^` and `chr(MAPQ+33)` before the base of a read that starts at the column. The parser removes only indel strings, then counts `.`/`T` or `,`/`a`, so the MAPQ character is counted as a base (`_bam_to_allc.py:207-250`, `:283-285`) |
| AC-4 | `bam-to-allc` | a `-` row for a covered G, as for a C on `+` | reference `CRG`: the `+` row at the C is written with context `CRG`, but the `-` row at the G is missing. Reference `CAG` gives both | the complement table has only `ACGTN`; the `KeyError` on `R` is caught by `except: continue` (`:178`, `:274-282`) |
| AC-10 | `extract-allc --strandness merge` | help: "merge the count on adjacent CpG in +/- strands"; every CG keeps one row, at its C on `+` | the last row of the file is lost (chrX 6); a lone `-` row that ends a chromosome stays `-` at the G (chrW 21) where other lone `-` rows move to `+` at the C | the loop holds one row back and writes it only when the next row arrives or the chromosome changes. There is no write after the loop, and the chromosome-change write skips the strand conversion (`_extract_allc.py:19-58`) |

## What this means for comparisons

- **AC-1** is a minor practical issue: silent and systematic, but its size depends on the
  aligner's MAPQ values and may be zero (not measured on a real BAM). It changes counts, not just rows. It needs a read to start exactly at the site
  with one of four MAPQ values, and all four pass the default `--min_mapq 10`. bwa's MAPQ
  (and so bwa-meth's) ranges over 0-60, which includes 11, 13 and 51. The source-reading
  inventory found that Bismark's single-end path can emit 11, but none of the 19,936
  records in our paired-end Bismark simulation has an affected value. hisat-3n emits only
  0, 1 and 60, so the YAP hisat-3n path avoids the bug. methylpy's `call-methylated-sites` uses the same
  parser, so it applies there too (not run here). A validation against ALLCools
  should expect disagreement at read-start positions with these MAPQs. alnbase counts
  aligned bases directly, never pileup text, and gives coverage 1 for every read.
- **AC-4** drops bottom-strand rows next to IUPAC codes. It matters only for references
  that keep IUPAC codes; alnbase's `index` keeps them.
- **AC-10** is aggregation, which alnbase leaves to its aggregation utility (strand merging
  is downstream). A comparison of merged ALLCs should check the last row of each
  chromosome.
