# ScaleMethyl bugs, reproduced

Five ScaleMethyl behaviours that give wrong calls, drop reads, or write exports their
documented readers cannot use. The candidates come from reading the source
(`docs/design/research/methylation-callers-hardcoded.md`, IDs SM-1, SM-2, SM-3, SM-4,
SM-14); this directory confirms each one by running the tool.

Tool: ScaleMethyl GitHub master `2a4973c`, `bin/met_extract.py`, `bin/write_allc.py` and
`bin/write_bismark.py`, run with the pipeline's own tool environment
(`envs/scaleMethylTools.conda.yml`, copied to `../../envs/scalemethyl.yaml`: polars
0.20.18, pysam 0.22.1, pyfaidx 0.8.1.4, duckdb 0.10.1, bwameth 0.2.7, bwa-mem2 2.2.1).
The exports are read by the tools they are meant for: ALLCools 1.1.1 (`c9f7be2`,
`../../envs/allcools.yaml`) and MethSCAn 1.1.0 (`72f38fb`, `../../envs/methscan.yaml`).
Checked on 2026-09-17 with alnbase 0.1.10.

```
SCALEMETHYL=... SCALEMETHYL_ENV=... ALLCOOLS_ENV=... METHSCAN_ENV=... \
  PY=python-with-duckdb ALNBASE=... ./run.sh [workdir]
```

`expected.txt` is the output of that run (reruns are identical). SM-1 uses simulated
paired-end reads aligned by bwa-meth exactly as the pipeline aligns them, so the FLAG and
`YD` values are bwa-meth's own; the other cases use hand-written reads. alnbase's calls for
the same reads come from the query file `cg_ch.toml` in this directory. `met_extract.py`
runs with `--aligner bwa-meth`, the pipeline default.

| ID | Command | Expected | Observed | Cause |
|---|---|---|---|---|
| SM-1 | `met_extract.py` | both mates of a directional fragment show the same conversion (bwa-meth tags both `YD:Z:f` for OT, `YD:Z:r` for OB). Cell with CG Cs retained and every other C converted: CG 100%, CH 0%. Cell with every C converted: CG 0%, CH 0% | `--threshold 1.0`: CH 50.03% in both cells, CG 49.92% in the unmethylated cell. At the default `--threshold 0.5`: `CH_high` = 50 in each cell (every second-in-pair mate dropped) and CG and CH coverage halved; the remaining calls match alnbase's on first-in-pair records at all 2534 positions | strand is `read.is_reverse` for every record (`met_extract.py:276`). Second-in-pair mates align in the opposite orientation to their fragment's strand, so they are tested at the other base of each pair (reference G for OT, C for OB), which is never converted and always reads as methylated (`:188-227`) |
| SM-2 | `met_extract.py --threshold 0.5` | docs: "Reads with greater than this percentage CH methylation will be discarded". Two reads with 1 of 10 CH methylated (10%) are both kept, with CG calls at 46 and 106 | the read whose methylated CH C comes first is dropped (`CH_high` 1) and its CG call at 46 is lost; the read whose methylated CH C comes last is kept | the ratio is tested after every call and the read is dropped at the first call where it exceeds the threshold (`:287-303`); the BSBolt path uses whole-read counts (`:102-106`) |
| SM-3 | `write_allc.py` | docs: ALLC columns "as described" by ALLCools, whose spec uses 1-based positions. The C at 0-based 106 is written at 107 | written at 106, a reference A. A tabix query for `chrC:107` (the pipeline indexes with `-b2 -e2`) returns nothing | parquet positions are 0-based (`read.reference_start + y`, `:225`) and the writer copies them unchanged (`write_allc.py:23-31`; also `write_bismark.py:23-29`) |
| SM-4 | `write_allc.py`, then `allcools extract-allc --mc_contexts CGN` | ALLC context is the reference 3-mer (`CGA`), so the CGN extract keeps the row | context is the literal `CG`; ALLCools writes 0 rows | `'{context}' as context` (`write_allc.py:27`) |
| SM-14 | `write_bismark.py`, then `methscan prepare` | docs: "bismark .cov format" with a `percent_methylated` column. Bismark's `.cov` is `chr start end percent M U`; MethSCAn, which `write_bismark.py` links to, reads M and U from columns 5 and 6 by default | `chrC 106 1.0 1 0`: no end column, a 0-1 fraction, 0-based. MethSCAn stops with `IndexError: list index out of range` | `ROUND(methylated / (methylated + unmethylated), 2)` in a 5-column `SELECT` (`write_bismark.py:23-29`) |

## What this means for comparisons

- **SM-1** affects every run through the default bwa-meth path (the parabricks path shares
  the code): the pipeline aligns the two reads as a pair (`modules/alignment.nf:44`). The
  `--aligner bsbolt` path takes strand from BSBolt's `YS` tag and is not affected.
  The default CH filter hides the wrong calls by discarding almost every second-in-pair mate,
  because each reference G (OT) or C (OB) it covers in CH context reads as methylated. What
  remains is roughly half the coverage, and a `CH_high` count that reports about half of all
  mates as having high CH methylation in every cell, whatever the cell's true CH level. Mates
  that pass the filter anyway, because they cover no CH position on the wrong base, would
  add 100%-methylated calls (follows from the code; not in this demo). With the filter relaxed (`chReadsThreshold` 100), CH methylation
  reads about 50% in every cell of the simulation. A comparison against ScaleMethyl should restrict
  alnbase to first-in-pair records; on the simulated reads that matches ScaleMethyl's default
  output at every position. alnbase takes each read's strand from its conversion tag, so both
  mates are called on their fragment's strand.
- **SM-2** biases what is kept by where methylated CH falls along the read: a read is
  dropped if the running ratio exceeds the threshold at any call, however low its whole-read
  ratio. Any read whose first CH call is methylated is dropped at the default 0.5. alnbase leaves per-read filters to
  queries over its hits: the whole-read CH count shown per read is a `group by`.
- **SM-3** shifts every exported ALLC and Bismark `.cov` row one base to the left, onto the
  base before the C for `+` rows. Region queries, joins with other tools' output and
  context lookups all miss by one. The parquet itself is 0-based without saying so.
- **SM-4** makes ScaleMethyl's ALLC files unusable with ALLCools' context patterns
  (`CGN`, `CHN`), which `generate-dataset` and `extract-allc` take.
- **SM-14** makes the `.cov` export unreadable by MethSCAn's default format; a custom
  `--input-format` (not tried here) could read it, but the positions are still 0-based (SM-3).
