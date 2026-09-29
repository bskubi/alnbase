# MethylDackel bugs, reproduced

Five MethylDackel behaviours that contradict its own help, usage text or README, each
shown on a few hand-written reads. The candidates come from reading the source
(`docs/design/research/methylation-callers-hardcoded.md`, IDs MD-8, MD-14, PR-3, PR-4,
PR-5); this directory confirms each one by running the tool.

Tool: MethylDackel 0.6.1 built from GitHub master `3c77bda` (`../../envs/build_methyldackel.sh`),
HTSlib 1.24. Checked on 2026-09-17 with alnbase 0.1.10.

```
METHYLDACKEL=... SAMTOOLS=... PY=python-with-duckdb ALNBASE=... ./run.sh [workdir]
```

`expected.txt` is the output of that run. Every case prints what MethylDackel promises,
then what it does. Where the behaviour is something alnbase computes, alnbase's answer
for the same reads follows, from the mirror query `../../../presets/methyldackel/queries.toml`.

| ID | Command | Promised | Observed | Cause |
|---|---|---|---|---|
| MD-8 | `extract --OT A,B,C,D` | calls at 1-based read positions A through B ("`--OT 5,0,0,0` would include all but the first 4 bases") | position A is excluded: `--OT 3,0,0,0` keeps Cs at 5, 7, 9, not 3, 5, 7, 9. The right bound is inclusive, as documented | `trimAlignment` masks SEQ indices `0..A-1` without converting A to 0-based (`common.c:137-172`) |
| MD-14 | `extract` bedGraph | column 4 is "the methylation percentage rounded to an integer" | 2 of 3 methylated gives `66` | `(int)(100.0*M/(M+U))` truncates (`extract.c:50`) |
| PR-3 | `perRead` | "if an NH tag is present and its value is >1 then an entry is ignored"; `--ignoreNH` turns this off | both `NH:i:2` records are written; `--ignoreNH` is `unrecognized option` | no NH check exists and the option is not registered (`perRead.c:185-193`, `:296-303`) |
| PR-4 | `perRead -p 5` | "Minimum Phred threshold to include a base" | a CG C at Phred 2 is counted when the base before it is also below Q5; a soft-clipped C is counted when the last aligned base is below Q5 | after skipping a low-quality base the loop does not `continue`, so the next base is scored with no quality or CIGAR check (`perRead.c:57-63`) |
| PR-5 | `perRead -l BED` | "A BED file listing regions for inclusion" | a read 900 kb outside the only interval is written | the BED is consulted per 1 Mb chunk, not per read (`perRead.c:158-172`) |

## What this means for comparisons

- **MD-8** matters most. MBias suggests the left bound as the first position to keep, the
  help's convention (`lthresh = i+2` after the last biased 0-based bin `i`, `svg.c:281`),
  so a user who passes its suggestion to `extract` loses one unbiased base on the 5' side. (MBias
  and `--OT` also count positions along SEQ as stored, not from the read's 5' end; that
  is MD-9, a documented-by-omission choice rather than a bug, and is not reproduced
  here.) alnbase's equivalent filter is written on `off_5p`, which counts from the read's
  5' end as sequenced: 1-based position A is `off_5p >= A-1`, and the run shows the four
  Cs it keeps.
- **PR-4** inflates `perRead` methylation from low-quality and soft-clipped bases. In
  alnbase a soft-clipped base is a clip column, which no base pattern matches, and
  quality is a per-hit column filtered downstream; both reads get 0 CG calls at Q ≥ 5.
- **PR-3** and **PR-5** are filters that belong upstream of alnbase (`samtools view -e
  '[NH]<=1'`, `samtools view -L`), so there is no alnbase counterpart. A comparison with
  `perRead` output must apply them to its input instead of relying on the options.
- **MD-14** changes only the rounded percentage column; the counts in columns 5 and 6
  are right.
