# Simulations

Every figure the validation suite reports is a property of two things: the tools it ran and
the data it ran them on. The second half used to be implicit — a BSReadSim command typed
once at a shell prompt, recoverable only from a manifest in a scratch directory. It is now
`profiles.toml`, one named profile per simulation, and `../bin/simulate.py`, which runs a
profile and writes a `provenance.json` beside the output saying which one produced it.

```
validation/bin/simulate.py lambda-illumina --out $SCRATCH/sims/lambda-illumina --ref-dir $SCRATCH/val/ref
```

Read a profile's `models` field before quoting any number derived from it. That field is
the part that cannot be inferred from the parameters.

The reference is not in the repository; `simulate.py` prints the command to fetch it when
it is missing, and that command has been checked to reproduce the digest pinned in
`profiles.toml`:

```
curl -sL '<the url in profiles.toml>' | sed '1s/.*/>chr/' > ref/lambda.fa
```

## The two profiles

`lambda-stress` is the simulation behind the validation suite's equivalence results.
`lambda-illumina` exists because several of the other
numbers there are rates a reader would take as descriptive of the tools, measured on input
that is not descriptive of real data; every rate is now re-reported from it.

| | `lambda-stress` | `lambda-illumina` |
|---|---:|---:|
| records | 29,102 | 38,802 |
| records with an indel in CIGAR | 16,678 (57.3%) | 374 (1.0%) |
| simulated SNVs / indels | 368 / 619 | 43 / 5 |
| cytosine observations | 721,204 | 1,451,568 |
| CG methylation | 0.501 | 0.728 |
| CH methylation | 0.171 | 0.009 |
| conversion efficiency | 0.998 | 0.995 |
| base error rate | 0.0050 | 0.0010 |
| Phred written | 40 | 30 |
| read length / insert | 100 / 300 | 150 / 350 |

Measured from BSReadSim's own `zt` truth via `../bin/truth_observations.py`, not from the
parameters, which is how the profiles' claims about themselves are checked. Conversion
efficiency is measured over unmethylated CGs that carry neither a sequencing error nor a
variant: in `lambda-stress` the raw figure is 0.996 rather than 0.998, because a 0.5% error
rate on bases whose Phred says 0.0001 is large enough to move it.

**Which profile supports which claim:**

- **Equivalence** — "alnbase reproduces Bismark's XM call for call", "alnbase and
  MethylDackel agree at all 24,056 sites" — belongs to `lambda-stress`. Two callers fed the
  same deliberately extreme input either agree record for record or they do not, and the
  extremity is what makes the comparison searching: 57% of records carry an indel.
- **Rates** — per-call disagreement percentages, what a default filter costs, precision and
  recall — belong to `lambda-illumina`. A rate measured on `lambda-stress` is a rate for a
  genome with one variant every 50 bases.

## What lambda cannot answer

**Nothing measured here can establish what a MAPQ filter costs.** Lambda is 48 kb with no
repeats, so no read is ever ambiguously placed and low MAPQ can only come from mismatch
load. That is why MethylDackel's default MAPQ 10 discards 36% of observations on
`lambda-stress`: at one variant per 50 bases, reads carry enough mismatches to look badly
aligned. On `lambda-illumina` the same filter discards **nothing** — 0 of 1,066,119 calls,
every read uniquely mapped at MAPQ 40 — and that is just as unrepresentative in the other
direction, because on a real genome the filter mostly removes reads in repeats, which lambda
does not have. The two profiles bracket the answer at 36% and 0%, which is another way of
saying neither measures it. The honest measurement needs a reference with repeats in it,
which is the reason the suite wants a few-Mb human profile; until then, the 36% figure
should be reported as a property of `lambda-stress` or not reported at all.

**Indel-adjacent behaviour cannot be given a realistic rate here either.** A human-like
variant rate on 48 kb produces 5 indels. That is enough to check that a caller still
handles them, and nowhere near enough to say how often the handling matters. The division
of labour is deliberate: `lambda-stress` says *what* each caller does around an indel,
`lambda-illumina` says how often the rest of the calls disagree, and neither says how often
indel handling matters on a real genome.

## Reproducibility

Re-running a profile reproduces its reads exactly: `lambda-stress` regenerated the
simulation behind the existing results with every record identical, on a different day in a
different directory. Two caveats, both discovered by doing it:

- **The BAMs are not byte-identical.** BSReadSim stamps each run with a fresh UUID read
  group, so `RG:Z:` differs and nothing else does. Compare with `RG` stripped.
- **BSReadSim's `configuration_sha256` covers the output directory path**, so the same
  simulation run in two places gets two hashes and the hash cannot answer "is this the
  input those numbers came from?". `provenance.json` therefore carries its own
  `simulation_sha256` over the reference sequence, the technology and the parameters, and
  nothing else. That is the field to compare across machines.

`simulate.py` also verifies the reference before running: it digests the FASTA's bases and
refuses to start unless they match the `sequence_sha256` pinned in `profiles.toml`. Pointing
a profile at the wrong copy of a reference would invalidate everything measured from it
without producing any error, which is the kind of failure the suite exists to catch in other
people's tools.

## Adding a profile

Add a `[profile.NAME]` table with `reference`, `technology`, `summary`, `models` and
`supports`, and a `[profile.NAME.args]` table of BSReadSim flags without their dashes.
Parameters left out take BSReadSim's defaults, which the written manifest records in full.

`models` should say what real thing the profile stands for and where it knowingly departs
from it; `supports` should say which claims may be quoted from it. A profile whose
parameters are self-inconsistent — `lambda-stress` writes Phred 40 on bases that are wrong
1 time in 200 — must say so, because that inconsistency silently disables every downstream
base-quality filter and makes those filters look free.
