# The cross-species re-test

Porath, Carmi & Levanon (2014) reported human A-to-I hyper-editing roughly
tenfold above mouse, rat and opossum, measured as hyper-edited sites per mapped
read. Two things have happened since.

The authors themselves revised it. Porath et al. (2017, *Genome Biol* 18:185)
screened 21 metazoans, switched the metric to events per million mapped bases,
standardised reads to 70–80 bp, and concluded that human is **not** exceptional.
So the tenfold claim is already historical, and re-testing it as stated would be
pushing on an open door.

What neither paper controls is **how well reads could be placed in the first
place**. Every species in both analyses went through the same
`bwa aln -n 2 -o 0 -N`: an edit-distance budget of two, no gaps, end to end.
After the A→G transform an editing event becomes a match, so that budget is
spent entirely on non-editing differences — SNPs, sequencing error, divergence
from the reference. Detection sensitivity is therefore a function of how far a
sample sits from its reference, which is exactly what differs across species.
`-N` is worth noting here: it is bwa's *non-iterative* mode, searching for all
hits at the limit, so the published pipeline is already searching as hard as
`bwa aln` can. The constraint is the budget, not the search.

And the references have improved by wildly unequal amounts:

| species | 2014-era contig N50 | today | gain |
|---|---:|---:|---:|
| rat (rn5 → GRCr8) | 52 kb | 64.3 Mb | 1,224× |
| opossum (MonDom5 → mMonDom1.pri) | 108 kb | 3.91 Mb | 36× |
| mouse (mm9 → GRCm39) | 25.6 Mb | 59.5 Mb | 2.3× |
| human (hg19 → GRCh38) | 38.4 Mb | 57.9 Mb | 1.5× |
| fly (dm3 → dm6) | 21,485,538 bp | 21,485,538 bp | 1.00× |

Fly is a free negative control: its assembly is, for practical purposes, the one
the 2014 paper used.

## What this runs

A 2×2 per species — reference quality crossed with aligner sensitivity:

|                | `bwa aln` (published) | HISAT2 |
|----------------|-----------------------|--------|
| old assembly   | **A** = Porath 2014   | **B**  |
| new assembly   | **C**                 | **D**  |

Cell A reproduces the 2014 result rather than resembling it: same assembly, same
pipeline at the same commit, and for opossum the same reads (Brawand 2011
SRR306743, the source Porath 2014 cites).

- **A→C** isolates sixteen years of assembly improvement, published pipeline untouched.
- **A→B** isolates the aligner, on the 2014 reference.
- **C→D** is the one that decides whether the aligner work still matters. If a
  modern assembly absorbs what the aligner swap buys, they are substitutes; if
  the gains are additive, both confounders are real and both survive into 2017.

## Running it

### On the OHSU ARC cluster

ARC has **no module for `bwa`, `hisat2` or `samtools`** — `module avail` offers
`bowtie2`, `bamtools`, `bedtools2` and `kallisto`, but not these — so the
toolchain is built with mamba rather than loaded.

`arc/site.sh` holds the lab's paths and slurm account, so they need not be
retyped; `submit.sh` reads it automatically. Edit it for a different lab, or
delete it and pass `-A` and `--env` explicitly.

```bash
source arc/site.sh

# 1. Build the environment once (~1.5 GB), on gscratch where every node sees it.
arc/setup_env.sh "$XSPEC_ENV"

# 2. Submit. Account and environment come from site.sh.
arc/submit.sh "$XSPEC_WORK" opossum
```

which resolves to `/home/exacloud/gscratch/YardimciLab/skubi/xspec` for the
work and `.../skubi/envs/hyperedit` for the tools, billed to `YardimciLab`.
Confirm the account name with `sshare -U -u $USER` — it is the ARC project,
which is usually but not always the lab's name.

The jobs download their own references and reads — ARC's compute nodes have
network access. `arc/fetch.sh "$XSPEC_WORK" opossum` pre-stages them, which is
optional but returns the core-hours a 16-core allocation would otherwise spend
idle during several GB of transfer.

Concurrent jobs sharing a work directory are safe: the git clone and the read
download are guarded by an atomic-`mkdir` lock, with `flock` ruled out because
gscratch does not support POSIX locking.

`--dry-run` prints the `sbatch` commands without submitting. `--modules` is
still there for anything the environment does not carry, but is empty by
default.

`environment.yml` pins bwa 0.7.18, hisat2 2.2.1 and samtools 1.21 — the exact
versions the chr21 and chr16–22 runs used, so a cluster result and a laptop
result are comparable rather than merely similar.

**Where the data lives.** `WORKDIR` must be on gscratch: it has to survive
between jobs and be visible from every compute node. Not `/home/users`, which
the ARC storage guide says outright is not for data processing; not
`/mnt/scratch`, which is node-local and wiped; `/arc/scratch1` works and is
unquota'd but is a technology preview that "will eventually be deleted", which
is a poor home for a week of indexing. `arc/fetch.sh` refuses the first two and
checks free space before starting.

**Why a chain of jobs rather than one.** ARC's `batch` partition caps at 36
hours and the opossum pair needs roughly forty. Splitting per assembly puts each
job at ~20 h, inside the default limit, so no QOS is needed — and short jobs
start sooner under ARC's backfill scheduler. The two assembly jobs are
independent and will run concurrently if the cluster has room; only the scoring
job is dependency-ordered, with `afterany` so it reports whatever finished.

If you would rather run it as one job, `--qos long_jobs --time=4-0` (10-day
cap) does that. `very_long_jobs` reaches 30 days but is capped at 24 CPUs.

**Node-local scratch.** Stage 2 rewrites the whole unmapped pool twelve times
per arm, as SAM and again as BAM, and only the beds are kept. That churn goes on
`/mnt/scratch`, requested via `--gres disk:400` and set up and torn down with
ACC's `mkdir-scratch.sh`/`rmdir-scratch.sh`; teardown is on a trap so it runs
whether or not the job succeeds. The indexes — the expensive, twenty-hour part —
stay on gscratch so a resubmission resumes instead of rebuilding.

**Resubmitting is the intended way to handle a timeout.** Every stage is guarded
by a marker file written only after its output is complete, so the same command
picks up where it stopped.

### Anywhere else

`run_xspec.sh` is plain bash and knows nothing about Slurm:

```bash
THREADS=16 ./run_xspec.sh /scratch/$USER/xspec opossum          # both assemblies
THREADS=16 ./run_xspec.sh /scratch/$USER/xspec opossum monDom5  # just one
THREADS=16 ./run_xspec.sh /scratch/$USER/xspec opossum score    # re-score
```

Species and assemblies come from `assemblies.tsv`; read accessions from
`reads.tsv`, whose FTP paths are resolved from ENA at fetch time. The rat, mouse
and fly rows are present but commented out — uncomment to extend the panel.

**Cost, per species pair:** ~170 GB persistent plus ~400 GB node-local scratch,
16 cores, 64 GB RAM, and on the order of 30–40 h. The fourteen indexes per
assembly are built up front so the six transforms can index concurrently; on a
disk-bound machine this is the wrong shape and the transforms should be
serialised instead (~25 GB peak, roughly six times the wall clock).

## Reading the output

`score_xspec.py` prints every term of the rate separately, and this is
deliberate. A better assembly places more reads in stage 1, raising the
*denominator*, so a modern reference can lower sites-per-read while finding
strictly more sites. A more sensitive aligner raises the *numerator* without
touching stage 1, which is held at `bwa aln`/`bwa mem` in all four cells. Report
the ratio alone and the two effects can cancel invisibly.

Specificity travels with every number, since a cell that finds more sites by
aligning more badly is not a better cell. Both checks are truth-set-free: the
control channels (of six mismatch classes only A2G can carry an A-to-I event, so
the other five are a noise floor by construction) and the ADAR motif (G depleted
at −1, enriched at +1, against ~21% background).

## Caveats to carry

- **The junction lists cannot be equal across the reference axis.** MonDom5's
  only gene set is a 7 MB Ensembl GTF; mMonDom1 has NCBI RefSeq RS_2026_06.
  Union with the sample's own first pass mitigates it — the first pass is
  read-derived — but opossum is the hardest case for the HISAT2 cells for
  exactly this reason. If cell D underperforms, check the list before believing
  the result.
- **Opossum yields little in absolute terms** in every published survey. The
  result is in the ratios between cells, not the absolute counts.
- **The caller is held fixed** across the aligner axis, unlike `run_real.sh`,
  which pairs `bwa aln` with the patched Perl caller and HISAT2 with the alnbase
  caller. That pairing is right when the caller is itself under test and would
  confound the only contrast this experiment is about.
- **Strandedness of the Cardoso-Moreira runs is inferred from the kit name, not
  measured.** Confirm before trusting A>G against T>C on those. The Brawand
  reads used by the opossum pilot are unstranded, as they were in 2014.
- **Take FASTA from UCSC, never NCBI.** NCBI's `*_genomic.fna.gz` is soft-masked
  with WindowMasker rather than RepeatMasker, and ships no RepeatMasker product
  at all for opossum or fly.
- **UCSC's `monDom5` browser and Ensembl's ASM229v1 are both still the 2006
  genome.** mMonDom1.pri is reachable only via GenArk or NCBI Datasets, so the
  path of least resistance silently re-runs the 2014 experiment on the 2014
  reference. `assemblies.tsv` points at GenArk for this reason, and
  `run_xspec.sh` hard-fails if a GTF and its FASTA share no sequence names.
