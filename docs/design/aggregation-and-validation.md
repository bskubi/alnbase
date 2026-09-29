# Design proposal: aggregation/export, and validation against existing tools

Status: **phase 0, revision 2.** Nothing here is built yet. Revision 2 applies the
decisions from the first review (summarised in §0) and the corrections that came from
documenting alnbase against its code (`docs/reference/`). The research behind
it, with source-code citations for every tool default, is in
[`research/`](research/):

- [`methylation-tools.md`](research/methylation-tools.md): Bismark, MethylDackel,
  BISCUIT, ALLCools/methylpy, CGmapTools, bwa-meth; formats, implicit queries,
  a cross-tool comparison table (§9), suspected bugs (§10), install notes and test data.
- [`singlecell-formats.md`](research/singlecell-formats.md): Amethyst v2.0.0 H5,
  premethyst cellInfo, ALLC/MCDS, scbs, epiBED, and filter conventions.
- [`other-assays.md`](research/other-assays.md): 22 non-methylation
  read-vs-reference assays, their comparators, formats, expressibility in
  alnbase, and candidate query extensions.

Unless marked "verified by execution", every comparator behaviour cited here
was found by **reading source**, and every suspected bug still has to be
reproduced before it goes anywhere near the paper. Reproducing them is part of phases 1–2.

- [0. Decisions from the first review](#0-decisions-from-the-first-review)
- [1. Goals and non-goals](#1-goals-and-non-goals)
- [2. Comparators, formats, and what they imply](#2-comparators-formats-and-what-they-imply)
- [3. The aggregation utility](#3-the-aggregation-utility)
- [4. Validation](#4-validation)
- [5. Changes to alnbase itself](#5-changes-to-alnbase-itself)
- [6. Phases and check-ins](#6-phases-and-check-ins)
- [7. Decisions needed at this check-in](#7-decisions-needed-at-this-check-in)

---

## 0. Decisions from the first review

- **Scope.** Directional libraries only. Non-directional and PBAT strand resolution is a
  known alnbase rigidity: no tests for it now, recorded as future work.
- **Mate overlap.** alnbase has its own `overlap`; other tools' overlap rules are not
  reimplemented. Comparisons use fragments whose mates do not overlap. alnbase's overlap
  handling is shown separately on Methyl-HiC, where supplementary alignments defeat
  interval-based overlap detection.
- **Naming.** CG, not CpG, except where a downstream format requires the latter.
- **Implicit-query inventories** list only hardcoded behaviour, not configurable
  thresholds, and call a behaviour "undocumented" only if it appears in no help text,
  manual or paper. Each tool is marked molecule-level (per-read output) or pileup-only,
  and molecule-level and single-cell callers (ScaleMethyl `met_extract.py`, premethyst,
  BISCUIT epiread, epialleleR, …) are in scope.
- **Formats.** Not every format: at least one per ecosystem, or an intermediate that
  existing converters turn into what downstream tools need. No windowing unless some
  tool genuinely lacks a way to window pileups (amethyst-facet already windows Amethyst
  H5 files, so the exporter writes base resolution only).
- **Cell barcodes** are read from a BAM tag; the user moves them out of read names upstream.
- **Mainstream tools over custom code**, legible scripts, pre-verified correctness
  (samtools/pyfaidx for FASTA context, existing simulators and converters).
- **Reports** include publication-quality figures with interpretive text, documentation
  for running and reading the suite, and runtime, memory, disk and less obvious resource
  use (threads, file handles, memory maps) with scaling across cores. This machine builds
  confidence; full runs happen on the cluster.
- **Bismark versions:** Perl 0.25.1 as the reference, with an equivalence check against Rust 3.1.0.
- **Candidate alnbase extensions** from revision 1 are resolved in §5: most were already
  expressible with queries.

---

## 1. Goals and non-goals

**Aggregation/export.** Turn the sharded hit rows from `query --parquet` or
`extract` into per-site, per-window, per-region, per-read and per-cell tables,
QC statistics, and the formats downstream tools read. It must keep alnbase's
paradigm (the user declares contexts, call states and filters; nothing is a fixed
CG/CH menu), exploit the sharding guarantee, produce Amethyst v2.0.0 H5 and
cellInfo files for any user-defined context, and explain its plan before running.

Non-goals: upstream filtering that `samtools view` does, statistical modelling
(DMRs, kinetics), plotting.

**Validation.** Reproducibly show that alnbase, given a query that mirrors a
tool's implicit query, reproduces that tool; that every residual difference has a
named, evidenced cause; which implicit choices are contestable or wrong; and
where such a choice changes a biological result.

---

## 2. Comparators, formats, and what they imply

### 2.1 Bulk methylation: the tools disagree with each other

The headline for the paper: the mainstream extractors run **different implicit
queries**. Abridged from `methylation-tools.md` §9, keeping only behaviour the tools do
not let the user change. This table is provisional: the hardcoded-behaviour inventory,
with documentation status and molecule-level/single-cell columns, is the first research
task of phase 1.

| Hardcoded choice | Bismark | MethylDackel | BISCUIT | ALLCools | CGmapTools |
|---|---|---|---|---|---|
| Context taken from | the reference bases aligned to the read, **skipping deleted reference bases**; inserted or clipped read bases make the context unknown; 2 reference bases past the read end | the reference at the position | a 5-base reference window centred on the C | the reference, 3 bases | the reference |
| A reference `N` in the context | call becomes `U`/`u`, never reported | treated as H (not G) | the whole site becomes context `CN` and drops out of CG output, even for an `N` 2 bases upstream | kept literally as `N` | context `--` |
| C within 2 bases of a contig end | whole read discarded at alignment | last base is labelled CHH | context `CN` | row skipped | context `--` |
| C adjacent to a deletion in the read | context read across the deletion, so a genomic CG can be called CHG/CHH | context unaffected (reference) | unaffected | unaffected | the base before an indel is not counted |
| Strand from | XR + XG tags | XG, else flags assuming directional | YD > ZS > XG > inference | alignment orientation | alignment orientation |
| Opposite-strand SNP veto | none | none by default | always on; cannot be disabled | none | none |

### 2.2 Formats the export utility must produce

Coordinates, strand symbols and context labels differ between tools, so each writer
owns those conventions.

| Family | Format | Key conventions |
|---|---|---|
| Bismark | per-call files (`CpG_OT_*.txt`) | read id, `+`/`-` = **methylation state** (not strand), chr, 1-based pos, call letter |
| | bedGraph | 0-based half-open, percentage |
| | `.bismark.cov` | 1-based, start = end, pct, m, u |
| | CpG / CX cytosine report | 1-based; **every** C in the genome incl. zero coverage; strand; trinucleotide context |
| | M-bias table | per context, per read, per position from 5' |
| MethylDackel | bedGraph (+ `--mergeContext`) | `track` header; 0-based; integer (truncated) percent |
| | methylKit | 1-based, `F`/`R`, `%6.2f` |
| | cytosine report | as Bismark, but includes last-base C |
| BISCUIT | pileup VCF → vcf2bed / mergecg, epiBED | VCF 1-based, BED 0-based, betas `%1.3f` |
| methylpy / ALLCools | ALLC | 7 columns, 1-based, k-mer context string |
| CGmapTools | CGmap / ATCGmap | 1-based, context + dinucleotide, level |
| Downstream readers | methylKit, bsseq (cytosine report), DSS (`chr pos N X`), ENCODE bedMethyl (bed9+2), bigWig | |
| Single cell | **Amethyst v2.0.0 H5** | `/metadata/version` = `"amethyst2.0.0"` (a dataset); `/<context>/<barcode>/1` compound `chr S10, pos i8, c i8, t i8`; windows `/<context>/<barcode>/<name>` with `chr, start, end, c, t, c_nz, t_nz` `[start,end)`; rows contiguous per chr (the R `indexChr` requires it); context is only a group name, so any context works; gzip; 1-based positions (premethyst convention) |
| | **premethyst cellInfo** | TSV, no header, 10 columns: `cell_id, cov, cg_cov (distinct sites), mcg_pct, ch_cov, mch_pct, n_frag, n_pair, n_single, xpct_m`; cell id = read name up to first `:` |
| | ALLCools MCDS (zarr), scbs/MethSCAn (per-chr CSR `.npz`), epiBED, scMET long table | see `singlecell-formats.md` §3 |
| Other assays | SLAMdunk `tcount.tsv`, bakR per-read counts CSV, REDItools TSV, JACUSA2 BED6+, ShapeMapper `profile.txt`, SEISMIC position table, rf-count `.rc`, mapDamage `misincorporation.txt`, htseq-clip / PureCLIP BED6, rastair BED/VCF, CRISPResso2 nucleotide tables, VCF / pysamstats TSV | see `other-assays.md` Part 1 |

### 2.3 Non-methylation panel

Ranked for a versatility-and-correctness paper, by signal class, not only by
popularity (`other-assays.md` §01.1):

| # | Assay | Signal | Comparator | Expressible today |
|---|---|---|---|---|
| 1 | SLAM-seq / TimeLapse | T>C per read | SLAMdunk, bakR counter, GRAND-SLAM | yes |
| 2 | A-to-I editing | A>G | JACUSA2, REDItools3 | yes |
| 3 | SHAPE-MaP / DMS-MaPseq | mismatches + indels | ShapeMapper2, SEISMIC-RNA, rf-count | substitutions yes; exact parity needs indel placement (X3) |
| 4 | iCLIP / eCLIP | truncation (read start −1) | htseq-clip, PureCLIP | with caveats; clean with X1+X2 |
| 5 | Ancient DNA damage | C>T by distance from read ends | mapDamage2, DamageProfiler | yes |
| 6 | TAPS or NOMe-seq | inverted logic / dual context (GCH vs HCG) | rastair / Bismark `--nome-seq`, BISCUIT `-N` | yes |
| 7 | SNV baseline | any substitution | bcftools mpileup, pysamstats | yes |

Alternates: GLORI/eTAM (m6A), TRIBE/STAMP, PAR-CLIP, DART-seq, 8-oxoG artifacts.

### 2.4 Requirements this puts on the infrastructure

1. **Coordinate, strand and context conventions per writer**, tested against
   files the real tool produced, not against our reading of its docs.
2. **A site universe.** Cytosine reports list every matching reference position,
   including zero-coverage ones. Generalised: every position where a query's
   *reference row* matches, on either strand. That is a reference-only scan, which
   alnbase can do far more cheaply than the aggregator (§5).
3. **Reference context per site** (trinucleotide, k-mer, dinucleotide): alnbase
   already reports it. Either one query per context (the query name is the context), or
   `^` marks on extra columns so their read and reference bases are written as capture
   columns. Synthetic test data takes context from the FASTA with `samtools faidx` or
   pyfaidx.
4. **Conversion strand per hit** (OT/OB/CTOT/CTOB) for per-strand files and
   M-bias: derivable in SQL from `is_reverse` and `is_last_in_template` for directional
   libraries; a record field is optional (§5).
5. **Cell keys come from a tag** (`-F CB:Z`, `--partition-by CB:Z`); moving a
   barcode out of the read name is the user's upstream step.
6. **Overlapping mates are excluded from comparisons**, not emulated: the validator
   selects fragments whose mates do not overlap and gives every tool that BAM.
7. **Per-read counts** (SLAM-seq conversions, incomplete-conversion filters,
   Amethyst's `c_nz`/`t_nz`, cellInfo's distinct-site counts).
8. **Positions from read ends by mate** (M-bias, damage profiles): `off_5p`,
   `off_3p` and the mate flags are already there.
9. **Tooling logistics.** Most comparators are on bioconda; exceptions need pip,
   tarballs or Docker (Docker here needs socket permission). Several tools only run
   on their own aligner's BAMs (Bismark's extractor needs XM/XR/XG; SLAMdunk needs
   NextGenMap `MP` tags), so validators need aligner stages, not only BAMs.

### 2.5 Early observations, unverified (to be confirmed or dropped in phases 1–2)

**Where a tool's choice can change a result:**

- A C/T SNP at a CG reads as loss of methylation in Bismark, MethylDackel and
  ALLCools; BISCUIT vetoes it. Allele-specific and population-level differential
  methylation are the exposed analyses.
- Deletions next to a C change Bismark's context call, and an insertion or clip
  hides the C.
- Overlap policy changes effective coverage roughly 2× in overlaps (MethylDackel
  issue #170 reports ~50% lower depth than Bismark), and the choice of mate changes
  M-bias.
- Incomplete-conversion thresholds differ roughly 4× between pipelines (premethyst
  0.4 vs 0.1 suggested; ScaleMethyl 0.5), moving neuronal mCH.
- Truncation offset conventions shift CLIP sites by 1 nt relative to motifs, and
  SHAPE/DMS tools place the same ambiguous deletion at different nucleotides.
- Several non-methylation tools ignore or misparse CIGARs. Two cases are
  **verified by execution**: wavClusteR's MD parser after deletions and clips, and
  SAILOR's ≥10-nt left soft clip.

**Limitations of the tools compared with alnbase (preliminary):**

- Context and state rules are fixed, and each tool fixes them differently.
- Contexts are limited to fixed CG/CHG/CHH menus (NOMe needs special modes).
- Filters are bundled into extraction and are often undocumented.
- Strand logic is tied to one library type.
- Several tools lack CIGAR-exact walks.
- There is no record of the definitions used in the output.

---

## 3. The aggregation utility

### 3.1 Shape of the problem

```
hit rows ─▶ hit filter ─▶ fragment policy (overlaps) ─▶ unit filters (read / fragment / cell)
         ─▶ classify into states ─▶ GROUP BY key ─▶ count per state ─▶ derived stats ─▶ writer
```

| Product | Group key | Examples |
|---|---|---|
| site | contig, pos, strand, measure | bedGraph, `.cov`, cytosine report, ALLC, CGmap, methylKit, REDItools, allele counts |
| window / region | contig, bin or BED/GTF feature | Amethyst windows, gene-body mCH, T>C per gene |
| read / fragment | record_id / qname | epiBED, per-read conversions, read-level filters |
| cell / sample | CB tag, RG, qname-derived key | cellInfo, Amethyst per-cell tables |
| position-in-read | measure, mate, off_5p / off_3p | M-bias, damage profiles |
| substitution | read_base × refr_base × context | mutation/structure profiles, damage matrices |

### 3.2 Concepts

**Dataset.** One run's shards plus the manifest (§5): partition key, file count,
query files, walk options, library, and contigs with their lengths. Incomplete
file sets are refused.

**Measure.** A named classification of query names into states. This is how
contexts stay user-defined:

```toml
[measure.CG]
m = ["CG"]            # protected
u = ["TG"]            # converted

[measure.GCH]         # NOMe accessibility, from the user's own queries
m = ["GCH_C"]
u = ["GCH_T"]

[measure.TtoC]        # SLAM-seq: states need not be methylation
conv  = ["T_as_C"]
total = ["T_as_C", "T_as_T"]
```

A query may belong to several measures. Derived values are expressions the user
writes over the measure's states, evaluated as SQL:

```toml
[measure.CG]
m = ["CG"]
u = ["TG"]
level = "m / (m + u)"
depth = "m + u"
```

**Stages and filters.** Each filter is a SQL boolean over the columns of its level:

| Stage | Columns | Example |
|---|---|---|
| `hit` | hit + record columns | `qual >= 20`; M-bias trim `off_5p >= 10 AND off_3p >= 2` |
| `read` | per-read measure counts | Bismark `filter_non_conversion`: `CH.m < 3`; ScaleMethyl: `CH.level <= 0.5` |
| `fragment` | per-qname measure counts | premethyst: `CH.m + CH.u > 0 AND CH.level <= 0.4` |
| `cell` | per-cell summaries | `reads >= 10000`; Amethyst vignette `cov > 1e6 AND CH.pct < 12` |
| `site` | per-site counts | `CG.depth >= 5`; exclude a BED; SNP veto from evidence queries (below) |

**SNP evidence** is a measure like any other. Queries whose anchor lands on the base
opposite a cytosine, read from the other conversion strand (where conversion does not
change it), count reference-matching and variant bases there; a site stage then joins that
evidence to the call site and drops sites whose variant fraction passes a threshold.

**Product.** A group key, measures, stages, and a writer with its options.

### 3.3 Using the sharding guarantee

A grouping `G` runs **inside each shard with no merge** exactly when the partition
key `P` is determined by `G`: equal `G` ⇒ equal `P` ⇒ same file (within a run).

| Grouping | Shard-local when partitioned by |
|---|---|
| record | always |
| fragment | `qname`, or anything constant across a fragment (`RG`, `CB`, a qname-derived cell key) |
| cell | the cell key |
| site | `ref_name` only; otherwise merged |

Execution is two-phase:

1. **Per shard, in parallel**, each with its own memory-limited DuckDB
   connection: hit filters, fragment policy, every shard-local unit filter, then
   partial aggregation to the finest group the products need (e.g. `(cell, contig,
   pos, strand, measure) → counts`). Read- and fragment-level work is where this
   pays off: a `GROUP BY record_id` over a whole run is a hash table the size of the
   read count, while per shard it is `1/N` of that.
2. **Merge**: concatenate when the group is shard-local (per-cell H5 needs no
   merge), otherwise `SUM` partial counts, which is exact because counts are additive.
   Statistics that don't compose (medians, distinct counts across cells) are
   computed after the merge or refused with a message.

If the dataset's partition key does not make a requested stage shard-local, the
planner runs that stage globally and says so in `explain`. Correctness never
depends on how the data was sharded, only speed does.

### 3.4 Interface

**1. Recipes (TOML),** stored inside or beside every output they produce, the same way
query files travel in the BAM header. The kind of a table is in its key, as in alnbase's
own `[tag.XM.bases]`: an `[output.NAME]` table is one aggregation, and each writer is a
subtable keyed by format, so one aggregation can feed several writers and each writer
has its own strictly checked options.

```toml
[input]
hits = "calls_*.parquet"

[measure.CG]
m = ["CG"]
u = ["TG"]
level = "m / (m + u)"

[measure.CH]
m = ["CHG", "CHH"]
u = ["THG", "THH"]
level = "m / (m + u)"

[stage.hit]
where = "qual >= 20"

[stage.read]                         # incomplete conversion
where = "CH.m + CH.u = 0 OR CH.level <= 0.1"

[output.cg_sites]                    # one aggregation ...
by = "site"
measures = ["CG"]

[output.cg_sites.bismark-cov]        # ... two writers
path = "out/sample.CG.cov.gz"

[output.cg_sites.bedgraph]
path = "out/sample.CG.bedGraph"
value = "CG.level"

[output.cells]
by = "cell"
cell = "CB:Z"
measures = ["CG", "CH"]

[output.cells.amethyst-h5]           # base resolution; group names are the measure names
path = "out/sample.h5"

[output.cells.premethyst-cellinfo]
path = "out/sample.cellInfo.txt"
columns = { cg = "CG", ch = "CH" }   # which measures fill the CG and CH columns
```

**Palindromic contexts.** alnbase reports each strand's cytosine at its own coordinate:
for a CG, the C on the plus strand and the C on the minus strand (the plus strand's G)
are two hits at adjacent positions. Some formats want them summed into one site. That
is aggregation, so it is an option of a site aggregation, not something alnbase does,
and it is only meaningful when the context reads the same on both strands (CG, CWG).
The option takes the offset explicitly (`merge_strands = { minus_offset = -1 }` for CG)
and the planner checks it against the measure's query definitions.

**2. CLI.**

```
alnbase-agg explain recipe.toml             # stages, shard-local or merged, generated SQL
alnbase-agg run recipe.toml [-j N] [--memory-limit 2GB]
alnbase-agg export --preset methyldackel-bedgraph --measure CG=CG/TG calls_*.parquet out.bedGraph
alnbase-agg presets                         # recipes reproducing each tool's defaults
alnbase-agg stats mbias|conversion|coverage|cells|strand-concordance  calls_*.parquet
```

**Presets** (`bismark`, `methyldackel`, `biscuit`, `allcools`, `premethyst`,
`slamdunk`, `mapdamage`, …) are ordinary recipes, and each one is validated against
the tool it names. That makes a preset a documented, tested statement of that
tool's implicit choices, which is useful in itself.

**3. Python API.** The same engine as a library, returning DuckDB relations so
users can continue in SQL, pandas or polars. This is also the basis for interactive
vignettes.

```python
import alnbase_agg as aa
ds  = aa.open("calls_*.parquet")
cpg = ds.measure("CG", m=["CG"], u=["TG"]).stage_hit("qual >= 20")
cpg.by_site(merge_strands=True).write("bedgraph", "out.bedGraph")
cpg.by_site().relation.df()
```

### 3.5 Implementation

**Recommended: a Python package, `alnbase-agg`, on DuckDB + pyarrow, with h5py,
pysam and pyBigWig for formats.**

- DuckDB does the heavy work in C++, out of core, directly on the parquet shards.
  Python plans, dispatches shards, and writes formats.
- The format libraries (h5py, pysam/tabix, pyBigWig, anndata, zarr) are in
  Python, and so are the target users and their notebooks.
- The validation suite uses the package directly.
- The contract is the recipe format plus the manifest, so a Rust port into the
  main binary stays open once the design settles.

Considered and not recommended now: building it into the Rust binary with
`duckdb-rs` and the `hdf5` crate. That bundles DuckDB's C++ build and needs
libhdf5 at build time, which complicates the cluster cross-build, and every format
change becomes a Rust release.

```
python/alnbase_agg/
  dataset.py     manifest, completeness and schema checks, .aref reader (contexts, site universe)
  recipe.py      TOML → typed plan; errors by table/key
  plan.py        stage placement, SQL generation, explain
  engine.py      per-shard process pool, DuckDB memory limits, merge
  policies.py    fragment observation policies
  stats/         mbias, conversion, coverage, cell QC, strand concordance
  writers/       one module per format, each tested against real tool output
  presets/       recipes reproducing tool defaults
```

Outputs follow the manifest's contig order (`@SQ`), not lexicographic order.

---

## 4. Validation

### 4.1 What makes it convincing

1. **Agreement where definitions agree**: per call when the tool exposes calls,
   per site otherwise.
2. **Every difference explained** by ablation, with the evidence kept.
3. **Ground truth** on simulated data, where every tool is scored against the truth
   rather than against each other.
4. **Breadth**: one engine across assay families, changing only query files and
   recipes.
5. **Consequence**: at least one case where a tool's choice changes a downstream
   result.

### 4.2 Mirroring each tool's implicit query

Each tool gets `validation/tools/<tool>/`. It holds the mirror query file, the
mirror recipe, the `samtools view` prefilter for its flag/MAPQ behaviour, and an
**inventory file**. The inventory has one row per implicit choice, with the source
citation and how alnbase expresses it (query, walk option, prefilter, recipe stage,
or *not expressible*). Examples of what mirroring involves:

- **MethylDackel's "N counts as H, contig end counts as CHH":** a boolean over
  patterns, e.g. `C~~@C~~ and not C~@CG and not C~~@C~G`. A reference `N` is
  not a subset of `G`, so it fails the exclusions, and `~` matches the pad past the
  contig end. Expressible today; to be confirmed by probe.
- **Bismark's read-projected context:** deleted reference bases are dropped from
  the context. alnbase emits deletion columns, so a fixed-width pattern can only
  enumerate bounded deletion lengths. A `--deletions skip` walk option (symmetric
  with `--insertions skip`) would make it exact (§5).
- **Overlap rules:** not mirrored; comparisons use fragments whose mates do not overlap (§0).

Where a choice is not expressible, the validator records it and the ablation
attributes the residual to it. That is how candidate extensions are justified with
numbers.

### 4.3 The strongest methylation test: per-read, per-base, in the same BAM

Bismark's calls live in its `XM` tag. Run alnbase on the Bismark BAM with the
mirror query writing a local tag (`xm`), and every read carries both call strings
over identical bases. Differences are then exact to the read and base, with no
coordinate conversion or aggregation in between. Feeding Bismark's own extractor
the alnbase-tagged BAM (with `XM` replaced) closes the loop through its
aggregation too.

For site-level tools (MethylDackel, BISCUIT, ALLCools), comparison is per site,
and each discrepant site is attributed by tracing the reads that cover it
(`--trace-records` on those reads) and by ablation.

### 4.4 Layout

```
validation/
  run.py                 central runner: select validators × datasets × scale; execute; aggregate; report
  datasets.toml          registry: URLs, checksums, sizes, which scales use them
  envs/                  pinned micromamba env per tool (lock files)
  lib/
    env.py               download micromamba, create envs, record exact tool versions
    data.py              download with checksum, cache, `clean`
    reference.py         FASTA → faidx → `alnbase index` (.aref) → tool genome preps
    simulate/            synthetic genomes, reads with per-base truth, oracle BAMs
    compare.py           uniform result schema
    ablate.py            ablation ladders and attribution
    report.py            single-file HTML + JSON across validators
  tools/<tool>/
    validator.py         prepare / run_tool / run_alnbase / compare
    query.toml           mirror query
    recipe.toml          mirror aggregation
    inventory.toml       implicit choices, citations, alnbase equivalents
    ladder.toml          ablation steps
```

### 4.5 Uniform result schema

- `calls.parquet`: `validator, dataset, measure, qname, mate, contig, pos, strand,
  tool_state, alnbase_state, truth_state, category, cause`.
- `sites.parquet`: `validator, dataset, measure, contig, pos, strand,
  tool_{states}, alnbase_{states}, truth_{states}, category, cause`, with `category`
  ∈ `agree | count_diff | tool_only | alnbase_only`.
- `metrics.parquet` (long): `validator, dataset, scale, level, measure, metric,
  value`. Metrics: sites in union/both/each only, exact agreement fraction,
  Σ|Δcount|, level RMSE and correlation, and precision/recall/F1 against truth.
- `ablation.parquet`: `validator, dataset, step, description, discrepant_calls,
  discrepant_sites`.
- `inventory.parquet`: all tools' inventory files, so the cross-tool
  implicit-query table is generated, not hand-maintained.

The report shows, per validator: tool versions and exact commands; agreement
metrics; an **ablation waterfall** (discrepancy removed by each implicit choice);
the top residual sites with read traces; and the generated cross-tool inventory.

### 4.6 Data, simulation, scale

| Scale | Data | Runtime target |
|---|---|---|
| `tiny` | simulated ~100 kb genome with designed hazards; ~10k pairs; oracle BAM + tool-aligned BAMs | seconds to a minute |
| `small` | simulated few-Mb genome; tool test sets (nf-core methylseq lambda/E. coli, Bismark and MethylDackel test files, NEB EM-seq fixtures) | minutes |
| `full` | public datasets per assay | hours; paper figures |

Simulation uses established simulators where they give per-read truth (the survey of
Mason, Sherman, BSBolt and others is a phase 1 task), and small hand-written SAM files
with a tiny FASTA for edge cases, in the style of `docs/reference/examples/`, which are
easier to read than any generator. Hazards to cover: CGs at read ends and contig ends;
`N` runs, including an `N` two bases upstream; deletions and insertions inside
contexts; soft clips; C/T SNPs at CGs; incompletely converted reads; duplicates,
secondary and supplementary alignments; low-quality tails. Directional libraries only;
overlapping mates are excluded from tool comparisons. Oracle BAMs with true CIGARs test
extraction apart from alignment; FASTQs feed the aligner arm (Bismark, bwa-meth).

This machine has 20 cores and 15 GB RAM (about 4 GB free now), so `full` runs of
public data will need to be subsampled here or run elsewhere.

### 4.7 Source reading vs. probing: the experiment

Run before any empirical probing. Both arms are run by investigators who have not
seen the source-derived research, and both are **barred from
`docs/design/research/`**:

- **Arm A, source first:** read the Bismark extractor's (and aligner's calling)
  source, write the mirror query, recipe and inventory, then design probes to
  confirm.
- **Arm B, probe first:** treat the same tool as a black box, derive the mirror
  query from synthetic probe BAMs, then read source to check.

Both are scored against an adjudicated answer key (the research inventory, with
each item confirmed by probe), on dimensions found / wrong / missed, wall time, tool
calls, tokens, and number of probes. A second pair on MethylDackel replicates it if
the first result is not clear-cut. The outcome sets the method for the remaining
tools.

---

## 5. Changes to alnbase itself

**Done in 0.1.2** (see `CHANGELOG.md`): offsets as sequenced for read 2; soft clips
ignored; `--batch-rows` as a total; `.aref` v2 with MD5 and `@SQ M5` checking; clean
errors and output removal under file and thread limits.

**Done in 0.1.8 for aggregation:**

1. **Run manifest.** In every file's footer (`alnbase_manifest`, beside
   `format_version` and `coordinate_base = 0`), plus `{stem}.manifest.json` written
   last as the completeness marker: alnbase version, `run_id`, command line, input and
   reference, library, partition key, worker and shard counts, file names by slot,
   fields, capture layout, query file text, walk options with defaults resolved, and
   contigs in `@SQ` order with lengths and MD5s; the file adds record counts,
   concordance totals and rows per file. See `docs/cli-reference.md`.
2. A `conv_strand` record field (OT/OB/CTOT/CTOB).

Hive-partitioned output was dropped: the sharded files plus the manifest are the
dataset.

**Off-by-one safeguards:** see `off-by-one-safeguards.md` (resolved in 0.1.8).

**Candidate extensions from revision 1, resolved:**

| Revision 1 proposal | Outcome |
|---|---|
| `--deletions skip` for Bismark's context | Not needed: a query enumerating deletions of up to k bases between C and G reproduces it (verified); deletions longer than k are counted in the ablation. |
| End-specific pad codes | Not needed: a pad's position in the pattern already distinguishes the ends. |
| Indel placement normalisation (SHAPE/DMS) | Not needed: enumerate the placements of an ambiguous indel and anchor on a fixed base (verified). |
| Library/strand modes for RNA assays | Out of scope (directional only). For stranded RNA libraries a pattern is written in walk orientation, reverse-complemented where needed. |
| Soft-clip code or projection | Dropped; soft clips are ignored (0.1.2). |
| Strand/mate predicate in `where` | Dropped; hit tables filter on flags downstream. |
| Per-read count tag `[tag.XX.count]` | Dropped. |
| Reference site scan for zero-coverage reports | Not in alnbase; if a writer needs a site universe, `seqkit locate` or a FASTA scan provides it. |

---

## 6. Phases and check-ins

| Phase | Work | Ends with |
|---|---|---|
| 0 | research; this design; reference documentation of alnbase; 0.1.2 fixes | **check-in (now)** |
| 1 | hardcoded-behaviour inventory incl. molecule-level and single-cell callers; minimal writer set and existing converters; simulators and workflow tooling; manifest; `alnbase-agg` MVP; MethylDackel and Bismark validators (incl. same-BAM XM diff) on non-overlapping fragments; resource and performance monitoring; central report with figures | **check-in**: first discrepancy magnitudes |
| 2 | source-vs-probe experiment; ablation ladders; BISCUIT, ALLCools, ScaleMethyl, premethyst; confirm or drop suspected bugs | **check-in**: causes, with evidence |
| 3 | public data; non-methylation panel; Methyl-HiC overlap demonstration | **check-in** |
| 4 | the demonstration where a tool's choice changes a result; codebase tour of the discrepancies | **check-in** |

---

## 7. Decisions needed at this check-in

1. The off-by-one safeguards in `off-by-one-safeguards.md`, especially per-query reach
   (a query's hits currently depend on the other queries in the run) and the default-on
   concordance check.
2. The non-methylation panel for phase 3 (§2.3), with NOMe-seq recommended over TAPS.
