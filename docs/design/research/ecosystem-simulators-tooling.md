# Export intermediates, ground-truth simulators, and validation tooling

Research date: 2026-09-16. This note answers four questions for the
`alnbase-agg` export utility and the validation suite
([`../aggregation-and-validation.md`](../aggregation-and-validation.md)):

1. Which formats should `alnbase-agg` write, so that maintained converters and readers do the rest?
2. Which downstream analyses need read- or molecule-level aggregation that a base-resolution pileup cannot give?
3. Which simulators give ground truth, per assay?
4. Which workflow, benchmarking and reporting tools should the suite use?

It builds on [`methylation-tools.md`](methylation-tools.md),
[`singlecell-formats.md`](singlecell-formats.md) and [`other-assays.md`](other-assays.md),
and does not repeat them. Citations of the form `repo/path:line` refer to shallow clones in
`clones/`, at these commits:

| Repo | Commit / version |
|---|---|
| lhqing/ALLCools | `c9f7be2` (1.1.1); DingWB fork `4a1f8a1` has the same CLI |
| amethyst-facet | PyPI sdist 1.2.4 |
| adeylab/premethyst | `fa5a47b` |
| ScaleBio/ScaleMethyl | `2a4973c` |
| anders-biostat/MethSCAn | `72f38fb` (PyPI `methscan` 1.1.0) |
| dohlee/metheor | `33248b1` |
| BioinfoUNIBA/REDItools3 | `8cfebea` (3.7) |
| dincarnato/RNAFramework | 2.9.7 (docs from RNAFramework-docs) |
| rouskinlab/seismic-rna | `7c5ac78` (bioconda 0.26.0) |
| ginolhac/mapDamage | `484bb30` (bioconda `mapdamage2` 2.2.3) |
| EMBL-Hentze-group/htseq-clip | 2.19.0b0 |
| EZbakR | GitHub clone (see `other-assays.md` §01.3) |
| snakemake | v9.27.0 release files, and `main` `91763d6` |
| nextflow | master `command-trace.txt`, `TraceRecord.groovy` (release 26.04.6) |

Bioconda and PyPI versions were queried on 2026-09-16 (`api.anaconda.org/package/bioconda/<pkg>`, `pypi.org/pypi/<pkg>/json`).

---

## TL;DR

- **Bulk and per-cell methylation need one text writer.** That writer is the Bismark
  genome-wide **cytosine report** (7 columns, 1-based, stranded, trinucleotide context).
  It feeds bsseq and everything built on BSseq (DSS, dmrseq, methylSig), and methylKit
  directly. ALLCools reaches it through `allcools table-to-allc`, which also gives
  bigWig and MCDS. MethSCAn reaches it through a custom `--input-format`.
- **For Amethyst**, the most legible route is a per-cell parquet table
  (`chr,pos,context,methylated,unmethylated`), fed to `facet calls2h5`. This is DuckDB's
  native output with `PARTITION_BY`.
- **ENCODE bedMethyl is the one methylation format with no converter from anything
  else.** Write it only if ENCODE-style bigBeds are actually wanted.
- **Epiallele and read-level methylation tools need a BAM, not a table.** alnbase
  itself already writes Bismark `XM/XR/XG`, and that BAM feeds epialleleR, Metheor, the
  Bismark extractor (M-bias, `filter_non_conversion`) and premethyst `bam-extract -B`
  (cellInfo). No epiBED writer is needed unless biscuiteer is a target.
- **Non-methylation formats:**
  - SLAM-seq: cB table → bakR/EZbakR.
  - A-to-I editing: REDItools table → `reditools annotate`/`index`.
  - SHAPE/DMS: RNA Framework **MM** file → `rf-mmtools toRC` → `rf-norm`, `seismic importmm` → SEISMIC, and DRACO.
  - aDNA: mapDamage `misincorporation.txt` + `dnacomp.txt` → `mapDamage --plot-only/--stats-only`.
  - CLIP: per-read crosslink BED6 → `htseq-clip count` → DEWSeq; aggregated BED6 → clippy / `iCount peaks`.
  - GRAND-SLAM, JACUSA2 and PureCLIP take only BAMs, so they are compared on their outputs.
- **Windowing:** every windowed consumer found already has an existing tool that
  builds windows from pileups. **No windowing feature is needed.**
- **Aggregation that cannot come from a pileup** is all per-read (or per-fragment):
  - SLAM-seq per-read conversion counts (cB);
  - per-read incomplete-conversion filters and cellInfo fragment counts;
  - SHAPE/DMS per-read mutation vectors;
  - hyper-edited reads;
  - epiallele metrics.

  The last are covered by the XM BAM. The rest are one `GROUP BY` read in DuckDB.
- **Simulators:**
  - Bisulfite/EM-seq: **BSReadSim** (golden BAM with per-base truth; not on bioconda, so build from source), with **Mason** (bioconda) as the site-level cross-check.
  - SLAM-seq: splash (in slamdunk).
  - aDNA: NGSNGS.
  - SHAPE/DMS: `seismic sim`.
  - SNVs: Mason.
  - Editing: BAMSurgeon `addsnv`.
  - Edge cases: SAM text plus a tiny FASTA, built with samtools.
- **Tooling:**
  - Workflow: **Nextflow**, because the cluster runs it (local executor).
  - Measurement: every command wrapped in `resource_monitor.py`, extended with VmHWM, rlimits, affinity and `du`, plus GNU `/usr/bin/time` and `taskset` for core-count sweeps. The Nextflow trace is a cross-check for CPU and I/O. Neither Nextflow's trace nor Snakemake's benchmark records threads, fds or mmaps, so neither replaces the wrapper.
  - Timing: hyperfine for timing only.
  - Environments: pixi or micromamba, rootless.
  - Reporting: Quarto.

---

## 1. Minimal writer set, via intermediates

### 1.1 Findings behind the table

**Bioconductor.** `bsseq::read.bismark()` reads Bismark coverage files (6 columns, no
strand) or genome-wide cytosine reports (7 columns, stranded). The docs "strongly
recommend" the cytosine report, and strand collapse needs it
([read.bismark](https://rdrr.io/bioc/bsseq/man/read.bismark.html)). Loci are the union
over files, so a report that lists only covered sites is fine. `rmZeroCov` handles zeros.
- **methylSig** (bioconda 1.22.0) reads its data through `bsseq::read.bismark()` and has
  `tile_by_windows()` / `tile_by_regions()`
  ([vignette](https://bioconductor.org/packages/release/bioc/vignettes/methylSig/inst/doc/using-methylSig.html)).
- **DSS** (`DMLtest`) and **dmrseq** take a `BSseq` object directly.
- **methylKit** reads the same file with
  `methRead(pipeline="bismarkCytosineReport", context="CpG")`. With
  `pipeline="bismarkCoverage"` it loses strand
  ([methRead](https://rdrr.io/bioc/methylKit/man/methRead-methods.html)).

So one file serves all five packages.

**ALLCools.** `generate-dataset` accepts **only** `--allc_table`, a headerless TSV of
`cell_id<TAB>allc_path` (`ALLCools/count_matrix/dataset.py:299-300`). Each ALLC must be
bgzipped and tabix-indexed (`pysam.TabixFile`, `:162`). Region counting sums rows whose
**context string exactly equals** one of the k-mers the IUPAC `mc_type` expands to
(`:44-66`). A literal `CG` context, as ScaleMethyl's `write_allc.py` writes, therefore
matches nothing for `CGN`: contexts must be reference k-mers.

The CLI has a converter, `allcools table-to-allc` (`__main__.py:712-760`,
`table_to_allc.py`). It takes 0-based column indices for chrom, pos, strand, context, and
any two of mc/uc/cov/mc_frac. When strand and context columns are given, it does no FASTA
lookups (`table_to_allc.py:97-99`). It then runs `sort -k1,1 -k2,2n -S 10G | bgzip` and
`tabix -b 2 -e 2 -s 1` (`:342-352`).

Caveats:
- `-S 10G` is a sort buffer request, which matters on this low-RAM host.
- `frac` is `round(mc/cov)`, so cov=0 rows produce a NaN cast. Do not write zero-coverage rows.
- **The subcommand does not run at all**, in 1.1.1 or on master. `@doc_params` is
  called with `input_path_doc=…` while the docstring interpolates
  `{table_to_allc_input_path}` (`:227`, `:274`), and `doc_params` does
  `obj.__doc__.format(**kwds)` at decoration time, so importing the module raises
  `KeyError: 'table_to_allc_input_path'` before any argument is read. Only this
  subcommand is affected; `allc-to-bigwig` imports fine. Verified by execution on
  2026-09-21. The function itself is correct — neutralising the decorator and
  calling it works — so the fix is one keyword name, but until it is released a
  user of this route needs the patch.
- `--header` is declared `type=str` and handed to `pandas.read_csv`, which raises
  `ValueError: header must be integer or list of integers` on a string. So there
  is no working way to skip a header line: the input table must be headerless,
  which is why every convertible writer in `alnbase_agg/text.py` uses
  `HEADER FALSE`.

Other commands:
- `allcools allc-to-bigwig --bin_size 1` gives base-resolution `frac` and `cov` bigWigs (default bin 50; `_allc_to_bigwig.py:97-180`).
- `extract-allc --strandness merge --output_format bed5` gives CpG-merged tables.

**Amethyst v2.0.0 H5.** Four existing writers were found:

| Writer | Input | Output | Notes |
|---|---|---|---|
| `facet calls2h5` (PyPI `amethyst-facet` 1.2.4, not on bioconda) | per-cell `.cov`/`.cov.gz` with column indices set by flags (`--cov-chr-col`, `--cov-pos-col`, `--cov-t-col`, `--cov-c-col`), or per-cell **ScaleMethyl-style parquet** with columns `chr,pos,context,methylated,unmethylated`, or other v2 H5 (`cli/commands/calls2h5.py:22-31, 181-218, 566-760`) | v2.0.0 with `/metadata/version` | Barcode and context are parsed from the **path** with `--parse` templates. For parquet, context comes from the column, so **any context label works**. Suffix must include `.cov`, `.parquet` or `.h5`. dtype `chr S10` silently truncates contig names longer than 10 bytes |
| premethyst `calls2h5` (`premethyst_commands/calls2h5.py`) | folder of `<cell>.CG.cov`/`<cell>.CH.cov` (5 cols `chr pos pct t c`) | **v1** (`/CG/<cell>` datasets); then `facet convert` → v2 | CG and CH only, hard-coded |
| ScaleMethyl `bin/write_amethyst.py` | multi-cell parquet with `barcode` column | `/CG/<bc>/1`, no `/metadata/version` | Pipeline script; imports the pipeline's `reporting` module; not packaged |
| premethyst `bam-extract -B` | name-sorted BAM with Bismark `XM` and barcode as qname prefix (`bam_extract.pm:40-42, 75`) | per-cell `.cov` **and cellInfo** | The only existing producer of cellInfo from a BAM. Applies premethyst's mCH read filter |

No tool converts Bismark `.cov` directly except through facet's column flags. Amethyst's R
package writes no H5 (no `h5write` in `amethyst/R`).

**scbs / MethSCAn** (PyPI `methscan` 1.1.0; not on bioconda). `prepare --input-format`
accepts `bismark`, `allc`/`methylpy`, `biscuit`, `biscuit_short`, or a custom
`chrcol:poscol:methcol:<n>{u|c}:sep:header` (`methscan/cli.py:186-208`).

**Gotcha:** the built-in `allc` format sets `header=True` (`prepare.py:322-331`). Real ALLC
files have no header, so **the first site of every file is silently dropped**. Use
`'1:2:5:6c:\t:0'` for ALLC, or `'1:2:4:5u:\t:0'` for a cytosine report.

**Epialleles.**
- **epialleleR** (bioconda 1.18.0) requires `XG` and `XM` tags. For paired-end data it
  requires a QNAME-sorted BAM. For bwa-meth/BSMAP BAMs it offers `callMethylation`
  ([vignette](https://bioconductor.org/packages/release/bioc/vignettes/epialleleR/inst/doc/epialleleR.html)).
- **Metheor** (only on the `dohlee` conda channel, not bioconda) computes PDR, LPMD, MHL,
  PM, ME, FDRP and qFDRP from a Bismark-style BAM. It ships `metheor tag` to add `XM`
  (`README.md:20, 266-274`).
- **epiBED** is read by biscuiteer `readEpibed`
  ([BISCUIT NAR 2024](https://academic.oup.com/nar/article/52/6/e32/7614859)). No other
  maintained consumer was found.

alnbase already writes `XM`, `XR` and `XG` (`docs/alnbase.md:256-290`,
`docs/bismark-xm.toml`), so these tools need no alnbase-agg writer at all.

**SLAM-seq and metabolic labelling.**
- **bakR** requires `XF, sample, TC, nT, n` (`bakR/R/bakRData.R:182`).
- **EZbakR** requires `sample`, `n`, and an `n<X>` column for every mutation column
  (`EZbakR/R/EZbakRData.R:49-74`). It also accepts an **Arrow dataset** (`EZbakRArrowData`),
  so a DuckDB `COPY … (FORMAT parquet, PARTITION_BY (sample))` is read natively.
- **Feature columns** (e.g. `XF`/`GF`) come from `featureCounts -R BAM` tags, which alnbase
  carries through as record columns.
- **GRAND-SLAM** (bioconda `gedi` 1.0.6a) reads BAMs or CIT only.
- **SLAMdunk's `tcount.tsv` columns can be derived from a cB.** `TcReadCount` =
  `SUM(n) FILTER (TC >= k)`, `ConversionsOnTs` = `SUM(TC*n)`, and `ReadCount` = `SUM(n)`.
  SLAMdunk's read-span `CoverageOnTs` is the exception (`other-assays.md` §01.1).

**A-to-I editing.** REDItools3 (bioconda 3.7) has `annotate`, which fills the `g*` DNA
columns from a second table, and `index`, which computes the editing index from `analyze`
tables (`REDItools3/README.md:39-80`; `reditools/tools/{annotate,index}`). JACUSA2
(bioconda 2.1.17) reads only BAMs. JACUSA2helper reads only JACUSA2 output.

**SHAPE/DMS.**
- **RC.** `rf-norm` needs RC files, which are binary. No text→RC converter exists;
  `rf-json2rc` is for DRACO JSON only.
- **MM.** `rf-mmtools` has **`toRC`**, which converts MM (per-read mutation maps) into RC
  (`RNAFramework-docs/docs/rf-mmtools.md`).
- **SEISMIC** (bioconda 0.26.0) has **`seismic importmm`**: "Import RNA Framework Mutation
  Map (MM) files as IDmut outputs" (`seismic-rna/src/seismicrna/importmm/main.py:47`).
  After import, `filter`, `cluster`, `table` and `fold` all apply.
- **DRACO** reads MM natively.

So **MM covers three consumers**. The MM layout per transcript is:
- ID length as **uint16** (`lib/RF/Data/IO/MM.pm:74, 135`, and SEISMIC `importmm/mm.py:34`; the rf-count docs say uint32, which is wrong);
- the ID;
- sequence length (uint32);
- the 4-bit packed sequence;
- the read count;
- for each read: `start, end, n_mut` (uint32), then the mutation indices;
- the file ends with an `[mmeof]` marker.

**Ancient DNA.** `mapDamage --plot-only -d DIR` and `--stats-only -d DIR` need only
`misincorporation.txt` and `dnacomp.txt` in `DIR` (`mapdamage/config.py:401-456`;
`rscript.py:12-37`). `lgdistribution.txt` is also read for the length plot. `--rescale`
needs the BAM.

**CLIP.**
- **htseq-clip `count`** adds **+1 per BED line** and ignores the score
  (`clip/countCLIP.py:216-244`). Its input must therefore be one line per read, as written
  by `extract`. It then runs `createMatrix` → DEWSeq.
- **iCount `peaks`** (bioconda 2.0.0) takes "cross-links in BED6 format" with counts in the
  score (`iCount/analysis/peaks.py:383-402`). **clippy** (bioconda 1.5.0) consumes the same
  kind of file.
- **PureCLIP** (bioconda 1.3.1) needs a BAM.

### 1.2 Table

| Ecosystem | Intermediate written by alnbase-agg | Existing converter / reader (command) | Lost or constrained |
|---|---|---|---|
| **bsseq, DSS, dmrseq, methylSig** | Bismark genome-wide **cytosine report** (`chr pos strand M U ctx trinuc`, 1-based, gz), covered sites only | `BS <- bsseq::read.bismark(files, colData=…, strandCollapse=TRUE)`, then `DSS::DMLtest(BS, …)`, `dmrseq::dmrseq(BS, …)`, `methylSig::filter_loci_by_coverage(BS)` / `diff_binomial` | Context labels must be CG/CHG/CHH for context-aware readers. Nothing per read |
| **methylKit** | same file | `methRead(files, sample.id, assembly, pipeline="bismarkCytosineReport", context="CpG", mincov=…)`; windows: `tileMethylCounts` | methylKit's context filter expects Bismark labels (to be verified). Default `mincov=10` |
| **Genome browsers: bigWig** | none extra (cytosine report) | `allcools table-to-allc --input_path r.CX_report.txt.gz --output_prefix s --chrom 0 --pos 1 --strand 2 --mc 3 --uc 4 --context 6` → `allcools allc-to-bigwig --allc_path s.allc.tsv.gz --bin_size 1 --mc_contexts CGN --chrom_size_path chrom.sizes --output_prefix s` | One float stream per file (frac or cov). If the two hops are too slow, add a 4-column bedGraph writer → UCSC `bedGraphToBigWig in.bg chrom.sizes out.bw` (bioconda `ucsc-bedgraphtobigwig` 482; input sorted `LC_ALL=C sort -k1,1 -k2,2n`) |
| **ENCODE bedMethyl / bigBed** | **bed9+2** (`chrom start end name score=min(cov,1000) strand thickStart thickEnd rgb coverage pct`), *only if needed* | `bedToBigBed -type=bed9+2 -as=bedMethyl.as in.bed chrom.sizes out.bb` (bioconda `ucsc-bedtobigbed` 482) | No converter into this format exists. ENCODE4 gemBS files add three genotype columns ([ENCODE4 WGBS](https://www.encodeproject.org/data-standards/wgbs-encode4/)). `modkit bedmethyl tobigwig` (bioconda `ont-modkit` 0.6.4) expects modkit's 18-column bedMethyl with a mod code in column 4, not ENCODE's |
| **ALLCools / MCDS** (ALLC only) | per-cell cytosine report (as above), **or** ALLC written directly if `table-to-allc` is too slow for many cells | `table-to-allc` (above), then `allcools generate-dataset --allc_table cells.tsv --output_path x.mcds --chrom_size_path chrom.sizes --regions chrom100k 100000 --quantifiers chrom100k count CGN,CHN` | Contexts must be reference k-mers matching IUPAC patterns. User context names do not survive |
| **Amethyst v2.0.0** (base resolution) | per-cell **parquet** `chr,pos,context,methylated,unmethylated` (1-based pos) via `COPY … TO 'cells' (FORMAT parquet, PARTITION_BY (cell))` | `facet calls2h5 --parse '{d}/cell={barcode}/{f}.parquet' cells.h5 cells/cell=*/*.parquet`; then `facet agg` for windows | `chr` truncated at 10 bytes. cellInfo is not produced (see §2). Alternative with no new writer: per-cell, per-context cytosine report named `*.cov.gz` with `--cov-pos-col 1 --cov-c-col 3 --cov-t-col 4` (works, but less legible) |
| **scbs / MethSCAn** | per-cell cytosine report (or per-cell ALLC) | `methscan prepare --input-format '1:2:4:5u:\t:0' cells/*.txt.gz data_dir` (ALLC: `'1:2:5:6c:\t:0'`, **not** `allc`); windows: `methscan matrix regions.bed` | One context per run. Strand not merged. Mixed per-cell sites dropped unless `--round-sites` |
| **Epialleles: epialleleR, Metheor, BISCUIT epiBED consumers** | **none**: alnbase `query` writes Bismark `XM/XR/XG` into the BAM | epialleleR `preprocessBam` → `generateVEFReport`/`generateMhlReport`/`extractPatterns` (QNAME-sorted PE); `metheor pdr|mhl|lpmd|me|pm|fdrp|qfdrp -i x.bam`. epiBED only for biscuiteer `readEpibed` (skip unless required) | `XM` has Bismark's 8 letters only, so a user context (e.g. GCH) needs its own run mapped onto `Z/z`. Metheor is not on bioconda |
| **SLAM-seq: bakR, EZbakR** | **cB**: `sample, <feature cols>, TC, nT, n` (CSV or parquet partitioned by sample) | `bakR::bakRData(cB, metadf)`; `EZbakR::EZbakRData(cB, metadf)` or `EZbakRArrowData(arrow::open_dataset(dir), metadf)` | Mutation positions (the cU / `--mutPos` table) are not kept. GRAND-SLAM: no intermediate; compare its NTR/half-life output |
| **SLAMdunk-style per-UTR counts** | same cB | a SQL `GROUP BY feature` over cB (not a tool) | Read-span coverage semantics |
| **A-to-I: REDItools3, JACUSA2** | **REDItools table** (14 columns, 1-based, Strand `+/-/*`, BaseCount `[A,C,G,T]`) | `reditools annotate rna.tsv dna.tsv`; `reditools index rna.tsv -r regions.bed`. JACUSA2: no intermediate (replicate DM test runs on BAMs) | JACUSA2's statistics. REDItools' Frequency definition must be reproduced exactly |
| **SHAPE/DMS: RNA Framework, SEISMIC, DRACO** | **RNA Framework MM** (binary, per read) | `rf-mmtools index x.mm && rf-mmtools toRC x.mm` → `rf-norm`; `seismic importmm x.mm` → `seismic filter/cluster/table/fold`; DRACO `--mm x.mm` | Mutation identity (substitution vs indel type) and quality. Indel placement is fixed by the writer |
| **Ancient DNA: mapDamage** | `misincorporation.txt` + `dnacomp.txt` (+ `lgdistribution.txt`) | `mapDamage --plot-only -d DIR`; `mapDamage --stats-only -d DIR` | Rescaling (needs BAM). DamageProfiler and pydamage take BAMs |
| **CLIP: htseq-clip / DEWSeq** | per-read crosslink BED6 (one line per read, htseq-clip `extract` layout) | `htseq-clip count -i sites.bed.gz -a windows.txt.gz -o counts.tsv.gz` → `htseq-clip createMatrix` → DEWSeq | None beyond `extract`'s name field (`qname|length`) |
| **CLIP: iCount / clippy** | aggregated crosslink BED6 (`chr pos pos+1 . count strand`) | `iCount peaks annot.gtf sites.bed peaks.bed`; `clippy -i sites.bed …` | PureCLIP: BAM only, compare its outputs |

**Net writer set:**
1. cytosine report;
2. per-cell site parquet (the native table);
3. cB;
4. REDItools table;
5. MM;
6. mapDamage tables;
7. crosslink BED6, per read and aggregated.

Optional writers: ALLC, bedGraph and ENCODE bedMethyl, each only if its converter route
proves too slow or a consumer demands it. The read-level methylation tools use the
alnbase-tagged BAM.

---

## 2. Analyses that need more than a base-resolution pileup

A pileup here means per (sample/cell, contig, pos, strand, context) counts. Anything that
links observations *within one read or fragment* cannot come from it. alnbase hit rows keep
`qname`, the mate and `off_5p/off_3p`, so every item below is one `GROUP BY qname` (or
fragment) away.

### 2.1 Needs read- or fragment-level aggregation

| Analysis | Why a pileup cannot give it | Existing tool that consumes alnbase output instead |
|---|---|---|
| SLAM-seq / TimeLapse **per-read conversion counts** (`TC`, `nT` per read → cB; reads with ≥k conversions per gene; per-read conversion-rate distributions, `alleyoop rates`/`tcperreadpos`) | Mixture models need the joint (TC, nT) per read | bakR/EZbakR read the cB, but **building the cB is ours** (read-level GROUP BY). No maintained tool builds cB from a non-NGM BAM except `mut_call.py` with its own implicit filters. SLAMdunk needs NGM `MP` tags |
| **Incomplete-conversion read filters** (premethyst mCH fraction pooled per fragment; ScaleMethyl per mate, zero-CH kept; Bismark `filter_non_conversion` absolute/percentage) | The per-read non-CpG retention must be known before summing sites | Bismark `filter_non_conversion` on the alnbase XM BAM covers Bismark's rule. premethyst `bam-extract -B` covers premethyst's rule (and writes cellInfo). ScaleMethyl's per-mate rule and arbitrary contexts need a semi-join in DuckDB |
| **cellInfo** columns `n_frag, n_pair, n_single, xpct_m` (fragment counts, fragments with zero CH) | Fragment counts | premethyst `bam-extract -B` only, with CG/CH contexts from XM letters and barcode from qname. Otherwise ours. `cov, cg_cov, mcg_pct, ch_cov, mch_pct` do come from the per-cell pileup |
| **Epiallele metrics**: PDR, MHL/lMHL, VEF, LPMD, entropy, epipolymorphism, FDRP/qFDRP; allele-specific methylation | Need each read's multi-CpG pattern | epialleleR, Metheor, biscuiteer, from the XM BAM. **No aggregation code needed** |
| **SHAPE/DMS** co-mutation clustering (DRACO, SEISMIC `cluster`); per-read filters (max mutations per read, mutation gaps) | Need per-read mutation vectors | MM writer (read-level), then SEISMIC `filter` applies the read filters |
| **Hyper-edited reads** (≥n A>G per read); editing tools' per-read mismatch-quality filters (SPRINT) | Per-read counts | Ours (GROUP BY). SPRINT's hyper-editing needs remapping, out of scope |
| **aDNA** PMD per-read scores (PMDtools); read-length distribution | Per read | PMDtools and pydamage take BAMs. `lgdistribution.txt` is a trivial GROUP BY |
| **Position-in-read profiles**: M-bias, misincorporation by distance from ends, CLIP truncation sites | Not per read, but need `off_5p/off_3p` or read starts, which a site pileup drops | Bismark extractor `--mbias_only` on the XM BAM. For the rest, a GROUP BY offset over hit rows |

### 2.2 Windowed or regional consumers: all have existing windowing

| Consumer | Existing windowing from base resolution |
|---|---|
| Amethyst | `facet agg -u/-v`; R `makeWindows` |
| ALLCools MCDS | `generate-dataset --regions name binsize|bed` |
| MethSCAn/scbs | `methscan matrix regions.bed` (per-cell `total_sites`/`methylated_sites`, which is exactly scMET's `total_reads`/`met_reads` long input after a melt) |
| scMET | no windowing of its own (input is a Feature×Cell long table, `scMET/vignettes/scMET_vignette.Rmd:90`), covered by `methscan matrix` |
| methylKit | `tileMethylCounts`, `regionCounts` |
| bsseq / DSS / dmrseq | `bsseq::getCoverage(BS, regions, type, what="perRegionTotal")`; DMR callers segment themselves |
| methylSig | `tile_by_windows`, `tile_by_regions` |
| bigWig / bedGraph tracks | `bigWigAverageOverBed`, `bedtools map -c 4,5 -o sum`; `allc-to-bigwig --bin_size` |
| SLAM-seq genes | feature tags from `featureCounts -R BAM` (per read, carried into cB) |
| CLIP windows | `htseq-clip createSlidingWindows` + `count` |
| Editing index over Alu or regions | `reditools index`; RNAEditingIndexer |

**Conclusion:** no consumer found lacks a windowing path. alnbase-agg does not need a
windowing feature.

---

## 3. Ground-truth simulators

Directional bisulfite and EM-seq only. Behaviour below comes from reading source and docs;
none of these simulators has been run here yet. `src/...` paths are relative to
`clones/`.

### 3.1 Bisulfite / EM-seq

| Tool (version or commit read) | Bioconda | Ground truth | Golden BAM | Variants / indels | PE | Legibility / status |
|---|---|---|---|---|---|---|
| **BSReadSim** 0.4.0 (`wbvguo/BSReadSim` `9833bd7`, 2026-09-02) | **No** (docs: "coming soon"). Source build: `git clone --recurse-submodules` + `pip install .`, needs CMake and C++17 | **Per read, per base.** `zt:Z` has one character per SEQ base encoding context (CG/CHG/CHH), methylated, conversion succeeded, variant-affected and sequencing error. `zr:B:S` holds per-read counts, including conversion failures. `--save-truth` writes a per-site MethDB and the variant VCF (`docs/outputs/index.md`) | **Yes.** `--format bam` places reads at the true origin with an "indel-aware, query-complete CIGAR", MAPQ 60, and Bismark `XG/XR/YS` | SNVs and indels ≤4 bp (random or VCF). Phased diploid, allele-specific methylation, context changes caused by variants | Yes | Directional by default (`--undirectional` optional). Conversion is a per-base Bernoulli draw (`--conversion-rate`, default 0.998). Sequencing errors are substitutions only. Active, MIT, preprint [doi:10.1101/2024.12.24.627620](https://doi.org/10.1101/2024.12.24.627620). One maintainer |
| **Mason** 2.0.13 (SeqAn 2.5.2, `apps/mason2`) | **Yes**, `mason 2.0.13` (built for x86-64-v3, which needs AVX2) | **Per site only.** The `mason_methylation` FASTA holds `/TOP` and `/BOT` levels per contig in 0.0125 steps. **No per-read conversion truth** (confirmed in [seqan#2546](https://github.com/seqan/seqan/issues/2546)) | **Yes.** `-oa out.bam` with origin tags `oR/oP/oH/oS` and `XE/XS/XI` (`mason_options.cpp:1150, 1195-1212`). The **CIGAR is recomputed** by global realignment to the origin interval (`mason_simulator.cpp:203-241`), so indel placement is optimal rather than the simulated history, and conversions count in NM/MD. Bisulfite positions were wrong before [PR #2547](https://github.com/seqan/seqan/pull/2547) (2025-02, in SeqAn ≥2.5.0) | SNPs, indels and SVs via `mason_variator` → `-iv` | Yes | `--enable-bs-seq --bs-seq-protocol directional --bs-seq-conversion-rate` (default 0.99; `mason_options.cpp:276-289`). Each C: methylated with p = level, else converted with p = rate (`simulate_base.cpp:141-168`). Sequencing errors include indels. Stable, low memory |
| **Sherman** (`FelixKrueger/Sherman` `c6f1c9b`) | No (single Perl script) | `--truth_set` gives `chrom pos converted_CG|CH` for each conversion, **with no read ID**. Reverse-read coordinates were wrong before an unreleased fix (issue #15). One global CG rate and one CH rate; context taken from the read; a last-base C counts as CHH (`Sherman:828`) | No (origin in the read name only) | `-s` adds per-read random substitutions; no indels | Yes | **No `--pbat`**, which corrects `methylation-tools.md` §12. Keep only because the nf-core E. coli fixtures were made with it |
| **BSBolt Simulate** (`NuttyLogic/BSBolt` `ea4870e`, 2023) | No (`pip install bsbolt`) | Per read, per base, but encoded in the FASTQ comment (`chr:start:end:<codes>:W|C`; `SimulateMethylatedReads.py:159-235`) | No; turning the comment into a CIGAR would be custom work | SNVs and indels via a forked wgsim | Yes | **100 % conversion** of unmethylated C (`:147`). Dormant |
| MethylFASTQ (2020), DNemulator (site returns 404), RRBSsim (Python 2, RRBS), BSSim (2012), pWGBSSimla, WGBSSuite | No | Site-level at best | BSSim writes SAM | varies | varies | Unmaintained |
| "SimMethyl" | — | No tool by that name found | | | | |

**Gaps shared by all of them.** None simulates:
- per-molecule incomplete conversion (whole unconverted reads, the target of the premethyst and ScaleMethyl filters);
- overlapping mates with conflicting bases;
- duplicates, or secondary and supplementary records;
- N-context and contig-edge hazards.

EM-seq-specific failures are not modelled either. EM-seq reads are simply directional
reads with high conversion. These cases come from hand-written fixtures (§3.4), or from
editing a BSReadSim golden BAM with a short pysam script (for example, restoring the
original C at every `zt`-converted base of chosen reads).

**Recommendation.**
- **BSReadSim** is the primary simulator. It is the only one with a golden BAM plus per-base methylation, conversion, variant and error truth, which fills `truth_state` in `calls.parquet` directly. Pin it by commit in a source-built environment until a bioconda recipe appears, and verify its CIGARs on a tiny run first.
- **Mason** (bioconda) is the cross-check for site-level "every tool vs truth" comparisons. Setting `--bs-seq-conversion-rate 1` and zero Illumina error probabilities makes the read base at each reference C equal the methylation state, which gives per-read truth too. Compare indel-adjacent calls with care, because of the realigned CIGARs.

### 3.2 Non-methylation panel

| Assay | Recommended tool (bioconda) | Ground truth | Golden alignment | Notes |
|---|---|---|---|---|
| SLAM-seq | **splash**, inside `slamdunk 0.4.3` (entry point `splash=slamdunk.splash:run`, `setup.py:164-167`). The bioconda package called `splash` is a different tool | Per read: `TC:i:<n>` and read name `<utr>_<idx>_<nTC>`. Per UTR: T count, T>C count, converted reads, half-life (`dunks/simulator.py:194-229, 262-284`). **No conversion positions** | No (unmapped BAM; origin is in an intermediate BED) | 3′-UTR reads only; T>C is a Bernoulli draw (`-tc 0.024`); `splash eval-counts`. bakR `Simulate_bakRData` and `EZSimulate` simulate count tables, not reads |
| Ancient DNA | **NGSNGS 0.9.2** | Per read: `_modV1V2V3V4` name suffix (deamination / misincorporation / indel / error); `-DumpIndel`, `-DumpVCF` | **Yes** (`-f bam`) | Briggs model with gargammel's parameters. With errors off, every C>T/G>A is damage. CIGAR and strand not yet checked. Alternative **gargammel 1.1.4** lists per-base `_DEAM:<positions>` (`deamSim.cpp:2086`) but has no golden BAM |
| SHAPE / DMS | **SEISMIC-RNA `seismic sim`** (0.26.0) | Per-position, per-cluster mutation-rate parameters, plus per-read relationship vectors (`idmut`). Substitutions and deletions, **no insertions** (`sim/muts.py:182-221`) | No (FASTQ; `sim/fastq.py`) | `sim fold` needs RNAstructure or ViennaRNA unless a CT file is given. Deletion placement matters for parity (X3) |
| SNVs | **Mason** `mason_variator` → `mason_simulator -iv` | Truth VCF plus origin tags | Yes (realigned CIGAR) | Same binaries as the methylation cross-check. Alternative NEAT 4.7.0 (golden BAM + VCF, heavier) |
| A-to-I editing | **BAMSurgeon `addsnv`** (1.4.1); no maintained editing simulator exists | VCF with per-site `VAF`; `--tagreads` marks edited reads `BS:i:1` (`markreads.py:13`) | No (edited reads are realigned; `--aligner STAR` supported) | Spike A>G (+ genes) or T>C (− genes) into a real or simulated RNA BAM |

**Per-base truth without a truth-tagging simulator.** Turn off sequencing errors, and
variants where they don't matter. Every read-vs-reference difference left in a golden BAM
is then signal by construction. This does not work for bisulfite data with a conversion
rate below 1, where unconverted and methylated C look the same. That is why BSReadSim is
preferred there.

### 3.3 Hand-written edge-case fixtures

The most legible form is **SAM text plus a tiny FASTA**, with binaries built at test time
by samtools. Commit only text, so a reviewer sees the whole case on one screen.

```
tests/fixtures/c_before_deletion/
  ref.fa         >chrT / ACGTTCGAACGT...   (tens of bp)
  reads.sam      @HD VN:1.6 SO:unsorted, @SQ SN:chrT LN:40, records with explicit CIGAR and XG/YD tags
  expected.tsv   qname  mate  pos  expected_call
  README         1-3 lines: the hazard and which tool behaviour it probes
```

Build (samtools 1.24):
1. `samtools faidx ref.fa`
2. `samtools view -b --no-PG -o reads.bam reads.sam`. This fails if CIGAR and SEQ lengths disagree, which is a free check.
3. `samtools sort` + `index`, or `samtools sort -n` for Bismark's extractor.
4. `samtools calmd -b` only when a tool needs MD/NM.
5. `samtools fixmate -m` if mate fields (RNEXT/PNEXT/TLEN) should be filled rather than typed.

Use pysam only for geometries that are painful to type (long reads, generated
overlap/dovetail grids), and keep it emitting **SAM text lines** from a small table
(`pysam.AlignedSegment.fromstring(line, header)`) rather than building binary records in
code.

---

## 4. Workflow, benchmarking and reporting

Context from the coordinator:
- This machine is for building confidence. Full runs happen on a university cluster that runs **Nextflow with only the local executor**.
- Benchmarks must cover runtime, memory, disk, and the less obvious OS limits (threads, open files, memory maps), scaled across core counts.
- `validation/tools/resource_monitor.py` already samples a command's process tree with psutil: processes, threads, open fds, memory maps, RSS, plus per-process maxima for the per-process limits.

### 4.1 Checked on this machine
- **Not installed:** snakemake, nextflow, conda, mamba, micromamba, pixi, quarto, hyperfine.
- **Present:** GNU `/usr/bin/time`, kernel 6.8, cgroup v2, ext4.
- **cgroups:** The development account has no user systemd session (`systemd-run --user` fails), so `memory.peak` cannot be used here.
- **Rootless env creation works (verified by running):**
  - micromamba 2.9.0: `micromamba create -p envs/samtools -c conda-forge -c bioconda --strict-channel-priority samtools=1.22` took 88 s. The env is 124 MB; the root prefix and cache grew to 1.6 GB.
  - pixi 0.81.0: one `pixi.toml` with channels `conda-forge, bioconda`, one feature and environment per tool (`no-default-feature = true`) and a single `pixi.lock`. `pixi install -e methyldackel` worked, and `pixi run -e methyldackel MethylDackel --version` returned 0.6.1.
- **Current versions:** snakemake 9.27.0 and nextflow 26.04.6 (bioconda); quarto 1.9.38, papermill 2.7.0, hyperfine 1.20.0, micromamba 2.9.0, pixi 0.80.0 (conda-forge).

### 4.2 What each engine records

**Snakemake 9.27.0 `benchmark:`** (`benchmark.py`, `executors/local.py:366-393`)
- **Columns:** `s, h:m:s, max_rss, max_vms, max_uss, max_pss, io_in, io_out, mean_load, cpu_time` (`:30-42`). `--benchmark-extended` adds `jobid, rule_name, wildcards, params, threads, cpu_usage, resources, input_size_mb` (`:44-54`).
- **How it measures:** psutil polls the job shell and `children(recursive=True)`.
  - Memory is summed over the tree at each sample, and the max over samples is kept, in MiB.
  - `io_in`/`io_out` are the last sample's `read_bytes/write_bytes` divided by 1024², so **MiB** (the docs say bytes).
  - `cpu_time` is the last-seen user+sys of each PID.
  - `mean_load` is a percentage (400 ≈ 4 cores).
- **Sampling schedule:** every **0.5 s for ~30 samples, then every 30 s**, with no final sample at exit (`:20-23, 282-285, 504-506`). The timer includes conda activation.
- **Caveats:**
  - Spikes are missed, and after 15 s the window is 30 s.
  - Summed RSS double-counts shared pages.
  - `cpu_time` and I/O can miss up to 30 s at the end.
  - Very short jobs give `NA`.
- **Other features:**
  - `repeat("b.tsv", n)` runs a job n times.
  - `temp()` deletes a file once all its consumers finish, so sizes must be measured inside the producing rule.
  - `--report` builds a self-contained HTML report.
- **micromamba:** since 8.20.6 Snakemake calls a hard-coded `conda` (`conda.py:48, 735`). `--conda-frontend mamba` is accepted but ignored. **Released Snakemake does not work with micromamba.** It needs `conda`, which pixi can install rootless.

**Nextflow trace** (`command-trace.txt`, `TraceRecord.groovy:70-114`; collected by the `.command.run` wrapper, so the local executor works too)
- **Fields:** `realtime, duration, %cpu, %mem, rss, vmem, peak_rss, peak_vmem, rchar, wchar, syscr, syscw, read_bytes, write_bytes, vol_ctxt, inv_ctxt, cpus, memory, cpu_model, hostname, exit, attempt, workdir, …`
- **`peak_rss`/`peak_vmem`:** per-process kernel high-water marks (`VmHWM`/`VmPeak`), summed over the tree at each sample. Samples come every 1 s ×10, then every 5 s, then every 30 s.
  - Because the kernel keeps each peak, a spike inside a process that is still alive is caught.
  - Children that start and exit between samples are missed.
  - Summing peaks reached at different times overestimates.
- **`%cpu` and I/O:** differences of `/proc/$$/stat` and `/proc/$$/io` before and after the task. These include reaped children, so they are exact.
- **Conda:** `conda.useMicromamba = true` is supported. The `conda` directive also accepts the path of an existing environment, so pixi-built prefixes (`.pixi/envs/<tool>`) can be used, if the directory form is confirmed with micromamba activation.

**Neither engine records threads, open fds, memory maps, output disk size or rlimits.** So:

### 4.3 resource_monitor.py: complement, not replace

- **Keep `resource_monitor.py` as the per-command wrapper.** Neither Snakemake's benchmark nor Nextflow's trace replaces it, because neither records threads, fds or mmaps, which are the numbers behind cluster failures (`RLIMIT_NPROC` counts threads, `RLIMIT_NOFILE` and `vm.max_map_count` apply per process). The engine's record is a cross-check for wall time, CPU and I/O. Nextflow's `rchar/wchar/read_bytes/write_bytes` and `%cpu` are exact because they are before/after deltas, which a sampler cannot match.
- **Small changes would make it robust** (each a few lines, still psutil):
  1. Also record **`VmHWM`** from `/proc/<pid>/status`, per process and as a max over processes. A 0.1 s RSS sampler misses allocation spikes; VmHWM does not, for processes alive at sampling.
  2. Wrap the command in **GNU `/usr/bin/time -f '%e %U %S %M %P %c %w %I %O'`**. `%M` is `wait4` `ru_maxrss` (KiB) of the **largest single process**, which is exact even for sub-second jobs and catches children reaped between samples.
     - The fork verified that two parallel 60 MB children report about 70 MB, not 120 MB: max, not sum.
     - Record both numbers: the sampled tree sum is an estimate of node pressure, `%M` is the exact per-process peak.
  3. Record the **limits in force** at start: `psutil.Process().rlimit(RLIMIT_NOFILE|NPROC|AS)`, `/proc/sys/vm/max_map_count`, `nproc`, and the CPU affinity (`os.sched_getaffinity(0)`).
  4. Record **disk** after exit: `du --apparent-size --block-size=1 -s OUT` (logical bytes; counts hard links once; sparse holes not counted as allocated) and `du --block-size=1 -s OUT` (allocated). For peak temporary disk, set `TMPDIR` to a per-task directory and sample its apparent size in the same loop. That is the only way to see scratch files a tool deletes itself.
  5. Sampling `/proc/<pid>/maps` line by line is slow for processes with tens of thousands of maps. Reading `/proc/<pid>/smaps_rollup` is not a substitute (no count), so keep the interval ≥0.1 s and report the sampler's own CPU time.
- **Core-count scaling.** Pin each run to its allocation, with `taskset -c 0-$((n-1))` (or `numactl --physcpubind`) around the tool, as well as passing `-t n`. The existing demo (`validation/demos/resource-limits/README.md`) shows that library thread pools size themselves from visible cores, not from the tool's `-t`. The local executor does not enforce `cpus`. Run benchmark tasks one at a time (Nextflow `maxForks 1` on benchmark processes; Snakemake `resources: bench=1` + `--resources bench=1`). Warm the page cache with one discarded run, because dropping caches needs root.
- **hyperfine** (1.20) is only for repeated wall-time with warm-up (`--warmup 1 --runs 5 -L threads 1,2,4,8,16 --prepare 'rm -rf out' --export-json`). Its `memory_usage_byte` is a cumulative `RUSAGE_CHILDREN` maximum across **all** runs (`unix_timer.rs:30-66`), so it is wrong across a parameter sweep. Do not use it for memory.

### 4.4 Orchestrator

| | Nextflow 26.04 (local executor) | Snakemake 9.27 | Plain Python runner |
|---|---|---|---|
| Runs on the cluster as provided | **Yes** | Needs a rootless install (pixi) and permission | Yes (Python) |
| Legibility | Groovy DSL; outputs in hashed `work/`, `publishDir` for readable paths | Most legible; readable output paths | Explicit, but all custom |
| Pinned envs | `conda` (+ `useMicromamba`), env path, containers; nf-core modules exist for Bismark, MethylDackel, bwa-meth, samtools | `conda:` + pin files (needs `conda`), apptainer | custom (`pixi run -e`) |
| Resume, cleanup | `-resume`; `cleanup = true` removes `work/` after success (no per-file `temp()`) | built-in; `temp()` | custom |
| Built-in metrics | trace (above), `-with-report`, `-with-timeline` | benchmark (above), `--report` | none |

**Recommendation: Nextflow DSL2**, because it is what the cluster runs and the same
pipeline must give the laptop's "tiny" results and the cluster's "full" results.
Keep it legible:
- one process per tool invocation, whose `script:` is a short shell command wrapped as `resource_monitor.py --out … -- /usr/bin/time -o time.txt -f … taskset -c … <tool>`;
- `publishDir` to `results/<validator>/<dataset>/<cores>/`;
- a `params` sweep channel for core counts;
- environments from one `pixi.toml`/`pixi.lock` (tools as separate environments), referenced by path, or `conda.useMicromamba` with explicit `name=version=build` specs.

Snakemake would be the choice if the cluster did not dictate the engine. Its `benchmark:` directive adds nothing that the wrapper and Nextflow trace do not already cover. A plain Python runner would reimplement caching and resume, against the "use existing tools" priority.

### 4.5 Reporting: Quarto

Use **Quarto** (conda-forge 1.9.38) with the Python/Jupyter engine:
- A `report.qmd` reads the parquet result tables (`metrics.parquet`, `sites.parquet`, the benchmark TSVs) with DuckDB or polars and draws figures with matplotlib.
- It renders with `quarto render report.qmd -P scale:full --to html` (`embed-resources: true` gives one self-contained file) or `--to typst` / `pdf` for publication figures.
- Parameters use a `#| tags: [parameters]` cell. `freeze: auto` avoids re-execution.
- It runs as the final Nextflow process.

papermill + nbconvert would also work, but `.ipynb` JSON with outputs is noisy in git and
PDF export needs a separate LaTeX or Chromium setup. The report is two tools there instead
of one.

### 4.6 Sources for §§3–4
- **Simulators:**
  - Mason: https://github.com/seqan/seqan/tree/main/apps/mason2
  - BSReadSim: https://github.com/wbvguo/BSReadSim, https://wbvguo.github.io/BSReadSim/
  - Sherman: https://github.com/FelixKrueger/Sherman
  - BSBolt: https://github.com/NuttyLogic/BSBolt
  - SLAMdunk/splash: https://github.com/t-neumann/slamdunk
  - NGSNGS: https://github.com/RAHenriksen/NGSNGS
  - gargammel: https://github.com/grenaud/gargammel
  - SEISMIC-RNA: https://github.com/rouskinlab/seismic-rna
  - BAMSurgeon: https://github.com/adamewing/bamsurgeon
  - NEAT: https://github.com/ncsa/NEAT
- **Snakemake:** v9.27.0 `benchmark.py`, `conda.py`; https://snakemake.readthedocs.io (Benchmark Rules, reporting, installation)
- **Nextflow:** `command-trace.txt`, `TraceRecord.groovy`; https://docs.seqera.io/nextflow/reports, https://docs.seqera.io/nextflow/reference/config (conda scope)
- **Measurement:**
  - hyperfine v1.20.0: `src/timer/unix_timer.rs`
  - getrusage(2), time(1)
  - https://docs.kernel.org/admin-guide/cgroup-v2.html
- **Quarto:** https://quarto.org/docs/computations/parameters.html, https://quarto.org/docs/projects/code-execution.html
