# Single-cell methylation output formats: Amethyst H5, premethyst cellInfo, and others

Research date: 2026-09-16. Everything below comes from source code cloned into
`clones/` unless marked otherwise. Paths are given relative to that directory.

| Source | Where it came from | Commit / version |
|---|---|---|
| Amethyst R package | `github.com/lrylaarsdam/amethyst` (the canonical repo; `adeylab/amethyst` returns 404) | `cff9f08` (2026-01-29, v1.0.5) |
| Facet (`amethyst-facet`) | PyPI sdist `amethyst_facet-1.2.4.tar.gz` (`adeylab/amethyst-facet` returns 404; no public repo found). Author Ben Skubi | 1.2.4 |
| premethyst | `github.com/adeylab/premethyst` | `fa5a47b` (2024-11-22) |
| sciMETv2 (older Adey scripts) | `github.com/adeylab/sciMETv2` | `705f850` (2024-03-12) |
| ScaleMethyl (Scale Bio Nextflow) | `github.com/ScaleBio/ScaleMethyl` | `2a4973c` (2025-07-22) |
| ALLCools | `github.com/lhqing/ALLCools` (dir `ALLCools-lhqing`) | `c9f7be2` (2025-02-21) |
| MethSCAn (successor to scbs) | `github.com/anders-biostat/MethSCAn`; scbs at `github.com/LKremer/scbs` | `72f38fb` / `a68e2a9` |
| BISCUIT | `github.com/huishenlab/biscuit` | `0a5ceae` (2026-08-18) |
| scMET | `github.com/andreaskapou/scMET` | `0aeea3f` |
| MethylStar | `github.com/jlab-code/MethylStar` | `e07d8a3` |
| Epiclomal | only via WebFetch (the clone was too large) | master |

---

## 1. Amethyst HDF5 files

### 1.1 The format has two versions

| | v1 (Amethyst 0.0.0.9000, old premethyst / sciMETv2) | v2.0.0 (Amethyst ≥ 1.0.0, Facet) |
|---|---|---|
| Base-level data path | `/<context>/<barcode>` (a dataset directly under the context group) | `/<context>/<barcode>/1` (`<barcode>` is a group; `1` is the dataset) |
| Base-level fields | `chr S10, pos int64, pct float64, t int64, c int64` | `chr S10, pos <i8, c <i8, t <i8` (`pct` is dropped) |
| Window aggregations | nowhere to put them | `/<context>/<barcode>/<name>` with `chr S10, start <i8, end <i8, c <i8, t <i8, c_nz <i8, t_nz <i8` |
| Version marker | none | dataset `/metadata/version` = `"amethyst2.0.0"` |

Citations:
- **v1 writer:** `premethyst/premethyst_commands/calls2h5.py:34-41`, which is identical to `sciMETv2/sciMET_cellCalls2h5.py`:
  ```python
  cov = np.genfromtxt(file_path, delimiter='\t', dtype=[('chr', 'S10'), ('pos', int), ('pct', float), ('t', int), ('c', int)])
  cov = np.sort(cov, order=['chr', 'pos'])
  cg_group.create_dataset(cell_id, data=cov, compression='gzip', compression_opts=9)
  ```
  The groups `CG` and `CH` are always created (lines 25-26 of the py file).
- **Facet dtypes:** `amethyst_facet/h5/dataset.py:14-17`:
  ```python
  observations_v1_dtype = [("chr", "S10"), ("pos", "<i8"), ("pct", "<f8"), ("c", "<i8"), ("t", "<i8")]
  observations_v2_dtype = [("chr", "S10"), ("pos", "<i8"), ("c", "<i8"), ("t", "<i8")]
  windows_dtype = [("chr", "S10"), ("start", "<i8"), ("end", "<i8"), ("c", "<i8"), ("t", "<i8"), ("c_nz", "<i8"), ("t_nz", "<i8")]
  ```
- **Version string:** `amethyst_facet/h5/__init__.py:7` (`version="amethyst2.0.0"`) and `h5/handles.py:43-56`. The version is read with `file["/metadata/version"][()].decode()` and written with `file.create_dataset("/metadata/version", data=fct.h5.version)`. h5py stores a Python `str` as a variable-length UTF-8 scalar string.
- **Official schema description:** `amethyst_facet/cli/commands/convert.py:44-55`:
  > "The V2 format stores base-pair observations as (chr, pos, c, t) in an HDF5 dataset at /context/barcode/1. 1 is the conventional name for the bp-level unaggregated observations. The V2 format stores window aggregations in a dataset at /context/barcode/[window_dataset_name] as (chr, start, end, c, t, c_nz, t_nz) … It also contains a dataset /metadata/version='amethyst2.0.0'."
- **Diagram in the Amethyst README:** `amethyst/images/h5structure.png`. It shows `project.h5 → /CAC,/CH,/CG → /<barcode> → /1 (chr,pos,t,c)` and `/100000 (chr,start,end,t,c,t_nz,c_nz)`. The README text is at `amethyst/README.md:62-75` and `:115-116`, and the change is logged in `NEWS.md` v1.0.0: "base information in /context/barcode/1; aggregated information in /context/barcode/name".

### 1.2 Field semantics
- `c` is the number of methylated calls at (chr, pos) and `t` is the number of unmethylated calls. The docstring in Facet `cli/commands/calls2h5.py:72-92` says: "t is the count of unmethylated calls … c is the count of methylated calls".
- In premethyst, counts are per **fragment**, not per read. The two mates of a pair are merged into one hash before counting, so a site covered by both mates counts once (`bam_extract.pm:244-269`).
- `pct` exists only in v1 and ScaleMethyl files, and its scale is not consistent:
  - premethyst writes **percent** with 2 decimals: `sprintf("%.2f", meth/cov*100)` (`bam_extract.pm:280`).
  - ScaleMethyl writes a **fraction**: `ROUND(methylated / (methylated + unmethylated), 2) as pct` (`ScaleMethyl/bin/write_amethyst.py:41`).
  - Facet v1 output computes a fraction, `c/(c+t)` (`dataset.py:106-108`).
  - No Amethyst reader uses `pct`. `R/diff.R:297` drops it with `intersect(c("chr","pos","pct"), names(data)) := NULL`.
- Window fields (Facet `windows/windows_aggregator.py:26-38`):
  - `c` and `t` are sums over the observed positions in `[start, end)`.
  - `c_nz` is the number of positions with `c > 0`; `t_nz` is the number with `t > 0`.
  - Windows with no observations are not written.
  - Rows are sorted by `chr, start, end`.

### 1.3 Contexts
The context is just the name of the first-level group. Nothing restricts it to a fixed list.
- The image shows `/CAC`, `/CH` and `/CG`.
- premethyst `context-extract` makes new contexts such as `GCH`, `HCH`, `GCG`, `HCG` and `CAC` (`context_extract.pm:29-35`).
- Facet `calls2h5` takes the context from the file name or from the ScaleMethyl parquet `context` column (`calls2h5.py:181-218`).
- Amethyst R functions take `type = "<context>"` and read `paste0(type, "/", barcode, "/1")`. For contexts other than CG/CH, `makeWindows` with `metric` "score" or "ratio" needs a metadata column named `m<context lowercased>_pct` (`R/helper.R:74-78`, `:148`, `:228`). For example, context `GCH` needs `mgch_pct`. (The message in the source says "cac_pct", but the code reads `paste0("m", tolower(type), "_pct")`.)

### 1.4 Cells per file, and how the R object finds them
- A file can hold any number of cells. The R object has a data.frame `obj@h5paths` with columns `barcode`, `path` and an optional `prefix`. Each function reads `h5read(path, paste0(type,"/",barcode,"/1"))`; see `R/index.R:37-47` and `:61`, and `R/facet.R:45-55` and `:88`.
- The cell ID used in metadata is `paste0(prefix, sub("\\..*$", "", barcode))`, so everything after the first `.` in the barcode is ignored (`R/helper.R:145`, `:220`; `R/facet.R:87`).
- How each pipeline groups cells into files:
  - premethyst writes one file per run or sample, holding all cells (`premethyst calls2h5 <folder> <prefix>`).
  - ScaleMethyl writes one file per sample × Tn5 well: `<sample>.<tgmt_well>_cov.h5` with many barcodes (`ScaleMethyl/bin/write_amethyst.py:3,24`; `amethyst/R/prepare.R:101`, which is `createScaleObject`).
  - The PBMC vignette uses a single file for all cells (`vignettes/pbmc_vignette/pbmc_vignette.Rmd:124`).
- A barcode that appears in several files, or in several experiments, is handled with `prefix` (NEWS v1.0.0).

### 1.5 `indexChr` and the sorting requirement
`R/index.R:58-78`:
```r
h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", bar, "/1")))
h5[, index := 1:.N]
chrs <- as.list(unique(h5$chr[!grepl("_|EBV|M", h5$chr)]))   # default when chrList is NULL
ind <- ind[, .(cell_id = unique_id, chr = x, start = min(index), count = .N)]
```
The output is a list keyed by chr of data.frames `(cell_id, chr, start, count)`. `makeWindows` stores it in `obj@index[["chr_cg"]]`, then reads **only that hyperslab**:
`rhdf5::h5read(path, name = paste0(type,"/",barcode,"/1"), start = sites$start[...], count = sites$count[...])` (`R/helper.R:151-153`, `:222-224`, `:327-329`; `R/diff.R:293-295`; `R/metacells.R:648`, `:942`).

**Hard requirement:** all rows for one chr must be **contiguous** in `/1`, because the index stores only `min(index)` and `count`. Every writer sorts by `(chr, pos)` (lexicographic bytes for chr, numeric for pos), so that is the safe convention. Sorting by pos within a chr is not needed for correctness: window joins are non-equi joins and uniform windows use `round_any`. Still, all writers sort that way.

The default chr filter drops any chr name containing `_`, `EBV` or `M`, so `chrM` is excluded, and so is any other name containing "M". `makeWindows(stepsize=…)` also drops window names with more than 3 `_` parts or matching `chrEBV|chrM|KI` (`R/helper.R:270`).

Deprecated v1 functions: `indexGenes` and `getGeneM` in `R/deprecated.R:134-275` read `paste0(type,"/",bar)` (v1 layout), plus `start`/`count`.

### 1.6 Fields the R reader actually uses
- `makeWindows` (`R/helper.R:155-157`, `:225-230`, `:331-333`) uses `chr`, `pos`, `c`, `t`. The value is **binarized per site**:
  ```r
  value = round(sum(c != 0) / (sum(c != 0) + sum(t != 0)), 3); n = sum(c + t, na.rm = TRUE)   # keep n >= nmin (default 2)
  h5[pos %% stepsize == 0, pos := pos + 1]    # windows are (k*step, (k+1)*step], 1-based-friendly
  window := paste0(chr, "_", round_any(pos, stepsize, floor), "_", round_any(pos, stepsize, ceiling))
  ```
  A bed-file join uses `pos >= start & pos <= end` (closed on both ends).
- `calcSmoothedWindows` (`R/diff.R:293-299`) uses raw sums `c = sum(c)` and `t = sum(t)`.
- `loadWindows` (`R/facet.R:88-96`) reads `/<type>/<barcode>/<name>`, deduplicates rows with `unique(h5)`, and computes `value = sum(c_nz)/(sum(c_nz)+sum(t_nz))` with `n = sum(c_nz + t_nz) >= nmin` (default 10). So **`c_nz`/`t_nz` drive the per-cell matrix**, and `c`/`t` are used only by `loadSmoothedWindows` (`R/facet.R:255`: `c = sum(c), t = sum(t), n = sum(c_nz + t_nz)`).
- Metrics "score" and "ratio" need `obj@metadata[cell, "mcg_pct"]` or `"mch_pct"`, a percent that comes from cellInfo (`R/facet.R:92`).

### 1.7 What Facet adds (`amethyst-facet` 1.2.4)
Commands: `agg`, `convert`, `delete`, `calls2h5`, `version` (see the README in PKG-INFO).
- **`facet convert`** turns v1 into v2 (`cli/commands/convert.py`). It rewrites every `/ctx/bc` into `/ctx/bc/1` and writes `/metadata/version`.
- **`facet calls2h5`** (`cli/commands/calls2h5.py`):
  - Inputs: premethyst `.cov`, ScaleMethyl per-cell `.parquet` (columns `chr,pos,context,methylated,unmethylated`), and other v2 `.h5` files.
  - Output dtype: `AMETHYST_H5_DTYPE = [('chr','S10'),('pos',int),('t',int),('c',int)]` (line 22). The field order `t,c` differs from `dataset.py`; readers look fields up by name.
  - Rows are sorted `chr,pos` (line 23).
  - The default `.cov` column map is `chr=0,pos=1,pct=2,t=3,c=4` (lines 25-31, 582-586).
  - Compression is gzip level 6 by default (lines 580-581).
  - It raises an error if the target exists without `/metadata/version` or with a different version (lines 336-355).
  - Dataset name conflicts can be set to ERROR, OVERWRITE or SKIP.
- **`facet agg`** (`cli/commands/agg.py`, `windows/*.py`) reads v2 `/ctx/bc/<obs>` (default: every dataset whose fields include `chr,pos`) and writes window datasets beside it:
  - Uniform windows `-u [name=]size[:step][+offset]`:
    - `start = (pos - offset) // size * size + offset`, `end = start + size` (`uniform_windows_aggregator.py:55-63`).
    - **The default offset is 1** (`cli/parse/uniform_windows_parser.py:111`), so windows are `[1+k*size, 1+(k+1)*size)` on 1-based positions.
    - **The default name is `"{size}:{step}+{offset}"`**, for example `100000:100000+1` (`uniform_windows_aggregator.py:42-43`).
    - The README's "stored in `/[context]/[barcode]/[window_size]`" and the Amethyst diagram (start 0, name `100000`) show older behaviour. Example in `agg --help`: `--uniform-windows 10000=10000+0` gives name `10000` with windows starting at 0.
    - The step must divide the size.
  - Variable windows `-v [name=]path.tsv`: a CSV/TSV with header `chr,start,end`, sniffed by DuckDB. Positions count when `start <= pos < end`, and overlaps are allowed. The default name is the file stem (`variable_windows_aggregator.py:119-149`, `variable_windows_parser.py:37-51`).
  - Compression defaults to gzip 6 (`cli/decorators/decorators.py:10-23`), with no explicit chunking, so h5py picks auto chunks.
- **`facet delete context|barcode|dataset`**.
- Facet itself depends on duckdb, polars, h5py and numpy<=1.26.4.

### 1.8 Compression and chunking
No writer sets `chunks=`; h5py auto-chunks whenever compression is on.
- premethyst, sciMETv2 and ScaleMethyl: `compression='gzip', compression_opts=9`.
- Facet: gzip 6.

No HDF5 attributes are used anywhere; the version is a dataset, not an attribute.

### 1.9 Coordinates
- **premethyst: 1-based.** `pos` starts at SAM column 4 (`$pos = $P[3]`, `bam_extract.pm:180`) and moves along the XB/XM string (`:330-384`).
  - For reads on the reverse strand, the call is at the reference coordinate of the G, which is the C on the − strand. **Strands are not merged**: a CpG shows up as two positions, pos (+) and pos+1 (−).
  - In Bismark mode (`-B`), pos is incremented once per XM character without looking at the CIGAR, so indels move later calls.
  - In BSBolt mode, digits in XB are skip counts.
- **ScaleMethyl: 0-based**, inferred from the code rather than tested.
  - The BSBolt path uses `aligned_pairs[pos_offset[idx]-1]`, which is a pysam 0-based reference position (`met_extract.py:108-118`).
  - The bwa-meth path uses `read.reference_start + y` (`:225`).
  - `pos` is written to parquet, Amethyst h5, `.allc` and bismark `.cov` without a +1 (`write_amethyst.py:39-44`, `write_allc.py:24-31`). The ALLC spec says 1-based.
  - ScaleMethyl keeps a `strand` column but Amethyst files do not.
- The Amethyst R code does not care much about the convention. `makeWindows` stepsize windows behave as `(k*step, (k+1)*step]`, which suits 1-based data. **Recommendation for alnbase: write 1-based `pos`**, to match premethyst, the ALLC/Bismark convention and Facet's default offset=1.

### 1.10 Gotchas for an exporter
1. `chr` is `S10`, so numpy **silently truncates** names longer than 10 bytes (e.g. `chrUn_KI270302v1`). This can merge contigs or break contiguity. Filter or rename such contigs, or check the length.
2. Keep each chr contiguous (sort by chr, pos).
3. The `/1` dataset name is fixed by convention and hard-coded in R.
4. Field order is free, but the names `chr,pos,c,t` (and `start,end,c,t,c_nz,t_nz`) must match exactly. `loadWindows` needs `c_nz`/`t_nz`.
5. Write `/metadata/version` = `"amethyst2.0.0"` (scalar string) or Facet will refuse to append.
6. Barcodes must not contain `/`. Avoid `.`, because R strips everything after the first `.` when it builds cell IDs.
7. ScaleMethyl's `write_amethyst.py` writes v2 paths (`/CG/<bc>/1`) but **keeps `pct` and has no `/metadata/version`**. So "v2 layout" files seen in practice are not always Facet-valid.
8. For alnbase's sharding: a partition key of cell barcode puts all of a cell's rows in one shard, so each shard can write its own H5, or its own `/ctx/bc/*` subtree. The R object then needs a (barcode, path) table. premethyst `bam-extract` has the same requirement: reads for a barcode must be contiguous (name-sorted with barcode prefix), because it reopens `>$barc.meth` on every barcode change (`bam_extract.pm:163-177`).

A minimal writer that was run and checked is in `write_amethyst_v2_minimal.py`, next to this file. It ran with h5py 3.6.0 / numpy 1.21.5 / HDF5 1.10.7, and `h5ls -rv` shows the expected compound types. It was not tested against R/rhdf5, which is not installed.

### 1.11 Metacell H5 (a separate format)
`R/metacells.R:1045-1360` writes a job JSON for an **external Python** tool, which is not in the repo. The tool writes `metacell_windows.h5` with `/<type>/chr`, `/<type>/start`, `/<type>/end` (1-D) and 2-D matrices `/<type>/{pct,c,t}` of windows × metacells, in either orientation (`buildMetacellH5Index`, `readH5WindowsSubset`). This is not the per-cell format.

---

## 2. premethyst: cellInfo files and intermediate files

### 2.1 Pipeline steps (`premethyst_commands/run_pipeline.pm:18-20, 176-288`; README)
1. **Demultiplex with `unidex`.** The read name becomes `@<CELL_BARCODE>:NNN#0/1`. The barcode is the qname up to the first `:`, and the fragment ID is the qname up to `#` (`bam_rmdup.pm:72-73`, `bam_extract.pm:155-156`).
2. **`fastq-trim`** with TrimGalore.
3. **`fastq-align`**: BSBolt `Align -F1 R2 -F2 R1` (the reads are swapped), then `samtools sort -n` (`fastq_align.pm:76, 95`).
4. **`bam-rmdup`** (`bam_rmdup.pm`):
   - Keeps reads with `samtools view -q 10` (MAPQ ≥ 10; `:14`, `:61`).
   - Within each barcode, keeps the **first read seen at each (chr, POS)**, ignoring strand, and also keeps its mate by fragment tag (`:91-105`).
   - Writes `<O>.bbrd.q10.nsrt.bam`.
   - Writes **`<O>.complexity.txt`** with columns `rank, barcode, total_reads(MAPQ≥10), unique_reads, pct_unique(2dp)`, sorted by unique reads descending (`:116-122`).
5. **`plot-complexity`**.
6. **`bam-extract`** (`bam_extract.pm`) writes a cellCall folder `<O>/` plus `<O>.cellInfo.txt`, described below.
7. Optional: `calls-rename`, `calls-filter`, `context-extract`.
8. **`calls2h5 <folder> <prefix>`** writes a v1 `<prefix>.h5`.

### 2.2 Per-cell intermediate files
- `<O>/<barcode>.meth` is temporary and deleted at the end. It holds one line per fragment: `chr pos XB [chr pos XB]` for a pair, or 3 columns for a single read (`:183-197`).
- **`<O>/<barcode>.CG.cov`** and **`<O>/<barcode>.CH.cov`** are tab-separated, have no header, and are sorted `-k1,1 -k2,2n` (`:236-237`). Columns (`:282`, `:289`):
  `chr  pos(1-based)  pct(percent, %.2f)  t(unmethylated count)  c(methylated count)`
  - CG contains BSBolt `x/X` calls (Bismark `z/Z`). CH merges CHG and CHH: BSBolt `y/Y/z/Z`, Bismark `x/X/h/H` (`:337-369`).
  - `-H` suppresses the CH file (the CH stats are still computed).
- `context-extract` output, `<outdir>/<barcode>.<NAME>.cov` (and `<barcode>.<XNAME>.cov` for the reciprocal), uses the same 5 columns.
  - It keeps rows whose `chr.pos` appears in a BED "sites" file. Matching uses the BED **column 2 as-is**, plus an optional `-o` offset (`context_extract.pm:45-50, 66-74, 115-134`).
  - The Adey site-BED generator writes 1-based positions in both columns 2 and 3 (`sciMETv2/sciMET_find_GpCs.pl:27-38`: `chr  pos  pos  strand`).
  - A GpC site gives two rows: C at pos (+) and the G at pos-1 (−). This is how GCH/NOMe contexts are defined: **by reference site list, not by reading the context from the read**.
- `premethyst fasta-context` (a general motif-to-BED tool) is not finished and dies on start (`fasta_context.pm:36`).
- `bam-extract -E` calls `premethyst stream-h5`, but no such command exists in the repo (`bam_extract.pm:126-133`).

### 2.3 `<O>.cellInfo.txt` columns
Written per cell at `bam_extract.pm:296-303`; tab-separated, **no header**, 10 columns:
```perl
$cellInfo = "$ARGV[1]\t$AllCov\t$CG_cov\t$CGpct\t$CH_cov\t$CHpct\t$frag_ct\t$pair_ct\t$single_ct\t$excluded_mCH";
```
| # | Amethyst name (`R/helper.R:482-488`) | Computation (`bam_extract.pm`) |
|---|---|---|
| 1 | `cell_id` | barcode (qname up to the first `:`) |
| 2 | `cov` | `CG_bases + CH_bases`: total number of CG and CH **calls** (the sum over sites of per-site coverage), after the read filter (`:283, 290, 296`) |
| 3 | `cg_cov` | `CG_cov`: number of **distinct CG positions** covered (`$CG_cov++` per coord) |
| 4 | `mcg_pct` | `CG_meth / CG_bases * 100`, `%.2f` ("0.00" if there are no bases) (`:298`) |
| 5 | `ch_cov` | number of distinct CH positions covered |
| 6 | `mch_pct` | `CH_meth / CH_bases * 100`, `%.2f` (`:299`) |
| 7 | `n_frag` | `single_ct + pair_ct`: fragments read from `.meth`, **before** the mCH filter (`:297`) |
| 8 | `n_pair` | paired fragments |
| 9 | `n_single` | unpaired reads |
| 10 | `xpct_m` | `$excluded_mCH`. Because of how the code branches, this actually counts **fragments with zero CH calls**, which are also thrown away. Fragments that fail the mCH-fraction test are dropped *without* being counted (see §4.1). |

Details:
- Amethyst `addCellInfo` also accepts the older 6-column form (`cell_id,cov,cg_cov,mcg_pct,ch_cov,mch_pct`), and names only the first N columns if the count is something else (`R/helper.R:472-497`).
- The vignette file `amethyst/vignettes/pbmc_vignette/pbmc_vignette_cellInfo.txt` has 6 columns and no header. Example row: `ACGCGACGGCACGAGAATCACTGTCATG 27258073 923300 82.21 21341679 0.54`.
- `context-extract` writes `<outdir>.<NAME>.cellInfo.txt` with 4 columns: `barcode, n_sites, total_calls(t+c), mC_pct(%.2f)` (`context_extract.pm:86, 136-139`).
- `sciMETv2/sciMET_cellCall2info.pl` rebuilds cellInfo from a `.cov` folder with header `#CellID Coverage CG_Cov CG_mC_Pct CH_Cov CH_mC_Pct [<ctx>_Cov <ctx>_mC_Pct …]`. There, `*_Cov` = number of sites and Coverage = the sum of t+c over all contexts.
- **The ScaleMethyl "cellInfo.txt" looks similar but is different** (`ScaleMethyl/bin/met_extract.py:438-467`):
  - Header: `CellID, Coverage, CG_Cov, CG_mC_Pct, CH_Cov, CH_mC_Pct, CH_high`.
  - `<ctx>_Cov = SUM(methylated)+SUM(unmethylated)`, which counts **calls, not sites**.
  - `CH_high` = number of reads dropped by the CH filter.
  - `generate_metrics.py:118-129` reads it as `cell_id,cov,cg_cov,mcg_pct,ch_cov,mch_pct,ch_high_reads` and merges it into `<sample>.allCells.csv`. That file adds `total_reads, unique_reads, passing_reads, pass, ch_high_reads_percent`, well columns and more. Amethyst's `createScaleObject` keeps rows with `pass == "pass"` (`R/prepare.R:95`).

---

## 3. Other single-cell methylation formats

### 3.1 ALLC (per cell) and ALLCools MCDS (zarr / xarray)
- **ALLC** (`ALLCools-lhqing/docs/allcools/start/input_files.md:45-55`):
  - 7 columns, no header, bgzip-compressed and tabix-indexed.
  - Columns: `chrom, pos (1-based), strand (+/-), context (e.g. CGT; can be longer than 3 bases), mc, cov, methylated (1 if no test)`.
  - One file per cell. The context is written per row as a reference k-mer.
- **MCDS**, built by `allcools generate-dataset` (`ALLCools/count_matrix/dataset.py`):
  - The input is `allc_table`: TSV `cell_id <tab> allc_path` with no header.
  - Output directory `<out>.mcds/`:
    - `chrom_sizes.txt` (`:318`).
    - `.ALLCools` YAML: `{region_dim, ds_region_dim, ds_sample_dim}` (`mcds/utilities.py:407-431`).
    - One **zarr group per region set**: `<out>/<region_name>/` (`:353-371`).
  - Each group holds the data variable `<region_name>_da` with dims `(cell, <region_name>, mc_type, count_type)`, where `count_type = ["mc","cov"]` (`:185-186`). The dtype is `uint32`, clipped at the max (`:211, 227-230`).
  - Region coordinates are `<region_name>_chrom`, `<region_name>_start`, `<region_name>_end` (`:361`).
  - Cell IDs are fixed-width `<U{maxlen}` (`:305`).
  - Region IDs:
    - For bins, `chrom_i` (e.g. `chr1_0`), made with `bedtools makewindows` (0-based BED).
    - For a BED file, column 4, or `<name>_<i>` (`:72-125`).
  - Region counting uses `tabix.fetch(chrom, start, end)` over the ALLC file and sums mc/cov by context (`:56-58, 166-181`).
- **Context handling:** each `mc_type` is an IUPAC pattern (e.g. `CGN`, `CHN`, `GCYN`, `HCHN`). `parse_mc_pattern` expands it into all concrete k-mers and sums them (`count_matrix/dataset.py:62-68`, `utilities.py:19-35, 68-79`). **This is the closest existing model for user-defined contexts.**
- **Quantifier types:**
  - `count`
  - `hypo-score` / `hyper-score`: `float16` binomial survival-function values, zeroed below a cutoff of 0.9, stored as `<region>_da_<mc_type>-hypo-score` with dims `(cell, region)` (`:193-206, 232-254`).
  - `hyper-score` is suggested for GpC/NOMe (`input_files.md`).

### 3.2 scbs / MethSCAn `prepare`
`MethSCAn/methscan/prepare.py`; `scbs/scbs/prepare.py` is essentially the same.
- **Input:** one coverage file per cell. The cell name is the file name without extension.
  - Formats (`:309-357`):
    - `bismark`: 0-based column indices chr=0, pos=1, meth=4, unmeth=5; no header.
    - `allc`/`methylpy`: meth=4, cov=5, **with a header line skipped**.
    - `biscuit`, `biscuit_short`.
    - A custom format `chr:pos:meth:<n>{c|u}:sep:header` (1-based columns).
- **Contexts:** none. The input is assumed to be already CpG-only, so you supply one context per run.
- **Per-cell site binarization** (`:139-149`):
  - Only `n_meth==0` or `n_unmeth==0` counts as a clean site; the value is `+1` if methylated, otherwise `-1`.
  - Sites with mixed calls are dropped, unless `--round-sites` is given; with it, the majority wins and ties are dropped.
- **Output directory:**
  - `<chrom>.npz`: scipy CSR `int8` of shape `(max_pos+1, n_cells)`. **The row index is the genomic position** as given (1-based for Bismark); values are 1, -1, or 0 for missing (`:41-51, 195-218`).
  - `column_header.txt`: cell names in column order.
  - `cell_stats.csv`: `cell_name, n_obs, n_meth, global_meth_frac` (`:240-252`).
  - `run_info.txt`.
  - Temporary COO text files `<chrom>_chunk%07d.coo` hold `pos,cell_idx,value`.
- **`methscan matrix`** writes `methylation_fractions.csv.gz`, `mean_shrunken_residuals.csv.gz`, `total_sites.csv.gz` and `methylated_sites.csv.gz`, or sparse `matrix.mtx.gz`/`features.tsv.gz`/`barcodes.tsv.gz` (`matrix.py:149-174`).

### 3.3 BISCUIT epiread / epiBED
`biscuit/src/epiread.c`
- The default output is epiBED (`:192-270`). Columns: `chr, start(0-based), end, read_name, read_number(1|2), bisulfite_strand(+|-), CG_RLE, GC_RLE ('.' unless NOMe -N), variant_RLE`.
- RLE symbols (`:29-41`):
  - `M`/`U`: methylated or unmethylated CpG.
  - `O`/`S`: open or shut GpC accessibility in NOMe mode.
  - `F`: filtered. `x`: ignored. `P`: soft clip. `D`: deletion. `-`, `i`, `d`: skips.
  - `R`/`Y`: ambiguous.
- `-O` gives the legacy epiread format: `chr, read, bsstrand, …, comma lists of 0-based CpG/GpC positions and call strings` (`:285-420`). `-P` gives pairwise output.
- **Contexts are fixed:** CpG, plus GpC with `-N`. In NOMe mode, GCG is ambiguous and HCG/GCH are separated.
- Read-level filters (default values at `bisc_utils.h:101-115`; flags in `epiread.c:1197-1211`):

  | Filter | Flag | Default |
  |---|---|---|
  | min base quality | `-b` | 20 |
  | min MAPQ | `-m` | 40 |
  | min AS | `-a` | 40 |
  | max cytosine retention per read | `-t` | 999999 (off) |
  | min read length | `-l` | 10 |
  | distance to 5′ end | `-5` | 3 |
  | distance to 3′ end | `-3` | 3 |
  | max NM | `-n` | |
  | proper pair, dup, secondary, qcfail | | filtered |
  | overlapping-mate double counting | | avoided |

  "Retention" is the count of retained Cs in all reference-C (or G) positions of a read (`bisc_utils.c:76-99`, applied at `epiread.c:646-647`, `pileup.c:777-778`).

### 3.4 scMET (R)
`scMET/R/scmet.R:13-15`: the input `Y` is a long `data.table` with **4 named columns: `Feature, Cell, total_reads, met_reads`**. "total_reads" and "met_reads" really mean the number of CpGs observed and the number methylated per cell × feature; binarizing per CpG first is expected. Features are user-defined regions. There is no context concept (CpG is assumed). `utils.R:33-87` converts to and from sparse matrices.

### 3.5 Epiclomal
Via WebFetch from `shahcompbio/Epiclomal/process_real_data/`.
- Intermediate CpG-level TSV: `chr, CpG_start, CpG_end, region_start, region_end, region_cpgNum, region_length, region_id, meth_frac, count_meth, count_unmeth, cell_id`.
- Final inputs from `get_data_ready_Epiclomal.R`:
  - `input_Epiclomal_<id>.tsv.gz`: a cells × CpGs matrix with first column `cell_id` and CpG columns named `chr:CpG_start`. Values are binarized: `0` if meth_frac < 0.5, `1` if > 0.5, empty/NA if exactly 0.5.
  - `regionIDs_input_Epiclomal_<id>.tsv.gz`: `region_id, start, end`, the 0-based **column-index** range of each region's CpGs.
- CpG only.

### 3.6 MethylStar
- Bismark-based and mostly bulk; it has a scBS mapping mode with `--pbat --se` (`src/bash/bismark-mapper-scBS-Seq.sh:64-98`).
- Outputs are per sample:
  - Bismark CX reports.
  - METHimpute output with columns `seqnames, start, strand, context, counts.methylated, counts.total, posteriorMax, posteriorMeth, posteriorUnmeth, status, rc.meth.lvl`, and `context.trinucleotide` (`src/bash/methimpute.R:33-36`).
  - methylKit format (`chrBase` etc., `methylkit.R:25`).
  - DMRcaller RData (`dmr-caller.R:45-55`).
  - bedGraph/BigWig per context.
- **Contexts are fixed to CG/CHG/CHH** (`meth-bedgraph.R:21`).
- Not really a single-cell format; low priority.

### 3.7 ScaleMethyl per-cell outputs (also worth supporting)
- `<sample>.met_<ctx>.parquet`: `barcode, chr, pos(UInt32, 0-based per code), strand, methylated, unmethylated`, ZSTD (`met_extract.py:414-432`).
- Intermediate per-chromosome parquet files also carry a `context` column.
- Barcode = `qname` regex `:([ACGT]+\+[ACGT]+\+[ACGT]+)$`.
- Per-read calls are deduplicated on `(qname,pos)`, so overlapping mates count once (`:33-61`).
- ScaleMethyl contexts are only CG or CH, decided per call from the XB letter.
- Optional per-cell exports: ALLC (`write_allc.py`: `chr,pos,strand,'CG'|'CH',mc,cov,1`, where the context column is a literal "CG"/"CH", not a trinucleotide), Bismark `.cov` (`write_bismark.py`: `chr,pos,frac(0-1, 2dp),meth,unmeth`; note this uses a fraction where Bismark uses a percentage), and Amethyst h5.
- Genome-bin matrices: `<sample>.CG.score.mtx.gz`, `<sample>.CH.mtx.gz`, and `features.tsv`/`barcodes.tsv` in MatrixMarket format.

---

## 4. Filtering conventions (for the export utility's filter language)

### 4.1 Read / fragment level (incomplete bisulfite conversion)
| Pipeline | Rule | Default | Citation |
|---|---|---|---|
| premethyst `bam-extract` | Keep a fragment (both mates pooled) only if `read_CH > 0` **and** `read_mCH/read_CH <= M`. Fragments with **zero CH calls are dropped** and counted in col 10; fragments above the threshold are dropped silently. | `-M 0.4` ("For brain, rec: 0.4; for non mCH cell types, rec: 0.1") | `bam_extract.pm:19, 38-39, 255-273` |
| ScaleMethyl `met_extract.py` | Drop a read (per mate, not per fragment) if `total_ch > 0 and mCH/total_ch > threshold`; count these as `CH_high`. Reads with no CH calls are kept. | `chReadsThreshold = 50` (%) → 0.5 | `nextflow.config:68`, `modules/dedup_and_extract.nf:58`, `met_extract.py:102-106, 300-303` |
| BISCUIT | Absolute count of retained Cs per read (`-t`); off by default | 999999 | `bisc_utils.c:76-99` |
| ALLCools | Nothing per read; per-cell `mCCCFrac` is used as a proxy for non-conversion | — | §4.3 |

CH as defined by premethyst and ScaleMethyl is every non-CpG C (CHG plus CHH) on the read's own converted strand. Both use the aligner's XB/XM context letters, not a reference lookup; ScaleMethyl's bwa-meth path does look at the reference.

### 4.2 Alignment / duplicate level
| Pipeline | Rule |
|---|---|
| premethyst `bam-rmdup` | MAPQ ≥ 10 (`-q`). Dedup key is (barcode, chr, leftmost POS), strand-agnostic, keeping the first read and its mate (`bam_rmdup.pm:14, 61, 91-105`). |
| ScaleMethyl | `sc_dedup --duplicate-key Leftmost --min-mapq 10` (`nextflow.config:84-85`, `dedup_and_extract.nf:31`) |
| BISCUIT | MAPQ ≥ 40, baseQ ≥ 20, 3 bp trimmed from each end, proper pairs, no dup/secondary/qcfail |
| ALLCools `bam-to-allc` | `min_mapq=10`, `min_base_quality=20` (`_bam_to_allc.py:355-356`) |

### 4.3 Cell level
- **premethyst `bam-extract` pre-filter** from `complexity.txt`: `unique_reads >= -N` (default 10000), `pct_unique <= -P` (100) and `>= -p` (0) (`bam_extract.pm:13-16, 90-101`).
- **premethyst `calls-filter`** (`calls_filter.pm:25-30, 74-83`) works on cellInfo columns:
  - `-G` minimum CG sites (`cg_cov`) and `-H` minimum CH sites (`ch_cov`).
  - `-g MIN,MAX` for `mcg_pct` and `-h MIN,MAX` for `mch_pct` (percent).
  - `-C` lists contexts to keep.
  - Failing `.cov` files are *moved* to the fail folder.
  - Bugs: it runs `opendir` on the cellInfo *file* instead of the calls folder (`:85`), so as written it probably moves nothing, and the pass/fail counters are never incremented.
- **Amethyst vignettes** filter in R on the metadata:
  - `cov > 1000000 & mch_pct < 12` (brain/doublet/batch vignettes, e.g. `vignettes/brain_vignette/brain_vignette.Rmd:131`).
  - `cov > 100000 & cov < 40000000` (PBMC, `pbmc_vignette.Rmd:155`).
  - Text: "We recommend cells have a minimum of 1M cytosines covered" (`:151`).
  - Here `cov` means total CG+CH calls.
- **ScaleMethyl cell calling** (`bin/generate_metrics.py:150-169`, `docs/analysisParameters.md:54-66`):
  - `pass` if `unique_reads >= threshold` and `minUniqTotal <= percent_unique <= maxUniqTotal` and `unique_reads <= max_uniq`.
  - `threshold = max(minUniqCount=1000, P99(unique_reads of barcodes ≥ minUniqCount) / minCellRatio=20)`, unless the user sets a fixed threshold.
- **ALLCools basic cell QC** (`docs/allcools/cell_level/step_by_step/100kb/01-CellBasicFiltering.ipynb`):
  - `MappingRate > 0.5`, `FinalReads > 500000`.
  - `mCCCFrac < 0.03`: mCCC is "the proxy of the upper bound of the non-conversion rate".
  - `mCHFrac < 0.2`, `mCGFrac > 0.5`.
- **MethSCAn `filter`** (`methscan/filter.py:38-57`; `cli.py:248-272`): `--min-sites/--max-sites` on `n_obs`; `--min-meth/--max-meth` (percent) on `global_meth_frac*100`.

### 4.4 Site / feature level
- scbs/MethSCAn: drop sites with mixed calls, or majority vote with `--round-sites` (§3.2).
- Epiclomal: binarize at 0.5, with ties set to NA.
- Amethyst `makeWindows`:
  - per-site binarization `c != 0`, `t != 0`;
  - per-cell window kept if `n = sum(c+t) >= nmin` (default 2);
  - `loadWindows` requires `sum(c_nz+t_nz) >= nmin` (default 10).
- ALLCools hypo/hyper score: cutoff 0.9.

### 4.5 Filter primitives an alnbase export DSL needs
From the above, the export filter language should be able to express:

1. **Read/fragment predicates computed over hits:**
   - `n_calls(ctx)`, `n_meth(ctx)`, `frac_meth(ctx)` for any context set (e.g. CH, CHH, CCC);
   - pooling scope: per mate (ScaleMethyl) or per fragment/qname (premethyst);
   - comparison choice: `>` vs `<=`;
   - what to do when there are zero calls: keep (ScaleMethyl) or drop (premethyst);
   - plus MAPQ, flag (dup/secondary/qcfail/proper pair), base quality, distance from 5′/3′ ends (BISCUIT `off_5p/off_3p`, which alnbase already emits), NM, AS.
2. **Dedup** by (cell, chr, leftmost pos [, strand]) before counting.
3. **Site aggregation** per (cell, context, chr, pos [, strand]):
   - `c`, `t` (with mates or fragments counted once);
   - optional strand collapse;
   - binarization (`c>0`, `t>0`, majority, drop ties/mixed).
4. **Cell-level metrics:**
   - number of sites covered per context;
   - total calls per context;
   - `%m` per context as `100*sum(c)/sum(c+t)`;
   - number of fragments, pairs and singles, and reads dropped by the CH filter;
   - thresholds with min/max bounds.
5. **Context definition:**
   - by reference k-mer/IUPAC pattern with an anchored C position and strand (ALLCools-style `GCH`, `HCG`, `CCC`); or
   - by an explicit site BED (premethyst `context-extract`).
   - A "reciprocal" (complement) context should be expressible too.
