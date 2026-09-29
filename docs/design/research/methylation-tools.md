# Bisulfite / EM-seq methylation extractors: formats, built-in "queries", differences, logistics

Research for alnbase: (a) a validation suite against existing extractors, (b) an export utility that writes their formats.

Date: 2026-09-16. All claims about defaults come from reading the source. Repos were cloned into `clones/` at these commits:

| Repo | Commit / version read | Notes |
|---|---|---|
| FelixKrueger/Bismark | `ec5cb58` (2026-08-01). Perl scripts are **v0.25.1**. Rust suite is `rust/VERSION` **3.1.0** | Since v3.0.0 (2026-07-07) Bismark is one Rust binary that claims to be "byte-identical to Perl Bismark v0.25.1 on the faithful default path" (CHANGELOG.md). The Perl code is the reference implementation, so it is cited below. |
| dpryan79/MethylDackel | `3c77bda` (2023-07-05), **0.6.1** | |
| huishenlab/biscuit | `0a5ceae` (2026-08-18), **1.10.3-dev** | |
| lhqing/ALLCools | `c9f7be2` (2025-02-21), **1.1.1** (PyPI) | Original repo |
| DingWB/ALLCools | `4a1f8a1` (2026-07-22) | Active fork. Some defaults differ (noted below). |
| yupenghe/methylpy | `5bb95e0`, **1.4.7** | |
| guoweilong/cgmaptools | `afa6e40` (2023-12-22), 0.1.x | |
| brentp/bwa-meth | `730eb07`, **0.2.10** | |
| nf-core/methylseq | `cc11bd9` (2026-07-27) | |
| samtools/htslib | `fe2a822`, **1.24** | Used for pileup and overlap semantics |
| nebiolabs/EM-seq | HEAD | NEB reference EM-seq pipeline |

Citations take the form `repo/path:line`, relative to those commits. On GitHub they map to `https://github.com/<owner>/<repo>/blob/<commit>/<path>#L<line>`.

---

## 0. TL;DR: the implicit query each tool runs

Notation: the **read symbol** is the base in SEQ (BAM orientation). The **refr context** is taken from the reference unless stated otherwise. A "C-strand read" is a C→T converted read, i.e. OT/CTOT, Bismark `XG:Z:CT`, bwa-meth `YD:Z:f`, BISCUIT "BSW". A "G-strand read" is a G→A converted read, i.e. OB/CTOB.

| Tool | Informative position | Context source | Methylated / unmethylated |
|---|---|---|---|
| **Bismark** (aligner writes XM; extractor only parses XM) | C-strand reads: ref base C with read base C or T. G-strand reads: ref G with read G or A | **Read-projected reference**: the reference bases that align to the read, **with deleted reference bases removed** and **inserted or soft-clipped read bases replaced by an `X` placeholder**, plus **2 extra reference bases past the read end**. Any `N` or `X` in the context gives `U/u`, and U/u calls are **never output** by the extractor | C (or G) = meth, T (or A) = unmeth |
| **MethylDackel** | Pileup column where ref=C (C-strand reads) or ref=G (G-strand reads); `is_del` and `is_refskip` skipped | **Pure reference** (faidx) at the column. N counts as H. Chromosome end counts as CHH | Same symbols. Base quality ≥5, MAPQ ≥10 |
| **BISCUIT pileup** | Same idea, via its own CIGAR walk | **Pure reference 5-mer** centered on the C. Any N in the 5-mer (including the 2 **upstream** bases, or chromosome-edge padding) gives context `CN` (NA) | Same symbols. BQ ≥20, MAPQ ≥40, 3 bp trimmed at both SEQ ends. A per-site SNP veto applies |
| **ALLCools / methylpy** | `samtools mpileup` column with ref C, counting `.` (fwd match) and `T`; with ref G, counting `,` (rev match) and `a` | Pure reference, `num_downstr_bases=2` (3-mer string, N kept literally) | Strand is taken from the **alignment orientation**, not a tag (optionally rewritten from XG/YZ first) |
| **CGmapTools bam2cgmap** | Old samtools 0.1.18 pileup: fwd reads at ref C, rev reads at ref G | Pure reference. N gives context `--`. Within 2 bp of a chromosome edge gives `--` | No BQ or MAPQ filter |

---

## 1. Bismark (v0.25.1 Perl = v3.x Rust)

### 1.1 Pipeline pieces and what each reads

- `bismark_genome_preparation` builds C→T and G→A converted genomes (`Bisulfite_Genome/CT_conversion`, `GA_conversion`) and indexes them with **bowtie2** (default), HISAT2, or (3.x) minimap2 or rammap. Bioconda `bismark=3.1.0` depends on `bowtie2>=2.5.5, hisat2, minimap2`. Legacy Perl: `bismark=0.25.1`. The Rust suite does all BAM I/O itself and needs **no samtools**. The Perl version shells out to samtools.
- `bismark` (aligner) makes the **methylation call at alignment time** and writes it into BAM tags `XM:Z:` (call string), `XR:Z:` (read conversion CT/GA) and `XG:Z:` (genome conversion CT/GA). It also writes NM and MD, and **does not write AS or MC tags** (`Bismark/bismark` SAM writer around lines 8689–8711: fields `NM, MD, XM, XR, XG[, XA][, RG]`).
- `bismark_methylation_extractor` **only parses XM/XR/XG plus CIGAR/POS**. It never looks at the reference, base qualities, MAPQ or FLAG (grep: no flag or qual logic; `Bismark/bismark_methylation_extractor:1576–1604` SE, `:1867–1919` PE).
- `deduplicate_bismark` **removes** duplicates (it does not mark them).
- `bismark2bedGraph` turns per-call files into a bedGraph and a `.bismark.cov`.
- `coverage2cytosine` turns `.cov` plus the genome FASTA into a genome-wide cytosine report, optionally with merged CpGs.
- `filter_non_conversion` drops reads with too many methylated non-CG calls, based on XM.

### 1.2 How the aligner assigns calls (the real "query")

`Bismark/bismark: extract_corresponding_genomic_sequence_single_end` (4273–4470) and `methylation_call` (4800–5020):

1. It builds `unmodified_genomic_sequence` by walking the CIGAR:
   - `M`: reference bases appended.
   - `I` and `S`: `'X' x len` appended (4339, 4352).
   - `D`: **nothing appended** (only `pos += len`, 4362).
   - `N`: skipped.
   - `H`, `P`, `=`, `X` in the CIGAR: die.
2. It **pads 2 reference bases**. For index 0 and 2 (OT, CTOB) they go after the aligned end. For index 1 and 3 (CTOT, OB) they go before the start. If the read sits within 2 bp of a chromosome edge, the padding cannot be extracted and **the read or pair is discarded** ("Sequences which were discarded because genomic sequence could not be extracted", check at 3127).
3. Calls for CT reads (`methylation_call`, 4822–4893):
   - The reference position must be `C` and the read base `C` (meth) or `T` (unmeth). Any other read base gives `.`.
   - Next genomic char `G` gives `Z/z`.
   - Next char `N` or `X` gives `U/u`.
   - Otherwise, the char after that: `G` gives `X/x`, `N` or `X` gives `U/u`, anything else gives `H/h`.
   - GA reads are mirrored: the genomic position must be `G`, read `G`/`A`, and context is read upstream.
   - **Consequences:**
     - A C next to a **deletion** takes its context from the first reference base *after* the deleted block. For example, ref `C G T` with the G deleted in the read is called as CHH/CHG, not CpG. MethylDackel and BISCUIT would call it CpG.
     - A C whose next base in the read is an **insertion or soft clip** becomes `U/u` and is dropped.
     - Past the read end, context comes from the reference, so it is "readable past read end".
4. The XM string is aligned to SEQ as stored in the BAM. For reverse-strand alignments the extractor reverses it back to 5'→3' read order (`:1640`, PE `:1939–1945`).
5. Strand classification (`:1611–1631`), using **read 1's** XR plus XG:
   - CT/CT = OT, index 0, "+"
   - GA/CT = CTOT, index 1, "−"
   - GA/GA = CTOB, index 2, "+"
   - CT/GA = OB, index 3, "−"
   - Library type is fixed at alignment time. Default is directional (OT and OB only). `--non_directional` gives all 4. `--pbat` gives CTOT and CTOB only.

### 1.3 Extractor defaults (`process_commandline`, 921–1347)

| Setting | Default | Source |
|---|---|---|
| SE vs PE | Auto-detected from the `@PG ID:Bismark` header (`-1` and `-2` present means PE) | 1149–1179 |
| Overlap (PE) | `--no_overlap` **ON** for PE ("only methylation calls of read 1 are used"). `--include_overlap` turns it off | 1214–1226, help 5823–5831 |
| `--ignore` / `--ignore_r2` / `--ignore_3prime` / `--ignore_3prime_r2` | 0 / 0 / 0 / 0. Bases are counted from the **5' end in read orientation** after the XM reversal. POS is recomputed from the CIGAR | 1190–1212, 1650–1800 |
| Base quality | **none** | |
| MAPQ | **none** (the aligner reports only unique best hits) | |
| FLAG filters | **none**: dup, secondary, supplementary, QC-fail and proper-pair are all ignored. Input must be Bismark output: PE records must be adjacent R1 then R2 (checks QNAME equality, `test_positional_sorting` 1371–1432). **The first record of each pair is treated as read 1 regardless of FLAG** | 1855–1893 |
| Contexts written | Z/z (CpG), X/x (CHG), H/h (CHH). **U/u silently ignored** | 2940–2943 |
| bedGraph `--cutoff` | 1 | 937, 1277–1284 |
| `--comprehensive` | off. The default is 12 strand-specific files `{CpG,CHG,CHH}_{OT,CTOT,CTOB,OB}_<name>.txt` | |
| `--merge_non_CpG` | off | |
| `--CX` | off (bedGraph and cytosine report are CpG only) | |
| `--zero_based` | off | |
| Header line | `Bismark methylation extractor version v0.25.1` as line 1 of each call file unless `--no_header` | 5077–5103 |
| nf-core/methylseq defaults | `--ignore_r2 2 --ignore_3prime_r2 2 --no_overlap` (`methylseq/nextflow.config`, `conf/modules/bismark_methylationextractor.config`) | |

**How PE overlap removal works** (`print_individual_C_methylation_states_paired_end_files`, 2813–4227):

- `end_read_1` is R1's 3'-most reference coordinate, computed **after** `--ignore`/`--ignore_3prime` trimming of R1 (2388–2417).
- R2 is walked **from its 5' end**. At the first R2 position that reaches `end_read_1` (`>=` for a "+" R2 walk, `<=` for a "−" walk, e.g. 2905 and 2987) the function **returns**, dropping the rest of R2.
- Effects:
  - R1 always wins.
  - All R2 bases in the overlap are dropped, whatever their quality.
  - Dovetailed R2 bases that extend past R1's 5' end are also dropped.
  - Overlap-dropped R2 positions are **not counted in M-bias either**, so the R2 M-bias is truncated.
- M-bias positions are indices into the call string **after** `--ignore` trimming (`$mbias_1{..}->{$index+1}`), so they shift when `--ignore` is used.

### 1.4 Output formats

**a) Per-call files** (`CpG_OT_*.txt`, `CpG_context_*.txt` with `--comprehensive`, `Non_CpG_*` with `--merge_non_CpG`), tab-separated. The first line is the version header unless `--no_header`. File naming is at 5085–5143.
```
<read_id> <'+' if methylated / '-' if unmethylated> <chrom> <1-based pos> <call char Z|z|X|x|H|h>
```
- **The `+`/`−` column is methylation state, not strand.**
- Position is 1-based. For G-strand reads it is the **G position**, i.e. the bottom-strand C.
- `--yacht` (SE only) adds `<read start> <read end> <read orientation>`.

**b) bedGraph** (`<name>.bedGraph.gz`), from `bismark2bedGraph:generate_output` (590–620):
```
track type=bedGraph
<chrom> <start 0-based> <end = start+1> <meth %>
```
**c) Coverage** (`<name>.bismark.cov.gz`), no header, **1-based closed** (`start == end`):
```
<chrom> <pos> <pos> <meth %> <count M> <count U>
```
- `--zero_based` adds `.bismark.zero.cov` with a 0-based start and end = start+1.
- `%` is written as Perl's **unformatted float** `($m/$t)*100`, e.g. `66.6666666666667`, `100`, `0`.
- Rows exist only where there is coverage ≥ `--cutoff`.
- **No strand column.** Top-strand C and bottom-strand C (the G of a CpG) are separate rows one position apart.
- Calls are pooled across OT/OB/CTOT/CTOB by position.
- Input lines are sorted with `sort -k3,3V -k4,4n` (434–449).
- Calls whose `+/−` does not match the letter case are warned about and skipped (`validate_methylation_call`, 553–588).

**d) Genome-wide cytosine report** (`coverage2cytosine`, or extractor `--cytosine_report`). File name `<name>.CpG_report.txt`, or `.CX_report.txt` with `--CX`. No header. From 168–460:
```
<chrom> <pos 1-based> <strand +|-> <count M> <count U> <context CG|CHG|CHH> <trinucleotide>
```
- Context comes purely from the reference FASTA, **uppercased** (`read_genome_into_memory`, `uc`).
  - `tri_nt` for C positions is `substr(pos-1,3)`. For G positions it is the reverse complement of `substr(pos-3,3)`.
  - Classification: `^CG` gives CG. `^C.G$` gives CHG, **so `CNG` counts as CHG**. `^C..` gives CHH, **so `CNN`/`CAN` count as CHH**. This is the opposite of the aligner, which drops N contexts as `U`.
  - Positions whose trinucleotide cannot be extracted (<3 chars) are skipped, and so is the **last base of each chromosome** (343–349).
- Coverage threshold default **0**, so **zero-coverage cytosines are included** (2185). `--coverage_threshold N` turns this off.
- `--ffs` adds tetra, penta and hexa-nucleotide columns.
- `--gc/--nome-seq` add GpC/NOMe outputs, and `--drach` does m6A.
- Coordinates are 1-based unless `--zero_based`.

**e) `--merge_CpG`** (`combine_CpGs_to_single_CG_entity`, 1753–1960). CpG-only mode is required. Output `*.merged_CpG_evidence.cov`:
```
<chrom> <C pos (1-based)> <G pos> <%.6f meth %> <M_top+M_bot> <U_top+U_bot>
```
- With `--zero_based`: `<C pos 0-based> <C pos+2>`.
- It reads the CpG report two lines at a time and requires `+` then `−` one bp apart.
- Pairs with zero total coverage are skipped.
- `--discordance_filter N` writes pairs whose top and bottom percentages differ by more than N into `*.discordant_CpG_evidence.cov` instead of merging them.
- `--merge_CpG` cannot be combined with `--threshold`.

**f) M-bias** (`<name>.M-bias.txt`): blocks titled `CpG context (R1)` then `================`, followed by a header `position\tcount methylated\tcount unmethylated\t% methylation\tcoverage`, with `%.2f`.

### 1.5 filter_non_conversion and deduplicate_bismark

- `filter_non_conversion` (XM-based):
  - Default `--threshold 3` methylated non-CG calls (H or X) per read. For PE, either read fails the whole pair (602).
  - `--consecutive` requires those calls to be consecutive (any z, h or x resets the counter).
  - `--percentage_cutoff P --minimum_count 5` (534) is the alternative criterion.
- `deduplicate_bismark`: the key is `FLAG:chr:start` for SE, where the reverse-read start is its 3'-most (5') coordinate. For PE the key is `FLAG_R1:chr:start_R1:end_R2` (292–500). The first occurrence is kept and the rest removed. `--barcode/--umi` add the UMI to the key.

### 1.6 Using the Bismark extractor on non-Bismark BAMs

- It needs `XM`, `XR`, `XG` on every record, and `XM` length equal to SEQ length.
- PE records must be adjacent (name-grouped, R1 first).
- CIGAR may only contain M/I/D/N/S (others die in `check_cigar_string`, 4228–4346).
- It needs `-s`/`-p` if the `@PG ID:Bismark` line is missing.
- **So alnbase can write Bismark-compatible XM/XR/XG tags and feed its BAM to the real extractor.** That tests alnbase's calls against Bismark's aggregation and formatting separately from its own aggregation.
- MethylDackel and BISCUIT can read Bismark BAMs (via XG) after coordinate sorting and indexing.

---

## 2. MethylDackel 0.6.1

### 2.1 Subcommands and outputs

**`MethylDackel extract <ref.fa> <sorted.bam>`**: needs a coordinate-sorted, **indexed** BAM or CRAM (it tries to build the index) plus a faidx-indexed FASTA. Outputs (`extract.c:1359–1437`):

- Default `<prefix>_CpG.bedGraph` (+ `_CHG`, `_CHH` with `--CHG`/`--CHH`). Header `track type="bedGraph" description="<prefix> CpG methylation levels"`, with ` merged` inserted under `--mergeContext` (`printHeader`, 562–569). Rows (`writeCall`, 39–99):
  ```
  <chrom> <start 0-based> <end> <int(100*M/(M+U))> <M> <U>
  ```
  - The percentage is **truncated to an integer, not rounded** (`(int)(100.0*m/(m+u))`).
  - Width is 1 per C. **No strand column.** A bottom-strand C is the G row.
  - Rows only where `M+U ≥ minDepth` (default 1).
- `--fraction` writes `_CpG.meth.bedGraph` as `chrom start end %f`.
- `--counts` writes `.counts.bedGraph` as `chrom start end M+U`.
- `--logit` writes `.logit.bedGraph` as `%f`.
- `--mergeContext` writes CpG rows `start=C pos (0-based)`, `end=start+2`, and CHG rows `end=start+3`, with the same 6 columns. A G-only covered CpG is reported at the C. `minDepth` applies to the merged site (`processLast`, 207–222; 483–491).
- `--methylKit` writes `<prefix>_CpG.methylKit` with header `chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT`. Rows use `%s.%i\t%s\t%i\t%c\t%i\t%6.2f\t%6.2f`: `chr.pos1`, chr, pos1 (1-based), `F` (ref C) or `R` (ref G), coverage, freqC and freqT as percentages. `%6.2f` **left-pads with spaces** (e.g. ` 50.00`).
- `--cytosine_report` writes a single `<prefix>.cytosine_report.txt` (all contexts requested) with no header: `chrom pos(1-based) +|- M U C{G|HG|HH} tnc`.
  - The trinucleotide comes from the reference, and unknown or off-chromosome bases become `N` (e.g. `CNN`, `CGN`).
  - **Zero-coverage Cs are included, and `--minDepth` is ignored** (42).
  - The last C of a chromosome **is** reported (as CHH/`CNN`), whereas Bismark skips it.

**`MethylDackel mbias <ref.fa> <bam> <prefix>`**:
- Writes SVGs `<prefix>_OT.svg` etc.
- `--txt` prints `Strand\tRead\tPosition\tnMethylated\tnUnmethylated` with 1-based positions (`svg.c:makeTXT`, 440–454).
- It prints `Suggested inclusion options: --OT a,b,c,d ...` to stderr.
- **Overlap removal is not applied in mbias** (`MBias.c` comment "bam_mplp_init_overlaps ... excluded here").

**`MethylDackel perRead`**:
- Output `readname chrom pos(0-based leftmost) %CpG-meth(%f) nInformative`.
- One row **per alignment**, so no mate merging or overlap handling.
- Defaults: `-q 10 -p 5`, and **`--ignoreFlags` default 0** (duplicates, secondary and supplementary all included). keepSingleton and discordant filters are not applied (`perRead.c:280–292`).

**`MethylDackel mergeContext <ref.fa> <bedGraph>`**: merges an existing per-C bedGraph (6-column form only).

### 2.2 Implicit query and defaults (`extract.c:724–752`, `common.c`)

| Aspect | Behaviour | Source |
|---|---|---|
| MAPQ | `-q 10` (≥10 kept) | extract.c:725, common.c:417 |
| Base quality | `-p 5` (bases with qual <5 dropped; must be ≥1) | 725, common.c:127 |
| FLAG | `-F 0xF00` (secondary, QC-fail, duplicate, supplementary ignored). Unmapped always ignored. `--keepDupes` clears 0x400. `-R/--requireFlags` 0 | 744, common.c:416–420, 1002–1004 |
| Multimappers | Alignment dropped if `NH:i:` >1, unless `--ignoreNH` | common.c:421–427 |
| Singletons | Dropped if `(flag&0x9)==0x9` (paired, mate unmapped), unless `--keepSingleton` | common.c:429 |
| Discordant / improper | Dropped if `(flag&0x3)==0x1` (paired, not proper), unless `--keepDiscordant` | common.c:430 |
| Strand of origin | If `XG:Z:` starts with C or G (Bismark style): XG plus R1/R2 plus reverse flag give **OT/OB/CTOT/CTOB**, so non-directional and PBAT are supported. **Otherwise flag-only, directional assumption**: R1 fwd or R2 rev = OT, R1 rev or R2 fwd = OB, SE fwd = OT, SE rev = OB. **bwa-meth's `YD` tag is ignored**. PBAT or non-directional BAMs without XG are therefore mis-stranded | common.c:84–116 |
| Informative base | OT/CTOT reads count only at ref C (read C = meth, T = unmeth). OB/CTOB only at ref G (G/A). Other read bases (e.g. N) are ignored. Reads on the "wrong" strand for the column only feed the variant filter | extract.c:422–441, common.c:118–134 |
| Context | From the **reference only**. `isCpG`: C followed by G, or G preceded by C. `isCHG`: C at pos+2 G. `isCHH`: any other C/G. Lower-case (soft-masked) bases are accepted. **N in the H slots counts as H (CNG = CHG, CNA = CHH). A C at the last base of a chromosome, or with pos+2 off the end, counts as CHH** (README "Methylation Context"). A reference `N` at the position itself is never a C or G, so it is skipped | common.c:49–82 |
| Indels / CIGAR | htslib pileup. `is_del` and `is_refskip` entries are skipped. Inserted and soft-clipped bases never appear in a pileup column. **The read's neighbouring bases and indels do not affect context** | extract.c:423–424 |
| Soft clips | Ignored for calling. Trimming bounds (`--OT` etc.) index into SEQ **including** soft-clipped bases | common.c:137–208 |
| Overlapping mates | **On by default** via a custom htslib pileup constructor (`overlaps.c:121–139`), per base in the overlap (`cust_tweak_overlap_quality`, 54–119). Skipped if the mates are on different BS strands. **Bases agree**: the higher-quality mate gets `+20%` (capped at 255) and the other is set to 0. On a **tie the second-pushed mate (larger POS) wins**. **Bases disagree**: the higher-quality mate's quality becomes `qa-qb` and the other 0. With equal quality, or the higher mate is N, both become 0. Low-quality base calls from a single mate can be removed by `-p` after this subtraction (e.g. 30 vs 28 gives 2 < 5, so **both are lost**) | overlaps.c |
| M-bias trimming | Default none. `--OT/--OB/--CTOT/--CTOB A,B,C,D` are 1-based inclusion positions for R1 (A..B) and R2 (C..D). `--nOT` etc. exclude N bases from each end. **Positions are indices into SEQ as stored in the BAM (reference left→right), NOT read 5'→3'**. For a reverse-strand alignment "position 1" is the read's 3' end. Evidence: `trimAlignment` indexes `qual[i]` directly, mbias bins by `plp->qpos`, and NEB's EM-seq pipeline trims the R2 5' end with `--nOT 0,0,0,5 --nOB 0,0,5,0` (`nebiolabs/EM-seq/modules/methyldackel_extract.nf:19`) | common.c:137–208, MBias.c:193–206 |
| Trimming mechanics | Trimmed bases get SEQ=N and qual=0 before overlap handling. So a trimmed R1 base can no longer beat R2 on quality, and R2's call is kept (the comment at common.c:446–457 describes this intent) | common.c:458–459 |
| Minimum depth | `-d/--minDepth 1` | 728 |
| SNPs | Off by default: `--minOppositeDepth 0`, `--maxVariantFrac 0.0`. When enabled, it counts non-G (or non-C) bases on the opposite-strand reads at C (or G) columns and excludes the site. With `--mergeContext` the whole CpG is excluded | extract.c:225–239, 443–459 |
| Conversion filter | `--minConversionEfficiency 0.0` (off). Counts CHG+CHH conversion per read | common.c:356–404 |
| Mappability | Off unless `-M bigWig` or `-B bbm` | |
| Region / threads | `-r`, `-l BED` (+ `--keepStrand`), `-@ 1`, `--chunkSize 1000000`. Chunks are extended so CpG/CHG are never split (`adjustBounds`) | |
| Pileup depth cap | `bam_mplp_set_maxcnt(INT_MAX)` (no cap) | extract.c:395 |

---

## 3. BISCUIT 1.10.x (`biscuit pileup`, `vcf2bed`, `mergecg`, `epiread`)

### 3.1 Formats

**`biscuit pileup [opts] ref.fa in.bam [in2.bam...]`** writes a **VCFv4.1** (`print_vcf_header`, pileup.c:923–989).
- One row per reference position with any signal (1-based `POS`). `REF` is the reference base, `ALT` is the top non-reference allele or `.`.
- `INFO` has `NS`, `CX=CG|CHG|CHH|CN` and `N5=<5-mer>` (strand-flipped for G), plus `SS/SC` in somatic mode and `AB`.
- `FORMAT` is `GT:GL1:GQ:DP:SP[:AC:AF1][:CV:BT]`.
  - **`CV`** is the strand-specific effective coverage (M+U after filters).
  - **`BT`** is the beta `M/(M+U)` written with **`%1.3f`**, so counts are **not stored exactly**.
  - `DP` is the raw depth. `SP` is allele support (e.g. `C12T3`).
- `-w` also writes `<out>_meth_average.tsv`.

**`biscuit vcf2bed [-t cg|ch|c|hcg|gch|snp] [-k mincov=1] [-e] [-c] in.vcf`** (vcf2bed.c):
- Default columns: `chrom start(0-based) end beta(%1.3f) cov`, repeated `beta cov` per sample. **No strand column.** Bottom-strand C rows sit at the G.
- `-e` inserts `refbase CX 2-base 5-base` before beta.
- `-c` writes `pct = round(beta*100)`, `M = round(cov*beta)`, `U = cov-M`. M is reconstructed from the 3-decimal beta, which is **lossy for cov ≳ 500**.
- `-t cg` keeps rows whose `CX` string equals `CG`, so CN-context CpGs are dropped.
- The default `-k` changed from 3 to 1 in 1.6.0.

**`biscuit mergecg [-N] [-c] [-k 0] ref.fa in.bed`** (mergecg.c):
- Columns: `chrom start end beta(%1.3f) cov C:beta:cov,G:beta:cov` (per sample).
- A CpG pair (C at i, G at i+1 in the reference) is merged to `start=C 0-based`, `end=start+2`. `M = rint(betaC*covC) + rint(betaG*covG)`.
- An unpaired C or G whose reference neighbour makes a CpG is still widened to 2 bp.
- `-c` gives `pct M U` followed by the `C:..,G:..` column.

**`biscuit epiread`** writes epiBED v2 (https://huishenlab.github.io/biscuit/epibed_format/):
- Columns: `chrom start(0-based, soft-clip-adjusted: read_start - n_softclip) end read_name read_number(1|2) bsstrand(+|-) CpG_RLE GpC_RLE(. unless -N) variant_RLE`.
- RLE codes:
  - `M/U`: methylated / unmethylated CpG
  - `O/S`: methylated / unmethylated GpC
  - `F`: filtered
  - `P`: soft clip
  - `D/d`: deletion
  - `x`: ignored
  - `A/C/G/T/R/Y`: SNP bases
  - `a/c/g/t/i`: insertion
- Old formats: `-O` and `-P` (pairwise).
- epiread CpG context uses only the **adjacent** reference base (epiread.c ~800–860), not the 5-mer, so it can disagree with pileup `CX`.

### 3.2 Implicit query and defaults (`bisc_utils.h:98–116` `meth_filter_init`, pileup.c)

| Aspect | Default | Source |
|---|---|---|
| Base quality | `-b 20` | bisc_utils.h:101, pileup.c:429 |
| MAPQ | `-m 40` (Bismark's MAPQs 42 and 40 pass; 24, 23 etc. are dropped) | :105, pileup.c:762 |
| Alignment score | `-a 40`, applied only if an `AS` tag exists (Bismark writes none; bwa and bwa-meth do) | :106, pileup.c:774–775 |
| NM / retention | `-n 999999` / `-t 999999` (off) | |
| Read length | `-l 10` | |
| End trimming | `-5 3 -3 3`: drops the first 3 and last 3 bases **in SEQ order including soft clips** (`qpos<=3 || rlen<qpos+3`). A TODO admits the 3' side should be the mapping end, and "5'" is really the left end of SEQ | :103–104, pileup.c:430–432 |
| FLAG | Secondary, duplicate, QC-fail and improper pairs (`paired && !proper`) are dropped. `-c`, `-u`, `-p` disable those. **Supplementary (0x800) is NOT filtered.** Unmapped reads are only removed indirectly (MAPQ) | pileup.c:764–769 |
| Strand (BS strand) | `YD:Z:f/r` (bwa-meth, BISCUIT align) first, then `ZS` (bsmap), then Bismark `XG` (CT=BSW, GA=BSC), then **inference** from the high-quality C>T vs G>A mismatch count (ties go to BSW). It only needs the conversion type, so directional, non-directional and PBAT all work | bisc_utils.c:163–238 |
| Informative base | BSW reads at ref C (C=retention, T=conversion). BSC reads at ref G (G/A) | pileup.c:841–858 |
| Context | Reference **5-mer** centred on the C (reverse complement at G). `fivenuc[3]=='G'` gives CG, `fivenuc[4]=='G'` gives CHG, else CHH. **An N anywhere in the 5-mer, including 2 bp upstream and chromosome-edge padding (first and last 2 bp), gives `CN` (NA)**. Read bases and indels are irrelevant | bisc_utils.c:33–72 |
| Reference N at the site | Row skipped | pileup.c:470–471 |
| Overlapping mates | **On** (`-d` disables). **Read 2 bases** at reference positions in `[max(R2start,R1start), min(R2end,R1end)]` are dropped. R1's end comes from the `MC` tag, or is **assumed equal to R2's reference length when MC is absent (always the case for Bismark BAMs)**. It is quality-blind and ignores whether R1's base survived filters | pileup.c:784–827 |
| SNP handling | Always on. Counts are redistributed (Y/R ambiguity, `-r` disables). The site is **not methylation-callable** (no CV/BT, so it disappears from vcf2bed) if, at ref C, BSC reads show T with `T/C ≥ 0.05` (G/A mirrored) | pileup.c:388–419, 517–534 |
| Min coverage | vcf2bed `-k 1`, mergecg `-k 0` | |
| Hard clips | `BAM_CHARD_CLIP` advances `qpos` (pileup.c:872–874) even though hard-clipped bases are not in SEQ, so indexing shifts after a leading `H`. `cinread.c:177–180` compensates, pileup does not. With BISCUIT/bwa hard-clipped supplementaries (default; `-Y` gives soft clips, `lib/aln/align.c:284`) and no 0x800 filter, this **is a bug, verified by execution** (BC-2, `validation/demos/biscuit-bugs`) | |
| Threads | `-@ 3`, window step 100000 | bisc_utils.h:59–60 |

---

## 4. methylpy / ALLCools (ALLC format)

### 4.1 ALLC format
Tab-separated, no header, bgzipped and tabix-indexed (`tabix -b 2 -e 2 -s 1`):
```
chrom  pos(1-based)  strand(+|-)  context(e.g. CGA; num_upstr_bases + C + num_downstr_bases)  mc  cov  methylated
```
- ALLCools always writes `methylated=1` (`_bam_to_allc.py:255–265`).
- methylpy writes `1`, then **replaces it with a binomial test result** (0/1): non-conversion rate from `--unmethylated-control`, BH-FDR `sig_cutoff 0.01`, `min_cov 2` (`call_mc_se.py:2064–2140`).
- **methylpy strips the `chr` prefix by default** (`remove_chr_prefix=True`, call_mc_se.py:1412, 1536–1539).
- Context is `seq[pos-up : pos+down+1]` (uppercased reference) for `+`, and the reverse complement for `−`. `N` is kept literally. Rows are skipped if the context is truncated (chromosome edge) or `cov==0`.

### 4.2 `allcools bam-to-allc` implicit query

`ALLCools/_bam_to_allc.py`:
- It runs `samtools mpileup -Q {min_base_quality} -q {min_mapq} -B -f ref bam` (149, 164) and parses column 5.
  - Ref C: `mc = count('.')`, `cov = count('.') + count('T')`.
  - Ref G: `mc = count(',')`, `cov = count(',') + count('a')` (246–290).
  - Indel `+N`/`-N` sequences are stripped first.
- Defaults: `min_mapq=10`, `min_base_quality=20`, `num_upstr_bases=0`, `num_downstr_bases=2` (lhqing 348–360; the DingWB fork is the same).
- **Strand = alignment orientation**, so only these inputs work unmodified: Bismark **SE directional** (forward reads = CT conversion) and methylpy BAMs.
  - `convert_bam_strandness` rewrites `is_forward` from `XG` (Bismark) or `YZ` (hisat-3n) **without** reverse-complementing SEQ. It is needed for PE.
  - **Default `False` in lhqing/ALLCools 1.1.1**, `True` in the DingWB fork.
  - With `False`, a Bismark PE OT read 2 (reverse-mapped, shows C/T at ref C) is counted as `,` at ref **G** positions, i.e. wrong-strand pseudo-methylation.
  - Non-directional CTOT/CTOB SE reads are also mis-assigned without conversion.
- samtools mpileup (1.x) defaults inherited (htslib.org mpileup manual):
  - `--excl-flags SECONDARY,QCFAIL,DUP` plus unmapped. **Supplementary is kept.**
  - **Anomalous pairs (paired, not proper) are skipped** unless `-A`.
  - **Overlap removal ON** (`-x` disables). htslib 1.24 `tweak_overlap_quality` (sam.c:5877–6010): when bases agree, **one mate chosen by a hash of QNAME** gets `min(qa+qb,200)` and the other 0. When they disagree, the higher-quality mate is kept ×0.8. A tie is chosen by hash ×0.8.
  - `-d 8000` max depth.
  - BAQ off (`-B`).
- **Suspected parsing bug (methylpy and ALLCools):** mpileup marks a read start as `^` followed by **MAPQ+33 as a character**, which is never stripped. MAPQ 13 prints as `.`, 11 as `,`, 51 as `T`, 64 as `a`. So every read starting at a C or G column with those MAPQs adds a spurious count. Bismark can emit MAPQ 11 (`bismark:calc_mapq` constants include 11). bwa-meth and bwa can emit 11, 13 and 51. `grep '\^'` finds nothing in either parser. **Verify with a crafted BAM.**
- ALLCools `extract-allc --strandness merge` (`_extract_allc.py:19–58`):
  - Merges `+` at i with `−` at i+1 and keeps the `+` row's context. A lone `−` becomes `+` at pos−1 with the `−` row's context.
  - **Suspected bug: the last buffered line of the file is never written** (no flush after the loop).
- methylpy's own pipeline (`call_mc_pe.py:1238–1291`) flips read-2 strand flags before mpileup. `call_methylated_sites` default `min_mapq=30`, `min_base_quality=1`.

---

## 5. CGmapTools (`cgmaptools convert bam2cgmap` = `CGmapFromBAM`)

**CGmap** (no header):
```
chrom  nuc(C|G)  pos(1-based)  context(CG|CHG|CHH|--)  dinuc(CA|CC|CG|CT|C?)  meth_level(%.2f)  mC_count  total_C_count(C+T)
```
**ATCGmap** (every covered position, including A/T):
```
chrom nuc pos context dinuc  W_A W_T W_C W_G W_N  C_A C_T C_C C_G C_N  meth_level(%.2f|na)
```
W = Watson (forward-mapped reads), C = Crick (reverse-mapped reads). A `.wig.gz` is also written.

Implicit query (`src/CGmapFromBAM.c:555–760`):
- Built on samtools 0.1.18 `sampileup`. Its default mask is `BAM_DEF_MASK` (unmapped, secondary, QC-fail, dup), with a max depth of 8000 (`include/samtools-0.1.18/bam_pileup.c:171`). **There is no MAPQ or base-quality filter.**
- Methylation at ref C is `fwd C / (fwd C + fwd T)`. At ref G it is `rev G / (rev G + rev A)`. **Strand comes from the mapping orientation only**, which the usage text says is correct for BS-Seeker2 and Bismark SE, but not for Bismark PE ≥0.8.3.
- Context comes from the reference and needs `pos-2 ≥ 0` and `pos+2 < len`, otherwise `--`. `N` in the H slots gives `--`.
- Pileup entries with `pl->indel != 0` (i.e. **the base immediately before an indel**) are **not counted**, and `is_del` is not checked.
- `-O/--rmOverlap` is off by default. When on, the first pileup entry per cleaned QNAME wins (the cleaner strips a 2-char `.1`/`#1`/`:1` suffix).
- `cgmaptools convert bismark2cgmap` converts a Bismark CX/CpG report into CGmap (ratio `%.2f`, rows with coverage > 0).
- Install: not on bioconda (only the `HCC/cgmaptools 0.1.2` conda channel). Otherwise build from source with `install.sh`, which vendors samtools 0.1.18 and zlib.

---

## 6. bwa-meth → MethylDackel handoff

- `bwameth.py c2t` converts R1 C→T and R2 G→A (`bwameth.py:190–211`) and keeps the original sequence in a `YS` comment. It aligns with `bwa mem -T 40 -B 2 -L 10 -CM [-U 100 -p]` to a doubled C→T/G→A reference with `f`/`r` contig prefixes.
- `-M` marks shorter split hits **secondary** (not supplementary).
- Post-processing (`handle_reads`, 461–510):
  - Strips the `f`/`r` prefix and adds **`YD:Z:f|r`**.
  - Restores the original SEQ.
  - **Chimera heuristic:** if the longest M is < 44% of the read length, it sets **0x200 QC-fail**, clears proper-pair and caps MAPQ at 1, for the whole pair. `--do-not-penalize-chimeras` disables this.
  - `--set-as-failed f|r` sets QC-fail for one strand (targeted libraries).
- There are no XR/XG/XM tags, and the design supports **directional libraries only** (README).
- The handoff (bwa-meth README; nf-core methylseq bwameth route; NEB EM-seq pipeline):
  1. `bwameth.py | samtools sort`
  2. **Picard MarkDuplicates** (sets 0x400)
  3. `samtools index`
  4. `MethylDackel extract ref.fa bam` (defaults drop 0x400, 0x200 chimeras and 0x100 secondary, MAPQ<10, improper pairs)
- MethylDackel **ignores `YD`** and infers OT/OB from R1/R2 plus the reverse flag. BISCUIT does use `YD`.
- NEB EM-seq pipeline: `MethylDackel extract --methylKit -q <thr> --nOT 0,0,0,5 --nOB 0,0,5,0 --CHH --CHG` plus `--cytosine_report` (`nebiolabs/EM-seq/modules/methyldackel_extract.nf:19`, `extract_cytosine_report.nf:18–22`).
- nf-core/methylseq MethylDackel defaults: no trimming, `--CHG --CHH` only with `--all_contexts`, `--mergeContext`/`--methylKit`/`--minDepth` optional (`conf/modules/methyldackel_extract.config`).

---

## 7. gemBS and Bis-SNP (brief)

- **gemBS** (bioconda `gembs 3.5.5_IHEC`; Rust rewrite at heathsc/gemBS-rs):
  - GEM3 aligner. `bs_call` does joint genotype and methylation calling into BCF. `mextr` extracts homozygous-reference C sites in CpG/CHG/CHH as txt, bed, bedGraph or bigWig.
  - The ENCODE4 WGBS pipeline uses it. Its bedMethyl is **bed9+2**: `chrom start end name score(min(cov,1000)) strand thickStart thickEnd itemRgb coverage percent_meth`, plus gemBS-specific genotype columns (https://www.encodeproject.org/data-standards/wgbs-encode4/).
  - Default filters (MAPQ, BQ, overlap) were not verified from source. Treat it as out of scope unless ENCODE bedMethyl export is wanted.
- **Bis-SNP** (bioconda `bis-snp 1.0.1`): GATK-1.x/Java era. It needs read groups plus GATK preprocessing and emits CpG and SNP VCFs. It is unmaintained. Not recommended for the validation matrix.
- Adjacent format worth knowing for export: **modkit bedMethyl** (ONT), with 18 columns: `chrom start end mod_code score strand start end color Nvalid_cov pct_mod Nmod Ncanonical Nother_mod Ndelete Nfail Ndiff Nnocall`. Verify against the modkit docs before implementing.

---

## 8. Downstream inputs (R/Bioconductor)

| Package (bioconda) | Expected input | Notes |
|---|---|---|
| **methylKit** (`bioconductor-methylkit 1.36.0`) | `methRead(pipeline="amp")` generic format: `chrBase chr base strand coverage freqC freqT`, 1-based, header, `strand` F/R, freq in 0–100 (MethylDackel `--methylKit` matches). Also `pipeline="bismarkCoverage"` (6-col .cov, no strand) and `"bismarkCytosineReport"`. `processBismarkAln` reads Bismark SAM/BAM (needs XM, sorted) with `minqual=20`, `mincov=10`, `nolap=FALSE` | Default `mincov=10` at read time |
| **bsseq** (`bioconductor-bsseq 1.46.0`) | `read.bismark()`: Bismark `.cov` (6 col, no strand) or **genome-wide cytosine report** (7 col, stranded; recommended). `strandCollapse=TRUE` (default) needs stranded loci. Coordinates must be 1-based (not `--zero_based`). bedGraph and call files are unsupported | Merged-CpG `.cov` from `--merge_CpG` can be read as a coverage file, but loses strand |
| **DSS** (`bioconductor-dss 2.58.0`) | Per sample, a data frame or text `chr pos N X` (N = total, X = methylated), 1-based. `makeBSseqData(list, names)` | Usually strand-merged CpGs (convention, not enforced) |
| **dmrseq** (`bioconductor-dmrseq 1.30.0`) | A `BSseq` object with integer M and Cov | Docs recommend removing loci with zero coverage in all samples of a condition |

---

## 9. KEY COMPARISON: implicit choices of Bismark vs MethylDackel vs BISCUIT (plus ALLCools and CGmapTools)

| Choice | Bismark extractor (on Bismark BAM) | MethylDackel extract | BISCUIT pileup+vcf2bed | ALLCools bam-to-allc | CGmapTools bam2cgmap |
|---|---|---|---|---|---|
| **Where the call is made** | Aligner (XM tag). Extractor only aggregates | At extraction (pileup) | At extraction (own CIGAR walk) | samtools mpileup text | samtools-0.1.18 pileup |
| **Context source** | Read-projected reference: **deleted bases skipped**, **I/S bases = unknown (U, dropped)**, +2 ref bases past the read end | Reference | Reference 5-mer | Reference 3-mer | Reference |
| **N in context** | `U/u`, not reported. (coverage2cytosine: CNG→CHG, CNN→CHH) | N counts as H (CNG→CHG, CNx→CHH) | Any N in the 5-mer (incl. upstream) → `CN`, dropped by `-t cg` | Literal `N` in the context string | `--`, row written with context `--` |
| **Chromosome edges** | Reads within 2 bp of an edge **discarded entirely**. c2c skips the last base | C at the end → CHH (`CNN`) | First/last 2 bp → CN | Truncated context → row skipped | ±2 bp → `--` |
| **Deletion next to C** | Context jumps over the deletion (CpG may become CHG/CHH) | No effect | No effect | No effect | Base before an indel **not counted** |
| **Insertion / soft clip next to C** | U/u → dropped | No effect | No effect (but soft clips count toward the 3 bp end trim) | No effect | Base before an insertion not counted |
| **Strand of origin** | XR+XG → OT/CTOT/CTOB/OB. Library type fixed at alignment | XG (if present) + flags → 4 strands. Without XG: flags only, **directional assumption** | Conversion type only: YD > ZS > XG > mismatch inference. Library-agnostic | **Alignment orientation** (optional XG/YZ rewrite; default off in lhqing 1.1.1) | Alignment orientation |
| **MAPQ** | none | ≥10 | ≥40 (+ AS ≥40 if AS present) | ≥10 | none |
| **Base quality** | none | ≥5 | ≥20 | ≥20 | none |
| **Duplicates (0x400)** | not filtered (run `deduplicate_bismark`, which removes) | dropped (`--keepDupes`) | dropped (`-u`) | dropped | dropped |
| **Secondary (0x100)** | not filtered (Bismark emits none) | dropped | dropped | dropped | dropped |
| **Supplementary (0x800)** | not filtered (Bismark emits none) | dropped | **kept** (+ hard-clip indexing issue) | **kept** | kept |
| **QC fail (0x200)** | not filtered | dropped | dropped | dropped | dropped |
| **Improper pairs** | not filtered (Bismark uses `--no-mixed --no-discordant`) | dropped (`--keepDiscordant`); singletons dropped | dropped (`-p`) | dropped (mpileup default, `-A` keeps) | not filtered |
| **NH>1** | n/a | dropped | n/a | n/a | n/a |
| **Overlapping mates** | ON: **R1 wins wholesale**, R2 truncated from its 5' end at R1's end, quality-blind | ON: **per-base quality arbitration** (winner +20%, tie → later mate; mismatch → subtract; can lose both) | ON: **R1 wins**, R2 bases in the overlap interval dropped (interval from MC or assumed length), quality-blind | ON (htslib): **hash(QNAME) picks the mate**, qualities summed (≤200) | OFF (`-O` → first pileup entry wins) |
| **End trimming default** | none (`--ignore*`, 5'→3' read coords) | none (`--OT` etc., **SEQ-left→right coords**, left bound excludes position A) | **3 bp both SEQ ends** (incl. soft clip) | none | none |
| **SNP handling** | none | off (opt-in opposite-strand filter) | **on**: site dropped if opposite-strand T/C ≥ 5% | none | none |
| **Min depth** | bedGraph/cov `--cutoff 1`; c2c threshold 0 (all Cs) | 1 (report: 0) | vcf2bed `-k 1` | cov>0 | cov>0 |
| **Unit of count** | Per read (R1/R2 separately, overlap trimmed) | Per alignment after overlap tweak | Per alignment after R2 trim | Per alignment after overlap tweak | Per alignment |
| **Default contexts** | Writes all three contexts to call files; bedGraph/cov/report CpG only (`--CX`) | CpG only (`--CHG --CHH`) | All (VCF); vcf2bed `-t cg` default | All (3-mer) | All |
| **Strand merge** | c2c `--merge_CpG` → `.merged_CpG_evidence.cov` (1-based C,G; %.6f) | `--mergeContext` → bedGraph start=C, end=C+2 (CHG +3) | `mergecg` → start=C 0-based, end+2, `C:..,G:..` column | `extract-allc --strandness merge` | `CGmapCombineStrands.pl` (not audited) |
| **Percent formatting** | Perl default float (bedGraph/cov); `%.6f` (merged) | `int()` truncation; `%f` fraction; `%6.2f` methylKit | beta `%1.3f`; `-c` rounded int | counts only | ratio `%.2f` |
| **Coordinates** | calls & cov & report 1-based; bedGraph 0-based half-open; `--zero_based` optional | bedGraph 0-based half-open; methylKit & cytosine_report 1-based | VCF 1-based; BED 0-based half-open | 1-based | 1-based |
| **Input order** | Unsorted/name-grouped (PE adjacent) | Coordinate-sorted + indexed | Coordinate-sorted + indexed | Coordinate-sorted (indexes itself) | Coordinate-sorted |

### Harmonisation recipes for a like-for-like comparison on the same Bismark BAM

- To make MethylDackel mimic Bismark:
  - Use `-q 0 -p 1`, and `--keepDupes` only if the BAM was deduplicated already. (Bismark BAMs have no dup flags; keeping the default is harmless there.)
  - `--keepDiscordant --keepSingleton` are no-ops for Bismark PE.
  - Translate `--ignore_r2 2` using SEQ orientation: OT R2 is reverse, so `--nOT 0,0,0,2`; OB R2 is forward, so `--nOB 0,0,2,0`. Also `--nCTOT 0,0,2,0 --nCTOB 0,0,0,2` for non-directional data (derive from getStrand: CTOT R2 forward, CTOB R2 reverse).
  - Residual differences to expect:
    - Overlap arbitration (R1-wins vs quality-based)
    - Deletion- and insertion-adjacent context
    - Reads near chromosome edges
    - N-context handling
    - Bismark ignoring U/u
- To make BISCUIT mimic: `-m 0 -b 0 -5 0 -3 0 -u -p -c`. BISCUIT's SNP veto and 5-mer N rule cannot be switched off (`-r` only disables redistribution).
- For ALLCools on Bismark PE: `--convert_bam_strandness` is mandatory. Watch the `^`-MAPQ artefact.
- In **alnbase terms**, each tool is a query plus a filter plus an aggregation. Suggested explicit parameters to expose so every tool can be emulated:
  1. `context_from = {reference, read_projected_reference(skip D, I/S→unknown)}`
  2. N/edge policy `= {unknown-drop, as-H, whole-kmer-NA}`
  3. Strand source `= {XG+flags, YD/ZS/XG/infer, alignment orientation}`
  4. Overlap policy `= {R1-wins-truncate-R2 (Bismark), R2-interval-drop (BISCUIT), quality-arbitration (MethylDackel), htslib hash-sum, none}`
  5. End trim coordinates `= {read 5'→3', SEQ left→right incl. soft clips}`
  6. FLAG/MAPQ/BQ filters
  7. SNP veto

---

## 10. Known correctness issues and contestable choices

**Bismark**
1. Context is built from the read-projected genome, so **deletions change the called context** and a C before an insertion or soft clip becomes U (dropped). This is contestable. alnbase with pure reference context will disagree at indel-adjacent Cs (`bismark:4273–4470`, `4800–5020`).
2. The aligner uses U/u for N contexts, while coverage2cytosine classifies CNG as CHG and CNN as CHH (`coverage2cytosine:365–378`). The two Bismark stages are internally inconsistent for N contexts.
3. The PE `--no_overlap` rule is R1-wins and quality-blind. It truncates all of R2 from the first overlapping position, and dovetailed R2 tails are lost too. R2 M-bias is computed only on non-overlapping bases. Discussion: https://github.com/FelixKrueger/Bismark/issues/233 ("ignore overlap or only R1, nothing in between").
4. The extractor applies no BQ, MAPQ or FLAG filtering and treats the first record of each pair as R1 regardless of FLAG. Any re-sorted or merged BAM breaks it: it dies if QNAMEs mismatch, but silently mislabels R1/R2 if the mates are swapped.
5. M-bias positions shift when `--ignore` is set (the index is taken after trimming).
6. Reads within 2 bp of a contig edge are discarded at alignment. This matters for small spike-in contigs (lambda, pUC19).
7. Cosmetic: `$second_read_conversion = s/\r//;` (PE, ~line 1905) assigns a match count instead of stripping. The value is unused.
8. The v0.25.0 CHANGELOG notes a fixed bismark2bedGraph CHH bug (#647, https://github.com/FelixKrueger/Bismark/issues/647).
9. bedGraph percentages are unformatted floats, and cov files are 1-based closed while bedGraphs are 0-based half-open. Both are frequent sources of off-by-one confusion when mixing Bismark files.

**MethylDackel**
1. **`--OT/--OB` left bound off by one:**
   - The help and README say `A` is the first included 1-based position ("--OT 5,0,0,0 would include all but the first 4 bases").
   - But `trimAlignment` masks `i < lb`, i.e. positions 1..A, so position A is **excluded** (`common.c:151–160`).
   - mbias's suggested `lthresh = i+2` (`svg.c:281–282`) therefore over-trims by one.
   - The right bound is inclusive as documented.
2. Trimming coordinates are in **BAM SEQ orientation**, not read 5'→3'. This is undocumented and a common source of confusion: https://github.com/dpryan79/MethylDackel/issues/163 and /issues/102 are open questions. NEB's pipeline relies on SEQ orientation (`--nOT 0,0,0,5`).
3. The overlap quality tweak can **drop both mates** on a base mismatch, or leave one below `-p`. The README says `-p` "defaults to 10" but the code default is 5.
4. Without XG it assumes a directional library. PBAT or non-directional bwa-meth-style BAMs get the wrong C/G column, and the `YD` tag is ignored.
5. `perRead`:
   - Low-quality skip: after skipping a base it does not `continue`, so the next base is evaluated without a quality or CIGAR-boundary check (`perRead.c:58–63`).
   - Defaults include duplicates, secondary and supplementary reads (`ignoreFlags=0`).
   - No overlap handling.
6. Percentages are truncated with `(int)`, and the methylKit output uses space-padded `%6.2f`.
7. Open issues:
   - #168: quality reads as zero under some compiler optimisation, losing counts. https://github.com/dpryan79/MethylDackel/issues/168
   - #170: ~50% lower depth than Bismark on the same data. Expected from the MAPQ, dup, improper-pair and overlap defaults. https://github.com/dpryan79/MethylDackel/issues/170
   - #157: different calls with different `-l` region sets.
8. The C at the last base of a chromosome is labelled CHH, and `--cytosine_report` includes it (Bismark does not).

**BISCUIT**
1. Supplementary alignments are not filtered, and hard-clip handling in `pileup.c:872–874` shifts `qpos`. Suspected wrong-base reads for hard-clipped records.
2. The `-5/-3` end filters use SEQ-left/right positions including soft clips, not 5'/3' of the aligned read (acknowledged TODO at pileup.c:430).
3. The overlap interval assumes equal mate length when there is no `MC` tag, which is the case for Bismark BAMs. That can drop too much or too little of R2.
4. An N anywhere in the 5-mer (including 2 bp **upstream**) makes the site `CN` and excludes it from `-t cg`, even for a real CpG.
5. Methylation counts in VCF/BED are reconstructed from 3-decimal betas (`vcf2bed -c`, `mergecg`). The rounding errors are lossy for exact count validation. Compare `CV` and `BT` with a tolerance, or use verbose `RN/CN` (`-v`).
6. The SNP veto silently removes CpGs at a ≥5% opposite-strand alt ratio. That is sensible, but differs from all the others.

**ALLCools / methylpy / CGmapTools**
1. `^`+MAPQ characters are not stripped from mpileup, giving spurious `.`, `,`, `T`, `a` counts for reads starting at a C or G with MAPQ 13, 11, 51 or 64 (verified by execution for ALLCools: `validation/demos/allcools-bugs`; methylpy shares the parser but was not run).
2. `convert_bam_strandness` default False in lhqing 1.1.1 means Bismark PE input silently gives wrong-strand counts.
3. `extract-allc --strandness merge` never flushes the final record (verified by execution: `validation/demos/allcools-bugs`).
4. methylpy strips `chr` prefixes by default.
5. CGmapTools has no quality filters, skips pileup entries adjacent to indels (`pl->indel`), and its orientation-based strand is wrong for Bismark PE (acknowledged in the tool's own usage text).

---

## 11. Installability and design of the comparison

| Tool | Bioconda (latest seen 2026-09) | Runtime deps | Genome prep | Runs on a foreign BAM? |
|---|---|---|---|---|
| Bismark | `bismark 3.1.0` (Rust; `bismark=0.25.1` = Perl legacy). Also `cargo install bismark`, `ghcr.io/felixkrueger/bismark:3.1.0` | bowtie2 ≥2.5.5 (or hisat2, minimap2); Perl version also needs samtools | `bismark_genome_preparation --bowtie2 <folder>` (C→T and G→A indexes); the extractor only needs the FASTA folder for `--cytosine_report` | **Only with XM/XR/XG tags** plus PE adjacency. alnbase could emit these tags to reuse Bismark's aggregation |
| MethylDackel | `methyldackel 0.6.1` (linux-64/aarch64, osx) | htslib 1.21, libBigWig, zlib | faidx only (`samtools faidx`) | **Yes** (any aligner; coordinate-sorted and indexed). 4-strand only if XG is present |
| BISCUIT | `biscuit 1.10.2.20260818` (+ `dupsifter 1.4.0` for dup marking) | libcurl, libdeflate, ncurses, zlib | `biscuit index` only for `biscuit align`; pileup needs FASTA + .fai | **Yes** (YD/ZS/XG or inference) |
| bwa-meth | `bwameth 0.2.10` (noarch) | bwa or bwa-mem2, samtools, toolshed | `bwameth.py index[-mem2] ref.fa` | n/a (aligner) |
| methylpy | `methylpy 1.4.7` (python 3.8 build) | bowtie/bowtie2/minimap2, samtools, picard, cutadapt, pysam, scipy | `methylpy build-reference` | `call-methylation-state` expects methylpy-style orientation |
| ALLCools | Not on conda. `pip install allcools` (1.1.1) or the DingWB fork from GitHub | samtools, tabix, pysam, pandas | faidx | Bismark SE directly. Bismark PE with `--convert_bam_strandness`. hisat-3n (YZ) |
| CGmapTools | Not on bioconda (`HCC/cgmaptools 0.1.2` channel) or `install.sh` source build | vendored samtools-0.1.18, zlib | faidx | Orientation-based. OK for BS-Seeker2 and Bismark SE only |
| gemBS | `gembs 3.5.5_IHEC` | GEM3, bcftools etc. | `gemBS index` | Needs a gemBS BAM |
| Bis-SNP | `bis-snp 1.0.1` | Java / GATK1 | none | Legacy |
| R pkgs | `bioconductor-bsseq 1.46.0`, `-methylkit 1.36.0`, `-dss 2.58.0`, `-dmrseq 1.30.0` | R | | |
| samtools / htslib | `samtools 1.24` | | | |

Docker is available on this host, and every bioconda package has a biocontainers image (`quay.io/biocontainers/<pkg>:<ver>--<build>`).

**Recommended validation design**
1. **Primary matrix: one Bismark BAM, many extractors.**
   - Align with `bismark --genome <prep>` (directional; also run `--non_directional` on the E. coli Sherman data).
   - Unsorted or name-sorted BAM → `bismark_methylation_extractor` (Perl 0.25.1 and Rust 3.1.0; they should match byte for byte).
   - `samtools sort` + `index` of the same BAM → `MethylDackel extract` (reads XG), `biscuit pileup` (reads XG), `allcools bam-to-allc --convert_bam_strandness`, and **alnbase**.
   - This isolates extraction differences from alignment differences.
2. **Native bwa-meth path:** `bwameth.py` → `samtools sort` → Picard MarkDuplicates → `MethylDackel extract`/`mbias`. Also run BISCUIT on it (YD) and alnbase (YD / flags).
3. **alnbase → XM round trip:** have alnbase write XM/XR/XG on the Bismark BAM (overwriting) and run `bismark_methylation_extractor`. The outputs should be byte-identical to the original if alnbase reproduces Bismark's read-projected-context rule. Diffs pinpoint indel/edge cases.
4. **Synthetic edge-case BAMs** (hand-written SAM), one per behaviour in §9/§10:
   - C followed by D, C followed by I/S
   - C at the chromosome start/end and within 2 bp of it
   - N in context (±1, ±2, upstream 2)
   - Overlapping mates with equal and unequal quality, agreeing and disagreeing
   - Dovetailed mates
   - Hard-clipped supplementary
   - Reads starting at a C with MAPQ 11/13/51
   - PBAT and non-directional strands with and without XG
   - `--OT A` boundary

---

## 12. Small public test datasets

| Dataset | Size | Why | URL / path |
|---|---|---|---|
| nf-core/test-datasets `methylseq` branch: lambda reference `reference/genome.fa` (one contig `chr`, 48,502 bp) + prebuilt `Bowtie2_Index.tar.gz`, `Bwameth_Index.tar.gz`, `Bwameth_mem2_Index.tar.gz`, `Hisat2_Index.tar.gz` | ~50 KB FASTA, indexes ≤375 KB | Tiny genome with prebuilt indexes, so alignment takes seconds | `https://raw.githubusercontent.com/nf-core/test-datasets/methylseq/reference/genome.fa` (etc.) |
| `testdata/SRR389222_sub{1,2,3}.fastq.gz` | 2.6–3.9 MB each, SE | The nf-core methylseq `test` profile samples (`methylseq/assets/samplesheet.csv`) | `https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub1.fastq.gz` |
| `testdata/Ecoli_10K_methylated_R{1,2}.fastq.gz` | ~470 KB each, 10k PE 100 bp | **Sherman-simulated ground truth**: 80% CpG methylation, 10% non-CG, **non-directional** (~25% per strand), so it exercises CTOT/CTOB handling. Needs the E. coli genome (below) | same base URL `/testdata/` (README on that branch) |
| Spike-in controls: `reference/Enterobacteria_phage_lambda/GCF_000840245.1_ViralProj14204_genomic.fa` (unmethylated), `reference/pUC19/pUC19.fa` (CpG-methylated) + indexes | 49 KB / 3 KB | Conversion-efficiency truth: lambda ≈0% and pUC19 CpG ≈96–98% (EM-seq, Vaisvila 2021). Note Bismark discards reads within 2 bp of contig ends | same branch `reference/` |
| Bismark repo `test_files/`: `test_R1.fastq.gz`, `test_R2.fastq.gz` (~220 KB), `NC_010473.fa.gz` (E. coli DH10B, 1.4 MB gz), `lambda_NC_001416.fa.gz`, `pUC19.fa.gz` | small | The official Bismark CI fixture; its Rust byte-identity gate uses these | `https://github.com/FelixKrueger/Bismark/tree/master/test_files` |
| Bismark `rust/bismark/test_files/tiny_pe_bismark.bam` | 22 KB | A ready-made Bismark v0.25.1 PE BAM (E. coli, bowtie2 2.4.5) for immediate extractor tests | same repo |
| MethylDackel `tests/`: `cg100.fa`, `ct100.fa`, `chgchh.fa`, `cg_aln.bam`, `cg_with_variants.bam`, `chgchh_aln.bam` (+ .bai, `cg_R1/R2.fq`) | <1 KB each | Unit fixtures for CpG, CHG/CHH and variant filtering | `https://github.com/dpryan79/MethylDackel/tree/master/tests` |
| NEB EM-seq pipeline `tests/fixtures/fastq/emseq-test{1,2}.ds.{1,2}.fastq.gz` + `methylation_controls.fa` (lambda, pUC19, T4) | 18–290 KB | Real **EM-seq** PE reads with spike-ins, run through bwa-meth + MethylDackel | `https://github.com/nebiolabs/EM-seq/tree/master/tests/fixtures/fastq` |
| BISCUIT test data: hg38 chr22 + 1000 Sherman PE 150 bp reads (`test/setup_tests.py`, `make_data.sh`) | chr22 ~12 MB gz | Human-sequence context (N runs, repeats) | `https://github.com/huishenlab/biscuit/tree/master/test` |
| Sherman simulator (FelixKrueger/Sherman) | tool | Generate controlled ground truth (`--CG`, `--CH` conversion %, `--non_dir`, `--pbat`, PE, error rate) at any size | `https://github.com/FelixKrueger/Sherman` |
| Larger real data: nf-core `test_full` human WGBS `SRR7961102…SRR7961150` (s3://ngi-igenomes/test-data/methylseq/); Vaisvila et al. 2021 EM-seq/WGBS NA12878 with lambda/pUC19 spike-ins (accessions in the Genome Research Data Access section, doi:10.1101/gr.266551.120) | GBs | Scale and realism | `test-datasets/methylseq/samplesheet/samplesheet_full.csv` |

---

### Sources (web)
- https://www.htslib.org/doc/samtools-mpileup.html (mpileup defaults, `^`/MAPQ encoding)
- https://huishenlab.github.io/biscuit/epibed_format/ ; https://huishenlab.github.io/biscuit/docs/methylextraction.html
- https://rdrr.io/bioc/bsseq/man/read.bismark.html ; https://rdrr.io/bioc/methylKit/man/methRead-methods.html ; https://rdrr.io/bioc/methylKit/man/processBismarkAln-methods.html ; https://rdrr.io/bioc/DSS/f/inst/doc/DSS.Rmd
- https://github.com/dpryan79/MethylDackel/issues/163, /102, /168, /170, /157
- https://github.com/FelixKrueger/Bismark/issues/233, /647
- https://www.encodeproject.org/data-standards/wgbs-encode4/
- https://github.com/nebiolabs/EM-seq ; https://github.com/heathsc/gemBS-rs
- Bioconda versions via `https://api.anaconda.org/package/bioconda/<pkg>` (queried 2026-09-16)
