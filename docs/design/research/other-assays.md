# alnbase beyond bisulfite: survey of read-vs-reference assays, comparators, and query-model gaps

Date: 2026-09-16. Scope: assays whose signal is a base-level comparison between an aligned read and the reference (substitution, deletion, insertion, truncation/read-end event, in sequence context, on a defined strand). For each assay: (a) signal, (b) comparators, (c) output formats, (d) implicit defaults, (e) alnbase expressibility and minimal query-syntax extension, (f) small validation data, (g) comparator correctness issues. Part 1 is the synthesis; Part 2 holds the per-family detail with citations (source file/line references are to shallow clones in `clones/`).

Method: tool source code was read directly for defaults (clones in `clones/`, fetched single files in `clones/raw/`); web docs, GitHub issues and benchmark papers are cited inline. **Unless marked "verified by execution", suspected comparator bugs were found by reading code and still need reproduction** before being claimed in a paper. Two were executed: the wavClusteR MD-tag parser on deletions (R) and SAILOR's soft-clip parser (Python).

---

## Part 0. alnbase semantics that decide expressibility (verified in `src/`)

1. **Walk orientation = original strand under `--library directional`** (`src/alignment.rs:117`, `src/tags.rs`): R1 forward (OT) and R2 reverse (CTOT) are walked top-strand; R1 reverse (OB) and R2 forward (CTOB) are walked reverse-complemented. Single-end reads are read 1. Consequences:
   - Forward-stranded RNA libraries (R1 sense: QuantSeq/SLAMseq kit, iCLIP, PAR-CLIP, eTAM R2-as-SE, 10x scDART R2): walk is RNA-sense; write patterns as RNA.
   - dUTP / reverse-stranded libraries (R1 antisense: TruSeq Stranded, most TimeLapse/TUC, DART bulk, STAMP, eCLIP R1): walk is **antisense** to the RNA; RNA T>C is `read G / refr A`, and multi-column motifs must be written reverse-complemented (column order flips).
   - Unstranded RNA and plain DNA: a fragment's two orientations are split across walk orientations; one pattern cannot collect "the same RNA strand".
   - For plain DNA-seq this is a feature: R1 G>T and R2 C>A (the 8-oxoG signature) are one pattern `read T / refr G`.
2. **Pads** (`_`) mark columns past either end of the aligned read (and off-contig), with reference sequence readable there. `_N` can only fire at the walk *start*, `N_` at the walk *end* — but walk start is the 5' end of R1 and the **3' end of R2**, and with `--soft-clips skip` (default) the pad abuts the first *aligned* base, not the sequenced 5' end. There is no code distinguishing "past the 5' end as sequenced" from "past the 3' end" or from "off contig".
3. **Soft clips and insertions emit identically** when enabled (`src/alignment.rs:277-320`): read base against a reference gap, `refr_pos` = nearest 5' coordinate. A clipped base therefore never carries the reference base it "would" cover, and a pattern cannot tell a clip from an insertion.
4. **Deletions** give one column per deleted reference base (read gap); placement is whatever the aligner chose (no normalization).
5. **No filters** (dup/secondary/supplementary/QC-fail/MAPQ/BQ); parquet rows carry `qual`, `mapq`, `flag`, aux tags. Per-read aggregation (conversions per read, mutations per read, read-level drops) and site aggregation are DuckDB jobs and are *not* counted as missing features below.
6. `alnbase overlap` merges mate overlaps with a quality policy; comparators each use a different overlap rule (see extension X7).

---

## Part 1. Synthesis

### 01.1 Overall ranking (value for a versatility + correctness paper)

Scoring considered popularity, comparator cleanliness (install, deterministic, diffable intermediate), ground truth / simulation ease, likelihood of a *meaningful* discrepancy, and expressibility today. "Signal class" is listed so a panel can cover distinct query shapes.

| # | Assay | Signal (class) | Comparator(s) (version; install) | Primary output format | Expressible today? |
|---|---|---|---|---|---|
| 1 | **SLAM-seq / TimeLapse / TUC** (metabolic labeling) | T>C in RNA (4sU; s6G G>A), per-read counts (substitution) | SLAMdunk 0.4.3 (bioconda); fastq2EZbakR/bam2bakR → bakR 1.0.1/EZbakR (Snakemake, CRAN); GRAND-SLAM 2.0.7 (Maven) | SLAMdunk `tcount.tsv` (Chromosome, Start, End, Name, Length, Strand, ConversionRate, ReadsCPM, Tcontent, CoverageOnTs, ConversionsOnTs, ReadCount, TcReadCount, multimapCount, ConversionRateLower/Upper); bam2bakR/fastq2EZbakR per-read `*_counts.csv` (qname, nA, nC, nT, nG, rname, FR, sj, TA, CA, GA, NA, AT, CT, GT, NT, AC, TC, GC, NC, AG, TG, CG, NG, AN, TN, CN, GN, NN); GRAND-SLAM `.tsv.gz` (per gene), `.binom.tsv` (n, d, Condition, Type, count), `.mismatches.tsv`, `.snpdata` | **Y** (`read C/refr T` fwd-stranded; `read G/refr A` dUTP) |
| 2 | **A-to-I RNA editing** | A>G in RNA (substitution) | REDItools2 / REDItools3 3.7 (bioconda); JACUSA2 2.1.17 (bioconda); Alu Editing Index (Docker) | REDItools TSV: Region, Position (1-based), Reference, Strand (1/0/2), Coverage-q30, MeanQ, BaseCount[A,C,G,T], AllSubs, Frequency (top alt/(alt+ref)), gCoverage-q30, gMeanQ, gBaseCount[A,C,G,T], gAllSubs, gFrequency; JACUSA2 BED6+ (0-based): contig, start, end, name, score, strand, bases11 (A,C,G,T), info, filter, ref; AEI CSV | **Y** (homopolymer/proximity filters verbose → X4; unstranded AEI → X2) |
| 3 | **SHAPE-MaP / DMS-MaPseq** | all-base mismatches + deletions + insertions; DMS A/C (substitution + indel) | ShapeMapper 2.3.2 (tarball; not bioconda); SEISMIC-RNA 0.26.0 (bioconda); RNA Framework rf-count 2.9.7 (bioconda) | ShapeMapper `<name>_<RNA>_profile.txt` (Nucleotide 1-based, Sequence, per sample `<S>_mutations`, `_read_depth`, `_effective_depth`, `_rate`, … , Reactivity_profile, Std_err, HQ_profile, HQ_stderr, Norm_profile, Norm_stderr) + per-read `parsed.mut`; SEISMIC `{step}-position-table.csv` (Position 1-based, Base, Covered, Informative, Matched, Mutated, Subbed, Subbed-A/C/G/T, Deleted, Inserted) + read table; rf-count `.rc` → `rf-rctools view` (0-based) | **Needs X3** (indel-placement normalization) for parity; substitutions/stops Y |
| 4 | **iCLIP / eCLIP truncation** | cDNA start −1 (read-end event) | htseq-clip 2.19 (bioconda); PureCLIP 1.3.1 (bioconda); iCount | htseq-clip `extract` BED6 (one line per read, 0-based: chrom, start, end, name=`qname|qlen`, score=YB or 1, strand); PureCLIP BED6 (chr, start, end, state, score=log posterior ratio, strand) | **Y with caveats** (`_N` / `N_` by mate); **clean with X1 + X2** |
| 5 | **Ancient DNA damage** | C>T at 5' end, G>A at 3' end (ds) / C>T both ends (ss) (position-dependent substitution) | mapDamage2 2.2.3 (bioconda); DamageProfiler 1.1 (bioconda); PyDamage 1.0 (bioconda) | `misincorporation.txt`: Chr, End (5p/3p), Std (+/-), Pos (1..70 from that end), A, C, G, T, Total, G>A, C>T, A>G, T>C, A>C, A>T, C>G, C>A, T>G, T>A, G>C, G>T, A>-, T>-, C>-, G>-, ->A, ->T, ->C, ->G, S (DamageProfiler writes the same header); `5pCtoT_freq.txt`, `3pGtoA_freq.txt` | **Y** (DuckDB on `off_5p/off_3p`); exact mapDamage clip convention → X6 |
| 6 | **TAPS / Illumina 5-Base** | C>T at *modified* CpG (inverted logic) | rastair 2.2.0 (bioconda); asTair 3.3.3 (pip, deprecated); MethylDackel with flipped calls | rastair `call` BED: chr, start, end, name, beta_est, strand, unmod, mod, no_snp, snp, coverage, genotype, gt_p_score, gt_conf_score, cpg (+ VCF, per-read BED); asTair `.mods`: CHROM, START, END, MOD_LEVEL, MOD, UNMOD, REF, ALT, SPECIFIC_CONTEXT, CONTEXT, SNV, TOTAL_DEPTH | **Y** (directional); tagmented TAPS → X2 (strand inference) |
| 7 | **GLORI / eTAM-seq** (m6A) | unconverted A (A retained vs A>G) (substitution, per-read filters) | GLORI-tools (scripts); hisat-3n-table (HISAT-3N) | GLORI site TSV (Chr, Sites, Strand, Gene, CR, AGcov, Acov, Genecov, Ratio, Pvalue, P_adjust); hisat-3n-table (ref, pos 1-based, strand, convertedBaseQualities, convertedBaseCount, unconvertedBaseQualities, unconvertedBaseCount) | **Y** (hisat-3n parity → X2 `tag:YZ`) |
| 8 | **NOMe-seq** | GCH conversion (accessibility) vs HCG/WCG (endogenous) (substitution in context) | Bismark `coverage2cytosine --nome-seq` 0.24.x (bioconda); BISCUIT 1.x `pileup -N` (bioconda) | Bismark `*.NOMe.GpC.cov`/`*.NOMe.CpG.cov` (chr, start, end, pct, meth, unmeth; 1-based) and `*.NOMe.{GpC,CpG}_report.txt` (chr, pos, strand, count_meth, count_unmeth, C-context, trinucleotide); BISCUIT VCF (INFO CX=HCG/GCH…) → `vcf2bed -t gch|hcg` | **Y** (directional); scNMT/scNOMe → X2 (`pbat`/non-directional) |
| 9 | **BID-seq / PRAISE** (Ψ) | deletion at reference T (U-tract ambiguity) (deletion) | BID-pipe (scripts + binaries); PRAISE pipeline | `call_sites/{genes,genome}.tsv.gz`: chr, pos (1-based), strand, `<sample>_depth`, `<sample>_gap` (depth includes gaps) | **Y** for aligner-placed deletions; **X3** for placement-invariant counting |
| 10 | **DART-seq** (m6A) | C>U next to A in RAC (substitution + motif) | Bullseye (bioconda) | Find_edit_sites.pl BED (0-based): chr, start, end, gene (`Gene|region|C2U|mut=n|score|rep=N`), score (DART edit ratio), strand, control_ratio, control_total, dart_ratio, dart_total, conversion; then RACfilter | **Y** (dUTP motif written RC → X2 convenience) |
| 11 | **8-oxoG / FFPE / error-by-cycle** | G>T (R1) / C>A (R2), C>T; by read position (substitution + position) | Picard 3.5.0 CollectOxoGMetrics / CollectSequencingArtifactMetrics (bioconda); fgbio 4.1.1 ErrorRateByReadPosition | Picard OxoG metrics: SAMPLE_ALIAS, LIBRARY, CONTEXT, TOTAL_SITES, TOTAL_BASES, REF_NONOXO_BASES, REF_OXO_BASES, REF_TOTAL_BASES, ALT_NONOXO_BASES, ALT_OXO_BASES, OXIDATION_ERROR_RATE, OXIDATION_Q, …; fgbio `.error_rate_by_read_position.txt` (read_number, position, …) | **Y** (single pattern captures both mates) |
| 12 | **SNV pileups (baseline)** | any substitution | bcftools/samtools 1.24 (bioconda); pysamstats 1.1.2 (bioconda); GATK CollectAllelicCounts | VCF (FORMAT PL,AD by default; INFO DP4/I16), pysamstats `variation`/`variation_strand` TSV (per-base counts, `_fwd`/`_rev`), GATK CollectAllelicCounts TSV (CONTIG, POSITION 1-based, REF_COUNT, ALT_COUNT, …) | **Y** |
| 13 | **CRISPR base editing** | C>T (CBE) / A>G (ABE) in amplicon window (substitution + window) | CRISPResso2 2.3.4 (bioconda) | `Nucleotide_frequency_table.txt` / `Nucleotide_percentage_table.txt`, `Quantification_window_nucleotide_frequency_table.txt`, `Substitution_frequency_table.txt`, `Alleles_frequency_table.zip` (base-by-position matrices; header = amplicon bases) | **Y** (base editing); prime-editing classification N (out of scope) |
| 14 | **PAR-CLIP** | T>C (4SU) / G>A (6SG) (substitution) | PARalyzer 1.5 (bioconda, closed source); wavClusteR 2.46 (Bioconductor) | PARalyzer clusters/groups CSV (1-based closed) + distributions file; wavClusteR GRanges → BED (writes 1-based starts: off-by-one) | **Y** |
| 15 | **TRIBE / HyperTRIBE / STAMP** (RBP-fusion editing) | A>G (ADAR) / C>U (APOBEC1) (substitution) | HyperTRIBE scripts (Perl+MySQL); SAILOR/FLARE (Snakemake+Singularity) | HyperTRIBE edit-site TSV (1-based); SAILOR ranked BED6 (0-based) | **Y** |
| 16 | **HITS-CLIP CIMS/CITS** | 1-nt deletions (CIMS), truncations (CITS) | CTK (own conda channel; Perl) | `parseAlignment.pl --mutation-file` BED6+5 (0-based); `CIMS.pl`: #chrom, chromStart, chromEnd, name[k=K][m=M], score, strand, tagNumber, mutationFreq, FDR, count | **Y** (CITS on PE → X1) |
| 17 | **Bisulfite SNP / ASM** | C>T SNP vs conversion via opposite-strand G/A | BISCUIT 1.10.2 (bioconda); Bis-SNP 1.0.1 (GATK3) | BISCUIT VCF (FORMAT GT, DP, SP, AC, AF1, CV, BT, GL1, GQ) → `vcf2bed -t snp|cg` (chr, start 0-based, end, beta, cov) | **Y** (cross-strand join in DuckDB) |
| 18 | **ac4C-seq / RedaC:T-seq** | C>T at reduced ac4C | lab scripts | lab-specific site tables (see Part 2 §03) | **Y** (dUTP → X2) |
| 19 | **mim-tRNAseq / m1A-m3C RT signatures** | mismatch spectrum + RT stops | mim-tRNAseq (bioconda) | `mods/mismatchTable.csv` (isodecoder, pos, type, proportion, condition, bam, cov), `mods/RTstopTable.csv` (isodecoder, pos, proportion, cov, condition, bam) | **Y** (stops via pad; PE → X1) |
| 20 | **dSMF** | GCH/HCG (DGCHN / NWCGW) conversion | SingleMoleculeFootprinting 2.7 + QuasR (Bioconductor) | `CallContextMethylation()` → MethGR (GRanges, 1-based; `<sample>_Coverage`, `<sample>_MethRate`, GenomicContext) + MethSM per-read matrices | **Needs X2** (QuasR non-standard BAM) or realignment |
| 21 | **RNA bisulfite m5C** | unconverted C (substitution) | meRanTK (scripts); UBS-seq pipeline | meRanCall TSV (#SeqID, refPos, refStrand, refBase, cov, C_count, methRate, mut_count, mutRate, CR, SNR, CalledBase, CB_count, state, 95_CI_lower/upper, p-values, score, seqContext, geneName, candidateName) | **Y** (dUTP → X2) |
| 22 | TAB-seq / ACE-seq / oxBS | unconverted C (5hmC) | as WGBS + mlml | as WGBS | **Y** (adds no versatility) |
| — | Dropped | biomodal 5/6-base (resolved from R1+R2 pre-alignment, MM tags), Fiber-seq/SMAC-seq (kinetics), antibody m6A/Ψ peaks (coverage), SPRINT hyper-editing (realignment), prime-editing classification, duplex/NanoSeq consensus construction, truncation lesion maps (no clean comparator; motivate X1 only) | | | N / out of scope |

**Recommended panel** (covers every query shape with the cleanest comparators): (1) SLAM-seq vs SLAMdunk + bakR per-read CSV [substitution, per-read counts, fwd and dUTP strandedness], (2) A-to-I editing vs JACUSA2 + REDItools3 [substitution, strand conventions, overlaps], (3) SHAPE/DMS-MaP vs ShapeMapper2 + SEISMIC + rf-count [indels; three tools disagree with each other], (4) iCLIP/eCLIP vs htseq-clip + PureCLIP [truncation / read-end], (5) aDNA vs mapDamage2 + DamageProfiler [position from read ends, DNA], (6) TAPS vs rastair or NOMe-seq vs Bismark/BISCUIT [inverted logic / dual context], plus (7) bcftools/pysamstats as the SNV baseline. GLORI and TRIBE/STAMP are strong alternates.

### 01.2 Candidate extensions (query syntax first; walk/CLI modes that change what a pattern means second)

| ID | Extension | Needed by (hard need) | Helps (convenience/clarity) | Notes |
|---|---|---|---|---|
| **X1** | **End-specific pad codes**: distinct symbols for "past the 5' end as sequenced", "past the 3' end as sequenced", "off contig"; decide whether clipped bases count as read (end = sequenced end) or not (end = aligned end), possibly two variants | iCLIP/eCLIP on PE and BAM-tag output; CITS (CTK) on PE; JACUSA2 rt-arrest/lrt-arrest; truncation lesion maps (CPD-seq, HydEn-seq); tRNA RT stops on PE | aDNA (first/last base patterns), SHAPE RT-stop mode, TRAC-seq | Today `_N` vs `N_` + a read-2 flag filter works in parquet only; soft-clip skip silently moves "the end" to the aligned end. |
| **X2** | **Library/strand modes** (CLI, but it changes pattern semantics): `rna-forward` (fr-secondstrand), `rna-reverse` (fr-firststrand/dUTP), `unstranded`/reference-orientation walk, `as-sequenced` (each read in its own orientation, no mate flip), `pbat`, `non-directional` with strand from tag (`tag:XG`/`XR` Bismark, `YD` bwa-meth, `XQ` QuasR, `YZ` hisat-3n), per-read strand inference from conversion counts, flag-only strand | scNMT/scNOMe (pbat/non-directional), dSMF/QuasR BAMs, tagmented TAPS/FOODIE/SMRF-seq (inference), hisat-3n-table parity for GLORI/eTAM/UBS (`tag:YZ`), Bullseye-exact and Alu Editing Index (unstranded) | DART, STAMP, RNA-BS, BID, ac4C, TimeLapse/TUC, editing on dUTP libraries (write RNA-sense patterns once); aDNA/fgbio parity (`as-sequenced`) | `src/tags.rs` doc comment already anticipates a non-directional variant. Strictly not query syntax, but without it several patterns cannot be written at all (unstranded data cannot be collected by one pattern). |
| **X3** | **Indel placement**: walk option `--indels aligner|left|right` (normalize within repeat), plus a code/predicate "column lies in the equivalence window of an ambiguous indel" | SHAPE-MaP/DMS-MaPseq parity (ShapeMapper shifts 5', rf-count shifts 3' within ±10, SEISMIC marks all positions ambiguous); BID-seq/PRAISE U-tract deletions | CRISPResso2 indel placement, CTK CIMS, CRISPR windows | The one gap that blocks exact parity for the highest-discrepancy family. |
| **X4** | **Proximity / run predicates**: e.g. `near(gap|junction|clip|end, k)`, homopolymer run `A{5,}` | — (expressible as OR over shifted patterns) | REDItools homopolymer and splice-distance filters, JACUSA2 D/I/S/Y feature filters, SPRINT "≥5 nt from any block edge", SAILOR junction overhang, BISCUIT/REDItools end-trim | Pure convenience; keeps pattern width manageable. |
| **X5** | **Quality predicate inside patterns** (e.g. column Q≥30 incl. neighbours) | — (`^` captures + DuckDB) | ShapeMapper neighbour-Q30 rule; BAM-tag output where no DuckDB step exists | Conflicts with the "nothing filtered for you" principle; only for tag output. |
| **X6** | **Soft-clip semantics**: distinct clip code (vs insertion) and/or `--soft-clips project` placing clipped bases against the reference they would cover | rastair `--rescue-soft-clip-cpg`; exact mapDamage2 position convention (counts from sequenced end incl. clips vs DamageProfiler/PyDamage from aligned end) | SAILOR/HyperTRIBE CIGAR-based read drops; CRISPResso2 | Today a clip is `base / gap`, identical to an insertion (`src/alignment.rs:277`). |
| **X7** | **`overlap` policies** (not syntax): max-mean-phred mate (ShapeMapper), bitwise AND (SEISMIC), both-mates-agree (rf-count, GRAND-SLAM), drop R2 bases (BISCUIT), samtools-style Q adjust (mpileup/AEI) | exact parity with those tools on overlapping PE | all PE comparisons | Each comparator uses a different rule; making the rule explicit is itself a paper point. |
| **X8** | **Strand/mate predicate in `where`** (`strand in (OT, CTOT)`, `mate == 2`) | — for parquet (flag column); needed for **BAM tag output** only | Revelio-style SNP masking tags, R1-only oxoG tags, eCLIP R2-start tags | |
| **X9** | Output conveniences (not syntax): per-read count tag (`TC:i` NGM/SLAMdunk-compatible), alignment-column offset column | SLAMdunk-compatible BAM export | mapDamage parity | |
| X10 | Cross-mate / third row (mate base at same reference column) | biomodal 5/6-base only | — | **Not recommended**: no comparable intermediate exists. |

Not extensions (explicitly DuckDB/samtools): conversions-per-read cutoffs (SLAM, GLORI ≤3 A, SAILOR non-target mismatches), read-end exclusion by `off_5p/off_3p`, BQ/MAPQ/dup filters, SNP masking by opposite-strand evidence or VCF join, coverage thresholds, mutation collapsing within N nt, TAPS polarity flip, strand collapse, oxBS/TAB subtraction.

### 01.3 Cross-cutting findings worth a paper figure

- **The same event, different numbers**: SHAPE/DMS tools place an ambiguous deletion 5', 3' or nowhere; insertions are credited 5' (ShapeMapper, rf-count) vs 3' (SEISMIC); collapse windows 6 vs 2 vs off; BQ Q30-with-neighbours vs Q25 vs Q20-mismatch-only. In repeats this moves reactivity between adjacent nucleotides and can flip paired/unpaired calls.
- **Crosslink offset conventions**: htseq-clip default offset 0 (first aligned base), Skipper same, JACUSA2 arrest = terminal aligned base, PureCLIP/iCount = start −1. A 1-nt shift moves sites relative to motifs.
- **CIGAR-naive counters** are common and consequential: HyperTRIBE drops any read with I/D/S; bam2bakR/fastq2EZbakR drop reads with indels; wavClusteR's MD parser mis-positions sites after deletions/clips (**verified by execution**); SAILOR mis-parses ≥10-nt left soft clips as right clips (**verified by execution**); QuasR's methylation counter ignores the CIGAR (safe only for its ungapped aligner); PARalyzer never parses CIGAR; SPRINT likely misparses H/=/X. alnbase's CIGAR-exact walk is a clean contrast.
- **Strand conventions** are the most frequent silent failure: SLAMdunk takes strand from the reverse flag alone (PE and reverse-stranded broken); bam2bakR treats unstranded as forward; SAILOR's per-read sense uses `is_reverse` without mate number; PARalyzer tests `FLAG == 16 || 272`; asTair/SAILOR route by exact FLAG values, dropping improper pairs; REDItools `-s` vs JACUSA2 `FR-SECONDSTRAND`/`RF-FIRSTSTRAND` naming mix-ups (JACUSA2 issue #78).
- **Hidden gates**: SAILOR requires `bcftools call -c` to call a variant before scoring (drops low-fraction edits); SLAMdunk SNP masking at VAF ≥ 0.8 leaves heterozygous T/C SNPs as "labelled"; GRAND-SLAM masks sites with >30 % mismatches (can hide true labelling hotspots); bcftools subsamples to 250 reads; GLORI silently excludes genes with conversion rate ≥ 0.2; Bismark NOMe GpC output is silently empty without `--CX`.
- **Overlapping mates**: REDItools, pysamstats (despite overlap detection flag), meRanCall (when the second mate has higher BQ), DART/Bullseye count both; hisat-3n-table may retract the wrong observation on disagreement; BISCUIT drops R2 bases. alnbase `overlap` makes the policy explicit.
- **Biology-flipping examples with literature**: ac4C in human mRNA (two 2024 Mol Cell papers disagree; turns on dedup, replicate pooling, mismatch-type normalization); BID vs BACS Ψ sites agree ~40 %; DMS G>A counting makes Gs look reactive (Mitchell 2023); SNPs at CpG masquerade as methylation in TAPS; aDNA authentication verdict depends on position-1 convention under soft clips; DART Bullseye vs CTK ~50 % overlap (Bullseye issue #4).

### 01.4 Logistics

- **Use one shared BAM per comparison** and feed it to both tools wherever the comparator accepts a BAM, so differences are counting policy, not alignment. Exceptions: CRISPResso2 re-aligns even with `--bam_input` (use its `--bam_output` as the shared BAM); ShapeMapper merges/trims before aligning (use `shapemapper_mutation_parser` on external SAM); GLORI-tools BAMs have `_AG_converted` contig suffixes (rename before `alnbase index`); SLAMdunk needs NextGenMap `MP` tags (a STAR/HISAT BAM silently yields zero conversions); GRAND-SLAM needs its CIT conversion.
- **MD tags** are required by ShapeMapper parser, PyDamage, CTK, SAILOR, wavClusteR → run `samtools calmd` on the shared BAM.
- **Coordinates**: 1-based — REDItools, HyperTRIBE, ShapeMapper `profile.txt`, SEISMIC CSV, hisat-3n-table, GATK counts, mapDamage `Pos` (from read end); 0-based — BED outputs (JACUSA2, htseq-clip, PureCLIP, SAILOR), rf-rctools arrays, ShapeMapper `parsed.mut`; **wavClusteR BED export writes 1-based starts (off-by-one)**; GRAND-SLAM `.snpdata` header has Coverage/Mismatches swapped; SLAMdunk `alleyoop dump` header/column order disagree.
- **Install**: bioconda — slamdunk, jacusa2, reditools3, pureclip, htseq-clip, icount, seismic-rna, rnaframework, mapdamage2, damageprofiler, pydamage, rastair, bismark, biscuit, picard, fgbio, bcftools, samtools, pysamstats, crispresso2, paralyzer (academic licence), mim-tRNAseq, bullseye. Not bioconda — ShapeMapper 2 (tarball; bioconda `shapemapper` is v1.2), CTK (own channel), REDItools2, SPRINT (Python 2), Alu Editing Index (Docker, ~12 GB resources), HyperTRIBE (MySQL), SAILOR/FLARE (Singularity), GRAND-SLAM current (Maven; bioconda `gedi` is 1.0.6a), fastq2EZbakR (Snakemake), asTair (pip, deprecated), GLORI-tools/BID-pipe (scripts; some BID-pipe steps binary-only), SingleMoleculeFootprinting/QuasR/wavClusteR (Bioconductor), BE-Analyzer (web only), Bis-SNP (GATK3; effectively unrunnable).
- **Version pinning**: SEISMIC renamed steps in 0.25/0.26; REDItools v2 vs v3 defaults differ (v3 strict `-me 1`, `-C` meaning changed); mapDamage master (2.3.0a0) output layout differs from bioconda 2.2.3; rastair counting changed 2.1→2.2; CRISPResso2 window off-by-one fix (PR #651) unreleased.
- **Fast validation data** (sizes as found): ShapeMapper TPP riboswitch 3×8k pairs (4.6 MB total) + E. coli rRNA 15k pairs with 16S/23S `.ct`; SEISMIC simulated HIV RRE (~1.4 MB) + `seismic sim`; SLAMdunk repo test set + `splash` simulator; bam2bakR/fastq2EZbakR `.test` PE dUTP chr21 BAMs (2–3.3 MB); nf-core/slamseq chr8 (17–59 MB FASTQ); wavClusteR `example.bam` (61 kB); REDItools2 chr21 BAM (1.4 MB); Alu Editing Index bundled BAMs (71 MB) with expected CSV; htseq-clip test BAMs (9 kB, 73 kB) with expected BEDs; ENCODE eCLIP ENCFF280ONP (229 MB); GSE99249 ADAR1-KO; HyperTRIBE `examples/` chr2L SAMs (33 + 30 MB); nf-core/eager mammoth mtDNA BAM (452 kB); PyDamage `tests/data/aligned.bam` (880 kB); gargammel simulator; rastair `test.bam` (1.6 MB) + `ug_tagmented.bam` (32 MB); asTair tests (~1 MB); bcftools `test/mpileup/mpileup.1.bam` (68 kB); CRISPResso2 tests (~139 kB) + 25k-read BE demo; Bullseye test BAMs (860 kB); mim-tRNAseq bundled FASTQs (~52 MB); meRanTK test data (118 MB).
- **Caveats of this survey**: GitHub API rate limit was hit during the CLIP/editing search (some issue threads read via web fetch; release tags not all checked); iCount example server unreachable; PRAISE GEO accession, ACE-seq accession, NRF1pair.bam and Zenodo lambda sizes not verified.

### 1.5 Other candidates considered

- **SMRF-seq (R-loop footprinting, non-denaturing bisulfite on displaced ssDNA)**: genuinely read-vs-reference (strand-specific C>T runs), long-read PacBio CCS aligned with Bismark; comparator footLoop (https://github.com/srhartono/footLoop; Bismark 0.20, window of 20 Cs with ≥55 % conversion, peaks ≥100 bp, strand assigned by C>T vs G>A ratio). Needs per-read strand inference (X2). Niche; mention as a generality example.
- **MAPit-FENGC, DAF-seq, FOODIE**: deaminase/methyltransferase footprinting variants; MAPit expressible, FOODIE/DAF-seq need per-read strand inference (X2); comparators custom or long-read.
- **Allele-specific expression counting** (GATK ASEReadCounter, phASER): a site-restricted SNV pileup; covered by the SNV baseline.
- **Ribo-seq P-site, PRO-seq, ChIP/ATAC**: positional only, no read-vs-reference base comparison — out of scope (though a read-start pattern could reproduce 5'-end counts, it adds nothing over bedtools).
- **Deep mutational scanning** (Enrich2, DiMSum): variant calling from full read sequences/codons, not per-column comparison — dropped.

---

## Part 2. Per-family detail

## Family 01 — Nucleotide-conversion RNA metabolic labeling and PAR-CLIP

Scope: SLAM-seq (SLAMdunk), TimeLapse-seq and TUC-seq (bam2bakR / fastq2EZbakR -> bakR / EZbakR; GRAND-SLAM -> grandR), GRAND-SLAM itself, s6G / TILAC G>A, and PAR-CLIP T>C (PARalyzer, wavClusteR; others noted). All these assays have the same read-level signal: **a U>C (4sU) or G>A (s6G/6SG) substitution in RNA-sense orientation**, counted per base and aggregated per read (n = mutable reference bases covered, k = conversions) or per site. None needs anything beyond a single-column substitution pattern, so all are **expressible in alnbase today**. The differences between tools sit in filters, strand logic, mate-overlap handling and SNP masking. Those differences are where validation is interesting.

Source code read for this section (shallow clones in `clones/`):
- `slamdunk` @ 14270e5 (v0.4.3, 2024-04-25)
- `gedi` @ cd2a244 (GRAND-SLAM 2.0.7, 2026-09-15)
- `fastq2EZbakR` (2026-06-30) and `bam2bakR` (2025-09-16)
- `bakR` 1.0.1 and `EZbakR` 0.2.1 (GitHub)
- `wavClusteR` 2.45.1 (Bioconductor git)
- `PARalyzer_v1_5.tar.gz`, from the URL in the bioconda recipe; only compiled classes ship, and they were inspected with `javap -c`

#### How alnbase orients these libraries (common to every assay below)

Under `--library directional`, the walk follows the strand implied by read 1. R1 forward and R2 reverse are walked top-strand. R1 reverse and R2 forward are walked reverse-complemented. A single-end read counts as R1.

| Library type (typical kit) | RNA-sense signal of 4sU | Pattern under `--library directional` | Pattern for s6G |
|---|---|---|---|
| Forward-stranded SE/PE: R1 is sense (Lexogen QuantSeq FWD, the SLAMseq kit, PAR-CLIP small-RNA libraries, 10x 3') | T>C | `read="C" refr="T"`; coverage `read="N" refr="T"` | `read="A" refr="G"` |
| Reverse-stranded dUTP: R1 is antisense (TruSeq stranded, NEB Ultra II Directional, KAPA stranded; most TimeLapse-seq and TUC-seq total-RNA libraries) | T>C in RNA, seen as A>G in the walk | `read="G" refr="A"`; coverage `read="N" refr="A"` | `read="T" refr="C"` |
| Unstranded | Unknown per read | Emit both `C/T` and `G/A`, then join to the annotated gene strand in DuckDB | same idea |

`N` matches any base but no gap, pad or junction. That makes `read="N" refr="T"` exactly the set of T positions that GRAND-SLAM, bakR and SLAMdunk treat as covered.

Everything else these tools do can be rebuilt from alnbase parquet (`qual`, `off_5p`, `off_3p`, `refr_pos`, `-f flag,mapq,NH,...`) plus DuckDB, or from samtools upstream:
- base-quality thresholds
- read-end trimming
- per-read counts of n and k
- the at-least-k-conversions threshold
- SNP masking by joining to a VCF
- multimapper weighting by 1/NH

Mate-overlap "double hits" can be rebuilt by a self-join on (name, refr_pos), or approximated with `alnbase overlap --mismatch-qual zero-both`.

---

### 01.1 SLAM-seq (thiol(SH)-linked alkylation), analysed with SLAMdunk

#### (a) Signal
- 4sU is incorporated into nascent RNA. Iodoacetamide carboxyamidomethylates 4sU, and RT then reads the base as C. The read shows a T>C mismatch in RNA sense. Context: any T; no motif. Typical conversion is about 2-5% of Ts in labelled RNA ([Herzog et al. 2017, Nat Methods](https://www.nature.com/articles/nmeth.4435)).
- Canonical library: Lexogen QuantSeq 3' mRNA-seq FWD, single-end, R1 sense, reads piling up in 3'UTRs.
- SLAMdunk models exactly that design: T>C on forward reads, A>G on reverse reads (`SlamSeqFile.py:137-143`, `isTCMismatch`). It is single-end only by design ([issue #165](https://github.com/t-neumann/slamdunk/issues/165), [#57](https://github.com/t-neumann/slamdunk/issues/57)). Reverse-stranded input is not supported; the maintainer's advice is to reverse-complement the FASTQ ([issue #14](https://github.com/t-neumann/slamdunk/issues/14), [#149](https://github.com/t-neumann/slamdunk/issues/149), [#175](https://github.com/t-neumann/slamdunk/issues/175)).

#### (b) Comparator
- **SLAMdunk 0.4.3**, on bioconda as `slamdunk` 0.4.3. It pins NextGenMap 0.5.5 (`version.py`) and needs VarScan 2.4.x, samtools and R. There is a Docker image, and a nf-core/slamseq 1.0.0 pipeline (archived, DSL1).
- Install complexity: low-medium through conda. NGM path problems are known ([#45](https://github.com/t-neumann/slamdunk/issues/45)).
- Its hard requirement is **NGM-specific BAM tags**, discussed under (g).

#### (c) Output formats
1. **`*_tcount.tsv`** (from `slamdunk count`). Two `#` header lines, then the columns (`SlamSeqFile.py:80`): `Chromosome Start End Name Length Strand ConversionRate ReadsCPM Tcontent CoverageOnTs ConversionsOnTs ReadCount TcReadCount multimapCount ConversionRateLower ConversionRateUpper`
   - `Start`/`End` are copied from the input BED: 0-based, half-open.
   - `Tcontent` = number of T (A when `Strand` is `-`) in the reference interval.
   - `CoverageOnTs` = sum over T positions of reads whose span covers the position.
   - `ConversionsOnTs` = T>C events on those positions.
   - `ConversionRate` = ConversionsOnTs / CoverageOnTs.
   - `ReadCount` = strand-matching reads overlapping the interval.
   - `TcReadCount` = reads with at least `-c` conversions.
   - `multimapCount` = reads with MAPQ 0.
   - `ConversionRateLower`/`Upper` are always `-1.0` (`tcounter.py:125-330`).
2. **`alleyoop collapse`** gives per gene: `gene_name length readsCPM conversionRate Tcontent coverageOnTs conversionsOnTs readCount tcReadCount multimapCount` (`tcounter.py:96`).
3. **`alleyoop positional-tracks`** writes bedGraphs: `_TC_rates_genomewide`, `_AG_rates_genomewide`, `_coverage_{plus,minus}`, `_TC_conversions`, `_AG_conversions`, `_coverage_T`, `_coverage_A`.
   - Plus-strand reads feed the T>C tracks and minus-strand reads the A>G tracks.
   - Coordinates: internally `referencePosition` is shifted by -1 (`startPosition=1`, `SlamSeqFile.py:451`), and the writer adds +1 back (`tcounter.py:417-470`). Net output is 0-based bedGraph.
   - The run-length writer never flushes the last run of each chromosome, so trailing values are dropped.
4. **`alleyoop dump`**, per read. The header is `Name Direction Sequence Mismatches tcCount ConversionRates`, but rows are written as name, direction, sequence, tcCount, conversionRates, mismatches (`SlamSeqFile.py:198-219`). **Header and data column order disagree** ([#139](https://github.com/t-neumann/slamdunk/issues/139) asks what the output means).
5. **BAM tags written by NGM:**
   - `TC:i` = number of T>C conversions, or A>G on reverse reads.
   - `RA:Z` = 25-cell ref x read conversion matrix.
   - `MP:Z` = `type:readPos:refPos,...` (1-based).

   Documented in the [SLAMdunk supplement](https://static-content.springer.com/esm/art%3A10.1186%2Fs12859-019-2849-7/MediaObjects/12859_2019_2849_MOESM1_ESM.pdf) and [#104](https://github.com/t-neumann/slamdunk/issues/104).

#### (d) Implicit filters and defaults (`slamdunk.py:342-430` unless noted)
- **Mapping:**
  - `-5/--trim-5p 12`: 12 nt clipped from the 5' end before NGM.
  - `-a/--max-polya 4`
  - `-n/--topn 1`
  - NGM run with `--slam-seq 2`, a T>C-tolerant scoring scheme (`mapper.py:77`).
- **Filter** (`filter.py:243-260`):
  - `-mq 2`: MAPQ >= 2.
  - `-mi 0.95`: `XI` identity >= 0.95.
  - `-nm -1`: no NM cap.
  - Supplying `-b/-fb` BED switches to multimapper-rescue mode. The MQ filter is overridden to 0, and a MAPQ-0 read is kept if all its hits fall in the same UTR, with one hit written at random and an `RD` tag (`filter.py:79-180`).
  - Duplicates are not removed; `alleyoop dedup` is optional.
- **SNPs** (`snps.py`):
  - `samtools mpileup -B -A -f ref` piped to `varscan mpileup2snp --strand-filter 0 --min-var-freq 0.8 --min-coverage 10 --variants 1`.
  - mpileup's default BQ >= 13 applies.
  - Only REF T/ALT C SNPs mask forward reads, and only A/G SNPs mask reverse reads (`utils/SNPtools.py`).
- **Counting** (`SlamSeqFile.py:313-400`):
  - `-q/--min-base-qual 27`. A mismatch below Q27 is dropped from the mismatch list, but its T position still counts in `CoverageOnTs`.
  - `-c/--conversion-threshold 1`. Help text elsewhere still mentions ">=2" ([#177](https://github.com/t-neumann/slamdunk/issues/177)).
  - Reads on the strand opposite the BED interval are skipped (`SlamSeqFile.py:369`).
  - `isMultimapper` = MAPQ == 0.
  - Mismatches come **only from the `MP` tag**, never from MD or the reference.
  - Per-read `n` for the MLE output = count of T in `query_sequence`, which includes soft clips and read Ts that are C>T mismatches, plus reference-T mismatches (`SlamSeqFile.py:181`, `tcounter.py:233-240`).
  - Coverage is `range(startRefPos, endRefPos)`, the whole reference span, so deleted and intronic positions count as covered (`tcounter.py:246-248`).

#### (e) Expressible in alnbase? **Yes.**
- Queries: `tc: read="C" refr="T"`; `t: read="N" refr="T"`.
- SLAMdunk's exact semantics are all downstream in DuckDB:
  - Q27 applies only to the conversion rows.
  - Coverage from read span needs `-f pos,end` or a CIGAR span.
  - 5' trimming is `off_5p >= 12` if the trimming was not done before alignment.
  - The strand-to-BED match uses `flag & 16`.
- To reproduce SLAMdunk's per-read `n` exactly, use `--soft-clips emit` and count read Ts, which reproduces the soft-clip inclusion.
- **No query-syntax extension needed.** One optional output-side (non-query) extension: a per-read count tag kind (e.g. `[tag.TC.count]`). It would reproduce NGM's `TC:i` and let `alleyoop read-separator`-style splitting run in samtools.

#### (f) Test data
- SLAMdunk repo `slamdunk/test/data/` (in git): `reads.fq` (1.2 kB, 40 nt reads), `ref.fa` (237 kB, mm10 chr5 excerpt), `actb.bed`. Expected output `reads_slamdunk_mapped_filtered_tcount.tsv`: 8 reads, 4 TC reads, `ConversionRate` 0.0222. Run with `test/test_sample.sh`.
- **nf-core/test-datasets branch `slamseq`** ([README](https://github.com/nf-core/test-datasets/tree/slamseq)):
  - 6 SE100 FASTQs from MOLM-13 DMSO and NVP-2, subset to chr8, 17-59 MB each.
  - Plus `hg38_chr8.fa.gz` (46 MB) and `hg38_refseq_3UTR.chr8.bed`.
  - Origin: Muhar et al. 2018 *Science*, GEO [GSE111463](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111463).
- Original SLAM-seq: Herzog 2017, GEO GSE99978 (mESC pulse-chase, QuantSeq FWD).
- **Ground truth:** SLAMdunk ships `splash`, a read-level simulator (`dunks/simulator.py`). It builds 3'UTR reads with set half-lives and conversion rates and writes the true conversion positions to a SAM (`convertRead`, `addTcConversionsToReads`). That makes it the easiest per-base ground truth in this family.

#### (g) Correctness issues and contestable choices
1. **Silent zero signal on non-NGM BAMs.** `fillMismatchesNGM` only reads `MP:Z`. A STAR or HISAT2 BAM without `MP` gives `mismatchList=[]` and a zero conversion rate for every read, with no error at `count`; `filter` crashes on missing `XI` ([#136](https://github.com/t-neumann/slamdunk/issues/136); maintainer: "Slamdunk really only supports NGM alignments"). Users feeding HISAT-3N BAMs were told to "mind the tags" ([#143](https://github.com/t-neumann/slamdunk/issues/143)). An alnbase comparison on the same BAM makes this failure mode visible at once.
2. **Heterozygous T/C SNPs are not masked.** `--min-var-freq 0.8` keeps only homozygous variants, so a het T/C SNP contributes about 50% "conversions" at that site in every sample, no-4sU controls included. For a lowly expressed 3'UTR with one het SNP, that is enough to push `ConversionRate` or `TcReadCount` above background, making a gene look "newly transcribed" or "fast-turnover". Comparing against alnbase T>C per-site rates in a no-4sU control shows these sites directly.
3. **Paired-end and reverse-stranded data are mis-assigned.** Strand is taken from `is_reverse` alone. Mate 2 of a forward-stranded pair is dropped as "antisense" or counted with the wrong conversion type. Users saw T>C only in R2 and A>G in R1 ([#175](https://github.com/t-neumann/slamdunk/issues/175)). alnbase's FLAG-based `directional` orientation handles both mates.
4. **Low-quality mismatches are counted as matches.** Below Q27 the conversion is discarded but the T stays in the denominator (`SlamSeqFile.py:330`; `CoverageOnTs` counts the read span). That deflates `ConversionRate` in low-quality read tails. It is contestable: GRAND-SLAM uses no base quality at all, and bakR drops the base from both numerator and denominator.
5. **Coverage counts read span, not aligned bases.** Deletions and `N` gaps count as covered Ts (`tcounter.py:246`). This is harmless for QuantSeq but wrong for spliced total-RNA reads (the HISAT-3N workaround in #143/#149).
6. **The MLE `n` is taken from the read, not the reference.** A C>T sequencing error adds a T to `n`, and soft clips add Ts (`SlamSeqFile.py:181`).
7. The `alleyoop dump` header is out of order (see c.4). The last bedGraph run per chromosome is not written. The BED parser is buggy ([#11](https://github.com/t-neumann/slamdunk/issues/11), open).

---

### 01.2 GRAND-SLAM (GEDI), for SLAM-seq, TimeLapse-seq, TUC-seq and scSLAM-seq

#### (a) Signal
- Same T>C, oriented by a `-strandness Sense|Antisense|Unspecific` setting (default auto-detect).
- For each gene it considers reads on the gene's strand (Sense), the opposite strand (Antisense), or both. For opposite-strand reads the reference is reverse-complemented and the conversion base becomes A (`SlamCollector.java:214-262`).
- Within a pair, mate 2 is expected to show the complement: `checkGenomic/checkRead` flip to A/G for the second read (`SlamCollector.java:301-313`).

#### (b) Comparator
- **GRAND-SLAM 2.0.7**, built into GEDI (`executables/Slam.java:20-25`, "2.0.7: geneData/extData/snpData files are removed in the end, tsv is gzipped").
- Install: Java 11 plus a Maven build from git (`mvn -f gedi/Gedi package`), then `gedi -e IndexGenome` on FASTA and GTF. Bioconda `gedi` is only 1.0.6a and PRICE-oriented. Complexity: **medium-high**.
- Input must carry `MD` (and ideally `NH`); the [wiki](https://github.com/erhard-lab/gedi/wiki/GRAND-SLAM) recommends STAR `--outSAMattributes MD NH`. BAMs are converted to CIT with `bamlist2cit`.
- Downstream: grandR (CRAN 0.2.7; [Rummel et al. 2023, Nat Commun](https://www.nature.com/articles/s41467-023-39163-4)).
- Paper: [Jürges, Dölken & Erhard 2018, Bioinformatics](https://academic.oup.com/bioinformatics/article/34/13/i218/5045735).

#### (c) Output formats (`SlamParameterSet.java`, `SlamInfer.java:237-252`, `SlamCollectStatistics.java:188-300`)
1. **`${prefix}.tsv.gz`**, one row per gene, in this column order:
   1. `Gene`, `Symbol`
   2. for each condition: `<cond> Readcount`
   3. for each non-no4sU condition: `<cond> 0.05 quantile`, `<cond> Mean`, `<cond> MAP`, `<cond> 0.95 quantile`
   4. for each non-no4sU condition: `<cond> alpha`, `<cond> beta` (beta posterior of NTR)
   5. with `-full`, for each condition: `<cond> Conversions`, `<cond> Coverage`, `<cond> Double-Hits`, `<cond> Double-Hit Coverage`
   6. with `-full`: `<cond> min2`
   7. `Length`

   Read counts are fractional under the default multimapper `Weight` mode.
2. **`${prefix}.mismatches.tsv`**: `Category Condition Orientation(First|Second) Genomic Read Coverage Mismatches`. Categories include `Exonic`, `ExonicAntisense`, `Intronic`, ...
3. **`${prefix}.mismatchdetails.tsv`**: `Category Genomic Read Position Overlap(0/1) Opposite(0/1) Coverage Mismatches`. `Position` is 0-based along the concatenated read1+read2 geometry. This is the file used to choose `-trim5p`/`-trim3p`.
4. **`${prefix}.binom.tsv`** / **`binomOverlap.tsv`**: `n d Condition Type count`. This is the histogram of reads by (number of Ts, number of conversions). **It is the file most directly comparable to an alnbase per-read aggregation.**
5. **`${prefix}.doublehit.tsv`**: `Category Condition Genomic Read Coverage Hits`.
6. **`${prefix}.ntrstat.tsv`**: `Mode Type Condition T->C p_new p_old ntr_lower ntr ntr_upper`. `${prefix}.rates.tsv` holds the estimated `single_old`/`single_new`/`double_*` rates.
7. **`${prefix}.snpdata`**: header `Location Coverage Mismatches P value`. The code appends `mismatches` then `coverage` (`SlamDetectSnps.java:83` vs `:142-146`), so **the header swaps Coverage and Mismatches**. Location is `chr:pos` with GEDI 0-based positions. The file is deleted at the end since 2.0.7 (`setRemoveFile(true)`).
8. **`${prefix}.strandness`**: the detected mode.

#### (d) Implicit filters and defaults
- **No base-quality filter at all.** CIT stores no qualities; when GEDI reconstitutes records it writes `I` for every base (`BamUtils.java:279`). A separate `gedi -e BamPhredFilter` (default `-min 28`) can edit MD before conversion.
- `-trim5p 0`, `-trim3p 0` (`SlamParameterSet.java:24-25`; the `SlamCollector` field default of 25 is overridden). Trimming applies only outside the mate overlap (`SlamCollector.java:294-298`).
- **SNPs** (`SlamDetectSnps.collect`):
  - Any mismatch type is considered; the T>C-only restriction is commented out.
  - Data are pooled over all conditions within a gene ±1 kb.
  - A position is called a SNP if `Beta.cumulative(0.3; k+0.3, n-k+1) < 0.001` (`-snpConv 0.3`, `-snppval 0.001`).
  - Mismatches inside the mate overlap count only if seen in both mates.
  - `-newsnp` restricts SNP calling to no4sU samples. The code comment says "SNPs *MUST* be called with pooled data!".
- **Strandness auto-detect:** Sense if sense reads > 2 x antisense reads, the reverse for Antisense, otherwise Unspecific.
- **Multimappers:** `-mode Weight` divides by NH. **Overlapping genes:** `-overlap All`.
- **Genes:** protein-coding only unless `-allGenes`. Reads must be transcript-compatible, or `-lenient`.
- **Mate overlap:** a conversion in the overlap counts only as a "double hit", and only if both mates show it (`SlamCollector.java:318-325`). Overlap Ts are excluded from the single `n`.
- Mismatches involving N or a non-ACGT genomic base are skipped. An inconsistency between MD and the genome throws `RuntimeException` (`SlamCollector.java:343`).
- **Error model:** no4sU conditions are recognised by `-no4sUpattern "no4sU|nos4U"`. Without one, the T>C error rate is `-errlm "TA*1.434"`, i.e. 1.434 x the T>A mismatch rate.

#### (e) Expressible in alnbase? **Yes.**
- Queries: the same two patterns (T>C, T coverage) in the orientation matching the library.
- Double-hit semantics: emit per-mate rows without `overlap`, then self-join on `name, refr_pos` in DuckDB.
- `binom.tsv` equivalence: `GROUP BY record/template -> (n, d)`.
- No extension needed.

#### (f) Test data
- fastq2EZbakR / bam2bakR `.test/data/bams/`: `WT_nos4U.bam` 2.0 MB, `WT_replicate_1.bam` 2.5 MB, `WT_replicate_2.bam` 3.3 MB. Also `genome/genome.fasta` (hg38 chr21, 47 MB) and `annotation/genome.gtf` (348 kB).
  - STAR 2.7.10b, PE, **reverse-stranded**, tags `NH HI AS NM MD` (checked from the BAM header).
  - The sample name `WT_nos4U` already matches GRAND-SLAM's default `-no4sUpattern`, so **one small BAM set runs GRAND-SLAM, bakR and alnbase**.
- Three-chemistry comparison in one cell line: [Boileau et al. 2021, Brief Bioinform](https://pmc.ncbi.nlm.nih.gov/articles/PMC8574959/), SRA BioProject **PRJNA726397**.
  - MCF-7, 4sU-tag, SLAM-seq, TimeLapse-seq and TUC-seq.
  - TruSeq stranded mRNA (antisense), PE100.
  - Their own counting used BQ > 20, >= 5 nt from read ends, and >= 1 T>C per pair for "labelled".
  - Too large for CI; subset one chromosome.

#### (g) Correctness issues and contestable choices
1. **No base quality.** In low-quality runs (e.g. NovaSeq Q2 tails), T>C errors inflate the old-RNA error rate. GRAND-SLAM absorbs this through its no4sU error estimate, but without a no4sU sample it falls back to `TA*1.434`. Tail artefacts are not T>A-proportional, so NTR can be biased up. An alnbase Q-stratified T>C rate in the no4sU control quantifies this.
2. **Any-mismatch SNP calling at 30% conversion.** Heterozygous SNPs (VAF ~ 0.5) are masked, unlike SLAMdunk. Genuine hotspots of 4sU conversion at high labelling (> 30% of reads at a site, e.g. long pulses) can be **masked as SNPs**, deflating NTR for highly labelled genes. This is the opposite failure to SLAMdunk's.
3. **Double hits only within the overlap.** PE libraries with short inserts lose single-mate conversions in the overlap. That is conservative, but `n` and `k` are no longer comparable to bam2bakR's rule of "the better-quality mate wins".
4. **Header swap in `.snpdata`** (c.7).
5. **`BamPhredFilter` indexes base quality by the MD offset** (`pos += genomic length`), ignoring soft clips and insertions (`BamPhredFilter.java` main loop). With soft clipping, the wrong base's quality decides whether a mismatch is removed.
6. Halfpipe authors report GRAND-SLAM's EM being outperformed across labelling efficiencies ([Halfpipe bioRxiv 2024](https://www.biorxiv.org/content/10.1101/2024.09.19.613510.full)). That comparison is model-level, not base-level.

---

### 01.3 TimeLapse-seq and TUC-seq, analysed with bam2bakR / fastq2EZbakR -> bakR / EZbakR

#### (a) Signal
- **TimeLapse-seq:** NaIO4 + 2,2,2-trifluoroethylamine converts 4sU to a C analogue (T>C) ([Schofield et al. 2018, Nat Methods](https://www.nature.com/articles/nmeth.4582)). The s6G variant gives G>A. **TILAC** uses 4sU plus s6G.
- **TUC-seq:** OsO4/NH4Cl converts 4sU to genuine C ([Riml et al. 2017, Angew Chem](https://onlinelibrary.wiley.com/doi/10.1002/anie.201707465)).
- Libraries are usually standard stranded total-RNA or poly(A) kits: **dUTP reverse-stranded PE**, where R2 is the RNA sense. The same signal as SLAM-seq, in whichever orientation the library dictates.
- `mut_call.py` describes itself as "python implementation of TimeLapse mutation calling" (header, line 23). It is the Simon lab's maintained successor of TimeLapse.R.

#### (b) Comparators
- **fastq2EZbakR** (Snakemake; docs v0.3.0; [readthedocs](https://fastq2ezbakr.readthedocs.io)) and its predecessor **bam2bakR**.
- Downstream: **bakR** 1.0.1 (CRAN; [Vock & Simon 2023, RNA](https://doi.org/10.1261/rna.079451.122)) and **EZbakR** (CRAN 0.1.0 / GitHub 0.2.1; [Vock et al. 2024 EZbakR suite](https://pmc.ncbi.nlm.nih.gov/articles/PMC11507695/)).
- Neither pipeline is on bioconda; they install through a Snakemake conda env (`workflow/envs/full.yaml`). Complexity: medium (Snakemake + conda + STAR/HISAT2).
- Other options: GRAND-SLAM (1.2); `dynast` (PyPI 1.0.1) for scNT-seq / scSLAM; `rnalib`; Halfpipe (GitHub `IMSBCompBio/Halfpipe`, not on bioconda).

#### (c) Output formats
1. **`*_counts.csv(.gz)`**, one row per read or read pair (`mut_call.py:62-63`): `qname nA nC nT nG rname FR sj TA CA GA NA AT CT GT NT AC TC GC NC AG TG CG NG AN TN CN GN NN`. With `--mutPos`, `gmutloc tp` is added.
   - `FR` = whether the pair was complemented to RNA sense.
   - `sj` = the pair contains a `N` CIGAR operation.
   - Mutation columns are `<ref><read>` in RNA-sense orientation.
   - `n{A,C,G,T}` = reference base counts over quality-passing aligned bases.
2. **`cB`** CSV or partitioned parquet ("arrow"), written by DuckDB SQL in `merge_features_and_muts.R:338-378`: `sample, rname, sj, <feature columns e.g. GF, XF, exonic_bins, ...>, <mut types e.g. TC>, <nT...>, n`, where `n` = number of reads with that combination.
   - bakR requires `XF, sample, TC, nT, n` (`bakR/R/bakRData.R:182`).
   - EZbakR requires `sample`, `n` and matching `n<base>` for every mutation column (`EZbakR/R/EZbakRData.R:49-74`).
   - `cUP` is the same with `nT` averaged per mutation-count group.
3. **`cU`** (`--mutPos`): `rname gloc tp trials n`, with `gloc` 1-based (`mut_call.py:260-270`, `+1`).

#### (d) Implicit filters and defaults
- **Pre-filter** (`sort_filter.sh`): `samtools view -q 2 -F 0x4 -F 0x8 -F 0x100 -F 0x200 -F 0x400 -F 0x800`. MAPQ >= 2; unmapped, mate-unmapped, secondary, QC-fail, duplicate-flagged and supplementary reads removed; name-sorted plus `fixmate`.
- **Base quality:** `minqual: 62` in config. This is a raw ASCII threshold, not Phred: `qual + 33 > minqual`, i.e. **Phred >= 30** (`config.yaml` note; `mut_call.py:196,208,217`). The script's own default `--minQual 40` means Phred >= 8.
- **Read-end distance:** `--minDist -1`, and the pipeline does not pass it, so no end trimming. The comparison is also `dist + 33 > minDist`, the same ASCII quirk.
- **Reads with any `I` or `D` in their CIGAR are ignored entirely for mutation counting** (`mut_call.py:118`, `if ('I' not in r.cigarstring) and ('D' not in r.cigarstring)`).
  - With paired-end data, the pair is built from the other mate alone.
  - If both mates have indels, or the SE read does, no row is written.
- Soft clips are excluded (`get_aligned_pairs(matches_only=True)`).
- **Mate overlap ("dovetail"):** asymmetric rules keep a mutation from R2 even at lower quality if R2's quality clears the threshold. The code calls it "a hack to simulate TimeLapse.R behaviour, but does not necessarily mean that it is a correct dovetail mutations handling" (`mut_call.py:140-148`).
- **SNPs** (`call_snps.sh`):
  - `bcftools mpileup -f genome -b <-s4U control BAMs> | bcftools call -mv` with bcftools defaults: min BQ 13, MAPQ 0, `--max-depth 250` per file.
  - Every alt allele at a position masks mutations of **any type** at that position (the key is `chrom:pos`).
  - No -s4U controls means no SNP masking.
  - SNP positions are excluded from mutation counts **but still counted in `nT`**; the script's own comment says this should be fixed (`mut_call.py:193`).
- **Strandedness:** `reverse` maps to `R`; `yes`/`forward` **and `no` (unstranded)** map to `F` (`rules/common.smk:701-706`).

#### (e) Expressible in alnbase? **Yes.**
- Queries: `read="G" refr="A"` and `read="N" refr="A"` under `directional` for dUTP; `read="T" refr="C"` for s6G.
- The indel exclusion is upstream (`samtools view -e 'cigar !~ "[ID]"'`, or DuckDB on `-f cigar`).
- The Phred-30 threshold is `qual >= 30`. The any-alt SNP mask is a DuckDB anti-join.
- The dovetail hack can be approximated only by a DuckDB self-join. It is not worth reproducing bit-for-bit; report the delta instead.
- No query extension needed.

#### (f) Test data
- As 1.2(f): the fastq2EZbakR/bam2bakR `.test` BAMs (about 8 MB total BAMs + 47 MB chr21 FASTA + GTF) and matching FASTQs `.test/data/WT{1,2,ctl}/*_R{1,2}.fastq.gz` (0.4-0.7 MB each).
- Also `.test/data/ThreePseq_{1,2}/Muhar_3pseq_K562_*.fastq` (1.6 MB each, SLAM-seq 3'-end, SE forward).
- Simulation without reads: `EZbakR::EZSimulate()` and `bakR::Simulate_bakRData()` produce cB directly, which suits the model, not alnbase. For read-level truth, use SLAMdunk `splash` (1.1f).

#### (g) Correctness issues and contestable choices
1. **Indel-containing reads are silently dropped.** A homozygous germline indel in a 3'UTR (common in cancer lines) removes every read spanning it from `counts.csv`. The gene's read count `n` and its conversion evidence then come from a positionally biased subset. With PE, one mate's indel halves the Ts inspected for that pair while `n` still counts the pair, so the per-read `nT` distribution shifts and the mixture-model NTR moves. alnbase's per-base output makes the dropped fraction measurable per gene.
2. **Unstranded libraries are treated as forward.** For an unstranded library, R1 is antisense for about half the fragments. Their genomic T>C is complemented and recorded as `AG`, not `TC`. If only `TC` is modelled (the default `mut_tracks: "TC"`), about half the labelled reads look unlabelled, so fraction-new is underestimated roughly by that share and half-lives are overestimated. This is inferred from reading the code and should be confirmed on the test data run as `strandedness: no`.
3. **The ASCII quality quirk** means a user who sets `minqual: 30` believing it to be Phred actually gets no filtering (Phred >= -3).
4. **SNP positions count in `nT` but not in mutations**, a small systematic deflation of the per-read mutation rate at het/hom SNP sites. **Mpileup max-depth 250** means SNP calls at deep sites rest on 250 reads, which is usually still fine.
5. **The dovetail rule is asymmetric between R1 and R2**, a self-described hack. Mutation calls in overlaps depend on which mate is R2.
6. Separate from the counter: 4sU-induced **read dropout and alignment bias** (reads carrying many T>C fail to align) biases NTR downward for highly labelled RNAs. EZbakR/grandR now model this ([Berg et al. 2024, NAR "Correcting 4sU induced quantification bias"](https://academic.oup.com/nar/article/52/7/e35/7612100)). It happens at alignment, before any tool here, but alnbase per-read k distributions from conversion-aware and naive aligners expose it.

---

### 01.4 PAR-CLIP (4SU T>C, 6SG G>A), analysed with PARalyzer and wavClusteR

#### (a) Signal
- 4SU-labelled RNA is UV-365 crosslinked to protein. The crosslinked 4SU reads as C, giving T>C at the crosslink site; 6SG gives G>A ([Hafner et al. 2010, Cell](https://www.cell.com/fulltext/S0092-8674(10)00245-X)).
- Libraries are small-RNA style: **single-end, R1 sense to the RNA**, short reads of about 20-50 nt, often collapsed.
- Context: any T. Clusters of T>C mark binding sites; motif context (e.g. the PUM2 UGUA-UA motif) is analysed downstream.
- Background T>C (SNPs, errors) is modelled from non-T>C mismatches (wavClusteR/BMix) or matched RNA-seq.

#### (b) Comparators
- **PARalyzer 1.5** ([Corcoran et al. 2011, Genome Biol](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-8-r79)): bioconda `paralyzer` 1.5, noarch Java, academic-only licence. Install is trivial, but it needs a UCSC `.2bit` genome and an INI file. Sources are not public; only `.class` files ship.
- **wavClusteR** ([Comoglio et al. 2015, BMC Bioinformatics](https://link.springer.com/article/10.1186/s12859-015-0470-y)): Bioconductor 2.46.0 (release 3.23). No bioconda package found (`bioconductor-wavclusterr` is absent). Install: R/Bioconductor, medium.
- Others:
  - **PCLIPtools** (2025/26, NAR; [PMC12667781](https://pmc.ncbi.nlm.nih.gov/articles/PMC12667781/), Bash + samtools; MAPQ 55, BQ 20, drops indel events)
  - **BMix** ([arXiv 1504.01619](https://arxiv.org/pdf/1504.01619))
  - **omniCLIP** (Drewe-Boss et al. 2022, Genome Biol; heavy Python install)
  - **PARpipe** (Ohler lab wrapper around PARalyzer)

#### (c) Output formats
- **PARalyzer clusters CSV** (README): `Chromosome, Strand, ClusterStart, ClusterEnd, ClusterID, ClusterSequence, ReadCount, ModeLocation, ModeScore, ConversionLocationCount, ConversionEventCount, NonConversionEventCount, FilterType`.
- **Groups CSV:** `Chromosome, Strand, GroupStart, GroupEnd, GroupID, ReadCount, FilterType`.
- **Distributions file:** 4 lines per group (Signal KDE, Background KDE, Conversion Percent, ReadCount; -1 where no conversion is possible), always in genomic left-to-right order even for `-` groups.
- **PARalyzer coordinates:** POS is parsed straight from the SAM (1-based), with end = start + len(SEQ) - 1 (bytecode, `PARCLIPsamParser.parseFile`). Clusters are therefore 1-based closed; confirm on test output.
- **wavClusteR:**
  - `getAllSub()` returns a GRanges count table with `substitutions` (e.g. `"TC"`), `coverage` and `count`, 1-based.
  - `getHighConfSub()` adds `rsf` = count/coverage.
  - `filterClusters()` gives clusters with `Ntransitions, MeanCov, NTotTransitions, WaveletPeak, ClusterSequence, SumLogOdds, RelLogOdds`.
  - `exportHighConfSub`/`exportClusters` write a UCSC BED with `start(GR)` **(1-based) as the BED start** (`R/export.R:21-47`). That is an off-by-one: a single-base site p comes out as `p p`, a zero-length interval in BED semantics. The strand colour vector also assumes GR is sorted `+` then `-`.

#### (d) Implicit filters and defaults
- **PARalyzer** (`Default_PARalyzer_Parameters.ini` shipped values, with README defaults in parentheses):
  - BANDWIDTH=3; CONVERSION=T>C
  - MINIMUM_READ_COUNT_PER_GROUP=5 (10); MINIMUM_READ_COUNT_PER_CLUSTER=5 (1)
  - MINIMUM_READ_COUNT_FOR_KDE=5 (1); MINIMUM_CLUSTER_SIZE=8 (1)
  - MINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER=1; MINIMUM_CONVERSION_COUNT_FOR_CLUSTER=1
  - MINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION=5 (1)
  - MINIMUM_READ_LENGTH=13; **MAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES=0** (5); EXTEND_BY_READ
  - No base-quality use. Recommended aligner: `bowtie -v 2 -m 10 --best --strata`, i.e. ungapped.
  - **SAM parser (from bytecode):**
    - Strand is `-` only if FLAG equals exactly 16 or exactly 272; any other FLAG, including duplicate-marked 1040 or paired reverse, is treated as `+`.
    - The CIGAR is never parsed: read length = SEQ length, and MD offsets index SEQ directly. Offsets are held in a `byte`, so reads over 127 nt overflow.
    - Consecutive lines with the same QNAME are multimappers: the alignment with the fewest mismatches is kept, and ties become non-unique and are discarded.
    - `=COLLAPSED` takes the copy number from `name-count`.
- **wavClusteR:**
  - `readSortedBam` keeps all mapped reads, with `simpleCigar=FALSE`, i.e. any CIGAR.
  - `getAllSub(minCov=20)`: substitutions come from MD + SEQ; minus-strand substitutions are complemented; coverage = strand-specific read coverage over `IRanges(pos, width=qwidth)`.
  - `getHighConfSub(substitution="TC")`, keeping sites within the RSF support from `getExpInterval`.
  - `filterClusters(minWidth=12)`.
  - No BQ, MAPQ or duplicate filters.

#### (e) Expressible in alnbase? **Yes.**
- Queries: `read="C" refr="T"`, `read="N" refr="T"`, plus an "any other mismatch" pattern `read="/" refr="N"` so DuckDB can count non-conversion mismatches per read (PARalyzer's `MAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES`, BMix/wavClusteR's background).
- 6SG: `read="A" refr="G"`.
- Clustering and KDE are downstream and out of scope. The validation target is the per-site count table (`substitutions/coverage/count`) and per-read conversion or mismatch counts.

#### (f) Test data
- wavClusteR ships `inst/extdata/example.bam` (61 kB + index). This is PUM2 PAR-CLIP from the vignette; STAR/Bowtie provenance is not stated.
- Hafner 2010, GEO [GSE21578](https://www.omicsdi.org/dataset/geo/GSE21578): PUM2 PAR-CLIP SRR048967/SRR048968, about 10 M reads of about 32 nt, SE.
- MATR3 PAR-CLIP GSE262647 (SRR28479449-51, 35 M reads each), DHX36 SRR6191036, HNRNPK SRR23199821. These are PCLIPtools' benchmark set.
- Ground truth by simulation is easy: inject T>C into reads at known sites. A useful stress test is **STAR-aligned reads with soft clips and small indels**, which modern pipelines produce and PARalyzer/wavClusteR assume away.

#### (g) Correctness issues and contestable choices
1. **wavClusteR mishandles any non-trivial CIGAR** (verified in R on its regex from `processMD.R`):
   - MD `10^AC5T3` tokenises to `[10][][A][C][5][T][3]`. The empty token becomes `NA`, all following positions are `NA`, and the **deleted reference bases A and C are reported as substitutions**.
   - The read base is taken from `qseq` at the MD offset, which ignores soft clips and insertions, so a soft-clipped read yields the wrong read base.
   - Genomic position = `start + offset - 1` ignores `N` gaps, so spliced reads are misplaced.
   - Coverage uses `qwidth` (includes soft clips) from POS.
   - This matches PCLIPtools' report that wavClusteR "occasionally fails with certain BAM files generated by STAR aligners". **On STAR BAMs it can create spurious "TC" sites next to deletions, or shift sites.** For an RBP whose motif sits near homopolymers, where indels are common, that can produce false binding-site calls. alnbase CIGAR-aware per-site counts give the corrected table.
2. **PARalyzer's exact FLAG equality:**
   - Duplicate-marked, paired, or supplementary reverse reads are assigned to the + strand with uncomplemented bases, putting their A>G on the wrong strand.
   - No CIGAR parsing: soft clips shift conversion positions, and MD `^` is handled only by nudging start/end.
   - `byte` offsets overflow for reads over 127 nt (modern 150 nt SE runs).
   - The shipped INI sets **0 non-conversion mismatches** (README default 5). Every read with one sequencing error or a SNP is discarded; for a SNP-dense 3'UTR this removes whole clusters.
   - PCLIPtools reports that PARalyzer does not enforce its minimum-read-depth settings: about 50% of PARalyzer-unique clusters have < 2 reads despite a minimum of 3 ([PMC12667781](https://pmc.ncbi.nlm.nih.gov/articles/PMC12667781/)).
3. **wavClusteR BED export off-by-one** (c).
4. **Neither tool uses base quality or masks SNPs directly.** wavClusteR relies on the RSF mixture model, and the FDR option uses RNA-seq. A het T/C SNP with RSF ~ 0.5 falls inside the typical high-confidence RSF support, so it becomes a "binding site" unless RNA-seq is supplied.

---

### Ranking inputs

- **SLAM-seq / SLAMdunk**
  - Popularity: high (Lexogen kit, about 1k citations); SLAMdunk is the kit-endorsed pipeline.
  - Comparator cleanliness: medium. Easy conda install, but it is tied to NGM `MP` tags, SE only, and forward-stranded only.
  - Ground truth: excellent (built-in read-level simulator `splash`, tiny test data, nf-core chr8 subset).
  - Discrepancy likelihood: high (het SNP not masked at VAF 0.8; Q27 numerator-only; span coverage; silent zero on non-NGM BAMs; PE strand mishandling).
  - Expressibility: **Y**.
- **GRAND-SLAM**
  - Popularity: high (SLAM, TLS, TUC and scSLAM; grandR ecosystem).
  - Comparator cleanliness: medium-low (Maven build, genome index, CIT conversion, model-level outputs).
  - Ground truth: good via `binom.tsv` (n, d histograms) and the test BAMs.
  - Discrepancy likelihood: medium-high (no BQ, any-mismatch SNP at 30% masking true hotspots, double-hit-only overlap, `.snpdata` header swap).
  - Expressibility: **Y**.
- **TimeLapse / TUC via bam2bakR / fastq2EZbakR**
  - Popularity: medium, rising (EZbakR).
  - Comparator cleanliness: medium (Snakemake, but per-read `counts.csv` is trivially diffable against alnbase + DuckDB, and the pipeline itself uses DuckDB/parquet).
  - Ground truth: small bundled PE reverse-stranded chr21 BAMs; PRJNA726397 for three chemistries.
  - Discrepancy likelihood: **high** (indel reads dropped, unstranded treated as forward, ASCII quality quirk, SNP `nT` inconsistency, asymmetric dovetail).
  - Expressibility: **Y**. It exercises reverse-stranded orientation under `directional`.
- **PAR-CLIP / PARalyzer, wavClusteR**
  - Popularity: medium (PAR-CLIP is older; eCLIP dominates today).
  - Comparator cleanliness: low (closed-source Java; aging R package; clustering outputs).
  - Ground truth: easy simulation; small example BAM.
  - Discrepancy likelihood: **very high on STAR BAMs** (CIGAR-naive MD parsing verified; FLAG equality; byte overflow; BED off-by-one).
  - Expressibility: **Y**. The per-site table is the comparable layer.
- **Candidate extensions from this family.** None is required.
  - (i) Optional RNA library modes (`--library fr-secondstrand|fr-firststrand`) so one query file is kit-agnostic, walking in RNA-sense orientation. `directional` plus complemented patterns already covers stranded libraries, and unstranded data needs gene-strand joins downstream in any case.
  - (ii) Optional non-query output kind: a per-read count tag (e.g. `TC:i`) for NGM/SLAMdunk-compatible BAM export.


---

## Family 02 — CLIP crosslink signatures and A-to-I RNA editing

Scope: (1) HITS-CLIP crosslink-induced mutation sites (CIMS), (2) iCLIP/eCLIP crosslink-induced truncation sites (CITS / "read start − 1"), (3) A-to-I RNA editing (REDItools2/3, JACUSA2, SPRINT, RNAEditingIndexer/AEI). PAR-CLIP T>C belongs to the metabolic-labeling/PAR-CLIP section and is only mentioned where CTK overlaps.

Source code read (shallow clones in `clones/`): `ctk` (commit 4ff5243, 2024-03), `htseq-clip` (5809d44, 2025-07, v2.19.0b0), `PureCLIP` (2021-08, v1.3.1), `iCount` (2020-03), `skipper` (8b6405b, 2026-09), `clippy`, `REDItools2` (1e9d396, 2025-07), `REDItools3` (8cfebea, 2026-07, v3.7), `JACUSA2` (6f4afc5, 2026-06, pom 2.1.15), `SPRINT` (a73c23b, 2025-01), `RNAEditingIndexer` (80ab14a, 2026-07).

Bioconda (queried via api.anaconda.org, 2026-09-16): jacusa2 2.1.17, reditools3 3.7, pureclip 1.3.1, htseq-clip 2.19.0b0, icount 2.0.0, clippy 1.5.0, umi_tools 1.1.6. **Not on bioconda:** CTK (own channel `chaolinzhanglab/ctk` 1.1.4, plus czplib Perl lib), REDItools2 (git + Python 2.7-era install), SPRINT (Python 2 binaries/Docker), RNAEditingIndexer (Docker `levanonlab/rna-editing-index`, ~12 GB resources), Skipper (Snakemake + conda envs).

#### alnbase semantics that matter for this family (verified in source)

- Walk orientation under `--library directional` (`src/tags.rs` `directional_strand`, `src/alignment.rs:117`): R1 fwd = OT and R2 rev = CTOT are walked along the top strand; R1 rev = OB and R2 fwd = CTOB are walked reverse-complemented. Single-end reads count as R1. So:
  - a **forward/sense-stranded** RNA library (R1 = RNA sense; iCLIP, seCLIP R1, QuantSeq FWD, REDItools `-s 1`, JACUSA2 `FR-SECONDSTRAND`) is walked RNA-sense for both mates;
  - a **dUTP/reverse-stranded** library (R1 antisense; TruSeq stranded, REDItools `-s 2`, JACUSA2 `RF-FIRSTSTRAND`, paired-end eCLIP where R2 is sense) is walked RNA-antisense for both mates, so A>G editing must be written `read C / refr T`, and motifs written reverse-complemented;
  - **unstranded** libraries cannot be oriented by FLAG; each read is walked in an arbitrary (R1-defined) orientation.
- `off_5p`/`off_3p` count from the ends of the read *as sequenced* regardless of walk direction (`src/hits.rs` header comment), so read-end trimming is a DuckDB `WHERE`.
- Pads (`_`) are emitted in the flank past each aligned end, with reference bases readable. In the pattern the **walk order** distinguishes the two ends: `_N` can only fire at the walk start and `N_` at the walk end. Walk start = 5' end of R1 (OT/OB) but the **3' end of R2** (CTOT/CTOB). There is no code meaning "5' end of this read as sequenced".
- With soft clips skipped (default) the pad column is the reference base adjacent to the first/last *aligned* base, i.e. the same "aligned start − 1" coordinate PureCLIP, iCount and htseq-clip use.
- Bases tags refuse anchors on gap/intron/pad (`src/bam_out.rs:18-30`); parquet output has no such refusal documented, but anchoring on a pad has "no position" in tags. Safest parquet idiom: anchor on the read base next to the pad and capture the pad with `^`.

---

### 1. HITS-CLIP crosslink-induced mutation sites (CIMS)

#### a. Signal
UV-crosslinked amino-acid adducts make reverse transcriptase skip or misincorporate at the crosslinked nucleotide. In standard HITS-CLIP the dominant signature is a **single-nucleotide deletion** at the crosslink site (deletions in ~8–20% of Nova/Ago mRNA tags; Zhang & Darnell 2011, [Nat Biotechnol 29:607](https://www.nature.com/articles/nbt.1873)). Substitutions are mostly sequencing error except for protein-specific cases and PAR-CLIP (4SU T>C), and insertions are rare. Deletions prefer U-rich sequence (TTT) for Nova, while iCLIP truncations find YCAY ([Sugimoto 2012](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-8-r67)). The reads are single-end and sense to the RNA (standard/BrdU-CLIP), so the walk under `directional` is RNA-sense. For eCLIP the CTK tutorial uses R2 only, which is sense ([CTK eCLIP tutorial](https://zhanglab.c2b2.columbia.edu/index.php/ECLIP_data_analysis_using_CTK)).

#### b. Comparator
- **CTK (CLIP Tool Kit) 1.1.3/1.1.4** (Shah et al. 2017 Bioinformatics). Perl, needs the `czplib` Perl library and Math::CDF. Installed from its own conda channel (`conda install -c chaolinzhanglab ctk`), **not bioconda**. Medium install complexity: Python 2/3 wrapper issues are open (#4). Pipeline: `parseAlignment.pl` → `tag2collapse.pl` → `getMutationType.pl` → `CIMS.pl`.
- Alternatives: the older standalone CIMS package ([SourceForge ngs-cims](https://sourceforge.net/projects/ngs-cims/)) and `htseq-clip extract -s d|i` (deletion/insertion sites; no statistics).

#### c. Output formats (from source)
- `parseAlignment.pl --mutation-file` (`parseAlignment.pl:115-140`): BED6+5, 0-based half-open (czplib `Bed.pm` stores chromEnd inclusive and `bedToLine` adds 1):
  `chrom, chromStart, chromEnd, name(read), score(=offset of the mutation in the read, SEQ/+ orientation, soft clips removed), strand(read strand), score(repeated), refBase(+ strand), type('>' sub | '-' del | '+' ins), altBase(+ strand; '.' for deletions), matchStart(1 if the read's 5' end is not clipped)`.
  Insertions are written with chromStart = chromEnd (width 1 on output) at the reference base *after* the insertion.
- `getMutationType.pl -t del|ins|sub [--subt t2c]` writes BED6 (cols 0–5). For substitutions on the `-` strand, from/to are reverse-complemented before matching `--subt` (`getMutationType.pl:112-121`).
- `CIMS.pl` output (`CIMS.pl`, header line): `#chrom, chromStart, chromEnd, name(+"[k=K][m=M]"), score, strand, tagNumber(k), mutationFreq(m), FDR, count(>=m,k)`. Here k is the number of unique tags overlapping the site (strand-specific, `tag2profile.pl -ss`) and m is the number of tags carrying the mutation. Coordinates are 0-based half-open.

#### d. Implicit filters and defaults
- `parseAlignment.pl` (`parseAlignment.pl:15-21`): `--map-qual 0` by default (tutorial uses `--map-qual 1` "to keep only unique mappings" with BWA), `--min-len 0` (tutorial 18), `--indel-to-end 5`. There is **no base-quality filter**; the tutorial instead filters reads by mean quality ≥20 over the first 30–39 cycles with `fastq_filter.pl -f mean:0-29:20`. The parser requires an MD tag (run `samtools calmd`). It croaks on `=`/`X` CIGAR ops and on ambiguity codes (issue [#7](https://github.com/chaolinzhanglab/ctk/issues/7)). Soft/hard clips are stripped before any offset is computed.
- Mutation "flag==2" rule (`parseMD`): a substitution counts only if the MD tag has a match run ≥1 somewhere on each side. An indel needs a match run ≥5 (`--indel-to-end`) **somewhere** on each side, and the run need not be adjacent: the flags accumulate over tokens. For insertions the check is the read offset ≥5 from both ends. In practice mismatches at the terminal aligned base are dropped.
- `tag2collapse.pl --random-barcode -EM 30 --seq-error-model alignment -weight --weight-in-name --keep-max-score --keep-tag-name` collapses PCR duplicates by 5' position plus a UMI EM model. Mutations are then restricted to the surviving representative tags, so **which duplicate is kept decides which mutations survive**.
- `CIMS.pl`: `-w 1` means only mutations of exactly size 1 are analysed; multi-nt deletions are silently ignored unless `-w 2`. `-n 5` permutation iterations (tutorial 10), `srand(0)`, `--FDR 1` and `-mkr 0` (all sites reported; the tutorial filters `FDR<=0.001`). The FDR depends on input tag order (docs Note 3).
- No SNP masking. Multi-mappers are removed only through MAPQ.

#### e. alnbase expressibility
**Yes, today.** A deletion: `read="N.N" refr="NNN"`, anchor col 1 (a gap anchor is fine for parquet, not for a `bases` tag), which reproduces `-w 1` (exactly one deleted base flanked by aligned bases). Size-2: `N..N`. Substitution subtype: `read="C" refr="T"` with walk sense orientation, or the relational `/` for any substitution. Insertion CIMS needs `--insertions emit` and `read="N" refr="."`. The per-site k/m counts and CTK's "indel-to-end" rule are DuckDB (`off_5p>=5 AND off_3p>=5`). CTK's non-adjacent MD-run quirk is not worth reproducing, only documenting. The FDR permutation is downstream statistics.
No syntax extension is needed. The only friction is eCLIP R2 in a paired BAM (walk is antisense; see §2).

#### f. Test data
- CTK tutorial (mouse brain Rbfox HITS-CLIP + BrdU-CLIP): SRA [SRP035321](http://www.ncbi.nlm.nih.gov/sra/?term=SRP035321). The tutorial offers downloads of the unique-tag files (26 MB), CIMS outputs (1.6 MB) and CITS outputs (641 KB) as reference results ([tutorial](https://zhanglab.c2b2.columbia.edu/index.php/Standard/BrdU-CLIP_data_analysis_using_CTK)). The expected answer is deletion CIMS enriched at G2/G5 of UGCAUG, a built-in biological sanity check.
- ENCODE RBFOX2 HepG2 eCLIP R2 FASTQs [ENCFF647KDW](https://www.encodeproject.org/files/ENCFF647KDW/) (351 MB) and [ENCFF289OFA](https://www.encodeproject.org/files/ENCFF289OFA/) (354 MB), experiment ENCSR987FTF.
- `htseq-clip/tests/testBamCLIP/test02.bam` (8.9 KB) and `test03.bam` (73 KB) ship with expected deletion/insertion BEDs for `extract -s d/i`: a tiny golden test.

#### g. Known issues and contestable choices
- **Aligner indel normalisation.** BWA/STAR left-align deletions in + reference coordinates. In a homopolymer (e.g. the TTT context Nova deletions prefer) the reported CIMS base is the leftmost, which is the 3'-most for minus-strand RNA, so the "crosslinked nucleotide" differs by strand. Both tools see the same CIGAR, but motif-position conclusions such as "which U is crosslinked" flip with strand. alnbase can show this explicitly by adding `refr` context to the pattern.
- `-w 1` drops 2-nt deletions silently. Sites supported mostly by multi-nt deletions get m=0.
- The CIMS FDR depends on tag order ([docs](https://zhanglab.c2b2.columbia.edu/index.php/Standard/BrdU-CLIP_data_analysis_using_CTK) Note 3), so pooled-replicate concatenation order changes site lists near the threshold.
- CITS removes every tag carrying a deletion before truncation analysis ("deletions are read-through"). That makes CIMS and CITS mutually exclusive by construction.
- With `--map-qual 1` on STAR BAMs, MAPQ 0/1/3 behave differently from BWA ("Other aligners might not use a positive MAPQ as an indication of unique mapping", tutorial Note 2).

---

### 2. iCLIP / eCLIP crosslink-induced truncation sites (CITS)

#### a. Signal
RT stops at the crosslinked adduct. The cDNA 3' end, which is the **first base of the sense read** (iCLIP R1, seCLIP R1, paired-end eCLIP R2), lies immediately downstream of the crosslink, so the crosslink = **reference position 1 nt upstream (5', RNA-sense) of the read's first aligned base**. There is no substitution; the signal is a read-boundary event. In iCLIP, 80–95% of cDNAs are truncated ([Sugimoto 2012](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-8-r67)). eCLIP amplifies read-through as well ([Van Nostrand 2016](https://yeolab.github.io/papers/2016/nmeth_eric_2016.pdf)). cDNA-start positions depend on cDNA length through constrained cDNA ends ([Haberman 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5240381/)).

#### b. Comparators
- **PureCLIP 1.3.1** (Krakau 2017 Genome Biol; HMM for crosslink sites; bioconda `pureclip`; C++/SeqAn, easy via conda). Its own ENCODE tutorial exists.
- **htseq-clip 2.19.0b0** (bioconda/PyPI; pure Python/pysam; easy). A counting front end for DEWSeq with no statistics, so it is the cleanest 1:1 comparator of site extraction.
- Also relevant: **CTK CITS.pl** (above); **iCount 2.0.0** `xlsites` (bioconda; UMI-aware per-site cDNA counts, BED6); **clippy 1.5.0** (peak caller consuming iCount-style crosslink BED, not BAM); **Skipper** (Yeo lab, ENCODE4 eCLIP; window counts, Snakemake).

#### c. Output formats and coordinate conventions (from source)
- **htseq-clip `extract`** (`clip/bamCLIP.py:336-404`): BED6 `chrom, start, end, name="<qname>|<query_length>", score=YB tag or 1, strand(read strand)`, 0-based half-open, one line per read. Site choices: `-s s` start, `e` end (**default e**), `m` middle, `i` insertion, `d` deletion. With `-s s` the start is `reference_start + offset` for + reads and `reference_end − offset − 1` for − reads, so **the default offset 0 is the read's first aligned base, not −1**. The docs say to use `--mate 2 --site s --offset -1` for eCLIP ([documentation.rst:139](https://github.com/EMBL-Hentze-group/htseq-clip)).
- **PureCLIP** (`docs/PureCLIPTutorial/output.rst`): BED6 `chr, start, end(start+1), state('3'), score(log posterior ratio), strand`. With `-or`, regions are merged within `-dm 8`. Read starts are counted at `beginPos` (+) or `beginPos + alignedRefLen − 1` (−), i.e. the aligned (soft-clip-excluded) 5' base (`src/parse_alignments.h:74-83`). By default sites are then assigned **1 nt upstream**; `-ctr` assigns them to the read start itself (`pureclip.cpp:86`).
- **CTK CITS.pl**: `bedExt.pl -n up -l -1 -r -1` gives the 1-nt upstream position, and `tag2peak.pl` outputs BED6 `name=CITS_<n>[gene=..][PH=..][PH0=..][P=..]`, `score=peak height`. `--gap 25` merges nearby sites (tutorial `-p 0.001 --gap 25`, no Bonferroni).
- **iCount xlsites** (`iCount/mapping/xlsites.py:469-477`): `xlink_pos = poss[0]-1` (+) or `poss[-1]+1` (−) from `get_reference_positions()`, i.e. aligned bases only. Output is BED6 `chrom, pos, pos+1, '.', score(cDNA or reads), strand`, written separately for unique (`sites_single`) and multi-mapped reads. Oddity: `xlink_pos = 1 if xlink_pos < 1` (a 0-based −1 becomes 1, not 0).
- **Skipper** (`rules/genome_windows.smk:69-76`): `bedtools bamtobed | bedtools flank -s -l 1 -r 0 | bedtools shift -p 1 -m -1` moves one base upstream and then back, so it counts the **read's first aligned base** (not −1) into windows. That is harmless for ≥100-nt windows but differs by 1 nt from every other tool.
- **JACUSA2 rt-arrest** (see §3) uses the terminal aligned base itself, not −1.

#### d. Implicit filters and defaults
- htseq-clip (`clip/clip.py:386-490`, `bamCLIP.py:268-307`): `-q/--minAlignmentQuality 10`; `--minReadLength 0`; `--maxReadLength 500`, so **longer reads are silently dropped** (the helper signature default is 100); `--maxReadInterval 10000` (reference span incl. introns); QC-fail and unmapped reads dropped; secondary reads dropped only with `--primary`; supplementary **not** dropped; duplicates dropped only with `--ignore_PCR_duplicates` and only if flagged; `-e/--mate` required (issue [#11](https://github.com/EMBL-Hentze-group/htseq-clip/issues/11): confusion for SE libraries; SE reads pass either mate setting).
- PureCLIP: **no FLAG, MAPQ or duplicate filtering in the parser**. It relies on preprocessing: STAR `--outFilterMultimapNmax 1 --alignEndsType EndToEnd`, `samtools view -f 2`/`-f 130` (R2 only), `umi_tools dedup --paired` ([preprocessing.rst](https://pureclip.readthedocs.io/)). `-ur` defaults to *all reads*: a paired BAM without `-ur 2`/`-f 130` counts R1 starts (the fragment's other end) as truncations. `-mtc 500` (intervals with >500 starts at one position are skipped for learning), `-mtc2 65000` count cap, `-dm 8`, `-bdw 50`.
- iCount xlsites (`run()` defaults): `mapq_th=0`, `multimax=50`, `mismatches=1` (UMI mismatches merged), `gap_th=4`, `ratio_th=0.1`, `group_by='start'`, `quant='cDNA'`. It requires the `NH` tag and a randomer in the read name (`:rbc:`).
- CTK CITS: tags with deletions are removed first; clusters must have >2 tags; `-p 0.01` default (tutorial 0.001); `--gap -1` default (tutorial 25).
- UMI dedup matters a lot: truncation sites are start-position identities, so position-only dedup (Picard) collapses genuine independent truncations at a crosslink. JACUSA2's manual explicitly says **do not deduplicate for arrest methods** without UMIs.

#### e. alnbase expressibility
**Expressible today with caveats; a small extension would make it clean.**
- Single-end iCLIP/seCLIP (R1 sense, walk sense): `read="_N" refr="NN"`, anchor col 1 and capture col 0 (`mark="^+"`), gives the aligned start −1 with the reference base there. It fires once per read, only at the 5' end, because `_N` cannot match at the walk end. Soft clips skipped by default match PureCLIP/iCount/htseq-clip's aligned-start convention; `--soft-clips emit` gives the unclipped-start variant as a sensitivity analysis.
- Paired-end eCLIP (R2 sense): under `directional`, R2 is CTOT/CTOB and walked antisense, so the R2 5' end is the **walk end**. Use `read="N_"`, then keep `flag & 128` (or `samtools view -f 130` upstream). `refr_pos` of the pad is still the correct genomic crosslink coordinate, but the reported `refr_base` is complemented relative to the RNA.
- Read-through or "coverage" denominators (PureCLIP binomial n, JACUSA2 through counts) need an every-base query plus DuckDB. That is legitimate but large.
- **Extension candidates:** (i) **read-end-anchored pad codes** that mean "past the 5' end as sequenced" vs "past the 3' end as sequenced" (e.g. `<` and `>`), so truncation queries do not depend on `--library` or mate number; (ii) a **`--library` RNA mode** where R2 is the sense strand (eCLIP PE; `reverse`) and one where R1 is sense (`forward`), so `refr_base`/motif context reads RNA-sense.

#### f. Test data
- htseq-clip `tests/testBamCLIP/test02.bam` (8.9 KB) and `test03.bam` (73 KB) with expected `*_SS_check.bed` (start sites), `_MS_` and `_ES_`: tiny golden tests for exact coordinate agreement.
- PureCLIP tutorial: ENCODE PUM2 K562 eCLIP processed BAM [ENCFF280ONP](https://www.encodeproject.org/files/ENCFF280ONP/) (**229 MB**, hg19), filtered with `samtools view -hb -f 130`.
- CTK eCLIP RBFOX2 HepG2 R2 FASTQs (above) and the CTK BrdU-CLIP CITS reference outputs (641 KB).
- clippy `tests/data/rbfox/HepG2_RBFOX2.xl.bed.gz` (37 MB crosslink BED) and `crosslinkcounts.bed` (14.6 MB): downstream-only, but a reference for iCount-format output.
- iCount hnRNPC iCLIP (König 2010; ArrayExpress E-MTAB-432). The example script's `icount.fri.uni-lj.si` URL may be dead (did not respond).
- Simulation: [sim_iCLIP](https://github.com/skrakau/sim_iCLIP) (PureCLIP authors' iCLIP/eCLIP simulator with ground-truth crosslink sites).

#### g. Known issues and contestable choices
- **Off-by-one conventions differ across tools.** CTK, PureCLIP, iCount and htseq-clip `-g -1` use −1. htseq-clip's *default* offset 0, Skipper's counting, and JACUSA2 arrest use the first aligned base. PureCLIP `-ctr` uses the start. A 1-nt shift moves crosslink peaks relative to motifs (e.g. Rbfox G2 vs G5 of UGCAUG), which is exactly the single-nucleotide conclusion these assays exist for.
- **Soft clipping.** Every tool uses the *aligned* start. RT adds untemplated nucleotides at cDNA 3' ends, and STAR local alignment soft-clips them (and sometimes true mismatching bases near the crosslink). The truncation position then moves by the clip length, which is why PureCLIP insists on `--alignEndsType EndToEnd`. alnbase can quantify this directly (`--soft-clips emit`, clipped-base identity at the 5' end).
- **Mate choice.** Any tool given a PE BAM without mate restriction (PureCLIP default `-ur` unset; iCount; CTK if R1 is not excluded) counts adapter-side ends as truncations.
- **Duplicates.** Position-only dedup destroys truncation signal (JACUSA2 manual). UMI dedup needs the UMI in the name (iCount `:rbc:`, umi_tools `_UMI`).
- htseq-clip silently drops reads >500 nt (merged/long reads) and keeps supplementary alignments by default.
- Constrained cDNA ends ([Haberman 2017](https://pubmed.ncbi.nlm.nih.gov/28093074/)) make cDNA-start distributions depend on read length filters (`--min-len 18` in CTK vs none in PureCLIP), which changes site rankings.

---

### 3. A-to-I RNA editing

#### a. Signal
ADAR deamination A→I is read as G, so the read shows **A>G on the RNA-sense strand**. On the reference + strand that is A>G for + transcripts and T>C for − transcripts. There is no required sequence context (a weak 5' neighbor U/A preference and 3' G preference, ADAR-dependent). Hyper-edited reads carry many A>G per read and often fail to map. Library strandedness decides whether a T>C in reads is A-to-I on a minus-strand gene or a U-to-C artifact. Common pitfalls: SNPs, RT/sequencing errors at read ends, mis-mapping near splice junctions and indels, homopolymers, paralog/Alu multi-mapping.

#### b. Comparators
- **REDItools2** (Flati et al. 2020; git 1e9d396, 2025-07). Python (README says 2.7-tested), MPI parallel version, **not on bioconda**; medium install. Nature Protocols workflow: [Lo Giudice 2020](https://www.nature.com/articles/s41596-019-0279-7).
- **REDItools3 v3.7** (PyPI `REDItools3`, bioconda `reditools3` 3.7; easy). A rewrite with **changed defaults and option semantics** (below).
- **JACUSA2 2.1.x** (Piechotta 2022 BMC Bioinf; bioconda `jacusa2` 2.1.17; single Java 17 jar; easy). `call-1` (vs reference/MD), `call-2` (DNA vs RNA or condition vs condition), `pileup`, `rt-arrest`, `lrt-arrest`.
- **SPRINT 0.1.8** (Zhang 2017 Bioinformatics; Python 2 binaries/Docker; not on bioconda; medium–hard). SNP-free cluster-based calling plus hyper-editing via A→G-masked remapping.
- **RNAEditingIndexer / Alu Editing Index (AEI)** (Roth, Levanon & Eisenberg 2019 Nat Methods; Docker `levanonlab/rna-editing-index`; needs samtools, bedtools, bamUtil, Java, ~12 GB resources; medium–hard). The dominant global-editing metric in large cohorts (GTEx, TCGA).
- BCFtools is used as a naive baseline in benchmarks ([Morales 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10527054/)).

#### c. Output formats (from source)
- **REDItools2** (`src/cineca/reditools.py:778, 1159-1170`): TSV with header
  `Region, Position, Reference, Strand, Coverage-q30, MeanQ, BaseCount[A,C,G,T], AllSubs, Frequency, gCoverage-q30, gMeanQ, gBaseCount[A,C,G,T], gAllSubs, gFrequency`.
  `Position` is **1-based**. `Strand` is **numeric**: 1 = +, 0 = −, 2 = unknown/unstranded. On a − site, `Reference` and `BaseCount` are complemented. `AllSubs` looks like `AG AT` (ref+alt for each non-ref base, sorted by count). `Frequency` = **(most frequent alt)/(most frequent alt + ref)**, not total non-ref and not A>G-specific. The `g*` columns are `-` until `annotate_with_DNA.py`. The column name says q30 regardless of `-bq`.
- **REDItools3** (`README.md`, `reditools/compiled_position.py`): the same 14 columns but `Coverage`/`gCoverage` without "-q30" and `Strand` as `+`/`-`/`*` (PR #110 "Strand symbols"). The README says Frequency is the "ratio of non-reference bases to reference bases", but the code docstring (`compiled_position.py:214`) says most-frequent-variant/(ref+variant), so docs and code disagree. `index` subcommand computes AEI-style indices.
- **JACUSA2** ([manual](https://github.com/dieterich-lab/JACUSA2/blob/master/manual/manual.pdf), `manual.tex:325-400`): BED6-extended, `##` header lines, **0-based [start,end)**. Columns: `contig, start, end, name, score, strand, bases11[,bases12..,bases21..], info, filter, ref`. `bases` = `A,C,G,T` counts, **complemented on − strand**. `score` is the Dirichlet-multinomial LRT statistic (`call-*`) or total coverage (`pileup`). `strand` is `.` for UNSTRANDED. `filter` holds a `;`-separated list of feature-filter IDs or `*`. `ref` is complemented on −. `-f V` gives VCF (unstranded). **rt-arrest** replaces each `basesIJ` with `arrestIJ, throughIJ` and uses a `pvalue` score with `arrest_score=` in `info`. **lrt-arrest** writes multi-line sites with an `arrest_pos` column. `-D/-I` add `deletion_pvalue`, `deletion_score`, `deletions11=del,total` to `info`. `-B A2G` adds read-stratified lines.
- **SPRINT** (`sprint/sprint_from_bam.py:303`, `sprint_main.py:640-652`): `SPRINT_identified_{regular,hyper,all}.res`, header `#Chrom, Start(0base), End(1base), Type(e.g. AG|TC|CT|GA on + reference), Supporting_reads, Strand, AD:DP`. A-to-I = `AG +` or `TC -`. The intermediate `tmp/*.zz` has per-read columns `Chr, SAM_Flag, MapQ, Loc(blocks), SNV(ref+alt:pos1;…), BaseQ, Read-loc, Seq, Read-name(_1/_2), Fragment-loc`.
- **AEI** (`TestResources/CompareTo/EditingIndex.csv`): CSV, one row per sample. `StrandDecidingMethod, Group, Sample, SamplePath, A2CEditingIndex, A2GEditingIndex, A2TEditingIndex, C2AEditingIndex, C2GEditingIndex, C2TEditingIndex, …` plus per-type canonical/mismatch counts split by exonic/intergenic/intronic and +/−/unknown strand. The index is in **percent** (e.g. A2G 2.39 vs C2T 0.03 in the test output).

#### d. Implicit filters and defaults
**REDItools2** (`reditools.py:1217-1245, 441-570`):
- MAPQ `-q 20`, base quality `-bq 30`, `-mrl 30` minimum read length, `-os 5` homopolymer span (only if `-m/-c` given), `-ss 4` splice span (only with `-sf` annotation), `-mbp 0`/`-Mbp 0` read-end trim, `-l 1` min coverage, `-men 1`, `-me 0`, `-Men 100`, `-S` strict off (prints *every covered position*), `-s 0` unstranded, `-T 1` (maxValue strand inference), `-Tv 0.7`, `-C` off.
- Read filter: drops QC-fail, secondary, supplementary, **duplicate-flagged**, and any read with an `SA` tag. **Paired reads must have FLAG exactly 99/147/83/163**, so improper pairs and pairs with an unmapped mate vanish.
- Column filter bug-like line: `if column["mean_quality"] < MIN_QUALITY` compares mean *base* quality to the *MAPQ* threshold (`reditools.py:543`). It is harmless at defaults and bites when `-bq < -q`.
- `-mbp/-Mbp` are measured on `alignment_index` in **SEQ (BAM) order including soft clips**, not in read orientation. On reverse reads "first X bases" trims the read's 3' end, and `length − pos < Mbp` trims **Mbp−1** bases from the right (apparent off-by-one; verify with a synthetic read).
- Strand (`-s 1|2`): each read gets +/− from FLAG, then the whole column is assigned the majority strand (ties go to `+`, `vstand()` at `reditools.py:1200-1211`). Without `-C`, **reads from the opposite strand are still counted, complemented** into the column.
- Overlapping mates are counted twice (no read-name dedup), and `N` bases are skipped.

**REDItools3 v3.7** (`reditools/tools/analyze/parse_args/parse_args.py`, `alignment_file.py:30,76`, `compiled_reads.py:138-215`):
- `-q 20`, `-bq 30`, `-mrl 30`, `-mbp 0`/`-Mpb 0`, `-l 1`. **`-me 1` by default** (only sites with ≥1 variant, equivalent to legacy `--strict`) and **`-men 0`**, `-Men 4`, `-s 0`, `-T 0.7`.
- **`-C` changed meaning**: in REDItools2 it restricted counting to the inferred strand; in v3 it only complements the reported bases on − sites.
- Keeps only FLAG ∈ {0, 16, 83, 99, 147, 163} and no `SA` tag, so duplicates, secondary, supplementary and improper pairs go (implicitly via FLAG).
- Stranded mode: forward flags {0, 99, 147} for `-s 1`, {16, 83, 163} for `-s 2`.
- Right trim `read_pos > query_length − Mbp` also keeps index `L−Mbp`, i.e. trims Mbp−1 bases (same apparent off-by-one), again in SEQ order.
- Regressions and bugs: [#102](https://github.com/BioinfoUNIBA/REDItools3/issues/102) IndexError on heavily soft-clipped reads in v3.7 (e.g. after `bamutil clipOverlap`); [#86](https://github.com/BioinfoUNIBA/REDItools3) stranded-SE fix.

**JACUSA2** (CLI tables `manual/CLI/*.tex`, `AbstractBCCfilterFactory.java:52-53`, `HomopolymerFilterFactory.java:46`):
- `-q 20` min base quality. `-m` min MAPQ is **−1 (off) in `call-1`** but **20 in `call-2`/`pileup`/`rt-arrest`**. `-c 5` min coverage, `-F 0` (so **duplicates are not filtered unless `-F 1024`**), `-P UNSTRANDED` default (arrest methods require stranded), `-filterNH`/`-filterNM` off.
- Feature filters are **opt-in** via `-a`: `D` (combined I+B+S), `B` read start/end, `I` indel, `S` splice site (distance ≤6 affecting ≥50% of reads by default), `Y` homopolymer ≥7, `M` max alleles, `H` homozygous in condition, `E` exclude-site file.
- Arrest location (`lib/data/storage/arrest/*LocInterpreter.java`): the arrest position is the **first/last aligned base of the informative read** (R2 5' end for `RF-FIRSTSTRAND` PE, R1 5' end for `FR-SECONDSTRAND`), not −1. The other mate is placed via mate start or TLEN, and improper pairs have no arrest site. For SE `RF-FIRSTSTRAND`, the arrest is the read's **3'** end.

**SPRINT** (`sprint/tools_zf/zz2snv.py`, `sam2zz.py`, `get_baseq_cutoff.py`, `dedup.py`):
- MAPQ **20 ≤ MAPQ < 200**, so **STAR unique reads (MAPQ 255) are all excluded** unless remapped with `changeSAMmapQ.py`, a classic "no output" trap.
- Base quality ≥25 (ASCII 58/89 auto-detected), and the **mean quality of all mismatches in the read** must be ≥25.
- "Fragment-loc" > 5: the mismatch must be ≥5 bases from **both edges of its CIGAR M block**, so any indel, intron, soft clip or read end excludes it.
- Dedup of reads with identical chrom + blocks + SEQ.
- Clustering: SNV duplets within `-cd 200` bp. Cluster size thresholds: Alu AD≥1: 3, Alu AD≥2: 2, non-Alu repeat: 5, non-repeat: 7, no-repeat-file: 5.
- `sprint_from_bam` renames all reads to `read1`, so the strand comes from FLAG alone (reverse-stranded data gets flipped labels unless handled).
- CIGAR parser handles only I/D/M/S/P/N. The number list includes H/=/X lengths, so **hard-clipped or `=`/`X` CIGARs misalign the op/length lists** (`sam2zz.py:doCG`); inferred from code, verify.

**AEI** (`Configs/DefaultsConfig.ini`, `FullConfig*.ini`):
- `bamUtil trimBam` **5 bases from each read end** (`bases_to_trim = 5`), then `samtools mpileup -B -d 1000000` with **samtools defaults: MAPQ ≥0 (multi-mappers included), BQ ≥13, `--ff UNMAP,SECONDARY,QCFAIL,DUP`, overlapping-mate quality adjustment on**. Mismatches count only at **BQ ≥30** (`PileupToCount_quality_threshold`).
- Common-SNP BED masking. Strand per Alu from RefSeq, then from mismatch excess (`StrandDecidingMethod`). The README says "alignment should be unique" but does not enforce it (issue [#27](https://github.com/a2iEditing/RNAEditingIndexer/issues/27)).

#### e. alnbase expressibility
**Yes for all per-read, per-base signals; convenience gaps for strandedness and proximity filters.**
- Site-level A-to-I: sense walk `read="G" refr="A"`; dUTP-stranded under `directional`, `read="C" refr="T"`. Coverage/BaseCount needs an every-base query (`read="N" refr="A"` or `refr="{AT}"`) grouped in DuckDB by `refr_pos` plus the strand derived from FLAG.
- Read-end trims (`-mbp/-Mbp`, AEI 5 bases, JACUSA2 `B`), BQ/MAPQ thresholds, duplicate/FLAG/SA-tag rules, SNP/annotation masks, splice-site annotation spans, cluster rules (SPRINT duplets), per-read mismatch-quality means and hyper-edited-read counts are all DuckDB/samtools and **not** extension needs.
- **SPRINT's "Fragment-loc > 5" is exactly expressible as one pattern:** with `--insertions emit`, `read="NNNNNGNNNNN" refr="NNNNNANNNNN"` anchor col 5. `N` rejects gap, pad and junction columns, so a deletion, intron, soft-clip boundary or read end within 5 bases kills the match. Pattern context doing real work, which makes it a good demo.
- **REDItools homopolymer (`-os 5`) and JACUSA2 `Y` (≥7), `I`/`S` distance ≤6:** expressible today only as an `or` over every shifted placement (e.g. 5 offsets × 4 bases for homopolymers; 12 shifted gap/junction patterns for ±6). Verbose but mechanical.
- SPRINT hyper-editing (A→G masked remapping) is **not expressible** (needs realignment). Out of scope.
- **Extension candidates:** (i) a **proximity/any-within predicate** (e.g. `near(gap|junction, 6)`) or a **run quantifier** on the reference row (`A{5,}`) for homopolymer and indel/splice proximity filters; (ii) **RNA library modes** (`--library forward|reverse|unstranded`) so REDItools `-s 1/-s 2/-s 0` and JACUSA2 `FR-SECONDSTRAND/RF-FIRSTSTRAND/UNSTRANDED` map 1:1 and `refr_base` is RNA-sense (unstranded = reference orientation for everyone, matching REDItools `-s 0` and AEI's pileup).

#### f. Test data
- **REDItools2 bundled** `test/SRR2135332.bam` (**1.4 MB**, hg19 chr21; `prepare_test.sh` downloads chr21.fa). The quickest REDItools2/3 vs alnbase check.
- **AEI bundled** `TestResources/BAMs/*/SRR59622xx_sampled_with_0.1…AluChr1Only.bam` (4 BAMs, **71 MB** total, stranded PE, hg38 chr1 Alu only) plus annotation BEDs and **expected `CompareTo/EditingIndex.csv`**: a golden end-to-end number (A2G index ≈1.6–2.4%).
- **Ground truth by genetics:** GEO [GSE99249](https://pmc.ncbi.nlm.nih.gov/articles/PMC10527054/), HEK293T WT (SRR5564274-6) vs ADAR1-KO (SRR5564268, -72, -73), paired-end. Used by the 2023 benchmark; editing should drop ~64–74% in KO.
- **JACUSA2 simulated**: [gDNA vs cDNA](https://data.dieterichlab.org/s/gDNA_VS_cDNA) and [cDNA vs cDNA](https://data.dieterichlab.org/s/cDNA_VS_cDNA) (human chr1, ART/Flux simulator, with `variants.txt` truth). Arrest events: [Zhou 2018 CMC-Ψ rRNA BAMs](https://data.dieterichlab.org/s/arrest_events).
- REDItools3 `test/sam_gen.py` + `test/aligner.py` generate synthetic SAM/FASTA, which can be reused to build adversarial reads (soft clips at ends, overlapping mates, strand ties).

#### g. Known issues and contestable choices
- **Strand-label trap.** REDItools2 `-s 1` = "secondstrand" (R1 sense) and `-s 2` = "firststrand" (dUTP); JACUSA2's names are the other way round in spirit (`FR-SECONDSTRAND` = R1 sense). [JACUSA2 #78](https://github.com/dieterich-lab/JACUSA2/issues/78): a user got reverse-complement results (A>G ↔ T>C) running REDItools2 `-s 2` and JACUSA2 `FR-SECONDSTRAND` on the same BAM. That flips "A-to-I" into "U-to-C", and in antisense-overlapping loci it flips gene assignment. SPRINT issues #19/#32 are the same confusion.
- **REDItools2 vs 3 defaults.** v3 `-me 1` (strict) and `-men 0` vs v2 non-strict / `-men 1`; `-C` semantics changed; FLAG whitelist. [#55](https://github.com/BioinfoUNIBA/REDItools3/issues/55) reports output of 200–300 MB (v3) vs 20–30 GB (v2). Strict-by-default is the likely explanation; the maintainers did not answer.
- **Frequency definition.** REDItools reports the most frequent alt over (ref + that alt). At a SNP-plus-editing or error-rich site the "Frequency" may be for A>C, not A>G, which biases editing-level distributions.
- **Opposite-strand reads** are pooled into the column in REDItools2 without `-C`, which diluted or inflated editing where sense/antisense transcripts overlap.
- **Double-counting overlapping mates** (REDItools2/3, likely JACUSA2) vs samtools-mpileup overlap handling (AEI). Supporting-read thresholds (e.g. "≥2 reads") can be met by one fragment, and pseudo-replication inflates JACUSA2 LRT scores. alnbase `overlap` makes this testable.
- **Multi-mappers in Alu.** AEI keeps MAPQ 0 by default, REDItools needs ≥20, SPRINT excludes 255. Alu editing index is dominated by repeats, so this choice moves AEI across cohorts.
- **Duplicates.** JACUSA2 does not filter FLAG 1024 by default; REDItools does. [Morales 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10527054/) finds BCFtools, then REDItools2, report the most SNP-derived "RES", and dedup plus replicate merging change counts severalfold.
- **Version drift.** [JACUSA2 #47](https://github.com/dieterich-lab/JACUSA2/issues/47): 2.0.1 found ~202k vs 537k sites in 1.3.0 on identical inputs, mostly low-ratio sites lost, with no changelog explanation.
- **Tool agreement is low in general.** [Diroma 2019](https://doi.org/10.1093/BIB/BBX129): REDItools/JACUSA recover more Alu sites, JACUSA/RES-Scanner more non-Alu sites.
- SPRINT's `zz2sam.py` hard-codes strand (FLAG forced to 0) and inflates reads ([#23](https://github.com/jumphone/SPRINT/issues/23), unresolved).
- JACUSA2 arrest positions for R1 derive from mate start/TLEN, which is inconsistent with R2-derived positions when soft clips differ.

---

### Ranking inputs

**HITS-CLIP CIMS (CTK)**
- Popularity: moderate and declining (HITS-CLIP largely superseded by eCLIP/iCLIP; CTK still the reference implementation).
- Comparator cleanliness: low–medium. Perl plus a custom conda channel; MD-tag-based; multi-step with collapse and permutation randomness. The mutation extraction step (`parseAlignment.pl --mutation-file`) is deterministic and cleanly comparable.
- Ground truth/simulation: CTK ships `simulateCIMS.pl`; deletions are trivial to simulate; motif enrichment (UGCAUG) is a biological check.
- Discrepancy likelihood: medium. Clip stripping, the MD "match run anywhere" rule, `-w 1` dropping multi-nt deletions, terminal-base exclusions.
- Expressibility: **Y**.

**iCLIP/eCLIP truncation (PureCLIP, htseq-clip, iCount)**
- Popularity: high (ENCODE eCLIP, >300 RBPs; iCLIP widely used).
- Comparator cleanliness: htseq-clip `extract` is very clean (1 line per read, golden test BEDs); PureCLIP is clean for site extraction but its HMM output is not 1:1.
- Ground truth/simulation: sim_iCLIP simulator; exact coordinates are trivially simulated.
- Discrepancy likelihood: **high and consequential**. −1 vs 0 offset conventions (htseq-clip default, Skipper, JACUSA2), soft-clip handling, mate selection, supplementary/long-read drops.
- Expressibility: **Y with caveats** (pad-order trick depends on walk orientation and mate); **needs 5'/3'-as-sequenced pad codes + RNA library modes** for clean queries.

**A-to-I editing, REDItools2/3**
- Popularity: very high (REDIportal ecosystem, Nature Protocols).
- Comparator cleanliness: medium. Per-position TSV is directly comparable to an alnbase pileup, but v2/v3 defaults differ and there are several code quirks (SEQ-order trimming off-by-one, MeanQ vs MAPQ threshold, majority-strand pooling).
- Ground truth/simulation: bundled 1.4 MB BAM; GSE99249 ADAR1-KO; REDItools3 synthetic SAM generator; JACUSA2 simulated sets.
- Discrepancy likelihood: high (strand pooling, overlapping mates, trim orientation, Frequency definition).
- Expressibility: **Y** (homopolymer filter verbose; **benefits from proximity/run predicate and RNA library modes**).

**A-to-I editing, JACUSA2**
- Popularity: high and growing (also used for m6A, Nanopore, arrest).
- Comparator cleanliness: good (single jar, bioconda, BED6-extended with explicit 0-based coordinates and A,C,G,T vectors); statistics are layered on counts that should match alnbase exactly.
- Ground truth/simulation: shipped simulated datasets with truth tables.
- Discrepancy likelihood: medium–high (MAPQ default −1 in call-1, no dup filtering, arrest = terminal base not −1, strand naming).
- Expressibility: **Y** for counts; feature filters D/I/S/Y **need proximity predicate** (or verbose OR-patterns); rt-arrest positions **need read-end pad codes** for clarity.

**A-to-I editing, SPRINT**
- Popularity: medium.
- Comparator cleanliness: low (Python 2, bwa-aln-centric, MAPQ window 20–199, renaming reads, CIGAR parser limits).
- Ground truth/simulation: easy to simulate; the regular (non-hyper) path is comparable.
- Discrepancy likelihood: high (MAPQ 255 exclusion, H/=/X CIGAR misparse, strand from FLAG only).
- Expressibility: regular-site signal and the "5 nt from any block edge" rule **Y (single pattern, a nice showcase)**; hyper-editing **N** (needs realignment, out of scope).

**Alu Editing Index (RNAEditingIndexer)**
- Popularity: very high for cohort-level editing (GTEx/TCGA papers).
- Comparator cleanliness: medium (Docker, big resources) but has a **golden expected CSV on 71 MB of bundled BAMs**.
- Ground truth/simulation: ADAR1-KO (GSE99249) as a biological control; the index is a simple ratio, easy to recompute from alnbase rows in DuckDB.
- Discrepancy likelihood: medium–high (MAPQ 0 included, BQ≥30 after mpileup's BQ≥13 and overlap adjustment, trimBam 5 nt in SEQ order including soft clips, strand-deciding heuristics).
- Expressibility: **Y** (A/T reference bases in Alu BED; unstranded **benefits from `--library unstranded`** reference-orientation mode).


---

## Family 03 — RNA modifications detected by chemical/enzymatic conversion or RT signatures

Scope: assays in which a modification (m6A, m5C, Ψ, m1A, m3C, ac4C, m7G) is read out as a
difference between read and reference: substitution, deletion, or RT truncation (read start).
Antibody-enrichment peak assays (MeRIP/m6A-seq, m1A-MeRIP peaks, acRIP) are dropped because their
signal is coverage, not read-vs-reference. miCLIP CIMS/CITS belongs with the CLIP section.

Source was read for GLORI-tools, Bullseye, pseudoU-BIDseq (BID-pipe), eTAM-seq_workflow,
hisat-3n-table, m5C-UBSseq, meRanTK, PRAISE and mim-tRNAseq. The clones are in
`clones/`. Line references are to those clones as of 2026-09-16.

#### How alnbase orients reads (applies to every section)

Under `--library directional`, R1 fwd and R2 rev are walked on the top strand. R1 rev and R2 fwd are
walked reverse-complemented. A single-end read counts as R1. So:

- **Sense-stranded SE, or R1-sense PE** (GLORI-tools input, eTAM R2-as-SE, 10x 3′ R2 in a CellRanger BAM,
  Lexogen/SMARTer-type sense libraries, meRanTK `--readDir fr`): the walk runs 5′→3′ in RNA sense, so
  patterns are written as RNA sequence (U→T).
- **dUTP / reverse-stranded** (NEBNext Ultra Directional, TruSeq Stranded): the walk is *antisense*.
  RNA C>U becomes read A / refr G with the context reverse-complemented, so `RAC`→`GTY` read
  right-to-left. You can still write this today by complementing the pattern by hand, but it is
  error-prone. → candidate extension **`--library rna-reverse` (and `rna-forward`)**.
- **Unstranded** (the Bullseye default parse, many older RNA-seq sets): no flag-based orientation is
  correct. Reproducing Bullseye needs either the union of both complemented patterns plus annotation
  strand joined in DuckDB, or a **reference-orientation (no-complement) walk mode**.
- 3-letter aligners (HISAT-3N, meRanGs/meRanT, Bismark-style) write the conversion strand to an aux tag
  (`YZ:A:+/-`, `XG`). hisat-3n-table takes strand only from `YZ`. The `Library` doc comment in
  alnbase already anticipates a tag-based variant. → candidate extension **`--library tag:YZ`**.

---

### 1. GLORI / GLORI 2.0 / GLORI 3.0 (m6A, A→G of unmethylated A)

#### (a) Signal
- Glyoxal plus nitrite deaminates unmodified A to I, which reads as G. m6A resists and stays A. Per site,
  m6A level = A/(A+G) at a reference A. This is "inverted" bisulfite logic
  ([Liu 2022 Nat Biotech](https://www.nature.com/articles/s41587-022-01487-9)). GLORI 2.0/3.0 are
  milder, faster chemistries for low input and use the same readout
  ([Sun 2025 Nat Methods](https://www.nature.com/articles/s41592-025-02680-9)).
- No sequence context is required. Sites are reported in DRACH and non-DRACH.
- Libraries are stranded. GLORI-tools accepts only **single-end, A→G-converted (sense) reads**
  (README "Annotations").

#### (b) Comparators
- **GLORI-tools v1.0** ([github.com/liucongcas/GLORI-tools](https://github.com/liucongcas/GLORI-tools)).
  Python scripts that drive STAR ≥2.7.5c, bowtie1, samtools and pysam. Not on bioconda; install by
  cloning. Building the index is heavy: the README warns about STAR RAM, and issues #24/#27 report
  STAR index failures.
- The GLORI 2.0 code is on [Zenodo 14233421](https://zenodo.org/records/14233421). It was not
  inspected here.
- **Neutral alternative: `hisat-3n` + `hisat-3n-table --base-change A,G`**, on bioconda as
  `hisat-3n` 0.0.3. It is simple, but see §2(g).

#### (c) Output formats (from source)
- `${prefix}_referbase.mpi` (pipelines/pileup_genome_multiprocessing.py + get_referbase.py): TSV
  `chr, pos(1-based), strand(+/-), ref_base (strand-complemented for -), comma-list of per-read bases,
  comma-list of per-read A-counts`.
- Formatter output (m6A_pileup_formatter.py): `chr pos dir gene name trans isoform biotype total AG A T C G`,
  followed by 12 `cutoff;total,AG,A` bins for A-cutoffs 1..10,15,20.
- `${prefix}.totalm6A.FDR.csv`, the final site list. README §5.2 gives `Chr, Sites, Strand, Gene, CR, AGcov, Acov, Genecov, Ratio,
  Pvalue, P_adjust`. `Sites` is 1-based. `Ratio` = A/(A+G) at the chosen A-cutoff.
- `${prefix}.totalCR.txt`: `SA, A-to-G_ratio` per chromosome and gene.

#### (d) Implicit filters and defaults (source)
- **Mapping:** reads are converted A→G in silico and aligned to an A→G genome. STAR
  `--outFilterMismatchNmax 2`, `--outFilterMultimapNmax 1` (unique only),
  `--outFilterScoreMinOverLread 0.5`, `--outSAMmultNmax 1`. `samtools view -F 20` drops reverse hits,
  which are re-mapped to a T→C "rvsCom" genome when `--rvs_fac` is set. Unmapped reads go to a bowtie
  transcriptome pass (`-k 1 -m 1 -v 2`). The original read SEQ is restored afterwards. Contig names
  carry the suffix `_AG_converted`, so they must be renamed before an alnbase index can match them.
- **Pileup:** pysam `pileup(min_base_quality=0, ignore_overlaps=False)` with a custom BQ ≥ **10**.
  The argparse help says "default=30", but the code default is 10. pysam's default stepper drops
  unmapped, secondary, QC-fail and duplicate reads. Deletions are skipped (query_position None). There
  is no read-end trimming (`--trim-head/--trim-tail` default 0 and run_GLORI.py does not pass them).
  When the same read name appears twice at a position (genome plus transcriptome alignment), the
  observation with the higher BQ wins. Ties go to the one with the lower A-count.
- **`get_referbase.py` drops every position with only one read** (`len(line[5].split(",")) > 1`).
- **Per-read A-cutoff (default 3):** a read counts only if it has ≤3 `A` characters (or `T` for reverse
  reads). The count covers the **entire SEQ** (`query_sequence.count('A')`): soft-clipped bases,
  mismatches and bases whose reference is not A all count. The "signal" filter then requires
  cov(≤cutoff)/cov(all) ≥ 0.8.
- **Site calls:** A+G ≥ 15, A ≥ 5, ratio ≥ 0.1, and (A+G)/total ≥ 0.8 at both all-read and cutoff
  levels. One-sided binomial against the **gene-level** non-conversion rate, p < 0.005 then BH FDR < 0.005.
  Genes whose non-conversion rate is ≥ 0.2 get background 1.0, so no site in them can ever be called.
  The gene CR only counts positions with A+G ≥ 15. Duplicates are removed upstream by UMI (seqkit
  rmdup -s, README).

#### (e) alnbase expressibility — **Yes, today**
- Queries (SE, sense walk): `m6A: read="A" refr="A"`, `unconv: read="G" refr="A"`.
- The per-read A-cutoff is a DuckDB group-by on `record_id`. To reproduce GLORI exactly, run with
  `--soft-clips emit` and a query `read="A" refr="~"`, because GLORI counts A over all of SEQ.
- Everything else (BQ ≥ 10, same-name conflicts, coverage ≥ 2, gene CR, binomial test) is downstream.
- Needed preprocessing only: rename the `_AG_converted` contigs.

#### (f) Test data
- GLORI 1.0: [GSE210563](https://www.omicsdi.org/dataset/geo/GSE210563). HEK293T mRNA
  SRR21356251 is about 450 M reads after dedup (GLORI-tools issue #31), and the published site file
  GSM6432590_293T-mRNA-1_35bp_m2.totalm6A.FDR.csv serves as a reference call set. Subset to one
  chromosome for speed.
- Nature Protocols data: [GSE233875](https://www.omicsdi.org/dataset/geo/GSE233875).
- GLORI 2.0/3.0: [GSE270643](https://www.omicsdi.org/dataset/geo/GSE270643).
- Ground truth: the GLORI papers use synthetic m6A/A probes at defined fractions. A→G conversion is
  also trivial to simulate.

#### (g) Correctness issues and contestable choices
- Issue [#31](https://github.com/liucongcas/GLORI-tools/issues/31): a re-implementation got about
  100 k sites against the published 214 k. A community reply says omitting or mis-running
  `--rvs_fac` loses reverse-strand sites. That makes the call set pipeline-fragile.
- The A-cutoff counts A over the whole SEQ, not over aligned reference-A positions. Reads from A-rich
  transcripts and longer reads are discarded more often at equal conversion efficiency. That biases
  coverage, and stoichiometry, against A-rich 3′UTRs. A reference-aware per-read count, which alnbase
  gives naturally, is the principled alternative, and the difference is directly testable.
- Gene-level background with a hard CR ≥ 0.2 exclusion makes whole genes invisible. That can be read
  as "gene X lacks m6A" when the gene simply converted poorly (structured, GC-rich).
- The BQ default is 10 in code but "30" in the help text.
- The coverage-1 drop happens before the formatter, which changes gene CR estimates at low depth.
- GLORI-tools is single-end only (issue [#28](https://github.com/liucongcas/GLORI-tools/issues/28)).
  PE users improvise.

#### Ranking inputs
Popularity high and rising (GLORI is the de facto m6A gold standard; [m6AConquer](https://www.biorxiv.org/content/10.1101/2024.09.10.612173v1.full.pdf) uses GLORI/eTAM as ground truth) · comparator messy (scripts, STAR index, suffixes) · simulation easy, spike-ins exist · discrepancy likely (whole-SEQ A-count, gene CR exclusion, cov≥2 drop) · expressible **Y**.

---

### 2. eTAM-seq (m6A, TadA-8.20 A→I, same logic as GLORI)

#### (a) Signal
- Evolved TadA deaminates unmodified A to I (read as G). m6A stays A. Stoichiometry = A/(A+G) at a
  reference A ([Xiao 2023 Nat Biotech](https://www.nature.com/articles/s41587-022-01587-6)). eTAM-seq-v2
  was preprinted in 2025 ([bioRxiv](https://www.biorxiv.org/content/biorxiv/early/2025/05/15/2025.05.11.653357.full.pdf)).
- The workflow aligns **R2 as single-end with `--rna-strandness F`**, i.e. R2 is sense.

#### (b) Comparators
- **eTAM-seq_workflow** ([github.com/shunliubio/eTAM-seq_workflow](https://github.com/shunliubio/eTAM-seq_workflow)).
  Bash scripts plus HISAT-3N, [pileup2var](https://github.com/shunliubio/pileup2var), and the R
  package `eTAMseq` shipped as a tarball. Moderate install. hisat-3n is on bioconda; pileup2var and the R package are not.
- **hisat-3n-table** is the documented alternative counter.

#### (c) Output formats
- **hisat-3n-table** (position_3n_table.h:334) is a TSV with header
  `ref  pos  strand  convertedBaseQualities  convertedBaseCount  unconvertedBaseQualities  unconvertedBaseCount`.
  `pos` is **1-based** (line 64), `strand` is `+`/`-`, and the quality strings are raw Phred+33.
- The eTAM R model input (3_run_model_ftop.R header) has columns
  `pos(chr_pos_strand) motif type(DRACH/nonDRACH) ftom_G ftom_A ftop_G ftop_A`.

#### (d) Implicit filters (source)
- HISAT-3N `--base-change A,G`. The rRNA pass uses `--no-softclip --norc`. The genome pass is followed
  by `samtools view -q 60` (unique).
- `umi_tools dedup --method=unique --spliced-is-unique`.
- **Per-read filter `[Yf]/([Yf]+[Zf]) >= 0.5`**: at least 50% of the read's reference A positions must be
  converted. Yf/Zf are HISAT-3N's per-read converted and unconverted counts.
- pileup2var `-f 524`, `-c 1`. The model requires coverage ≥ 10 in both FTO− and FTO+ (or IVT).
- **hisat-3n-table** (alignment_3n_table.h):
  - Strand comes **only from the YZ tag**. The MD tag is required.
  - There is **no base-quality or MAPQ filter**; `-u` means unique by NH.
  - Soft clips and insertions are skipped.
  - In a matched segment, a base counts only if the read shows the "from" base. At a mismatch it
    counts only if ref = from and read = to. Every other base is dropped.
  - Reads spanning more than 500 kb are ignored.
  - Overlapping mates are merged by read name. If they disagree, the position is removed
    (position_3n_table.h:131-166).

#### (e) alnbase — **Yes, today**
Same queries as GLORI, on the SE BAM from HISAT-3N. Yf/Zf come through `-f` aux fields, or can be
recomputed per `record_id` in DuckDB. Filtering is `samtools -q 60` upstream.

#### (f) Test data
[GSE201063](https://www.omicsdi.org/dataset/geo/GSE201063) (site-specific quantification subset; the full
HeLa/mESC sets are in the same SuperSeries). IVT controls give false-positive calibration.

#### (g) Issues
- **Possible hisat-3n-table miscount (from reading the source, not yet confirmed empirically).** When mates disagree,
  `appendReadNameID` removes the earlier observation by erasing the first character in the converted
  (or unconverted) quality string that equals the *incoming* base's quality. If the earlier
  observation had a different quality, the wrong entry, or none, is removed. The disagreeing
  observation is then only partly retracted, so counts and qualities drift. A fixture with overlapping
  mates of unequal quality would show it immediately. This is a good alnbase-vs-comparator
  demonstration.
- The per-read ≥50% conversion filter depends on read length and base composition. So does GLORI's
  cutoff, but in a different way, and m6AConquer notes that heterogeneous filters drive most cross-method
  discordance.
- GLORI-tools issue #29 asks whether it is eTAM-compatible, a sign that the two pipelines are
  not interchangeable.

#### Ranking inputs
Popularity medium-high · comparator moderate (hisat-3n-table is simple and bioconda) · simulation easy, IVT controls · discrepancy likely (mate-conflict code, 50% read filter) · **Y**.

---

### 3. DART-seq / scDART-seq (m6A-adjacent C→U by APOBEC1-YTH)

#### (a) Signal
- APOBEC1-YTH deaminates the C **immediately 3′ of the m6A** in RAC (m6A-C), giving C→U, read as T
  ([Meyer 2019 Nat Methods](https://www.nature.com/articles/s41592-019-0570-0);
  [Tegowski 2022 Mol Cell](https://www.cell.com/molecular-cell/fulltext/S1097-2765(21)01143-6)).
- Sites are compared with a YTH-mutant (YTHmut) or APOBEC1-only control. The original thresholds were
  coverage ≥ 10, edit ratio 5–95%, ≥1.2× over mutant, and ≥2 C→U reads.
- Bulk libraries are typically dUTP reverse-stranded. scDART uses 10x 3′, where R2 is sense.

#### (b) Comparators
- **Bullseye**, on bioconda as `bullseye` 1.0.0
  ([github.com/mflamand/Bullseye](https://github.com/mflamand/Bullseye)). Perl plus samtools, tabix and
  bedtools; easy install.
- **JACUSA2 + JACUSA2helper**, on bioconda as `jacusa2` 2.1.17, has a DART-seq vignette
  ([JACUSA2helper DART-seq](https://dieterich-lab.github.io/JACUSA2helper/articles/JACUSA2helper-dart-seq.html)).
  Strand-aware; covered in the editing section.

#### (c) Outputs (source)
- **parseBAM.pl matrix** (bgzip plus tabix `-b 2 -e 2`): `chr, pos (1-based SAM POS arithmetic), A, T, C, G, N,
  total[, strand]`. Strand is only present with `--stranded`.
- **Find_edit_sites.pl BED**: header
  `#chr start end gene score strand control_ratio control_total dart_ratio dart_total conversion` (plus
  `score` with `--score`). BED 0-based half-open. `gene` = `Gene|region|C2U|mut=<n C→T reads>|<score>|rep=N`
  (example_data/expected_output). Col5 = dart edit ratio.
- **RACfilter.sh** keeps rows whose 3-nt window `[pos-2,pos]` (strand-aware bedtools getfasta) matches `[AG]AC`.

#### (d) Implicit filters (source)
- **parseBAM:**
  - Unmapped reads are dropped.
  - Duplicates are removed only with `--removeDuplicates`; secondary alignments only with `--removeMultiMapped`.
  - **No BQ or MAPQ filter.**
  - Soft clips and insertions are removed and deletions are not counted.
  - **Library strand is ignored unless `--stranded`**. The maintainer confirms both strands are pooled
    ([issue #4](https://github.com/mflamand/Bullseye/issues/4)).
  - Overlapping mates are counted twice.
  - `--minCoverage` default 1; the example script uses 10.
- **Find_edit_sites:**
  - Only positions inside refFlat features are considered; strand comes from the annotation.
  - For "−" genes, G→A is counted.
  - Denominator = total − N, which includes all other substitutions.
  - `--minEdit 10`, `--maxEdit 80`, `--EditedMinCoverage 10`, `--ControlMinCoverage 10`, `--MinEditSites 2`.
    `--editFoldThreshold` code default 1.5; the help text says 1.
  - The **control must be ≥60% non-edited** (drops heterozygous SNPs).
  - SNPs are masked only via `--filterBed`, and that mask ignores strand.

#### (e) alnbase — **Yes today with caveats; better with a library mode**
- 10x or sense SE: `m6A_adj: read="~~T" refr="RAC"` with the anchor in column 2, plus `unedited: read="~~C"`.
- dUTP PE: complement by hand, `read="A~~" refr="GTY"`, anchor column 0.
- To reproduce Bullseye's unstranded pooling, query both orientations and assign strand by annotation
  in DuckDB. Cell barcodes come through `-f CB`.
- **Needs (for ergonomics and exact reproduction):** `--library rna-reverse` and an unstranded walk mode.

#### (f) Test data
- **Bullseye `example_data/`** (repo, 860 KB): four BAMs (WT/Mettl3KO soma, Apc locus, mm10), a refFlat,
  and an expected-output BED. It runs in seconds and is ideal for CI.
- Bulk and single-cell data: GSE180954 (per issue #3; the RAC-sites BED is the YTHmut-only file).

#### (g) Issues
- The unstranded parse assigns antisense-transcript reads and opposite-strand overlapping genes to the
  annotated strand. A G→A artefact on the minus strand can be called as C→U on a plus-strand gene.
- Mate overlap is double-counted, which inflates both coverage and edit counts. With no BQ filter,
  low-quality tails contribute false C→T.
- The edit ratio is computed over *total*, not C+T. Positions with other mismatches, such as an SNP or
  an RT error hotspot, have their ratio deflated.
- Issue #3: published bulk results could not be reproduced (partly a chr naming mismatch).
- Issue #4: about 50% overlap with the older CTK-based pipeline on the same data. That is a concrete
  "tool choice flips site set" case.
- Help text and code disagree on the fold-threshold default (1 vs 1.5).

#### Ranking inputs
Popularity medium (m6A without antibody; scDART) · comparator clean-ish, bioconda, tiny test data · simulation easy · discrepancy very likely (unstranded pooling, no BQ, mate double-count) · **Y** (needs `rna-reverse`/unstranded mode for ergonomics).

---

### 4. m6A-SAC-seq (m6A → allyl-m6A → cyclized adduct, RT misincorporation)

- **(a)** MjDim1 allyl-labels m6A, and iodine cyclization makes RT misincorporate at the m6A position.
  The signal is mixed mutations at a reference A (mostly A→T/A→G/A→C) plus deletions, compared with an
  untreated control ([Hu 2022 Nat Biotech; Nat Protoc 2023](https://www.nature.com/articles/s41596-022-00765-9)).
- **(b)** [y9c/m6A-SACseq](https://github.com/y9c/m6A-SACseq) (Snakemake/Docker style, same author as
  BID-pipe/UBS-seq) and [shunliubio/m6A-SAC-seq](https://github.com/shunliubio/m6A-SAC-seq). Not inspected in depth.
- **(e)** **Y:** `read="{CGT.}" refr="A"` (mutation or deletion). Stoichiometry is calibrated per motif
  downstream.
- **(g)** Same family of issues as BID-seq: motif-specific calibration curves and background from
  the untreated sample.
- **Ranking:** medium popularity · comparator Docker/Snakemake · **Y**. Lower priority than GLORI/eTAM.

---

### 5. BID-seq (Ψ → bisulfite adduct → RT deletion)

#### (a) Signal
- Neutral bisulfite forms a stable Ψ-monobisulfite adduct, which RT skips. The result is a **1-nt
  deletion at Ψ** with no C→U background
  ([Dai 2023 Nat Biotech](https://www.nature.com/articles/s41587-022-01505-w);
  [Nat Protoc 2023](https://www.nature.com/articles/s41596-023-00917-5)).
- Stoichiometry comes from per-5-mer calibration curves (NNΨNN spike-ins).
- The pipeline default is `forward_stranded: true` (R1 = RNA). If false, reads are reverse-complemented
  before mapping (`rcFastq`).

#### (b) Comparator
**BID-pipe v2.0** ([github.com/y9c/pseudoU-BIDseq](https://github.com/y9c/pseudoU-BIDseq); Zenodo
10.5281/zenodo.8158036). Ships only as a Docker/apptainer image, `y9ch/bidseq`. Several steps are **prebuilt
binaries with no source in the repo** (`cpup`, `samFilter`, `deletionFilter`, `joinFastq`, `rcFastq`;
issues #17 and #19 ask for the source). Easy to run, hard to audit. Not on bioconda.

#### (c) Outputs
- `call_sites/{genes,genome}.tsv.gz`: `chr, pos, strand, <sample>_depth, <sample>_gap …`. pos is 1-based
  (mpileup). **depth includes gaps**, per the 2022-07-26 change log in adjustGap.
- `filter_sites/*.tsv.gz`: the same columns plus `<group>_ratio`, `<group>_fraction` (calibrated) and
  `<group>_passed` (bin/pickSites.py).

#### (d) Implicit filters (Snakefile, config.yaml)
- **Alignment:**
  - Contamination, then genes (bowtie2 end-to-end `--norc -a`, `--rdg 1,2`), then genome (STAR local,
    `--scoreDelOpen -1 --scoreDelBase -1`, `--outFilterMultimapNmax 10`, `--outSAMmultNmax -1`, so every
    best-scoring multimapper is written).
  - **realignGap** re-aligns reads that have D, S or a non-numeric MD with parasail SW (gap open 3,
    extend 2, dnafull). Spliced reads are split and re-joined. Reads whose re-alignment changes the
    query length are dropped.
  - `samtools view -e '[NM]<=5 && [NM]/(qlen-sclen)<=0.1'`.
  - Dedup is UMIcollapse (`--data naive --merge avgqual --two-pass`) when UMIs are present.
- **Pileup:** `samtools mpileup -aa -B -d 0 -Q 5 --reverse-del`. Forward uses `--ff 3608`; reverse uses
  `--rf 16 --ff 3592`. These flag masks **do not exclude 256 (secondary)**. Combined with
  `--outSAMmultNmax -1`, that suggests multimappers can be counted at every locus. The maintainer says
  they should not be (issue [#15](https://github.com/y9c/pseudoU-BIDseq/issues/15)), so an
  empirical check is needed. In mpileup, a deletion's `-Q` check uses the quality at `p->qpos`, which is a
  neighbouring read base (samtools bam_plcmd.c:752-755).
- **adjustGap:** a deletion inside a repeat (the same k-mer shifted up to the run length) is
  **reassigned to the putative position with the highest upstream-weighted gap signal summed across
  all samples** (UPSTREAM_NUM = 5, weight 0.5^i). Depth at the other positions is adjusted.
- **Prefilter:** group gap ≥ 5, depth ≥ 10, ratio ≥ 0.01.
- **Final:** treated depth ≥ 20, input depth ≥ 20, treated gaps ≥ 5, ratio ≥ 0.02, calibrated
  fraction > 0.02, ratio ≥ 2× input, Fisher p < 1e-4. The paper's cutoffs were stricter (issue
  [#5](https://github.com/y9c/pseudoU-BIDseq/issues/5)).

#### (e) alnbase — **Yes for the raw signal**
- Query: `psi_del: read="~.~" refr="~T~"` with the anchor on the gap column (parquet only, since
  `bases` tags refuse gap anchors). Denominator: `read="~{T.}~"`.
- The motif calibration needs the 5-mer: `refr="NNTNN"` with `^` marks, or use refr_base plus
  neighbours.
- Upstream: realignment (it is a BAM rewrite). Downstream: gap redistribution across a U-run (record
  ±5 refr columns with `^` so the run extent is known), calibration, Fisher test.
- No query-syntax extension is required. A walk option to left- or right-normalize gaps in homopolymers
  would be a convenience, but GATK LeftAlignIndels does this upstream.

#### (f) Test data
- GSE179798 (original BID-seq; needs the barcode `NNNNNXXX-XXXNNNNNATCACG`, issue #3).
- Nature Protocols data and calibration probes: `calibration_curves.tsv` (256 motifs) ships in the repo.
- The `test/` yaml files reference FASTQs that are not in the repo.

#### (g) Issues
- **U-tract ambiguity.** A deletion in `UUU` can sit at any of three positions. Aligner placement (STAR
  left-aligns within its scoring), parasail re-placement and adjustGap's cross-sample reallocation can
  each move the site and change its stoichiometry. BACS authors report that BS/deletion methods "cannot
  determine the exact position of Ψ in consecutive uridine sequences" and that local realignment "can
  generate artifacts including overestimation of Ψ modification level"
  ([BACS, Nat Methods 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11541003/)). BID and BACS overlap
  only 40% (230/575). **Flip scenario:** a PUS7-dependent site called at U2 in `UUU` in one pipeline
  and at U3 in another looks like two different sites with half the stoichiometry. alnbase
  "as-aligned" plus explicit run-aware aggregation makes the choice visible.
- Deletions near read ends are soft-clipped by local aligners. That was the motivation for realignGap
  (docs/Algorithms.md example 2), and it means the signal depends on read-end handling.
- Closed binaries and the secondary-flag question above.

#### Ranking inputs
Popularity medium-high (Ψ quantification standard with PRAISE) · comparator Docker-only, partly closed · calibration spike-ins give ground truth; deletion simulation easy · discrepancy **very likely** (gap placement, secondary alignments) · **Y**.

---

### 6. PRAISE (Ψ bisulfite deletion; Yi lab)

- **(a)** Same chemistry as BID-seq (bisulfite at near-neutral pH gives a Ψ deletion)
  ([Zhang 2023 Nat Chem Biol](https://www.nature.com/articles/s41589-023-01304-7)). Takara, KAPA and
  eCLIP library variants exist; **R2 only** is used for Takara/eCLIP.
- **(b)** [github.com/Zhe-jiang/PRAISE](https://github.com/Zhe-jiang/PRAISE): Python scripts plus a
  **vendored, modified biopython** (custom substitution matrices). Not on bioconda; install friction is
  medium (the patched biopython build).
- **(c)** `samtools mpileup -d 15000000 -BQ0 --ff UNMAP,QCFAIL -aa` goes through `parse-mpileup.py`
  into a `.bmat` per-base matrix.
- **(d)**
  - Pre-map dedup with `seqkit rmdup -s` (sequence identity, before UMI extraction).
  - Trim 14 nt at 5′ (8 UMI + 6 template switch) and 6 nt at 3′ (random-primer indels).
  - hisat2 `--no-spliced-alignment --very-sensitive` against a **transcriptome** (realignment cannot
    handle introns).
  - Custom realignment (`realignment_forward.py`/`reverse.py -ms 4.8`).
  - `remove_end_signal.py` converts end-proximal gaps to soft clips.
  - `remove_multi_mapping.py`.
  - The pileup keeps **secondary and duplicate** flags (`--ff UNMAP,QCFAIL` only) with **BQ 0**.
- **(e)** **Y** (same as BID-seq). The read-end gap exclusion is a DuckDB filter on `off_5p`/`off_3p`.
- **(f)** Data accession not confirmed in this pass (check the article's data availability statement).
  HEK293T reported 2,209 sites.
- **(g)**
  - Pre-alignment sequence dedup collapses identical reads from highly expressed short RNAs. That
    understates depth for tRNA/snoRNA and can change their ratio estimates.
  - With BQ 0, all bases count.
  - The end-signal-to-softclip rule decides whether a Ψ near a read end is counted.
- **Ranking:** medium popularity · scripts plus patched biopython · **Y** · good companion to BID-seq
  for showing how two pipelines disagree on the same chemistry.

---

### 7. RBS-seq (m5C + Ψ + m1A from one bisulfite library)

- **(a)**
  - m5C: C retained against C→T.
  - Ψ: a modest bisulfite-induced deletion.
  - m1A: mismatch in the non-bisulfite library, lost after Dimroth rearrangement to m6A.
  - Source: [Khoddami 2019 PNAS](https://www.pnas.org/doi/10.1073/pnas.1817334116).
- **(b)** No maintained public tool; custom analysis.
- **(e)** **Y**. All three are pattern unions at single columns:
  - `read="C" refr="C"` (m5C)
  - `read="." refr="T"` (Ψ)
  - `read="/" refr="A"` on the untreated library (m1A)
- **(f)** [GSE90963](https://www.pnas.org/doi/10.1073/pnas.1817334116), HeLa.
- **Ranking:** low (no comparator). Useful only as a "one query file, three marks" illustration.

### 8. BACS (Ψ → C substitution)

- **(a)** 2-bromoacrylamide cyclization makes RT read Ψ as C (U→C in cDNA). It also reports A-to-I
  and m1A in the same data
  ([Xu 2024 Nat Methods](https://www.nature.com/articles/s41592-024-02439-8)).
- **(e)** **Y**: `read="C" refr="T"` (sense walk).
- **(g)** The authors present substitution signatures as resolving U-tract ambiguity. A BACS-vs-BID
  comparison through the same alnbase queries is a neat internal consistency demonstration.
- **Ranking:** medium-low (code availability unclear).

---

### 9. RNA bisulfite m5C: RNA-BS-seq (meRanTK) and UBS-seq

#### (a) Signal
- Unmethylated C → U (read T). m5C stays C, so the level is C/(C+T) at a reference C.
- Libraries are stranded. meRanTK `--readDir fr` (default) means R1 sense; `rf` handles dUTP.
- UBS-seq uses ultrafast high-concentration bisulfite, which lowers background
  ([Dai 2024 Nat Biotech](https://www.nature.com/articles/s41587-023-02034-w)).

#### (b) Comparators
- **meRanTK 1.3.0** ([github.com/icbi-lab/meRanTK](https://github.com/icbi-lab/meRanTK); the docs PDF is in
  the repo). Perl with bundled STAR/HISAT2/bowtie2 wrappers and a conda env yml. Not on bioconda;
  medium install.
- **m5C-UBSseq v0.1** ([github.com/y9c/m5C-UBSseq](https://github.com/y9c/m5C-UBSseq), Zenodo
  10.5281/zenodo.11046885). Snakemake plus hisat-3n / hisat-3n-table plus polars. Medium install.

#### (c) Outputs
- **meRanCall** (src/meRanCall.pl:503) is a TSV with header
  `#SeqID refPos refStrand refBase cov C_count methRate mut_count mutRate CR SNR CalledBase CB_count state
  95_CI_lower 95_CI_upper p-value_mState p-value_mRate score seqContext geneName candidateName`. With
  `--fdr` it adds `p-value_mState_adj` or `p-value_mRate_adj`. Optional BED6+3 (`-bed63`) and narrowPeak.
  refPos follows Bio::DB::Sam pileup, 1-based.
- **UBS-seq** writes hisat-3n-table columns (1-based) for four read subsets
  (unfiltered/filtered × unique/multi), which are joined. `detected_sites/filtered/*.tsv` has columns
  `ref pos strand u d ur pval passed`.

#### (d) Implicit filters
- **meRanCall defaults** (src/meRanCall.pl:63-84):
  - Base quality and mate overlap:
    - `minBaseQ 30`.
    - Overlapping mates are resolved per position by BQ. **If the second mate seen has the higher BQ,
      both mates are counted** (lines 807-812 only skip the later mate when it is not better).
  - Read-level filters:
    - **`C_cutoff 3`**: a read with more than 3 C (or G on the minus strand) across its **whole SEQ**
      is excluded. The CHANGES file says a `<=` fix is "not yet released", but HEAD already has `<=`.
    - `maxDup 0` (off; when on, duplicates are keyed by start+CIGAR+flag).
    - Read-end skipping `fs5/fs3/rs5/rs3` (default 0) is corrected using `readLength 100`. Note the
      off-by-one `>` rather than `>=` at 3′.
  - Site-level thresholds:
    - `minCov 20`, `minC 3`, `minMethR 0.2`.
    - `SNR 0.9`: the fraction of reads passing C_cutoff must exceed 0.9.
    - `minMutR 0.8` (only relevant with `--reportUP`).
  - Positions with N are skipped. So are indels and ref-skips.
  - In `--transcriptDBref` mode, a wrongly oriented mate is skipped with a warning.
- **UBS-seq:**
  - hisat-3n `--directional-mapping --base-change C,T`. PE requires proper pairs.
  - Per-read filter `[XM]*20 <= (qlen-sclen) && [Zf] <= 3 && 3*[Zf] <= [Zf]+[Yf]`: at most 3 unconverted
    C, at most 1/3 unconverted, mismatches ≤ 5%.
  - Prefilter: depth ≥ 20, support ≥ 3, unconverted ratio ≥ 0.02, clustered ratio < 0.5, multi ratio < 0.2.
  - Final: binomial against the mean background (non-candidate sites) with p < 0.001, u ≥ 2, d ≥ 10, ur > 0.02.

#### (e) alnbase — **Yes**
- Queries: `m5C: read="C" refr="C"`, `conv: read="T" refr="C"`, sense walk.
- The C_cutoff and Zf ≤ 3 filters are per-read aggregations in DuckDB. Emit clipped bases for
  meRanTK-exact counts.
- `alnbase overlap` handles mate overlap with an explicit policy.
- Needs `rna-reverse` for dUTP libraries, or hand-complemented patterns.

#### (f) Test data
- **meRanTK `testdata/`** (118 MB in repo): mm10 chr19 FASTA+GTF, refSeq subset, tRNAs, and FASTQs
  (`clean_KHOD_400k_test.fastq`, r1/r2). This makes a quick self-contained run.
- UBS-seq: [GSE225614](https://pmc.ncbi.nlm.nih.gov/articles/PMC11217147/). The repo `data/` FASTQs are
  placeholders (4 KB).

#### (g) Issues
- **meRanCall mate-overlap asymmetry** double-counts overlaps whenever the later-iterated mate has the
  higher BQ, so overlapping fragments contribute 1 or 2 votes depending on iteration order. That is
  directly demonstrable with a two-read fixture.
- The whole-SEQ C_cutoff has the same length and composition bias as GLORI. RNA m5C in mRNA is famously
  contentious: most early mRNA m5C calls were incomplete conversion in structured regions. The C_cutoff and
  SNR filters decide how many of those survive, so this is a real "choice flips conclusion" case
  (mRNA m5C abundant vs rare).
- hisat-3n-table inherits the mate-conflict quality-erase issue (§2g).

#### Ranking inputs
Popularity medium (m5C RNA-BS long-standing; UBS-seq newer) · comparator: meRanTK has built-in test data but is old Perl; UBS-seq is Snakemake · simulation easy (bisulfite) · discrepancy likely (mate overlap, C_cutoff) · **Y** — also closest to alnbase's native bisulfite use case, so least "versatility" value.

---

### 10. RT-signature marks without chemistry: m1A, m3C, m1G, m2,2G, I (tRNA-seq; HAMR)

#### (a) Signal
- Watson-Crick-face methylations (m1A, m3C, m1G, m2,2G, and m1I at tRNA 37) cause **mixed
  misincorporation plus RT stops**. The pattern of mismatch types at a reference base is diagnostic.
- Truncation shows up as read 5′ ends (cDNA 3′) at position +1 of the modified base.
- Demethylase (AlkB, DM-tRNA-seq, ARM-seq) or Dimroth (m1A→m6A, m1A-seq/m1A-MAP) comparisons remove the signal.
- TGIRT, MarathonRT and SuperScript IV differ strongly in read-through.
- tRNA-seq reads are generally sense.

#### (b) Comparators
- **mim-tRNAseq v1.3.11**, on bioconda as `mimseq`
  ([github.com/nedialkova-lab/mim-tRNAseq](https://github.com/nedialkova-lab/mim-tRNAseq)). GSNAP
  SNP-tolerant alignment to clustered tRNA references (modification sites supplied as "SNPs").
  Easy via conda; heavy runtime.
- **HAMR** ([github.com/wanglab-upenn/HAMR](https://github.com/wanglab-upenn/HAMR); newer HAMRLNC 2025).
  Python; not bioconda. Parameters: `min_read_qual`, `min_read_cov`, `seq_error_rate`, `max_p`/`max_fdr`,
  `refpercent`.

#### (c) Outputs (mim-tRNAseq, mmQuant.py)
- `mods/mismatchTable.csv`: `isodecoder, pos, type (A/C/G/T), proportion, condition, bam, cov`.
- `mods/RTstopTable.csv`: `isodecoder, pos, proportion, cov, condition, bam`.
- `mods/readthroughTable.csv`.
- `mods/predictedMods.csv` (unannotated sites above `--misinc-thresh 0.1`).
- `single_read_data/<sample>/<isodecoder>.tsv.gz`: per-read mismatch at each canonical position, plus a
  `Charged` column from CCA analysis.
- Positions are in Sprinzl-style canonical tRNA numbering (via ssAlign) and 1-based.

#### (d) Implicit filters (mim-tRNAseq)
- **No base-quality filter.** Mismatches are read from the MD tag.
- Soft clips at both ends are removed. RT non-templated additions are clipped, and the stop is recorded
  at alignment start + 1.
- GSNAP `--max-mismatches` (default: GSNAP's automatic "ultrafast" value; a value in 0–1 is taken as a
  fraction of read length); `--remap` pass with `--remap-mismatches`.
- `--min-cov 0.0005` (fraction of total reads per cluster).
- Isodecoder deconvolution rules: `--cluster-id 0.97`, `--deconv-cov-ratio 0.5`.
- Reads whose 5′ extends beyond the assigned member are treated as full-length.

#### (e) alnbase — **Yes**
- Mismatch at reference A: `read="/" refr="A"`, or explicit per type `read="G" refr="A"`.
- **RT stop:** `read="_~" refr="~A"` with the anchor on column 1. That is the first aligned base, with
  the modified base one column 5′ read from the reference past the read end. It works because pads
  carry reference context and skipped soft clips sit outside the walk.
  - Caveat: **`off_5p` is not 0 when the read has a 5′ soft clip**, so `off_5p==0` in DuckDB is
    *not* equivalent. Use the pad pattern.
  - For PE, restrict to R1 via flag, because a CTOT/R2 leading pad is a 3′ end.
- Isodecoder deconvolution and canonical numbering are out of scope; run on mim-tRNAseq's own BAMs.

#### (f) Test data
- **mim-tRNAseq repo:** `mimseq_hek_1/2.fastq.gz` and `mimseq_k562_1/2.fastq.gz`, about 52 MB each,
  with the `sampleData_HEKvsK562.txt` sample sheet. Self-contained.
- m1A-MAP data (Li 2017 Mol Cell) and DM-tRNA-seq (Zheng 2015) exist on GEO; accessions not verified here.

#### (g) Issues
- No BQ filter plus MD-based counting means sequencing errors in low-quality tails inflate
  misincorporation, and the effect differs by read length. Rerunning the comparison with an alnbase BQ
  predicate shows the impact.
- RT-stop placement depends on soft-clipping of non-templated additions. An aligner that does not clip
  them shifts stops by 1–3 nt. That matters for m1A58 vs m1A57 calls.
- GSNAP SNP-tolerance at known modification positions hides mismatches from the alignment score, so
  alignments differ from vanilla aligners and mismatch rates are not portable between pipelines.
- m1A in mRNA 5′UTRs was largely attributed to antibody cross-reactivity
  ([Grozhik 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6851129/)). Only RT-signature and
  demethylase-difference analyses survive, which makes read-level mismatch evidence the arbiter.

#### Ranking inputs
Popularity medium (tRNA-seq field) · comparator bioconda with bundled data, but tRNA-specific reference and deconvolution confound site comparison · ground truth: MODOMICS-annotated tRNA modifications · discrepancy likely (no BQ; stop offset) · **Y** (read-start via `_` pad).

---

### 11. ac4C-seq and RedaC:T-seq (N4-acetylcytidine → reduced base → C→T)

- **(a)** NaCNBH3 (ac4C-seq) or NaBH4 (RedaC:T) reduction makes ac4C misread as T. Controls are mock or
  deacetylated samples, or NAT10−/−. The motif is **CCG with the middle C modified** in yeast and
  archaea.
- **(b)** ac4C-seq: [SchwartzLab/ac4c-seq](https://github.com/SchwartzLab/ac4c-seq). STAR 2.5.3a local,
  **JACUSA (v1) pileup**, R/lme4
  ([Thalalla Gamage 2021 Nat Protoc](https://pmc.ncbi.nlm.nih.gov/articles/PMC9103714)).
- **(c)** `ac4c_significantSites.txt` holds misincorporation rates, p-values, a ±10 bp context and the
  motif flag.
- **(d)**
  - Filters: C→T must be the most frequent non-C, C→T reads ≥ 3, treated rate ≥ 2–3%, control ≤ 1–5%,
    difference ≥ 2%.
  - Exclude SNPs, repeats, NUMTs, paralogs and rRNA/tRNA repeats.
  - RedaC:T (original authors): MAPQ 60, BQ 20, depth ≥ 10, strand-aware mismatch calling on concordant
    pairs.
- **(e)** **Y**: `read="~T~" refr="CCG"` with the anchor on column 1, plus all mismatch types for the
  "C→T must dominate" rule.
- **(f)** [GSE135826](https://pmc.ncbi.nlm.nih.gov/articles/PMC9103714) (HeLa, *T. kodakarensis*).
- **(g) The key flip case.**
  - Georgeson & Schwartz reanalysed RedaC:T-seq and found no evidence for ac4C in human mRNA
    ([Mol Cell 2024](https://www.cell.com/molecular-cell/fulltext/S1097-2765(24)00227-2)):
    - C→T excess appeared in essentially one of two WT replicates.
    - All mismatch types, not only C→T, were elevated.
    - The excess disappeared after **removing reads with identical 5′ and 3′ ends (PCR duplicates)**
      and was attributed to low library complexity.
    - One "site" was SNP rs1065711.
  - Arango/Oberdoerffer countered
    ([Mol Cell 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11353019/)):
    - 70–95% of duplicates are natural, not PCR.
    - About 2,000 C:T sites are reproducible and absent in NAT10−/−.
    - They used mismatch-type-specific error from input (75th percentile).
  - **Duplicate handling, replicate pooling and mismatch-type normalization decide whether ac4C exists in
    mRNA.** alnbase's "no filters; filter in SQL" stance lets a paper show both conclusions from one
    table.
- **Ranking:** popularity low-medium but **high narrative value** · comparator JACUSA v1-based, dated ·
  **Y**.

---

### 12. m7G-MaP-seq / m7G-quant-seq (m7G → abasic site → mixed mutations and deletions); TRAC-seq

- **(a)** NaBH4 (MaP) or KBH4 plus mild depurination (quant-seq) turns internal m7G into an abasic
  site. RT then gives **mixed G→A/C/T mutations plus deletions**, and "total variation ratio" is the
  level ([Enroth 2019 NAR](https://pmc.ncbi.nlm.nih.gov/articles/PMC6847341);
  [Zhang 2022 ACS Chem Biol](https://pubs.acs.org/doi/10.1021/acschembio.2c00792)).
  **TRAC-seq** instead cleaves at m7G (aniline), so the signal is a read 5′ end.
- **(b)** [jeppevinther/m7g_map_seq](https://github.com/jeppevinther/m7g_map_seq): R `getFreq2000.R`
  over `samtools mpileup` with BAQ-style error modelling. bowtie2 `--local -N 1 -D 20 -R 3 -L 15`.
- **(d)** Coverage ≥ 500 (rRNA) or ≥ 1500 combined (tRNA/mRNA). Control mutation < 1%.
  Log-likelihood ratio test with p < 1e-5. Signal = treated − control.
- **(e)** **Y**: `read="{ACT.}" refr="G"`. Adjacent insertions need `--insertions emit` plus a 2-column
  pattern. TRAC-seq read-start uses the `_` pad pattern.
- **(f)** [GSE121927](https://pmc.ncbi.nlm.nih.gov/articles/PMC6847341).
- **(g)** Whether insertions and deletions count toward "variation" differs between MaP and quant-seq.
  Local alignment clips terminal mutations.
- **Ranking:** low popularity · small script comparator · **Y**.

---

### Cross-cutting candidate extensions (this family)

| Extension | Needed by | Why not DuckDB/samtools |
|---|---|---|
| `--library rna-forward` / `rna-reverse` (R1 sense / R1 antisense) | DART (dUTP bulk), RNA-BS (meRanTK `rf`), BID/PRAISE on reverse-stranded libs, ac4C/RedaC:T | Walk orientation decides the context and read order of multi-column patterns (RAC, CCG, NNΨNN). Hand-complementing works but is fragile. |
| Unstranded / reference-orientation walk | Bullseye-exact reproduction, unstranded RNA-seq | Under `directional`, reads of one RNA strand are split between two walk orientations, so one pattern cannot collect them. |
| `--library tag:YZ` (strand from aux tag) | hisat-3n-table parity (GLORI-alt, eTAM, UBS-seq) | hisat-3n-table orients by YZ, not flags. The alnbase doc comment already anticipates this. |
| (Not needed) per-read conversion counts, read-end trimming, BQ, dedup, gap redistribution, calibration | GLORI, eTAM, meRanTK, UBS, BID, PRAISE | Per-`record_id` aggregation, `off_5p/off_3p`, `qual`, `^`-recorded context and samtools filters cover them. Emit soft clips to match whole-SEQ counts. |
| (Not needed) RT stop / cleavage site | mim-tRNAseq, TRAC-seq, HAC-seq | Leading `_` pad pattern. Document that `off_5p==0` is wrong when the 5′ end is soft-clipped. |

### Summary ranking for this family (versatility-paper value)

1. **GLORI (+eTAM via hisat-3n-table)**: top m6A method, trivial queries, A-count/CR filters likely to
   disagree, spike-in ground truth.
2. **BID-seq (+PRAISE)**: deletion signal (not substitution), U-tract placement is a real,
   literature-acknowledged ambiguity, calibration probes in the repo.
3. **DART-seq / Bullseye**: tiny in-repo test BAMs, context motif query, clear comparator weaknesses
   (unstranded, no BQ, mate double-count). Motivates `rna-reverse`.
4. **ac4C/RedaC:T**: best "analysis choice flips biology" story (duplicates and replicate handling), using
   published data.
5. **mim-tRNAseq (m1A/m3C RT signatures)**: demonstrates mismatch-spectrum plus RT-stop (pad) queries;
   bundled data.
6. **RNA-BS m5C (meRanTK/UBS)**: real mate-overlap bug, but closest to bisulfite, so least versatility.
7. m6A-SAC-seq, BACS, m7G-MaP, RBS-seq: optional breadth.


---

## Family 04 — Structure probing by mutational profiling (SHAPE-MaP, DMS-MaPseq)

Sources were read from shallow clones made on 2026-09-16 under `clones/`: `shapemapper2` (v2.3.2 tree), `seismic-rna` (0.26.0), `dreem` (archived), `RNAFramework` (2.9.7), and `RNAFramework-docs`. File and line references point into those clones.

### How the signal works (all tools)

A chemical probe leaves an adduct on flexible nucleotides. SHAPE reagents (1M7, NAI, 2A3) acylate the 2'-OH of all four nucleotides. DMS methylates A-N1 and C-N3 strongly, and G-N1 and U-N3 when the buffer is right (Bicine, pH about 8). DMS also methylates G-N7, which says nothing about structure. Under "MaP" conditions (Mn2+ with SuperScript II, or TGIRT-III, or MarathonRT) the reverse transcriptase reads through the adduct and writes a **mismatch, deletion or insertion** into the cDNA at or next to that nucleotide.

- **What alnbase has to compare:** the read against the reference at every aligned position. Each position is a match, a substitution by type, a deletion of a reference base, or an insertion next to it. Reactivity is (mutations / effective depth) in the treated sample, minus the untreated rate, optionally divided by the denatured rate.
- **Strand:** the RNA is almost always the + strand of a small target such as an amplicon, rRNA or transcript. "3'" and "5'" in every tool mean higher and lower reference coordinate. For truncation methods (RT-stop), the stop is "read start - 1" on + strand reads.
- **Context:** none for SHAPE. For DMS the relevant context is base identity (A/C, plus G/U with filtered substitution types).

---

### ShapeMapper 2 (SHAPE-MaP; eDMS-MaP via `--dms`)

#### (a) Signal
- **What counts:** every mismatch, deletion and insertion, each classified as `AG`, `A-`, `-A`, `multinuc_deletion`, `complex_insertion`, and so on (`docs/file_formats.md`).
- **Adduct site:** the "reference position immediately 5' of the last unchanged reference nucleotide before a mutation, scanning 3'->5'" (`docs/analysis_steps.md`). In code this is `m.right - 1`, the 3'-most position the (possibly collapsed) mutation covers (`internals/cpp-src/src/MutationProcessing.h` about lines 730-790). An insertion (`right = left + 1`) is therefore credited to the reference base on its **5'** side.
- **Ambiguous indels:** during parsing, an indel with alternative placements is merged into a mutation spanning the whole ambiguous window. It is then **shifted 5' (left-aligned) by default** (`shiftAmbigIndels`, MutationProcessing.h:46). `--right-align-ambig[-dels|-ins]` reverses this. Adducts in a repeat are thus credited to the 3' edge of the left-aligned indel. The documented example is a 5-nt deletion in `...GACGTCAAGTCATC...`.
- **Multinucleotide collapse:** mutations separated by **fewer than 6 unchanged nucleotides** (`--min-mutation-separation 6`, passed to C++ as `--max_internal_match 5`) are merged into one event and credited at the 3'-most position. The positions the merged event covers are removed from effective depth, except the adduct site (`collapseMutations`, MutationProcessing.h:325; analysis_steps "Effective read depth").
- **`--dms` mode (v2.2+, Mitchell et al. 2023):**
  - It counts only mismatches. Indels are dropped at all nucleotides.
  - G->A mismatches (N7-G) are routed to a separate `*.txtga`/`*.mutga` channel.
  - `G_multinuc_mismatch` is excluded. The whitelist is `dms_tags` at MutationProcessing.h:621-627.
  - Normalization is per nucleotide.

#### (b) Version and install
- **Version:** v2.3.2 (2026-09-02; the tarball fix). v2.3.0/2.3 (2024-11 and 2026-03) added `--N7` and primer-adjacent masking. See https://github.com/Weeks-UNC/shapemapper2/releases.
- **Bioconda: none.** The `bioconda/shapemapper` package is ShapeMapper **1.2**, a different program.
- **Install:** the release tarball bundles Linux x86-64 binaries (bowtie2 2.3.4.3, STAR 2.5.2a, BBmerge 37.78, and the C++ parser and counter). Building from source needs boost and cmake (issues #68, #69). Moderate difficulty, Linux only.

#### (c) Output formats
- **`<name>_<RNA>_profile.txt`** (tab-separated, header row):
  - Columns: `Nucleotide` (**1-based**), `Sequence` (AUGC, lowercase = masked), then per sample `<S>_mutations`, `<S>_read_depth`, `<S>_effective_depth`, `<S>_rate`, `<S>_off_target_mapped_depth`, `<S>_low_mapq_mapped_depth`, `<S>_mapped_depth` (or `_primer_pair_<n>_mapped_depth`).
  - Then `Reactivity_profile`, `Std_err`, `HQ_profile`, `HQ_stderr`, `Norm_profile`, `Norm_stderr`.
  - `*_profile.txtga` holds the same for N7-G.
- **`.shape`:** 2 columns (1-based position, normalized reactivity; -999 = no data). **`.map`** adds stderr and sequence. There are also `_varna_colors.txt` and `_ribosketch_colors.txt`.
- **`--output-parsed-mutations` → `<name>_<sample>_<RNA>_parsed.mut`:** one line per read (**the best per-read comparison target**). Tab-separated fields:
  1. read type (`PAIRED_R1`, `MERGED`, `PAIRED`, ...)
  2. read name
  3. 0-based leftmost mapped position (inclusive)
  4. 0-based rightmost mapped position (inclusive)
  5. mapping category (`INCLUDED`/`LOW_MAPQ`/`OFF_TARGET`)
  6. primer pair index or -999
  7. mapped-depth 0/1 string
  8. effective-depth 0/1 string
  9. mutation-count 0/1 string
  10. mutations: space-separated groups of 5 fields, namely the 0-based nearest unchanged position on the left, the same on the right, the quoted read sequence replacing the target between them, the quoted quals, and the quoted class, with `_ambig` appended when relevant.
- **`--output-counted-mutations` → `_mutation_counts.txt`:** one column per mutation class, plus `read_depth`, `effective_depth`, `off_target_mapped_depth`, `low_mapq_mapped_depth`, `mapped_depth`. Rows run 5'->3'.
- **Downstream consumers:** DanceMapper (ensembles), RingMapper and PairMapper read `parsed.mut` plus `profile.txt` ("ShapeMapper2 should be run with the --output-parsed-mutations option", https://github.com/MustoeLab/DanceMapper). SuperFold and RNAstructure `Fold -sh` read `.shape`.

#### (d) Implicit filters and defaults (`internals/python/pyshapemap/pipeline_arg_parser.py:186-251`)
- **Read trimming before alignment:** a window of 5 (`--window-to-trim`) with mean Q < 20 (`--min-qual-to-trim`) cuts the rest of the read. Reads shorter than 25 after trimming are dropped (`--min-length-to-trim`).
- **Mates:** BBmerge `vstrict=t`. Unmerged pairs go through post-alignment merging when the fragment is 800 or less (`--max-paired-fragment-length`). In overlaps the higher quality is kept; for conflicting mutations, the mutation group with the higher mean phred over the mutation and its adjacent bases wins.
- **Aligner:** bowtie2 `--local --sensitive-local --ignore-quals --mp 3,1 --rdg 5,1 --rfg 5,1 --dpad 30 --maxins 800`. The docs warn that with bowtie2, heavily mutated reads get lower MAPQ.
- **MAPQ:** `--min-mapq 10`. The Python layer passes this; the C++ standalone default is 30, at MutationParserExe.cpp:60.
- **Primers and read ends:**
  - **Amplicon mode:** read ends must lie within ±10 of a primer pair (`--max-primer-offset 10`). Mutations overlapping primers are dropped and primer sites removed from depth.
  - **Changelog vs code:** the v2.3 changelog says the 3 nt next to each primer are masked. The code applies that **only when `dms` is true** (`stripPrimers`, MutationProcessing.h:220-315). The same function also prints the debug strings "Engaging dms logic…" and "This should never ever print" for every amplicon read.
  - **Non-amplicon mode:** `--exclude_3prime` is `random_primer_len + 1`, so **1 nt is excluded even at the default `--random-primer-len 0`** (components.py:826-827). It applies only to the reference-right end of reverse-strand R1, reverse-strand R2, and merged/paired reads (`trimRightEnd`, MutationProcessing.h:1180-1300). Forward R1 is never trimmed.
- **Basecall quality (`--min-qual-to-count 30`):**
  - An unmutated position counts toward depth only if **it and its immediate 5' and 3' neighbours** have Q30 or better.
  - A mutation is dropped (and its span removed from depth) if any basecall inside it is below Q30, or if the nearest basecall on either side is. Across a gap, the next real basecall is used (`filterQscoresCountDepths`, MutationProcessing.h:485-790).
  - Collapsed mutations include the internal *matched* bases' quals.
- **Reactivity:**
  - `--min-depth 5000`: positions below it are excluded from HQ/Norm.
  - `--max-bg 0.05` (0.02 with `--dms`): an untreated rate above this is excluded.
  - Lowercase positions are masked.
  - Normalization is boxplot style (drop outliers above max(1.5×IQR, 90th/95th pct), divide by the mean of the top 10%), done jointly across RNAs unless `--indiv-norm`.
- **Other defaults:** no duplicate removal. Multimappers are handled only through MAPQ (STAR runs with `--outSAMmultNmax 1`). `--correct-seq` (off by default) rewrites the reference at variants above 60%.

#### (e) Expressibility in alnbase
- **Raw per-position events: yes, today.**
  - Mismatch: `read="/" refr="N"`, or by type, e.g. `read="G" refr="A"`.
  - Deletion: `read="." refr="N"`, which gives one column per deleted base.
  - Insertion: `--insertions emit`, then `read="N" refr="."`.
  - Coverage: `read="N" refr="N"`.
  - Neighbour qualities: capture with `^` marks, e.g. `mark="^+^"`, then filter in DuckDB.
  - Collapse (<6 unchanged nt, credit to the 3'-most event), removal of the collapsed span from depth, and the primer/3'-end exclusions (`off_3p`, `refr_pos`) are all per-read aggregations. They can be done in DuckDB with window functions over `record_id` ordered by `refr_pos`, so they do not count as needed features.
- **DMS mode:** A/C/U mismatches with `read="/" refr="H"`. G->C/T with `read="{CT}" refr="G"`. G->A as a separate query. Yes, today.
- **Gaps that need an extension:**
  1. **Ambiguous indel placement (needed).** The aligner's placement is arbitrary. ShapeMapper shifts 5' and credits the 3' edge. A fixed-width IUPAC pattern cannot find the left-aligned equivalent at a variable distance, and there is no column-to-column equality on the reference row. Minimal fix: a walk option `--indels aligner|left|right` that normalizes each gap (and insertion) to its leftmost or rightmost equivalent placement **in walk orientation** before patterns match. Better still, add a code such as `?` or a column flag meaning "this column is inside the equivalence window of an indel", so SEISMIC's "mark all equivalent positions ambiguous" can be written too. An upstream approximation is GATK `LeftAlignIndels`, which is left-only, needs GATK, and is not samtools.
  2. **Neighbour-quality predicate inside patterns (convenience).** ShapeMapper's depth rule, "this base and both neighbours Q30 or better", can be done with `^` captures, but it triples the coverage table. A per-column quality row in the pattern (e.g. `qual=">>>"`, or `min_qual=30` over marked columns) would express it at the source. Across a deletion, the "neighbour" is the next *non-gap* column, which a fixed pattern handles only for bounded deletion lengths.
  3. **Orientation mode.** Under `--library directional`, R1-reverse and R2-forward reads are walked reverse-complemented. For an amplicon whose adapter-tailed forward primer carries Read 1, every pair is R1-fwd/R2-rev and walks + strand, so this works today. Nextera or random-primer libraries, and anything with pairs in both orientations, would produce complemented bases and a flipped "3'-most" rule. That needs a `--library` value that walks every read in reference (+) orientation, e.g. `forward`/`unstranded-top`, or a strand predicate in `where`.
  4. **Mate merging policy.** `alnbase overlap` exists, but ShapeMapper's rule (take the mutation group from the mate with the higher mean phred over mutation plus flanks) is not one of its `--mismatch-qual` options. This is an overlap-command extension, not query syntax.

#### (f) Test data
- **In repo, `example_data/`** (4.6 MB): TPP riboswitch amplicon (`TPP.fa`, primers lowercase). Three samples (TPPplus, TPPminus, TPPdenat), each with 4 chunks × R1/R2 of 2,000 reads, i.e. **8,000 read pairs per sample**. Run with `run_example.sh` (`--amplicon`).
- **In repo, `internals/test/data/`** (28 MB total):
  - `ribosome_plus`/`ribosome_minus`/`ribosome_denat`: **15,000 read pairs each**, about 2.7 MB each, from E. coli 16S/23S 1M7.
  - `16S.ct` and `23S.ct` hold accepted secondary structures, used by `ROC_tests.sh` with `--random-primer-len 9`. This is **structural ground truth** for ROC/AUC.
  - `TPP_with_gap.fa` and `TPP_with_insert.fa` are variant-correction tests.
- **C++ unit tests:** `internals/cpp-src/test/testMutationParser.cpp` has 79 tests with hand-built SAM lines and expected mutations, including `IdentifyAmbiguousMutations.*` (left/right-aligned ambiguous gaps, gaps next to mismatches, ambiguous gaps near read ends). These are ready-made golden cases for alnbase indel normalization.
- **Public data:** Busan & Weeks 2018 E. coli rRNA 1M7 data, distributed via the Weeks lab site (no SRA accession found). eDMS-MaP (Mitchell et al. 2023, NAR 51:8744) compares SSII, MarathonRT, TGIRT-III and eHIV on the same RNAs; the accession is in the supplement, not located here.

#### (g) Known issues and contestable choices
- **5' shift of ambiguous indels.** It was chosen empirically ("empirical improvements", Busan & Weeks 2018 SFig 2). SEISMIC and rf-count do otherwise (see below). In a homopolymer or repeat at a helix/loop boundary, the three conventions credit the signal to different nucleotides.
- **Collapse window of 6.** Optimized on one SSII/1M7 rRNA dataset; the docs say other RTs "may require re-optimization". Issue #4: `--min-mutation-separation 0` still collapsed adjacent mutations. Two genuinely reactive neighbours within 5 nt become one, **under-reporting the 5' one**, which can make a flexible loop look partly protected.
- **Quality-filter asymmetry.** A Q29 neighbour removes a real mutation *and* the depth there, so the rate at positions next to low-quality cycles is estimated from a biased subset.
- **`--dms` removes all indels.** Mitchell 2023: SSII makes 26% indels, Marathon 2.8%, TGIRT 7.4%. For an SSII DMS dataset this throws away about a quarter of the signal. Keeping G->A in the main profile makes Gs look reactive regardless of pairing (N7-G), which **flips G loop/stem calls**. That is exactly why ShapeMapper `--dms` filters G->A and SEISMIC masks G by default, while rf-count counts G->A unless `-om` is given.
- **The documented 3-nt primer-adjacent masking is applied only in `--dms` code paths**, plus a stdout debug print for every read in non-DMS amplicon mode.
- **External alignments.** The parser requires an MD tag and rejects `=`/`X` CIGARs (issue #56: minimap2 input fails with "CIGAR string incorrectly formatted" or "MD tag does not match CIGAR"). alnbase reads the reference directly and handles both, which is a concrete versatility point.
- **Other issues:** #2, the histogram "median" rate was actually about the 5th percentile. The docs note bowtie2 MAPQ depends on mutation count, so `--min-mapq` interacts with modification level.

---

### SEISMIC-RNA (successor to DREEM; DMS-MaPseq, SHAPE-MaP, ETC)

#### (a) Signal
- **Relation vectors:** each read gets one byte per reference position. The bits are MATCH=1, DELET=2, INS_5=4, INS_3=8, SUB_A=16, SUB_C=32, SUB_G=64, SUB_T=128, NOCOV=255, IRREC=0 (`core/rel/code.py`).
- **Low-quality bases** (Phred below `--min-phred 25`) and N become the OR of all possible relations, e.g. `ANY_N ^ SUB_ref` (`idmut/py/encode.py`).
- **Ambiguous indels (`--ambindel`, default on):** every position an indel could slide to gets `DELET|MATCH`, or the insertion bits (`idmut/py/ambindel.py`; docs `algos/ambindel.rst`). A position counts toward `Mutated` or `Matched` only if its byte **definitively** fits a yes or no pattern (`RelPattern.fits`, `core/rel/pattern.py:337-342`). **An ambiguous deletion therefore contributes to neither numerator nor denominator anywhere in its window.**
- **Insertions** are credited to the base on their **3'** side by default (`--insert3`).
- **Probe defaults** (`--probe`, default **DMS**; `core/arg/cli.py:30-100`, `filter/main.py:102-110`):
  - DMS: G and U masked, poly(A) runs of 5 or more masked, reads with two mutations closer than 4 nt **dropped** (`--min-mut-gap 4`, `--mut-collisions drop`), and observer-bias correction on (`--quick-unbias`).
  - SHAPE/ETC: mutations closer than 2 **merged**, keeping the 3'-most.
  - ETC: masks A and C.

#### (b) Version and install
- **Version:** 0.26.0 (2026-08-26). On PyPI and **bioconda (`seismic-rna` 0.26.0)**. Needs Python 3.13, bowtie2, fastp, samtools and RNAstructure (for fold). The conda install is easy.
- **Citation:** Allan et al. 2024, "Discovery and Quantification of Long-Range RNA Base Pairs in Coronavirus Genomes with SEARCH-MaP and SEISMIC-RNA", bioRxiv https://doi.org/10.1101/2024.04.29.591762 (PMC11326378).
- **Churn warning:** `relate`→`idmut` (0.25), `mask`→`filter`, and `ensembles`→`filterscan`/`clusterscan` (0.26) were renamed recently. Pin the version.

#### (c) Output formats
- **`idmut/{ref}/idmut-batch-N.brickle`:** Brotli-compressed pickled relation vectors, plus `idmut-report.json`.
- **`seismic table` → `{step}-position-table.csv`:** indexed by `Position` (**1-based**) and `Base`. Relationship columns: `Covered`, `Informative`, `Matched`, `Mutated`, `Subbed`, `Subbed-A`, `Subbed-C`, `Subbed-G`, `Subbed-T`, `Deleted`, `Inserted` (`core/table/base.py:29-52`). `Informative = Matched + Mutated` (`core/table/write.py:80`).
- **`{step}-read-table.csv`:** indexed by `Read Name`, with the same relationship columns counted per read. Clustered runs also write `cluster-abundance-table.csv`.
- **Other exports:** `seismic export` writes JSON for the "SEISMICgraph" web app. Fold writes `.ct`/`.db` plus a VARNA color file.

#### (d) Implicit filters and defaults
- **`align`:**
  - bowtie2 **local** mode, `--bt2-L 20`, `--bt2-gbar 4`, `--bt2-dpad 2`, `--bt2-X 600`, `--bt2-score-min-loc L,1,0.8`.
  - No discordant, mixed or dovetail alignments.
  - `--min-mapq 25`.
  - `--min-reads 1000` per reference (below that, the BAM is skipped).
  - fastp trimming.
- **`idmut`:**
  - `--min-phred 25`.
  - `--clip-end5 4`/`--clip-end3 4`: every read end ignores 4 positions, tied to the bowtie2 gbar.
  - `--overhangs` on.
  - Mates are merged by **bitwise AND** (`idmut/py/idmut.py:494-556`). Match AND substitution gives 0, "irreconcilable".
  - Separate BAM per reference, named after it.
  - `--sep-strands` off, so reverse-strand reads are mixed in as the same RNA.
- **`filter`:**
  - `--min-finfo-read 0.95`: drop reads with less than 95% informative positions.
  - `--max-fmut-read 1.0`: effectively off.
  - `--min-ninfo-pos 1000`: mask positions with fewer informative reads.
  - `--max-fmut-pos 1.0`.
  - `--drop-discontig`: drop pairs whose mates don't overlap or abut.
  - Probe masks as in (a); iterative until stable.
- **No duplicate removal.**

#### (e) Expressibility
- **Matched, Subbed-X and Deleted for unambiguous events: yes, today.** Use `read="/"` with a refr base set, `read="."`, and `read="="` for matched. Also expressible today:
  - Quality below 25 treated as "not informative": DuckDB on `qual`.
  - 4-nt end clipping: DuckDB on `off_5p`/`off_3p`, except these count SEQ including soft clips, so pull `cigar` with `-f` or add aligned-end offsets.
  - Read-level drops (min-finfo, min-mut-gap drop, discontiguous mates): DuckDB aggregation.
  - Observer-bias correction is downstream math.
- **Insertion 3' anchoring:** `--insertions emit` plus a capture of the next column. Yes.
- **Needs extension:**
  - The **ambiguous-indel window** ("mark all equivalent placements as neither match nor mutation") needs the same `--indels`/equivalence-window code proposed above. This case needs the *window*, not a single normalized placement.
  - **Mixed strands** need the reference-orientation library mode.
  - **Bitwise-AND mate merging** is an `overlap` policy (discard both calls on conflict and keep the high-quality mate's call when the other is low quality).

#### (f) Test data
- `src/userdocs/tutorials/amplicon/data/`: `hiv-rre.fa` plus simulated FASTQs `dms1` (about 23.9k pairs, 1.4 MB), `dms2` (1.3 MB) and `nodms` (0.17 MB). Read names are `batch-0_read-0`, quality is `I`/`!`, so the data is **simulated**.
- **`seismic sim total`** (`ref`, `fold`, `params`, `ends`, `muts`, `fastq`) generates FASTQs with *known per-position mutation rates* and cluster mixtures, including a DMS `min_mut_gap` bias model and SHAPE "injected" 5' bleed-through mutations. **This is directly usable as ground truth.**
- **Real data:** Zubradt et al. 2017 DMS-MaPseq, GEO **GSE84537** (large). Tomezsko et al. 2020 (HIV RRE, DREEM; Nature 582:438) is the classic ensemble dataset.

#### (g) Known issues and contestable choices
- **Ambiguous indels are dropped rather than placed.** This is the opposite of ShapeMapper (5' shift, counted) and rf-count (3' shift, counted). In G/C-rich repeats with SSII-type deletion signal, SEISMIC reports *lower* reactivity and *lower* informative depth there. SEISMIC can call a repeat loop unreactive where ShapeMapper calls it reactive.
- **Drop-read-on-close-mutations for DMS.** It removes exactly the molecules with adjacent reactive A/Cs (open loops), biasing rates downward in loops. SEISMIC compensates with its observer-bias correction; DREEM did not originally, and other tools don't. Cluster proportions (e.g. the 6% minor RRE structure detection) depend on this choice.
- **G/U masked by default for DMS.** Any comparison with ShapeMapper `--dms` (four-base) must switch it off. A user who leaves the defaults loses all G/U structure information.
- **Rapid API/CLI renames** (0.25/0.26) make version pinning essential for a reproducible benchmark.

---

### DREEM (obsolete; for reference only)
- https://github.com/rouskinlab/dreem. README: "DREEM is obsolete and has been replaced with SEISMIC-RNA". Last commit 2024-10.
- **Defaults:** `DEFAULT_MIN_PHRED = 25` (`dreem/util/cli.py:15`). Bit-vector clustering preprocessing: `signal_thresh=0.005`, `include_gu=False`, `include_del=False` (**deletions excluded from clustering**), `min_mut_dist=4`, `min_reads=1000` (`dreem/cluster/bitvector.py:25`).
- **Recommendation:** compare against SEISMIC-RNA instead. Cite DREEM only for history (Tomezsko et al. 2020, https://www.nature.com/articles/s41586-020-2253-5).

---

### RNA Framework `rf-count -m` (SHAPE-MaP / DMS-MaPseq; also RT-stop mode)

#### (a) Signal
- **Mismatches** are parsed from the **MD tag**. A mismatch counts if its base is Q20 or better (`-q 20`), and neighbours ±1 are also checked only with `-es` (`rf-count` lines 1826-1832).
- **Deletions** (up to `-md 10` nt):
  - A ±10-nt window is slid to test ambiguity (lines 1750-1800).
  - Ambiguous deletions are re-aligned to the **right-most (3')** valid placement by default; `-la` gives left, `-na` drops them.
  - **All deleted bases** are marked by default; `-rd`/`-ld` mark only the 3'/5' base.
  - **No quality check on deletions**: the `-es` block is commented out at lines 1781-1784 and 1800-1803.
- **Insertions** come from the CIGAR and are credited to the reference base 5' of the insertion (`push(@ins, $last - 1)`, line 1886). There is no ambiguity handling and no quality check.
- **Collapse (`-cc`, off by default):** consecutive events within `-mc 2` nt are collapsed to the 3'-most (`collapsemutations`, line 1939). `-dc N` discards them instead.
- **`-om`** restricts to substitution types given in IUPAC form (e.g. `A2T;C:N`) and disables indels.
- **RT-stop mode (no `-m`):** counts at read start - 1, **excludes reverse-strand reads** (line 252), and excludes soft-clipped reads unless `-ic`. The stop is added to coverage at that base.

#### (b) Version and install
- **Version:** 2.9.7 (2026-06-03). **bioconda `rnaframework` 2.9.7.** Perl plus samtools; R for plots; bowtie/bowtie2 for rf-map. Easy.
- **Citation:** Incarnato et al. 2018 NAR 46:e97, https://academic.oup.com/nar/article/46/16/e97/5035169. Docs: https://rnaframework-docs.readthedocs.io/en/latest/rf-count/.

#### (c) Output formats
- **RC binary**, one entry per transcript: `uint32 len_id`, `char[] id\0`, `uint32 len_seq`, a 4-bit packed sequence, `uint32[len] counts`, `uint32[len] coverage`, `uint64 mapped_reads`. It ends with a `uint64 total_reads`, a `uint16 version=1` and the marker `[eofrc]`, all little-endian. There is an `.rci` index.
- **`rf-rctools view`:** 4 lines per transcript (ID; sequence; comma-separated per-base counts; comma-separated per-base coverage). The arrays are indexed from 0. `-t` gives tabular output.
- **`-orc`:** raw counts by class (AC ... TG, ins, del). **`-mm`:** MM binary (per read: start, end, mutation indices) for DRACO.
- **`rf-norm`:** an XML per transcript (`<transcript id length><sequence/><reactivity>csv</reactivity>`) with scoring (Ding/Rouskin/Siegfried/Zubradt) and normalization (2-8%, 90% winsor, box-plot) attributes.
- **Mask file:** ranges are **0-based inclusive**.

#### (d) Implicit filters and defaults (lines 110-125, 248-256, 1183-1190)
- **samtools `-F UNMAP,QCFAIL,DUP`:** **duplicates are discarded by default** (`-ndd` keeps them). `-q 0` means no MAPQ filter. Secondary alignments are kept unless `-pn`.
- **Read-level drops:**
  - `-me 0.15`: drop a read if edit distance / aligned length exceeds 0.15. Consecutive deleted or inserted bases count once.
  - `-eq 20`: median read quality below 20 over aligned bases.
  - `-ds 1`: minimum aligned length.
- **Paired end (since 2.9.5):** mates are joined. In the overlap, **only mutations present in both mates are kept** (`-pam` keeps either; `-fsr` counts both mates separately).
- **No read-end trimming** in `-m` mode (`-t5` only affects RT-stop mode). No primer masking except via `-mf`. Collapse is off.

#### (e) Expressibility
- **Mismatch by type with Q20 or better:** yes.
- **Deletions, all bases marked:** yes, via `read="."`.
- **Right-most or left-most base only (`-rd`/`-ld`):** yes with fixed-width patterns for bounded deletion length, e.g. gap followed by non-gap. This holds only once the placement is normalized.
- **Insertion credited to the 5' base:** yes (`--insertions emit` with a capture).
- **Edit-distance and median-Q read drops, both-mates-agree:** DuckDB. Both-mates-agree can also be an `overlap` policy.
- **Needs extension:** ambiguous-deletion **right** alignment (the `--indels right` walk option). The same options cover `-la` and `-na`. rf-count's ±10-nt window limit is itself a quirk alnbase need not copy.
- **RT-stop mode:** expressible today in parquet. The walk's leading flank is a pad column carrying the reference base 1 nt before the first aligned base, so `read="_N" refr="NN"` anchored on column 0 hits exactly "start - 1". Anchors on pads are refused only for `bases` tags. rf-count's exclusion of reverse reads and of clipped reads is DuckDB (`flag`, `cigar`).

#### (f) Test data
- No bundled test reads (`tests/test-deps.sh` only).
- rf-count accepts any BAM, so the ShapeMapper example/rRNA data and SEISMIC's simulated data (once aligned) serve as well.

#### (g) Known issues and contestable choices
- **3' placement of ambiguous deletions is the opposite of ShapeMapper's 5' default**, and all deleted bases are marked. For a 2-nt deletion in `CAGAGU`, ShapeMapper credits one nucleotide while rf-count credits two nucleotides further 3'.
- **No quality filter on deletions or insertions,** while mismatches need Q20. Indel-rich SSII data on a poor run inflates indel signal relative to substitutions.
- **`-me 0.15` read drop** discards the most heavily modified molecules (high-dose in-cell probing, long reads). Reactivity is biased low exactly where it is highest. The same concern applies as with SEISMIC's drop policy, but uncorrected.
- **Duplicates dropped by default,** while the other tools keep them. With UMI-less amplicons, "duplicate" flags are mostly false and remove real depth.
- **Mate-agreement requirement** halves effective sensitivity in overlaps when one mate is low quality. Meanwhile, coverage counts the overlap once.

---

### Downstream tools (brief)
- **RNAstructure `Fold -sh file.shape`, SuperFold (Weeks lab), ViennaRNA:** consume normalized reactivities. No read-level logic, so they are not direct comparators.
- **DanceMapper / RingMapper / PairMapper (Mustoe/Weeks labs):** consume ShapeMapper `parsed.mut`. An alnbase exporter to `parsed.mut` would plug alnbase into ensemble and correlation analysis.
- **DRACO (Incarnato lab):** consumes rf-count `-mm` MM files. An MM exporter is simple (start, end, mutation indices per read).
- **SEISMIC `cluster`/`clusterscan`:** consume its own brickle batches, so there is no easy injection point short of `seismic importmm`, which exists (`src/seismicrna/importmm`, not inspected in detail).

---

### Logistics for a validation experiment
1. **Isolate counting from alignment.**
   - Run ShapeMapper with `--output-aligned-reads --output-parsed-mutations --output-counted-mutations` and feed the *same* SAM/BAM to `shapemapper_mutation_parser` (documented as possible for external SAM, with an MD tag), SEISMIC `idmut` (one BAM per reference, named after it), `rf-count -m`, and alnbase.
   - Differences then come only from counting policy.
   - Keep in mind ShapeMapper's own pre-alignment trimming and BBmerge step; for a pure comparison use its parser on unmerged BAMs.
2. **Per-read diff:** ShapeMapper `parsed.mut` field 10 (mutation groups with 0-based flank coordinates) against alnbase hit rows grouped by read name. Per-position diff: `profile.txt` (1-based) against SEISMIC `position-table.csv` (1-based) against `rf-rctools view` (0-based arrays) against alnbase `refr_pos` (0-based).
3. **Expected headline discrepancies:**
   - Ambiguous indel attribution: 5' vs 3' vs dropped.
   - Insertion side: 5' vs 3'.
   - Collapse window: 6 vs 2 vs off.
   - Neighbour-Q30 vs Q25 vs Q20 mismatch-only.
   - Mate overlap policy: max-mean-phred vs AND vs both-mates-agree.
   - Duplicate removal (rf-count only).
   - G->A handling in DMS.
   - The 1-nt `exclude_3prime` default in ShapeMapper.
4. **Ground truth:**
   - `seismic sim total` for known rates.
   - ShapeMapper rRNA test data plus 16S/23S `.ct` for ROC/AUC (accessible = unpaired).
   - Spike-in via mutating reads on a known alignment is also trivial.
   - Test ambiguous indels by simulating deletions in homopolymers with known true positions, then checking which convention recovers the right nucleotide. By construction none can, which is exactly the point to show.

### Ranking inputs
- **Popularity:** high. SHAPE-MaP and DMS-MaPseq are the dominant RNA structure-probing readouts. ShapeMapper 2 (RNA 2018) and RNA Framework (NAR 2018) are each widely cited, and SEISMIC-RNA/DREEM anchor the ensemble literature.
- **Comparator cleanliness:** medium-high.
  - ShapeMapper's `parsed.mut` is per read, a rare and excellent diff target, but the pipeline is monolithic, not on bioconda, and has interlocking heuristics.
  - SEISMIC is on bioconda with clean CSVs, but its CLI churns.
  - rf-count is on bioconda and simple, with binary RC plus a text view.
- **Ground truth / simulation ease:** high. SEISMIC has a built-in simulator, ShapeMapper ships rRNA reads with reference structures, and the targets are tiny, so runs take seconds.
- **Discrepancy likelihood:** very high. Three mainstream tools make three different, documented choices on the same event (ambiguous deletion, insertion side, collapse, mate conflicts), and the ambiguous-deletion choice moves reactivity between adjacent nucleotides in repeats. Two of the choices (G->A in DMS, drop-reads-with-close-mutations) have documented biological consequences (Mitchell 2023; SEISMIC observer bias).
- **Expressibility:**
  - Mismatches, deletions, insertions, coverage, DMS base/type filtering and RT-stop: **Y today**. Collapse, end exclusion and read-level drops go to DuckDB, and neighbour Q can use `^` captures.
  - Ambiguous-indel placement or window: **needs `--indels left|right` normalization plus an "in indel-equivalence window" code**.
  - Mixed-orientation libraries: **need a reference-orientation `--library` mode** (fixed-orientation amplicons work today).
  - ShapeMapper, SEISMIC and rf-count mate-conflict rules: **need `overlap` policy options**.
  - A per-column quality predicate is a convenience, not a necessity.

### References
- ShapeMapper 2: https://github.com/Weeks-UNC/shapemapper2. Docs: `docs/analysis_steps.md`, `docs/file_formats.md`, `docs/dmsmode.md`. Issues #2, #4, #56: https://github.com/Weeks-UNC/shapemapper2/issues
- Busan S, Weeks KM (2018) RNA 24:143, https://rnajournal.cshlp.org/content/24/2/143.full
- Mitchell D et al. (2023) NAR 51:8744, https://academic.oup.com/nar/article/51/16/8744/7201944
- SEISMIC-RNA: https://github.com/rouskinlab/seismic-rna and https://rouskinlab.github.io/seismic-rna/. Allan et al. 2024 bioRxiv https://doi.org/10.1101/2024.04.29.591762 and https://pmc.ncbi.nlm.nih.gov/articles/PMC11326378/
- DREEM: https://github.com/rouskinlab/dreem. Tomezsko et al. 2020 Nature 582:438, https://www.nature.com/articles/s41586-020-2253-5
- Zubradt M et al. 2017 Nat Methods 14:75, https://www.nature.com/articles/nmeth.4057 (GEO GSE84537)
- RNA Framework: https://github.com/dincarnato/RNAFramework and https://rnaframework-docs.readthedocs.io/en/latest/rf-count/. Incarnato et al. 2018 NAR 46:e97, https://academic.oup.com/nar/article/46/16/e97/5035169
- Smola MJ et al. 2015 Nat Protoc (SHAPE-MaP, deletion realignment), PMID 26426499
- DanceMapper: https://github.com/MustoeLab/DanceMapper. DRACO: https://github.com/dincarnato/draco


---

## Family 05 — mC→T chemistries (TAPS family, Illumina 5-Base), 5hmC assays, and methyltransferase footprinting (NOMe-seq, dSMF, MAPit/SMAC)

Scope note. Source was read from local clones under `clones/` (as of 2026-09-16):
rastair 2.2.0 (github.com/bsbludwig/rastair, HEAD b4d44fa), asTair 3.3.3 (bitbucket.org/bsblabludwig/astair),
taps-foundry (github.com/watchmaker-genomics/taps-foundry), MethylDackel 0.6.1, Bismark 3.1.0 (the Rust suite; its NOMe code is byte-identical to Perl v0.25.1),
BISCUIT 1.10.3 (bioconda has 1.10.2.20260818), methylpy 1.4.7, SingleMoleculeFootprinting 2.7.0 and QuasR 1.53.1 (git.bioconductor.org).
Bioconda versions were checked through api.anaconda.org on the same day.

Relevant alnbase facts (from `src/tags.rs` and `src/alignment.rs`):
- With `--library directional`, R1 forward is OT, R1 reverse is OB, R2 reverse is CTOT and R2 forward is CTOB. OB and CTOB reads are walked reverse-complemented, so each read is seen in the orientation of the original strand, and the reference context is that strand's context.
- `directional` is the only library mode.
- Offsets (`off_5p`, `off_3p`) count SEQ positions, including soft-clipped bases.

---

### 1. TAPS / TAPSβ / CAPS / PS / hmC-CATCH, and Illumina 5-Base (mC→T chemistries)

#### (a) Signal
- **TAPS** (Liu et al. 2019, *Nat Biotechnol* 37:424): TET oxidises 5mC and 5hmC to 5caC, pyridine borane reduces that to DHU, and PCR reads it as T. A **modified C reads as T**, and an unmodified C stays C. This is the inverse of bisulfite, but the reference-vs-read substitution is the same C>T on the converted strand (G>A in bottom-strand SEQ).
- The variants differ only in which modification converts:
  - **TAPSβ** (βGT blocks 5hmC): 5mC only.
  - **CAPS** (KRuO4 turns 5hmC into 5fC): 5hmC only.
  - **PS** (pyridine borane alone): 5fC and 5caC.
  - These are in Liu et al. 2021, *Nat Commun* 12:618. **hmC-CATCH** is analogous for 5hmC.
- **Illumina 5-Base** also turns 5mC into T and leaves unmodified C intact. Bismark 3.x documents it as "the chemical inverse of bisulfite … paired-end and directional" (`Bismark/docs/src/content/docs/rust/illumina-5-base.md`).
- **Context.**
  - CpG is the default for rastair; CHG and CHH are allowed by asTair `--context all`.
  - rastair 2.2 counts evidence **pairwise**: an OT read counts as modified only as `TG` over a reference `CG`, and as unmodified only as `CG`. OB reads are the mirror, `CA`/`CG` in top coordinates (`rastair/src/metrics/methylation.rs:151-163`: `read_counts` keys on `adj` = the neighbouring read base).
  - rastair 2.2's changelog says: "The reported beta values might change … only reads containing TG/CA are counted as methylated."
- **Strand.** Libraries are directional, because conversion happens after adapter ligation.
  - rastair takes OT/OB from the FLAG in default paired mode (include flags 3, exclude 3852). `--unpaired` treats forward as OT. `--guess-read-orientation` infers OT/OB **per read** from TG-vs-CA mismatch motifs, for tagmented (non-directional) libraries (`docs/src/calling/strand-guesser.md`).
- **Library artefact.** End repair fills in with unmodified dC, so fragment ends look *unmodified*, which is the opposite polarity to bisulfite M-bias. Watchmaker's pipeline masks R2 5′ bases with `--nOT/--nOB 0,0,20,0` and asTair `--start_clip 7 --end_clip 7` (taps-foundry README).

#### (b) Comparators
| Tool | Version | Install |
|---|---|---|
| **rastair** | 2.2.0 (2026-08-24) | bioconda `rastair` 2.2.0, a Rust binary with an ML model bundled. Non-commercial licence. nf-core/methylseq uses it for TAPS (subworkflow `bam_taps_conversion`: `rastair mbias` → `rastair call` → `rastair methylkit`). Easy to install. |
| **asTair** | 3.3.3 | pip only, not on bioconda. Deprecated in favour of rastair. Python with pysam. Easy to install. |
| **MethylDackel + flip** | 0.6.1 | bioconda. Used in jknightlab/TAPS-pipeline: bwa-mem → MethylDackel extract `-q 10 -p 13 --mergeContext --OT/--OB` → swap meth/unmeth. Easy to install. |
| Bismark `--illumina_5base` | 3.1.0 | bioconda. For 5-Base only, concordance-gated against DRAGEN 4.4.6 (r≈0.99 over 55M CpGs). DRAGEN itself is proprietary. |

#### (c) Output formats
- **rastair `call` VCF/BCF.** INFO has `CPG`, `CPGnovo`, `AS_SB_OT`, `AS_SB_OB`, `M5mC_Strands` (4 values: unmod, mod, no-SNP, SNP), `PIR`, and so on. FORMAT has `GT:GL:GC:DP:M5mC:DPM5mC:ADM5mC:ML`.
  - `M5mC` can hold **two betas** at a site that is both a reference CpG and a de-novo CpG.
  - `ALT=.` at CpGs with no variant.
  - FILTERs include `m_vaf`, `m_bq_ratio`, `m_pos`, `m_highDp`, `low_ml_score` (`docs/src/formats/vcf*.md`).
- **rastair `call -c` BED ("mods").** Columns: `#chr start end name beta_est strand unmod mod no_snp snp coverage genotype gt_p_score gt_conf_score cpg`.
  - `start` is 0-based and `end` is 1-based (half-open).
  - `cpg` is `REF` or `NEW` (de novo).
- **rastair `per-read` BED.** Columns: `chr start end read_id mapq orientation insert_size read_length flag num_cpg num_mod mod_cpgs unmod_cpgs snp_cpgs mod_denovos unmod_denovos`.
  - Positions are in read coordinates relative to the first *aligned* base. `--count-clipped` counts from the first base of SEQ instead.
- **rastair `bam`.** Writes an MM/ML modBAM, rewriting SEQ so that T becomes C at methylated sites, or legacy `XM/XR/XG`.
- **asTair `.mods`.** Columns: `#CHROM START END MOD_LEVEL MOD UNMOD REF ALT SPECIFIC_CONTEXT CONTEXT SNV TOTAL_DEPTH`.
  - 0-based START, END = START+1.
  - One row per strand: the C on + gives `REF=C ALT=T`, the G on − gives `REF=G ALT=A`.
  - `SNV` is `No` or `homozygous` (a heuristic), or `WGS_known` when `--known_snp` is given.
  - `TOTAL_DEPTH` includes reads of the uninformative orientation.
  - There is also a `.stats` file per context.
- **MethylDackel** bedGraph: `chrom start end pct nMeth nUnmeth`, 0-based. For TAPS these must be flipped: pct′ = 100 − pct, and meth and unmeth swapped.
- **taps-foundry `save_as_methylkit.py`**: `chrBase chr base strand coverage freqC freqT`, where freqC means *modified*, the reverse of its bisulfite meaning.

#### (d) Implicit filters and defaults
- **rastair `call`** (`docs/src/cli.md`):
  - Reads and bases: `--min-mapq 1`, `--min-baseq 10`, `-f 3 -F 3852` (proper pairs only; unmapped, mate-unmapped, secondary, QC-fail, duplicate and supplementary reads dropped).
  - Overlapping mates are deduplicated unless `--keep-overlapping-reads`.
  - `--max-coverage 1000`.
  - Methylation thresholds: `--m-min-depth 3`, `--m-vaf-min 0.2`, `--m-bq-ratio-min 0.27`, and `--m-read-position-min/max 0.2/0.8` (the `m_pos` filter on alt evidence position).
  - `--nOT/--nOB 0,0,0,0`, with distances in read length including soft clips.
  - The ML filter is on by default (`--ml 0.5`); `--no-ml` gives "very similar output to other callers".
  - Genotype correction: at a heterozygous C/T, β = max(M − (M+U)/2, 0)/(M_excess + U).
  - `--rescue-soft-clip-cpg` is off by default.
- **asTair `call`** (`astair/caller.py:42-70, 235-268, 400`):
  - `--minimum_base_quality 20`, `--minimum_mapping_quality 0`.
  - **`compute_baq=True`**: htslib BAQ recomputation is on by default and lowers qualities near indels before the BQ20 cut.
  - `--max_depth 250`: the pysam pileup cap, which silently truncates amplicon, mtDNA and spike-in depth.
  - `ignore_orphans=True`, and pysam stepper `samtools` drops unmapped, secondary, QC-fail and duplicate reads.
  - Flags are matched *exactly*: C sites use 99/147 (unexpected 83/163), G sites the reverse. Single-end uses 0/16. Reads with any other bit set, e.g. 0x200 or 0x800, are silently ignored.
  - `ignore_overlaps=True`, which pysam implements as random removal of overlapping mate bases.
  - Clip `0/0`.
  - A `LowCov` label below depth 3.
  - `--library directional|reverse`.
- **MethylDackel** (`extract.c:726-746`, `common.c:84-116`):
  - `minMapq 10`, `minPhred 5`, `ignoreFlags 0xF00`, duplicates excluded.
  - Strand comes from FLAG unless there is a Bismark-style XG tag.
  - Counts C versus T only, **without checking the neighbouring read base**.
  - `--minOppositeDepth 0` means no SNP masking.

#### (e) alnbase expressibility — **Y** (directional). The polarity flip happens downstream.
```toml
[query.T_CG]      # MethylDackel/asTair-style single-position call (modified in TAPS)
read = "T~"
refr = "CG"
[query.C_CG]
read = "C~"
refr = "CG"
[query.TG_CG]     # rastair >=2.2 pairwise evidence
read = "TG"
refr = "CG"
[query.CG_CG]
read = "CG"
refr = "CG"
[query.snp_opp]   # opposite-strand genotype evidence at the CpG C (read on the other strand's walk shows A over its G)
read = "~A"
refr = "CG"
mark = ".+"
```
- **Opposite-strand SNP evidence.** asTair's `SNV` heuristic, MethylDackel `--minOppositeDepth/--maxVariantFrac`, and rastair's genotype adjustment are joins in DuckDB on `refr_pos` over rows from the other strand. That is aggregation, not a query feature.
- **Read-end masks** (`--nOT/--nOB`, `--start_clip`) are `WHERE off_5p >= k`. alnbase's `off_5p` counts soft clips like rastair does.
- **De-novo CpGs** (ref `CA`, SNP A>G): a per-read pattern `read="{TC}G" refr="C{AT}"` gives the evidence. Whether the site counts depends on a genotype call, which is aggregation.
- **Gaps:**
  1. **Non-directional/tagmented TAPS** needs a per-read strand choice. rastair infers it from TG-vs-CA mismatch counts; alnbase has no mode that either emits both orientations or infers strand. A `--library infer` mode is needed, or a tag-driven mode (below).
  2. `--rescue-soft-clip-cpg` needs `--soft-clips emit`. alnbase would also have to provide the reference base adjacent to a clipped fringe base, which should be checked in `walk_alignment`.
  3. asTair's BAQ is an upstream `samtools calmd -Ar` equivalent, so it does not count.

#### (f) Datasets
- **rastair repo `tests/data/`**, all in git:
  - `test.bam`: 1.6 MB, paired TAPS.
  - `test.fasta.gz`: 19.5 MB.
  - `ug_tagmented.bam`: 32 MB, tagmented and non-directional. Good for testing strand inference.
  - `rastair1.bed`: 0.5 MB, rastair v1 expected output.
- **asTair tutorial**: 1:1 methylated/unmethylated lambda TAPS, `zenodo.org/record/2582855` (`lambda.phage_test_sample_{1,2}.fq.gz`, `lambda_phage.fa`); expected about 48.2 % CpG. Needs alignment with bwa mem; sizes not verified because the Zenodo API timed out. There is also `astair/tests/test_data/` (hg38 chr10:20000-60000 plus a common-SNP VCF, about 1 MB).
- **Bismark 3.x CI**: a lambda/pUC19 5-Base spike-in gate (in repo `validation/`).
- **Illumina 5-Base NA12878 / HG002 demo data**: on BaseSpace, login required.
- **Ground truth.** The lambda 0 %/100 % spike-ins and the NA12878 GIAB genotype truth set test SNP confounding directly.

#### (g) Known issues and contestable choices
1. **C>T SNPs at CpG look methylated in TAPS** (the reverse of bisulfite).
   - MethylDackel-flip and asTair without `--known_snp` report a heterozygous CpG→TpG as about 50 % modified on OT. A homozygous one is 100 % on OT, while OB (G unchanged) reads 0 %.
   - Strand-merged β lands in between.
   - Comparing genotypically different samples (tumour/normal, individuals) produces spurious hypermethylated DMRs. rastair's excess-M correction and pairwise counting target exactly this.
   - The rastair preprint (Etzioni et al. 2026, bioRxiv 10.64898/2026.03.19.712983) reports about 500k de-novo CpGs in NA12878.
2. **Pairwise versus single-position counting changed rastair's own betas between 2.1 and 2.2.** A CpG whose G carries a G>A SNP (hence `TA` on OT) contributes to MethylDackel but not to rastair 2.2. alnbase expresses both (`T~` and `TG`), which makes this a clean "same data, two definitions" demonstration.
3. **asTair `max_depth=250`** silently caps depth. On high-coverage spike-ins or amplicons, reads beyond 250 are dropped in file order. It also has BAQ on by default, and `MAPQ 0` versus rastair's `MAPQ 1`.
4. **Flag-exact matching in asTair** (99/147/83/163) silently discards reads carrying 0x200 or 0x800 bits and mate-unmapped reads, even when they are otherwise good.
5. **End-repair unmodified fill-in** means that without an R2 5′ mask, fragment-end CpGs are biased toward "unmodified". Short cfDNA fragments suffer most, which lowers global cfDNA methylation estimates.
6. **Polarity inversion is silent.** Running a bisulfite caller on TAPS without flipping yields β′ = 1 − β with no warning. The bwa-mem3 docs report r = −0.956 against truth for exactly this mistake.

---

### 2. 5hmC: TAB-seq, ACE-seq, oxBS-seq (subtractive); 5-letter/6-base (biomodal duet)

#### (a) Signal
- **TAB-seq** (Yu et al. 2012, *Cell* 149:1368): βGT protects 5hmC, TET turns 5mC into 5caC, and bisulfite converts C, 5mC (now 5caC) and U to T. **Retained C = 5hmC.** The signal is the same C>T/retention readout as WGBS, in all contexts, on a directional library.
- **ACE-seq** (Schutsky et al. 2018, *Nat Biotechnol* 36:1083): βGT protects 5hmC, and APOBEC3A deaminates C and 5mC to U/T with no bisulfite. **Retained C = 5hmC**; analysis is bisulfite-like. A3A has sequence preferences (it favours TpC) and incomplete deamination of 5mC produces false 5hmC ("Failure to convert", *Nat Chem Biol* 2018).
- **oxBS-seq** (Booth et al. 2012, *Science* 336:934): oxidation turns 5hmC into 5fC, which bisulfite converts to U. **Retained C = 5mC only.** 5hmC = BS − oxBS, which is subtractive across two libraries.
- **CAPS / hmC-CATCH** are the mC→T equivalents (section 1).
- **biomodal duet +modC / evoC / 6-base.**
  - Hairpin ligation joins each original strand to a copy strand. Deamination follows, and 5mC/5hmC are copied or protected differently on the two strands.
  - The modification is resolved **from the pair of bases (R1 original, R2 copy) at each position, before alignment**: "the nth base of read 1 is paired with the nth base of read 2", with a pairwise-alignment fallback.
  - The resolved single-end reads carry true genotype bases, and the modification sits in `MM` tags (`C+m`, `C+h`, `C+C?`), `XR:i` (resolution method) and `XL:i` (fragment length). Source: biomodal Data Interpretation Guide 1.4.1, software-docs.biomodal.com.
  - **After resolution there is no read-vs-reference mismatch to query.**

#### (b) Comparators
- TAB-seq and ACE-seq use **Bismark** or **MethylDackel** unchanged (covered in the bisulfite section).
- **dnmtools `mlml`** (bioconda `dnmtools` 1.5.1, easy) combines BS, oxBS and TAB into 5mC/5hmC/neither estimates.
- biomodal: the **duet pipeline** (Nextflow via the biomodal CLI 2.x). It is licensed and not on bioconda, so installation is heavy.

#### (c) Output formats
- Bismark and MethylDackel as for bisulfite.
- **mlml** output: `chrom pos strand label pm ph pu conflicts`.
  - `pm`, `ph`, `pu` are 5mC, 5hmC and unmodified fractions.
  - `conflicts` counts inputs whose binomial CI excludes the estimate (`-alpha 0.05`).
  - Input is dnmtools counts format (`chrom pos strand context level n_reads`), identically sorted.
- **duet quantification**:
  - a cytosine report with 1-based `contig, position, strand, modC count, C count, context`;
  - bedMethyl, bedGraph and Bismark-style formats, 0-based;
  - separate mC, hmC and modC files for evoC/6-base.

#### (d) Defaults
- mlml: `-t 1e-10`, `-a 0.05`, no coverage filter.
- duet: reads with more than 5 Ns or under 15 nt are discarded; poly-G tails over 8 nt are trimmed; secondary, supplementary and MAPQ 0 alignments are removed; duplicates are marked by 5′ position and orientation, optionally with XL; controls are downsampled to 200×.

#### (e) Expressibility
- **TAB-seq, ACE-seq, oxBS: Y.** They are the same queries as WGBS (`C~`/`T~` over `CG`). Interpretation and mlml-style combination happen downstream.
- **biomodal: N, and not a read-vs-reference assay after resolution.** Before resolution it would need a **three-row pattern (read, mate/copy-strand base, reference)** over hairpin mates. Mates would have to be aligned as an overlapping pair, and alnbase's `overlap` command already pairs overlapping mate bases per reference position, so a `mate` row is a natural extension. But the vendor pipeline resolves before alignment, so there is no comparable intermediate. **Recommend dropping** it from the validation set, or at most running on the resolved BAM's MM tags, which is not alnbase's domain.

#### (f) Datasets
- TAB-seq: GEO **GSE36173** (mESC and H1 WG TAB-seq; large, so subsample one chromosome).
- ACE-seq: GEO series from Schutsky 2018 (accession not verified here).
- TAB-seq spike-ins: M.SssI-methylated lambda (5mC → should read T) and a 5hmC-containing amplicon (should read C). These give internal ground truth for conversion and protection rates.
- dnmtools ships small test data for mlml.

#### (g) Issues
- **Non-conversion is the whole signal.** In TAB-seq, incomplete TET oxidation of 5mC (about 3–5 %) looks like 5hmC. Genome-wide 5hmC is often only a few percent, so the false-positive rate is comparable to the signal, and conclusions about 5hmC at low-hmC loci flip.
- The subtractive BS − oxBS produces negative estimates, which mlml's `conflicts` column flags.
- **None of this discriminates between alnbase and Bismark**: the discrepancy potential is in downstream modelling, not in read-vs-reference extraction.

---

### 3. NOMe-seq (M.CviPI GpC methyltransferase + bisulfite), scNOMe/scNMT-seq

#### (a) Signal
- Accessible DNA is methylated by M.CviPI at **GpC**, so bisulfite leaves **C retained = accessible** at GCH.
- Endogenous methylation is read at **HCG**, where C retained = methylated.
- **GCG** is ambiguous and excluded by every tool.
- Contexts differ between tools:
  - **Bismark `coverage2cytosine --nome-seq`**: CpG report restricted to **ACG and TCG only** (WCG, which excludes **CCG** as well), and GpC report restricted to **GCA, GCC, GCT** (`coverage2cytosine:387-393, 781-783, 919-934`). The docs give the rationale as "the influence of G-CG (intended) and C-CG (off-target)". M.CviPI has slight off-target activity at CCG (NOMe-HiC, PMC10018866).
  - **BISCUIT `pileup -N`**: context from the 5-mer reference (`bisc_utils.c:62-72`) as HCG (A/C/T before CG, **so CCG is included**), GCG, GCHG and GCHH. `vcf2bed -t gch` takes GCHG and GCHH (`cytosine_context_nome[]` maps both to "GCH"), and `-t hcg` takes HCG. `mergecg -N` merges only when both C and G are in HCGD.
  - **methylpy**: `--num-upstream-bases 1` (default 0) writes a 4-mer context; GCH and HCG filtering is left to the user.
- Strand: each strand's C is read against its own strand's context.
- The library is directional for WGBS-style NOMe. It is **non-directional/PBAT for scNOMe and scNMT-seq** (scBS-seq). Clark et al. 2018 ran Bismark then `coverage2cytosine --nome-seq`.

#### (b) Comparators
| Tool | Version | Install |
|---|---|---|
| Bismark `extract --CX` → `bedGraph --CX` → `cov2cyt --nome-seq` | 3.1.0 | bioconda, easy |
| BISCUIT `pileup -N` → `vcf2bed -t gch/hcg` | 1.10.x | bioconda, easy |
| methylpy `call-methylation-state --num-upstream-bases 1` | 1.4.7 | bioconda, moderate |

#### (c) Output formats
- **Bismark NOMe outputs:**
  - `*.NOMe.CpG_report.txt` and `*.NOMe.GpC_report.txt`: `chr pos strand count_meth count_unmeth C-context trinucleotide`. 1-based by default; `--zero_based` subtracts 1. GpC bottom-strand rows are at pos−1.
  - `*.NOMe.CpG.cov` and `*.NOMe.GpC.cov`: `chr start end pct meth unmeth` with start = end = pos.
  - A context summary (`print_context_summary`, issue #321).
  - **Only covered positions** are written in NOMe mode.
  - The GpC file is **empty unless the input `.cov` was made with `--CX`**.
- **BISCUIT pileup VCF:**
  - `INFO: CX=` holds HCG, HCHG, HCHH, GCG or GCH in NOMe mode; `N5=` holds the 5-mer.
  - `FORMAT` includes `BT` (beta) and `CV` (coverage), among others.
  - `vcf2bed` writes `chr beg end beta cov` (0-based half-open); `-c` gives `beta% M U`; `-e` adds context columns. Default `-k 1` (min coverage 1 since v1.6).
  - The pileup stats file has `sample chrm HCGn HCGb HCHGn HCHGb HCHHn HCHHb HCHn HCHb GCn GCb`.
- **methylpy allc:** `chr pos(1-based) strand context mc_count total methylated`.

#### (d) Defaults
- **Bismark extractor**: no base-quality or MAPQ filter; Bismark's own alignment uniqueness applies; paired-end overlap is removed (`--no_overlap` default for PE); M-bias ignore is 0 unless set.
- **cov2cyt `--nome-seq`**: coverage threshold 1 (`--coverage_threshold`); keeps ACG/TCG and GCA/GCC/GCT only.
- **BISCUIT pileup** (`bisc_utils.h:101-109`, `pileup.c:1033-1046`):
  - `-b 20` minimum BQ, `-m 40` minimum MAPQ, `-a 40` minimum AS.
  - **`-5 3 -3 3`**: bases within 3 nt of either read end are ignored by default.
  - `-l 10` minimum read length.
  - Duplicates, secondary alignments and improper pairs are filtered (`-u`, `-c`, `-p` disable).
  - Overlapping mates are not double-counted (`-d` enables).
  - `-t` max retention per read is off (999999).
- **methylpy**: `--min-mapq 30`, `--min-base-quality 1`, `--trim-reads True`, `--remove-clonal False`, `--min-cov 0`.

#### (e) Expressibility
**Y for directional NOMe; needs a library mode for scNOMe/scNMT (PBAT/non-directional).**
```toml
[query.GCH_acc]     # accessible (C retained)
read = "~C~"
refr = "GCH"
mark = ".+."
[query.GCH_inacc]
read = "~T~"
refr = "GCH"
mark = ".+."
[query.WCG_m]       # Bismark --nome-seq endogenous (ACG/TCG)
read = "~C~"
refr = "WCG"
mark = ".+."
[query.HCG_m]       # BISCUIT HCG (includes CCG)
read = "~C~"
refr = "HCG"
mark = ".+."
[query.HCGD_m]      # BISCUIT mergecg -N symmetric context
read = "~C~~"
refr = "HCGD"
mark = ".+.."
```
(and the corresponding `T` rows.)
- The upstream reference base is readable past the read's 5′ end. The read-side `~` matches pad, so a C at the first read base still gets its GCH context. Tools that take context from the reference behave the same way.
- **Gap:** scNMT/scNOMe reads are non-directional or PBAT. alnbase needs `--library pbat` (R1 on CTOT/CTOB) and a non-directional mode taking strand from the aligner's tag (Bismark `XR/XG`, bwa-meth `YD:Z:f/r`).

#### (f) Datasets
- IMR90 whole-genome NOMe-seq, GEO **GSE21823** (Kelly et al. 2012, *Genome Res* 22:2497): large, so subsample one chromosome or a 10 Mb region.
- scNMT-seq, GEO **GSE109262** (Clark et al. 2018): single cells are small (a few million reads each) and processed CpG/GpC reports are deposited, which makes a direct comparator output available.
- Ground truth: CTCF/TSS nucleosome-depleted-region footprints (an accessibility dip in GCH at about ±150 bp).
- Simulation: easy, since GCH and HCG contexts are independent conversion rates.

#### (g) Issues and contestable choices
1. **CCG inclusion** (BISCUIT HCG) versus exclusion (Bismark WCG). M.CviPI off-target methylation at CCG inflates "endogenous" methylation in accessible regions. Accessible promoters and CpG islands are exactly where endogenous methylation is expected to be low and anti-correlated with accessibility, so this can **weaken or flip the reported methylation–accessibility anticorrelation** at open chromatin.
2. **BISCUIT's default 3-bp end trim, MAPQ 40 and BQ 20** versus Bismark's none shifts per-cell GpC coverage substantially for short scNMT reads.
3. **Bismark's empty GpC output** when `--CX` was not used upstream is a silent failure mode.
4. **Bottom-strand coordinates**: Bismark writes GpC bottom-strand rows at pos−1 and CpG rows at pos. Merging strands differs between tools (SMF shifts ±1; BISCUIT `mergecg -N` requires HCGD on both).
5. Endogenous GpC methylation is negligible in mammals, but not in plants or *Neurospora*; the default assumptions then break.

---

### 4. dSMF — SingleMoleculeFootprinting (Krebs lab; M.CviPI ± M.SssI + bisulfite, amplicon/capture)

#### (a) Signal
- Retained C = methylated by the exogenous MTase = accessible.
- Experiment type is **detected from the sample name** (`DetectExperimentType`, `methylation_calling.r:184-199`):
  - `_NO_` (M.CviPI only): **DGCHN** and **NWCGW** (strict contexts) for GpC and endogenous CpG.
  - `_SS_` (M.SssI): `CG`.
  - `_DE_` (double enzyme): `GCH`, `GCG` (kept as its own context) and `HCG`.
- A permissive first pass takes GC and HCG, collapses strands (minus-strand C shifted −1 for CG-like contexts and +1 for GC), and applies the coverage filter; then comes the strict context subset (`methylation_calling.r:492-545`).
- Alignment is QuasR `qAlign(bisulfite="undir", paired="fr", aligner="Rbowtie")`, i.e. **non-directional**.

#### (b) Comparator
- Bioconductor **SingleMoleculeFootprinting** 2.7.0 (devel), bioconda `bioconductor-singlemoleculefootprinting` 2.0.0 (lagging).
- **QuasR** 1.53.1 (bioconda 1.50.0).
- Install is moderate to heavy: R, BSgenome, and the ExperimentHub cache (about 1 GB in the vignette).

#### (c) Output
`CallContextMethylation()` returns `list(MethGR, MethSM)`.
- **MethGR** is a GRanges of Cs (strand-collapsed), with per-sample `<sample>_Coverage` and `<sample>_MethRate` metadata and `GenomicContext` (for example `DGCHN`/`NWCGW`). Coordinates are 1-based, from BSgenome.
- **MethSM** is a per-sample list of per-context sparse matrices, reads (`aid` names) × C positions (column names = genomic start). Values: **2 = methylated, 1 = unmethylated, 0/absent = not covered**, from `qMeth(reportLevel="alignment")` `meth+1`.
- `qMeth(mode="allC")` bulk output adds `_T` (total) and `_M` (methylated) columns.

#### (d) Defaults
- `coverage = 20`.
- `ConvRate.thr = NULL`: the read-level conversion filter is off now; it was 0.2 in 0.0.1.
- `qMeth` has `mapqMin = 0`.
- **No base-quality filter.**
- **No CIGAR parsing.** The counting loop steps reference and read in lock-step (`for(i=pos, j=0; i<iend; i++, j++)`, `QuasR/src/quantify_methylation.cpp:114-140`). That is correct only for ungapped Rbowtie alignments; a BAM from Bismark/bowtie2, bwa-meth or BISCUIT with indels or soft clips is **silently mis-registered** downstream of the first indel.
- The left mate is trimmed so it does not overlap the right mate (`iend = mpos`).
- The conversion strand is **FREVERSE only**. QuasR aligns paired bisulfite reads as `--ff` (R2 reverse-complemented before alignment, `merge_reorder_sam.cpp:255-263`, adding `XQ:i:<bisQueue 0..3>`), so **QuasR BAMs are not standard `fr` SAM** for R2.
- No duplicate filter in qMeth. QuasR alignment keeps unique hits (`-k 2 --best --strata` then filters).

#### (e) Expressibility
**Contexts: Y.** `refr="DGCHN" read="~~C~~" mark="..+.."` and `refr="NWCGW"` similarly.

**Library/strand: needs extension.**
1. `undir` means reads come from all four strands, so alnbase needs a strand-from-tag mode:
   - QuasR `XQ:i`, which encodes the converted-read/converted-genome queue;
   - Bismark `XR/XG`;
   - bwa-meth `YD`.
2. QuasR's `--ff` mate convention means R2's FLAG strand is not the standard one. alnbase would need a "strand = own FREVERSE, ignore R1/R2" mode, or the data should simply be re-aligned with Bismark `--non_directional` or bwa-meth for the comparison (preferable).
3. The per-read conversion-rate filter and strand collapse are DuckDB work (`GROUP BY record_id`, shifting ±1).

#### (f) Datasets
- `SingleMoleculeFootprintingData::NRF1pair.bam()` (ExperimentHub): mm10 dSMF `_DE_` sample around the NRF1 locus on chr6. It is the vignette dataset; the whole cache is about 1 GB, BAM size not verified.
- Package `inst/extdata/*.qs`: precomputed Methylation and TFBS objects, directly usable as expected output.
- Published: Sönmezer et al. 2021 *Mol Cell*; Kleinendorst & Barzaghi et al. 2021 *Nat Protoc*.

#### (g) Issues
1. **The CIGAR-free counting loop** is a real correctness hazard with any gapped aligner. Rbowtie is ungapped, so the default path is safe; users feeding external BAMs through `sampleFile` are not, and an indel shifts every downstream C call on that read by one base. alnbase's CIGAR-aware walk is the natural fix, and this makes a strong demonstration.
2. **Context chosen by sample-name substring** (`_NO_`, `_SS_`, `_DE_`). A misnamed sample silently gets the wrong contexts: DE contexts on an NO sample keep GCG (conflating endogenous and exogenous signal). In M.CviPI-only experiments that **inflates accessibility at methylated CpG-rich loci**.
3. **No base-quality filter** and MAPQ ≥ 0.
4. Strict DGCHN/NWCGW discards far more sites than Bismark or BISCUIT GCH/HCG, which changes single-molecule TF-occupancy state counts.

---

### 5. MAPit / MAPit-FENGC, SMAC-seq, Fiber-seq, DAF-seq, FOODIE

- **MAPit(-FENGC)** (Darst lab; bioRxiv 10.1101/2022.11.08.515732; *Methods Mol Biol* 2026):
  - Signal and contexts: M.CviPI GCH plus endogenous HCG, **GCG omitted**, on bisulfite amplicons.
  - Sequencing is **PacBio CCS** of long amplicons, not Illumina.
  - Processing: custom "reAminator" pipeline and **methylscaper** (Bioconductor 1.18.0, which does its own read-to-amplicon alignment from FASTA rather than from BAM). Filters: ≥95 % conversion and ≥95 % reference length.
  - Data: BioProject **PRJNA752452**.
  - Expressible (Y) as a read-vs-reference query if CCS reads are aligned into a BAM, but there is no clean BAM-based comparator. **Low priority.**
- **SMAC-seq** (Shipony et al. 2020, *Nat Methods*; M.CviPI + M.SssI + EcoGII) is nanopore, and **Fiber-seq** (6mA) is PacBio or nanopore. Modifications are read from signal or kinetics, not as read-vs-reference mismatches, and 6mA does not convert. **Drop both.**
- **DAF-seq** (Stergachis lab, *Nat Biotechnol* 2025) is a genuine read-vs-reference signal: an SsDddA deaminase gives C>T or G>A at accessible C, per molecule on either strand. It is primarily PacBio HiFi, though, and strand must be inferred per read. **Mention only.**
- **FOODIE** (Xie lab, *PNAS* 2024, 10.1073/pnas.2423270121) is **short-read Illumina**. The DddB deaminase converts accessible C to T in all contexts (DddA only TpC), with strand inferred per read from C>T versus G>A. Processing is Bismark alignment and dedup plus **custom Python scripts**, with no standard comparator. It needs the same per-read strand-inference mode as tagmented TAPS. **Low priority**, but a nice non-bisulfite demonstration of generality.

---

### Candidate alnbase extensions surfaced by this family
1. **More library modes** (strictly a CLI/library mode rather than query syntax, but required):
   - `pbat` (R1 on CTOT/CTOB), needed for scNMT and scNOMe;
   - `non-directional` with **strand from tag** (Bismark XR/XG, bwa-meth YD, QuasR XQ), needed for scNMT, scNOMe and dSMF;
   - **per-read inference** from conversion-motif counts, rastair-style TG>CA, needed for tagmented TAPS and FOODIE;
   - `flag-strand-only` (QuasR `--ff` BAMs).
2. **Cross-mate / copy-strand row** (a third pattern row for the mate's base at the same reference column), needed for biomodal 5/6-base on unresolved reads. There is no comparator intermediate, so it is not recommended now.
3. **Soft-clip rescue with a reference base**: confirm that `--soft-clips emit` columns carry the adjacent reference base, needed for rastair `--rescue-soft-clip-cpg`.
4. Nothing else is needed. Opposite-strand SNP evidence, read-end masks, conversion-rate filters, strand collapse, mlml combination and TAPS polarity flips are all DuckDB work.

---

### Ranking inputs
- **TAPS / 5-Base (rastair, asTair, MethylDackel-flip)**
  - Popularity: rising. TAPS+ is commercial (Watchmaker), Illumina 5-Base is on DRAGEN, and nf-core/methylseq supports it.
  - Comparator cleanliness: high. rastair is a bioconda Rust binary with BED and VCF output; asTair is easy.
  - Ground truth: lambda 0/100 % spike-ins, NA12878 GIAB genotypes, and a small TAPS BAM in the rastair repo.
  - Discrepancy likelihood: **high.** Pairwise TG versus single-position counting (a change inside rastair 2.1→2.2), SNP confounding, asTair BAQ and depth-250 cap, and exact-flag matching.
  - Expressibility: **Y** (directional). Tagmented TAPS needs strand inference.
  - Overall: **top pick of this family.**
- **NOMe-seq (Bismark `--nome-seq`, BISCUIT `-N`)**
  - Popularity: moderate. scNMT-seq is widely used.
  - Comparator cleanliness: high; both are on bioconda with explicit context rules.
  - Ground truth: TSS/CTCF footprints, easy simulation, and GSE109262 with processed reports.
  - Discrepancy likelihood: **high.** Bismark WCG versus BISCUIT HCG (CCG), BISCUIT's 3-bp end trim, MAPQ 40 versus none, and the silent empty-GpC failure.
  - Expressibility: **Y** directional; scNMT needs `pbat`/non-directional.
  - Overall: **second pick.**
- **dSMF (SingleMoleculeFootprinting/QuasR)**
  - Popularity: niche (Krebs lab and followers).
  - Comparator cleanliness: moderate. R/Bioconductor, a non-standard BAM, and output as R objects.
  - Ground truth: vendor example BAM plus precomputed objects.
  - Discrepancy likelihood: **high when a gapped aligner is used** (the CIGAR-free counting loop), plus sample-name-driven contexts.
  - Expressibility: contexts Y; needs a tag- or flag-based non-directional mode, or re-alignment.
  - Overall: a good **correctness showcase**, but moderate effort.
- **TAB-seq / ACE-seq / oxBS**
  - Popularity: declining; being replaced by TAPS/CAPS, 6-base and SIMPLE-seq.
  - Comparator: identical to WGBS, so it adds no versatility; mlml is downstream.
  - Ground truth: spike-ins.
  - Discrepancy likelihood: low at the extraction level.
  - Expressibility: **Y.**
  - Overall: low value.
- **biomodal 5-/6-base**
  - Popularity: rising.
  - Comparator: proprietary pipeline that resolves before alignment.
  - Expressibility: **N.** Needs a cross-mate row and there is no comparable intermediate.
  - Overall: **drop.**
- **MAPit-FENGC / FOODIE / DAF-seq**
  - Popularity: niche.
  - Comparators: custom or non-BAM.
  - Expressibility: MAPit Y; FOODIE and DAF-seq need per-read strand inference.
  - Overall: mention as generality examples only. **SMAC-seq and Fiber-seq: drop** (long-read signal or kinetics, not read-vs-reference).

### Sources
- rastair: https://github.com/bsbludwig/rastair (docs/src/calling/methylation.md, strand-guesser.md, formats/*.md, cli.md, CHANGELOG.md); https://www.rastair.com/calling/methylation.html; preprint https://www.biorxiv.org/content/10.64898/2026.03.19.712983v1
- asTair: https://bitbucket.org/bsblabludwig/astair (README, astair/caller.py); lambda test data https://zenodo.org/record/2582855
- taps-foundry: https://github.com/watchmaker-genomics/taps-foundry
- nf-core TAPS subworkflow: https://nf-co.re/subworkflows/bam_taps_conversion/ ; https://github.com/nf-core/methylseq
- TAPS pipeline with MethylDackel: https://github.com/jknightlab/TAPS-pipeline/ ; bwa-mem3 methylation tags: https://bwa-mem3.readthedocs.io/en/latest/methylation/tags.html
- TAPS: Liu et al. 2019 https://www.nature.com/articles/s41587-019-0041-2 ; TAPSβ/CAPS/PS: https://www.nature.com/articles/s41467-021-20920-2
- Bismark 5-Base and NOMe: Bismark repo docs/src/content/docs/rust/illumina-5-base.md, docs/.../methylation-extraction.md, coverage2cytosine; https://github.com/FelixKrueger/Bismark/blob/master/coverage2cytosine
- MethylDackel: https://github.com/dpryan79/MethylDackel (extract.c, common.c)
- BISCUIT: https://github.com/huishenlab/biscuit ; https://huishenlab.github.io/biscuit/docs/methylextraction.html ; https://huishenlab.github.io/biscuit/biscuit_vcf2bed/
- methylpy: https://github.com/yupenghe/methylpy (parser.py)
- dnmtools mlml: https://dnmtools.readthedocs.io/en/latest/mlml/
- biomodal: https://software-docs.biomodal.com/projects/data-interpretation-guide/en/latest/reference.html
- TAB-seq: Yu et al. 2012 https://www.cell.com/fulltext/S0092-8674(12)00534-X (GSE36173); ACE-seq: Schutsky et al. 2018 (PubMed 32822044 protocol); "Failure to convert" https://www.nature.com/articles/s41589-018-0172-7
- NOMe-seq: Kelly et al. 2012 https://genome.cshlp.org/content/22/12/2497 (GSE21823); scNMT-seq: https://www.nature.com/articles/s41467-018-03149-4 (GSE109262); CCG off-target: https://pmc.ncbi.nlm.nih.gov/articles/PMC10018866/ ; NOMe/ATAC/DNase features: https://www.biorxiv.org/content/10.1101/547596v2.full.pdf
- SingleMoleculeFootprinting: https://bioconductor.org/packages/release/bioc/html/SingleMoleculeFootprinting.html ; vignette https://www.bioconductor.org/packages/release/bioc/vignettes/SingleMoleculeFootprinting/inst/doc/methylation_calling_and_QCs.html ; data http://www.bioconductor.org/packages/release/data/experiment/html/SingleMoleculeFootprintingData.html ; QuasR https://git.bioconductor.org/packages/QuasR (src/quantify_methylation.cpp, src/merge_reorder_sam.cpp)
- MAPit-FENGC: https://www.biorxiv.org/content/10.1101/2022.11.08.515732v2.full ; https://link.springer.com/protocol/10.1007/978-1-0716-5072-1_2
- FOODIE: https://www.pnas.org/doi/10.1073/pnas.2423270121 ; DAF-seq: https://www.nature.com/articles/s41587-025-02914-3 ; SMAC-seq: https://experiments.springernature.com/articles/10.1038/s41592-019-0730-2


---

## Family 06 — Ancient DNA damage, CRISPR editing outcomes, SNV pileups, bisulfite SNP/ASM, DNA-damage mismatch signatures

Scope note. Everything below was checked against comparator source code cloned or fetched on 2026-09-16
(clones in `clones/{mapDamage,DamageProfiler,pydamage,CRISPResso2,pysamstats,biscuit}`, single files in
`clones/raw/`). Line numbers refer to those snapshots. "Aggregation" means a DuckDB step over alnbase
parquet rows, which per the brief does not count as a missing feature.

alnbase facts this section relies on (verified in `src/`):

- Walk orientation under `--library directional` (`src/alignment.rs:117`, `src/tags.rs:139-146`): R1 forward and
  R2 reverse are walked top strand; R1 reverse and R2 forward are walked reverse-complemented. So **R1 is always
  walked in its own sequencing orientation, and R2 is always walked complemented relative to its own sequencing
  orientation** (it is shown in its mate's frame).
- `off_5p`/`off_3p` are SEQ offsets **as sequenced**, and they **count soft-clipped bases** even when clips are
  skipped (`docs/cli-reference.md:138`, `src/alignment.rs:16,33`, `src/hits.rs:17-31`).
- Pad `_` covers both flanks and contig ends. In parquet mode an anchor may sit on a gap or pad (the "anchor on a read
  base" rule applies only to `bases` tags, `src/tags.rs:229`, `src/query_toml.rs:108`). A pad column carries the
  reference coordinate, `qual = -1` and "nearest 5'" `read_off` (`src/column.rs:10-22`).
- `alnbase overlap` merges mate overlaps with `--match-qual sum-capped` (cap 40) and `--mismatch-qual subtract` by
  default (`docs/cli-reference.md:320-323`).

---

### 1. Ancient DNA post-mortem damage (mapDamage2, DamageProfiler, PyDamage)

#### (a) Signal
- Cytosine deamination in single-stranded overhangs: **C>T rising toward the 5' end of each read**, and for
  double-stranded (ds) libraries **G>A rising toward the 3' end**, because the 3' end is the fill-in copy of the
  complementary strand's 5' overhang. Single-stranded (ss) libraries show **C>T at both ends** and no G>A excess
  (Briggs et al. 2007 PNAS; Meyer et al. 2012 Science).
- Context: any C. CpG C>T is elevated where 5mC deaminates straight to T, which UDG cannot remove, so "UDG-half"
  libraries keep damage at CpG and at the terminal base. That is a context-dependent query.
- Strand: all in **read-sequencing orientation**, per read, relative to that read's own ends. For paired-end reads
  that were not collapsed, R1's 3' end is usually not the molecule's 3' end, and R2's 5' end is the molecule's
  other terminus. The tools treat mates as independent reads, which is contestable (see g).
- Also used: base composition of the **reference just outside the read** (positions -10..-1 and +1..+10). This shows
  purine enrichment at strand breaks from depurination, and it needs reference context past the read ends.

#### (b) Comparators
| Tool | Version (bioconda) | Install |
|---|---|---|
| mapDamage2 | 2.2.3 (bioconda `mapdamage2`); GitHub master is 2.3.0a0 with a different output layout | Python plus R (Rcpp, RcppGSL, gam, inline, ggplot2). The Bayesian step compiles C++ at runtime and fails on incomplete conda toolchains (mapDamage issue #27). Moderate. |
| DamageProfiler | 1.1 (bioconda `damageprofiler`) | Java/JavaFX jar. Easy. Used by nf-core/eager. |
| PyDamage | 1.0 (bioconda `pydamage`) | Pure Python. Easy. Per-contig damage test for metagenomes. |

#### (c) Output formats
**mapDamage 2.2.3 `misincorporation.txt`** (`mapdamage/tables.py` at tag 2.2.3): `#` comment lines, then the TSV header
`Chr End Std Pos A C G T Total G>A C>T A>G T>C A>C A>T C>G C>A T>G T>A G>C G>T A>- T>- C>- G>- ->A ->T ->C ->G S`
(`mapdamage/seq.py` `HEADER`).
- `End` is `5p` or `3p`. `Std` is the read's mapping strand (`+`/`-`).
- `Pos` is 1-based and runs 1..`--length` (default 70) from that end.
- `A C G T` count the reference base at that position, and `Total` is their sum.
- `X>Y` columns count reference X read as Y, with both taken in read orientation (reverse reads are
  reverse-complemented, `main.py`).
- `X>-` counts deletions, `->X` counts insertions, and `S` counts soft-clipped reads (clip length ≥ Pos).
- On GitHub master (2.3.0a0) `Chr` is replaced by `Sample Library` (`statistics.py:_write_freq_table`, issue #41).

**`5pCtoT_freq.txt` / `3pGtoA_freq.txt`** (R, `Rscripts/mapDamage.R:104-111,170-171` at 2.2.3): the header is
`pos\t5pC>T` (or `pos\t3pG>A`). The value is `sum(C>T)/sum(C)`, taken over Chr and Std at each Pos. These files do not
exist in master.

**`dnacomp.txt`**: `Chr End Std Pos A C G T Total`. For `5p`, Pos runs -10..-1 (reference before the read) and then
1..70 (read bases). For `3p` it is mirrored.

**`Stats_out_MCMC_*.csv`**: posterior summaries for the Briggs-model parameters (δs, δd, λ, θ, ...).

**`lgdistribution.txt`**: `Std Length Occurences`.

**DamageProfiler** (`io/OutputGenerator.java:321-360`) writes `misincorporation.txt` with the mapDamage 2.2.x header
(`Chr End Std Pos A C G T Total G>A C>T ... S`), plus `5pCtoT_freq.txt`, `3pGtoA_freq.txt` (lines 535-538),
`dnacomp.txt`, `DNA_comp_genome.txt`, `lgdistribution.txt`, `dmgprof.json`, `edit_distance.txt` and PDFs.

**PyDamage** writes `pydamage_results.csv`, one row per reference contig. Its columns are `reference`,
`null_model_p0`, `null_model_p0_stdev`, `damage_model_p`, `damage_model_p_stdev`, `damage_model_pmin`, `..._stdev`,
`damage_model_pmax`, `..._stdev`, `pvalue`, `qvalue`, `RMSE`, `nb_reads_aligned`, `coverage`, `reflen`,
`predicted_accuracy`, then `CtoT-0..CtoT-(wlen-1)` and `GtoA-0..` (`damage.py:240-250`). Positions are 0-based.

#### (d) Implicit filters and defaults (source-verified)
| | mapDamage 2.2.3 | DamageProfiler 1.1 | PyDamage 1.0 |
|---|---|---|---|
| Flags dropped | UNMAP, SECONDARY, QCFAIL, DUP, SUPPLEMENTARY (`main.py:_filter_reads`) | **only unmapped** (`DamageProfiler.java:137-146`), so duplicates, secondary and supplementary alignments are all counted | **only unmapped** (`damage.py:70`) |
| MAPQ | none | none | none |
| Base qual | `-Q 0`. When set, low-quality bases are masked to N in both read and reference, so they leave the denominator too (`align.py:align_with_qual`) | none | none |
| Position index | **alignment-column index from the first aligned base**. Soft clips are excluded (`read.query`) and deletion/insertion columns occupy positions (`align.py`, `statistics.py:update`) | **SEQ offset, including soft clips** (read string vs `makeReferenceFromAlignment`, where clipped positions hold ref `'0'` and are not counted). Deletions take no position. | **alignment-string index including soft-clip columns** (ref `' '`), deletions and insertions (`parse_damage.py:37-60`) |
| Ends / strands | both ends, both strands | both ends, both strands. `-sslib` only changes plotting/labels | **5' only**: forward reads give C>T from the left end. Reverse reads give ref-orientation G>A from the right end, which is the same 5' C>T event, but it is reported as "GtoA" (`parse_damage.py:62-80`). |
| Window | `--length 70`, plots `-m 25` | `-l 100` (frequencies), `-t 25` (plot/output) | `--wlen 13` |
| MD tag | not needed | used if present, else computed | **required** (`get_reference_sequence`) |
| Paired-end | mates are independent reads. The authors recommend collapsing (issues #27, #38). | independent; `-only_merged` keeps `M_`-prefixed names | independent |

#### (e) Expressibility in alnbase
**Yes, today, with aggregation.**
- Queries: `read=N refr=A`, `N/C`, `N/G`, `N/T` (or 4 queries `N@X`) give read_base, refr_base, `off_5p`, `off_3p` and
  `flag`. Deletions come from `read=. refr=N`. Insertions need `--insertions emit`. Soft-clip counts come from
  `--soft-clips emit` or CIGAR.
- DuckDB groups by off_5p (or off_3p) and refr_base/read_base. **R2 rows must be complemented** to get
  read-orientation substitutions, because of the directional walk.
- `dnacomp` outside-read context: anchor a query on the pad column (`read=_ refr=A`), then derive the signed distance
  from `pos`/`cigar` record fields.
- Exact **mapDamage parity** needs its alignment-column position. That means off_5p minus the leading soft clip plus
  the preceding indel columns. It can be derived from `cigar` in DuckDB, but the derivation is painful. A per-hit
  **alignment-column offset output column** would make it trivial. That is an output column, not a syntax change.
- CpG-specific damage for UDG-half libraries (`T~ / CG`), and "the terminal base only", are already expressible.

Nice-to-have, not required: `--library as-sequenced` (walk every read in its own orientation), which matches all
three comparators without the R2 complement step.

#### (f) Test data
- nf-core/eager test-datasets (branch `eager`):
  - `testdata/Mammoth/bam/JK2782_TGGCCGATCAACGA_L008_R1_001.fastq.gz.tengrand.fq.combined.fq.mapped.bam` (452 kB,
    collapsed reads vs mammoth mtDNA)
  - `reference/Mammoth/Mammoth_MT_Krause.fasta` (17 kB)
  - `testdata/Human/bam/JK2067_downsampled_s0.1.bam` (680 kB)
  - Raw URL prefix: `https://raw.githubusercontent.com/nf-core/test-datasets/eager/`
- PyDamage `tests/data/aligned.bam` (880 kB) with `NZ_JHCB02000018.1.fa` (910 kB).
- mapDamage `mapdamage/tests/test.bam` (588 bytes) with `fake1.fasta` (tag 2.2.3).
- **Ground truth by simulation:** gargammel (bioconda `gargammel`; Renaud et al. 2017 Bioinformatics) with Briggs
  parameters for ds or ss libraries. Or add damage to ART reads with a known per-position C>T profile.

#### (g) Known issues and contestable choices
1. **Three incompatible position conventions** (source-verified, table above). Under soft clipping these disagree by
   the clip length. bwa-mem soft-clips damaged terminal mismatches (Oliva et al. 2021 Brief Bioinf
   doi:10.1093/bib/bbab076). The consequences differ by tool:
   - mapDamage re-bases to the first aligned base.
   - DamageProfiler and PyDamage keep the clipped positions in the index but never count them.
   - Result: on bwa-mem BAMs, DamageProfiler/PyDamage position-1 C>T can be near zero while mapDamage shows a peak.

   This can flip **authentication** (PyDamage's `predicted_accuracy`/`qvalue` contig classification, and nf-core/eager
   damage-based filtering). alnbase makes the convention explicit.
2. **DamageProfiler counts duplicates, secondary and supplementary alignments.**
   - Supplementary pieces have artificial ends at chimeric junctions, which carry no damage, so they dilute terminal
     C>T.
   - Duplicates weight molecules unevenly.
   - mapDamage drops all of these.
3. **DamageProfiler window-shrink bug** (source reading, not yet reproduced). `Frequencies.comparePos` runs
   `if (seq.length < this.length) this.length = seq.length;` on a **field**. After the first read shorter than
   `-l`, every later read is counted only up to that length, so counts at positions beyond the shortest read seen are
   order-dependent. On a coordinate-sorted BAM they come from early contigs only. It is visible when `-t` exceeds the
   minimum read length.
4. **PyDamage never looks at 3' ends** and needs MD tags. Its "GtoA" columns are 5' C>T on reverse reads, so a ss
   library's 3' C>T is invisible.
5. **Paired, non-collapsed reads.** All three tools treat R1's 3' end as a molecule end, which dilutes 3' G>A. mapDamage
   authors recommend collapsing, but collapsing loses the non-overlapping pairs (issue #27 thread).
   `--single-stranded` in mapDamage affects only the Bayesian model and only 5' (issue #38).
6. mapDamage does not handle `N` CIGAR ops (not relevant for DNA). mapDamage 2.3.0a0 changed the column set, so
   pipelines parsing `Chr` break (issue #41).

---

### 2. CRISPR base editing (and prime editing) outcomes — CRISPResso2, BE-Analyzer

#### (a) Signal
- Base editors: CBE gives **C>T** (G>A when the protospacer is on the minus strand). ABE gives **A>G**.
- Position is restricted to an **editing window at fixed coordinates relative to the protospacer**. The usual
  CBE/ABE window is protospacer positions ~4-8 counted from the PAM-distal 5' end, i.e. 13-17 bp 5' of the PAM.
  Sequence context matters (e.g. TC preference of APOBEC-based CBEs).
- Readouts:
  - per-position conversion fraction
  - per-read allele classes and bystander co-editing (linkage within a read)
  - indel byproducts at the nick site
  - for nuclease experiments, indels overlapping the cut-site window
- Strand: amplicon (genomic) orientation. Protospacer-strand–relative numbering is presentation only.

#### (b) Comparators
- **CRISPResso2 2.3.4** (bioconda `crispresso2`, released 2026-04-23; Clement et al. 2019 Nat Biotech).
  - Python/Cython with fastp, bowtie2 and samtools; Docker available. Moderate install.
  - It **aligns reads itself** with a Needleman-Wunsch aligner that has a gap incentive at the cut site
    (`--needleman_wunsch_gap_open -20`, `_extend -2`, `_gap_incentive 1`, EDNAFULL).
  - `--bam_input` uses only SEQ (SAM column 10) from the BAM and **re-aligns** it (`CRISPRessoCORE.py:process_bam`,
    ~l.2035-2050). The BAM alignment is therefore ignored.
  - `--bam_output` writes CRISPResso's own alignments as BAM against `CRISPResso_output.fa` (CHANGELOG, v2.2.x). That
    is the clean way to hand identical alignments to alnbase.
- **BE-Analyzer** (Hwang et al. 2018 BMC Bioinf 19:542) is web-only: client-side WebAssembly at rgenome.net, with no
  CLI or source release. **Unsuitable as a scripted comparator.**
- Others: CRISPRessoBatch/Compare; the Broad base editor validation pipeline (wraps CRISPResso2); CRISPR-SURF (not
  base-level).

#### (c) Output formats (CRISPResso2)
- **`Nucleotide_frequency_table.txt`** and **`Nucleotide_percentage_table.txt`**
  (`plots/data_prep.py:3553-3570`)
  - TSV. The header row is the amplicon sequence, one base per column with an empty first cell (so column names
    repeat).
  - Row labels are `A C G T N -` and cells are read counts at each 0-based amplicon position.
  - The percentage table divides by total aligned reads.
- **`Quantification_window_nucleotide_frequency_table.txt`**: the same layout, restricted to the window.
- **`Substitution_frequency_table.txt`**, `Quantification_window_substitution_frequency_table.txt` (with
  `--base_editor_output`): rows `A C G T N`.
- **`Selected_nucleotide_{frequency,percentage}_table_around_sgRNA_<seq>.txt`**: columns are the targeted base
  numbered within the plot window (e.g. `C5`), and rows are nucleotides.
- **`Alleles_frequency_table.zip`**: `Aligned_Sequence Reference_Sequence Unedited n_deleted n_inserted n_mutated
  #Reads %Reads`. `--write_detailed_allele_table` adds ref_positions, all_/substitution_positions,
  substitution_values, and more.
- **`CRISPResso_quantification_of_editing_frequency.txt`**: `Amplicon Unmodified% Modified% Reads_in_input
  Reads_aligned_all_amplicons Reads_aligned Unmodified Modified Discarded Insertions Deletions Substitutions Only
  Insertions ...`
- **`Modification_count_vectors.txt`** and `Quantification_window_modification_count_vectors.txt`: the first row is
  the amplicon, then rows Insertions, Insertions_Left, Deletions, Substitutions, All_modifications, Total. An insertion
  increments **both** flanking positions in the "Insertions" row.
- `CRISPResso_output.vcf` with `--vcf_output --amplicon_coordinates` (v2.3.2+).

#### (d) Implicit filters and defaults (`CRISPResso2/args.json`)
- Reads:
  - `--min_average_read_quality 0`, `--min_single_bp_quality 0`
  - `--min_bp_quality_or_N 0`: bases below this are **turned into N**, not dropped
  - `--default_min_aln_score 60`: homology % needed to assign a read to an amplicon
  - Ambiguous reads (equal score to several amplicons) are **excluded** unless `--expand_ambiguous_alignments` or
    `--assign_ambiguous_alignments_to_first_reference` is set.
  - Duplicates are not removed. `--samtools_exclude_flags 4` applies for BAM input.
- Windows:
  - `--quantification_window_size 1` and `--quantification_window_center -3` (Cas9 cut site): a 2-bp window that
    decides modified vs unmodified.
  - Substitutions outside the window do not make a read "modified".
  - `--exclude_bp_from_left/right 15` for indel quantification.
  - Insertions count only if fully inside the window (≥2.1.0; `--use_legacy_insertion_quantification`).
- Paired-end FASTQ is merged with fastp (`--min_paired_end_reads_overlap 10`). Non-overlapping pairs are lost, and BAM
  input is **not** merged.
- Nucleotide tables: positions a read does not cover are counted as `-` (deletion), because of global alignment end
  gaps (`CRISPRessoCORE.py:~4075-4081` loops over `aln_seq` with `ref_pos>=0`). The denominator is all aligned reads.

#### (e) Expressibility in alnbase
**Per-base tables: Yes, today.** Run on CRISPResso's `--bam_output` BAM against `CRISPResso_output.fa`, or on a
bwa/minimap2 amplicon BAM.
- Queries `read=N refr=N` plus `read=. refr=N` give per-position A/C/G/T/N/- counts. The window is `refr_pos BETWEEN`,
  applied in DuckDB. `--insertions emit` gives insertion counts.
- Per-read "edited in window", bystander linkage and allele classes: `GROUP BY record_id`.
- Context-restricted conversion (e.g. `TC` for APOBEC): the pattern `T~ / TC` anchored at col 1, or `~T/TC`.
- **Prime editing: No (out of scope).** Classification means competitive alignment to WT versus edited amplicon
  references plus scaffold-incorporation detection. That is not a single-reference column comparison.
- No syntax extension is needed.

#### (f) Test data
- CRISPResso2 repo `tests/`: `FANC.Cas9.fastq` (139 kB), `FANC.Untreated.fastq` (139 kB), `HEK3.Cas9.fastq`
  (130 kB), `Both.Cas9.fastq.smallGenome.bam` (60 kB), with expected outputs in `tests/expectedResults/`
  (`testRelease.sh` diffs `Nucleotide_frequency_table.txt`).
- Base editor demo: `http://crispresso.pinellolab.partners.org/static/demo/base_editor.fastq.gz` (25,000 reads, EMX1).
  Command:
  `--guide_seq GAGTCCGAGCAGAAGAAGAA --quantification_window_size 10 --quantification_window_center -10 --base_editor_output`
  (README of the eric-erki mirror of CRISPResso2).
- Batch demo: SRR3305543-SRR3305546, first 25,000 reads each (BE1/BE2/BE3/untreated; Komor et al. 2016).
- Ground truth: simulate amplicon reads with known per-position C>T fractions plus seq errors (trivial).

#### (g) Known issues and contestable choices
1. The **default 2-bp window** plus counting **any substitution in the window** as "Modified" means sequencing errors
   inflate Modified% with no untreated-control subtraction (CRISPRessoCompare is needed).
2. **Gap incentive at the cut site**: indels in repeats are placed at the cut site rather than left-normalized.
   - Left-normalized BAM aligners put the same event just outside a 2-bp window, so the read counts as unmodified.
   - Result: different editing % on identical reads. That is an alignment convention, not biology. Keep alignments
     identical (`--bam_output`) when validating.
3. **Uncovered positions count as `-`**, so partial-length reads depress per-position base-editing percentages near
   amplicon ends.
4. **`--bam_input` ignores the input alignment and does not merge mates**. PE BAMs double-count fragments relative to
   FASTQ mode.
5. Recent bugs:
   - Off-by-one clipping of the quantification window when it reaches the reference end (PR #651, unreleased after
     2.3.4).
   - `--bam_output` truncated after the first unaligned read (PR #602).
   - Quantification-window coordinate inference across multiple amplicons regressed (PR #598).
6. **`--min_bp_quality_or_N`** turns low-quality bases into N. N is not a substitution ("N's don't count as
   substitutions", CHANGELOG), so raising the threshold silently lowers apparent editing.

---

### 3. SNV pileups (baseline): bcftools mpileup/call, samtools mpileup, pysamstats, GATK CollectAllelicCounts

#### (a) Signal
- Non-reference base at a position (read≠ref), on both strands, counted per position with strand split (ADF/ADR, DP4).
- Indels: `+n`/`-n` markers in pileup text, deletion placeholders `*`.
- Read-position bias (RPBZ) uses the offset in the read.

#### (b) Comparators
- bcftools 1.24 and samtools 1.24 (bioconda; trivial install).
- pysamstats 1.1.2 (bioconda; Cython on pysam ≥0.15; unmaintained since ~2019).
- GATK4 4.6.2.0 `CollectAllelicCounts` (bioconda `gatk4`; Java).

#### (c) Output formats
- **samtools mpileup**: text columns `chrom pos(1-based) ref depth read_bases base_quals`.
  - Optional columns: `-s` MAPQ, `-O`/`--output-BP` position in read, `--output-BP-5`, `--output-QNAME`, `-M` mods.
  - Base symbols: `.` or `,` for ref match on forward or reverse strand, uppercase or lowercase for mismatches, `*`
    for a deletion, `^Q` for read start plus MAPQ, `$` for read end.
- **bcftools mpileup**: VCF/BCF, POS 1-based.
  - Default INFO: `DP, I16, QS, SGB, RPBZ, MQBZ, MQSBZ, BQBZ, SCBZ, MQ0F, VDB, IDV, IMF` (`mpileup.c:1402`).
  - Default FORMAT: `PL, AD` ("high-quality bases").
  - Optional annotations: `ADF, ADR, DP, DP4 (deprecated), SP, DPR`.
  - `I16` is 16 sums: ref-fwd, ref-rev, alt-fwd, alt-rev counts, then sums of BQ, BQ², MQ, MQ², distance-to-end, etc.
- **bcftools call -m**: GT, PL and QUAL.
- **pysamstats `-t variation`** (`config.py:93-115`): `chrom pos ref reads_all reads_pp matches matches_pp mismatches
  mismatches_pp deletions deletions_pp insertions insertions_pp A A_pp C C_pp T T_pp G G_pp N N_pp`.
  - Positions are 1-based unless `--zero-based`.
  - `variation_strand` adds `_fwd`/`_rev` splits of each.
- **GATK CollectAllelicCounts**: TSV with a SAM-style `@` header, columns `CONTIG POSITION REF_COUNT ALT_COUNT
  REF_NUCLEOTIDE ALT_NUCLEOTIDE` (1-based; `AllelicCountCollection.java`).

#### (d) Implicit filters and defaults (source-verified)
| | bcftools mpileup ≥1.13 | samtools mpileup 1.24 | pysamstats 1.1.2 (via pysam pileup) | GATK CollectAllelicCounts |
|---|---|---|---|---|
| Flags skipped | UNMAP, SECONDARY, QCFAIL, DUP (`mpileup.c:1395`) | same (`bam_plcmd.c:1177`) | same (pysam `flag_filter`, `libcalignmentfile.pyx:2470`) | unmapped, duplicate, zero-ref-length; MQ filter |
| Orphans (paired, not proper) | **skipped** (`MPLP_NO_ORPHAN`; `-A` keeps) | **skipped** | **skipped** (pysam `ignore_orphans=True`, l.2474), even though the output has separate `reads_all` / `reads_pp` columns, so `reads_all` ≈ `reads_pp` for PE data | not skipped |
| Min MAPQ | 0 | 0 | 0 | **30** |
| Min BQ | **1** (13 until 1.12; `-X 1.12` restores) | **13** | **0** at C level. pysam's `min_base_quality=13` applies only to Python accessors, and pysamstats reads `col.plp` directly (`opt.pyx:1866-1905`) | **20** |
| BAQ | **partial** BAQ (only in problematic regions, `MPLP_REALN_PARTIAL`); `-D` full, `-B` off | **full** BAQ on (`MPLP_REALN`) | off (no fasta passed to pileup) | none |
| Other quality transforms | `--max-BQ 60`, `--delta-BQ 30` (BQ capped at neighbour BQ + 30) | — | — | — |
| Max depth | **250 per file**, subsampled (`-d 0` unlimited since 1.23) | 8000 | 8000 | no downsampling |
| Mate overlap | detected: htslib `tweak_overlap_quality` (`sam.c:5877`). If bases agree, one mate gets qual sum capped at 200 and the other 0. If they disagree, the higher-quality base is ×0.8 and the lower set to 0. The mate picked is deterministic by read-name hash. | same | pysam `ignore_overlaps=True` modifies quals, but pysamstats has min_baseq 0, so **both mates are still counted** | none (double counted) |

#### (e) Expressibility in alnbase
**Yes, today.**
- `read=N refr=N` gives a row per aligned base. `read=. refr=N` gives deletions. Insertions need `--insertions emit`.
- Strand and position-in-read come from `flag`/off_5p. Pileup counts are an aggregation.
- To match: run `samtools view -F 0x904 [-f 2]` upstream and `alnbase overlap` for mates. Compare with comparators
  run `-B` (no BAQ), `-d 0`, `-A` where appropriate. Note that alnbase's overlap qualities (cap 40, subtract) differ
  from htslib's (cap 200, ×0.8), so thresholded counts can differ at overlap-mismatch sites.
- No syntax extension needed.

#### (f) Test data
- bcftools repo `test/mpileup/mpileup.{1,2,3}.bam` (68 kB, 28 kB, 28 kB), `mpileup.ref.fa` (4 kB), with expected
  outputs `test/mpileup/mpileup.*.out`.
- nf-core/test-datasets branch `modules`:
  `data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam` plus the chr22 subset genome and VCFs (small,
  MB-scale; size not re-verified).
- GIAB HG002 truth on a small region for calling.
- Simulation: ART or wgsim reads with BAMSurgeon spikes at known VAF.

#### (g) Known issues and contestable choices
1. **Default BQ differs 1 / 13 / 0 / 20** across bcftools, samtools, pysamstats and GATK. AD or VAF at the same site
   differs by tool.
2. **bcftools' 250-read cap** subsamples amplicon or deep panel data. AD is then a random subsample, which matters for
   low-VAF somatic calls and for CRISPR amplicons.
3. **Orphans silently skipped** by all htslib-based tools: discordant pairs near SVs vanish.
4. **pysamstats counts both overlapping mates** despite overlap detection, because the qual-zeroing is bypassed. It
   also drops orphans, so `reads_pp == reads_all` for most PE BAMs.
5. **BAQ** can drop true SNVs near indels, and partial versus full BAQ changed between versions (bcftools PR #1474).
   BAQ is a realignment-derived quality transform that alnbase does not implement. It should be disabled in
   comparisons.
6. **Overlap handling picks the mate by read-name hash** when qualities tie. Mate-specific artifacts (R1 oxoG) are
   therefore removed at random from the counts.

---

### 4. Bisulfite SNP-aware calling and allele-specific methylation — BISCUIT, Bis-SNP, Revelio (+ CGmapTools, MethylExtract)

#### (a) Signal
- Under bisulfite (or EM-seq), a T at a reference C on a C-strand read is ambiguous: unmethylated C or a C>T SNP.
- The **opposite strand disambiguates**: at the same coordinate, reads from the complementary strand see G versus A,
  which conversion does not affect.
  - In BAM SEQ terms, BSW reads (XG=CT: OT, CTOT) make C/T ambiguous, while their G/A calls are informative.
  - BSC reads (XG=GA: OB, CTOB) make G/A ambiguous, while their C/T calls are informative.
- In alnbase's directional walk every conversion reads as C>T. So at a top-strand C position, bottom-walked reads
  show `read A / refr G` (walk frame) if a C/T SNP exists.
- ASM (allele-specific methylation): link the SNP allele and CpG states **on the same read or fragment**.

#### (b) Comparators
- **BISCUIT 1.10.2** (bioconda `biscuit`, 2026-08; Zhou et al. 2024 NAR). C; easy.
  - `biscuit pileup` produces VCF; `vcf2bed -t snp|cg`; `epiread` plus `asm` for ASM.
- **Bis-SNP 1.0.1** (bioconda `bis-snp`; Liu et al. 2012 Genome Biol). GATK-3-era Java: it needs read groups, a
  sequence dictionary and old Java. Hard. Unmaintained.
- **Revelio** (Nunn et al. 2022 BMC Genomics 23:477, github.com/bio15anu/revelio) masks ambiguous base qualities in the
  BAM. Callers such as GATK or freebayes then call SNPs normally.
- CGmapTools `snv` / `asm` (not on bioconda). MethylExtract. EpiDiverse/SNP.

#### (c) Output formats
- **BISCUIT pileup VCF** (`pileup.c:942-968`):
  - INFO: `NS, CX` (context CG/CHG/CHH, or HCG/GCH in NOMe mode), `N5` (5-mer context), `AB` (ambiguous alt).
  - FORMAT: `GT, DP, SP` (allele support after bisulfite inference, e.g. `C11,T3`), `AC, AF1, CV` (strand-specific
    cytosine coverage), `BT` (beta), `GL1, GQ`, and RN/CN with some options.
  - Verbose `-v` adds per-read `Bs0/Bs1, Sta, Bq, Str, Pos, Rret`.
- **`biscuit vcf2bed -t cg`**: `chr start(0-based) end beta cov`. With `-e`, context columns come before beta; with
  `-c`, beta M U. `-k 1` is the default min coverage since 1.6.0.
- **`vcf2bed -t snp`**: chr, start, end, ref, alt and allele info.
- **Bis-SNP**: `cpg.vcf` and `snp.vcf` with FORMAT `BRC6` (C, T, other on the cytosine strand; G, A, other on the
  guanine strand), plus `CM, CU, DP, GT, GQ`.
- **Revelio**: a BAM with modified qualities.

#### (d) Implicit filters and defaults
**BISCUIT pileup** (`bisc_utils.h:98-115`, `pileup.c:758-790,505-520`):
- `-b` min BQ 20, `-m` MAPQ **40**, `-a` min AS **40**, `-l` min read length 10.
- **`-5 3` / `-3 3`**: bases within 3 of either end are dropped. The check is done in **BAM SEQ orientation**
  (`d->qpos` is a SEQ index), so asymmetric settings are wrong for reverse reads. The code comment admits that the 3'
  end should use the aligned end because of soft clips.
- Improper pairs, secondary, duplicate and QC-fail reads are filtered.
- **Overlap handling drops R2 bases** in the mate-overlap interval, using coordinates from the MC tag (or assuming
  equal lengths). It is not quality-aware.
- Strand comes from the YD, ZS or XG tag, otherwise from inferred C>T vs G>A counts.
- SNP genotyping: error 0.001, mu 0.001, contam 0.01, priors 1/3.
- `redistribute_cnts`: ambiguous Y/R counts are pooled with T/A **if there is any T/A evidence across samples and no
  C/G** (cross-sample inference).
- A site is "methcallable" if T/C < 0.05 (`pileup.c:525-535`).

**Bis-SNP**: `-mmq 30`, `-mbq 5`, `minConv 1` (reads must show at least 1 conversion near the 5' end).

**Revelio**: sets BQ to 0 for ambiguous calls (T on C-strand reads, A on G-strand reads) and optionally clips read
ends.

#### (e) Expressibility in alnbase
**Yes, today, with aggregation.**
- Per-read rows at C positions (`read=N refr=C`) and at walk-G positions (`read=N refr=G`) carry refr_pos and flag.
  Pivot on refr_pos in DuckDB: T/C on walk-top plus A/G on walk-bottom at the same refr_pos is SNP evidence.
- ASM: `GROUP BY qname` joins the SNP-allele row and CpG call rows from the same record or fragment (`-f qname`).
- alnbase can also do something BISCUIT and Bis-SNP cannot state: **read-aware context**. The CpG pattern `YG/CG`
  requires the read itself to carry G at the G column, which excludes reads carrying a CpG-destroying allele.
- For the **tag output path** (e.g. a Revelio-like per-base mask, or strand-specific codes), a query-level **strand
  or mate predicate** would help, since tags cannot be post-filtered. This is an optional extension.

#### (f) Test data
- BISCUIT `test/` (Sherman-simulated 1,000 PE reads on hg38 chr22; fetched by `setup_tests.py`).
- nf-core/methylseq test-datasets (small bwa-meth and Bismark BAMs).
- Ground truth: simulate with Sherman on a genome carrying known SNPs (a VCF applied with `bcftools consensus`) at set
  conversion rates.
- Real: GIAB HG002 EM-seq or WGBS with GIAB SNP truth. NEB's EM-seq germline calling poster used HG002; a small chr20
  slice is enough.

#### (g) Known issues and contestable choices
1. **Where end-trimming is measured.** BISCUIT's end exclusion uses SEQ orientation and includes soft clips, while
   Bismark and MethylDackel use read orientation or M-bias. Methylation near read ends shifts between tools.
2. **R2 drop versus quality merge** in overlaps. When mates disagree (a PCR error on R1), BISCUIT keeps R1
   unconditionally, htslib keeps the higher-quality base, and alnbase `overlap` subtracts.
3. **Cross-sample redistribution** in BISCUIT: an ambiguous T in one sample is assigned by other samples' evidence.
   A joint call can therefore flip one sample's SNP or methylation status.
4. **T/C < 5% "methcallable" cutoff**: heterozygous C/T SNPs at low coverage are called as partial methylation. That
   is the classic false DMR/ASM source.
5. Bis-SNP depends on GATK3 and is effectively unrunnable on modern stacks. Excluding it is defensible.

---

### 5. DNA damage and library-artifact mismatch signatures (8-oxoG G>T, FFPE C>T), error-by-cycle, duplex error profiles, truncation-based lesion maps

#### 5a. Oxidative (8-oxoG) and deamination artifacts — Picard CollectOxoGMetrics / CollectSequencingArtifactMetrics, GATK CollectF1R2Counts, fgbio ErrorRateByReadPosition

**(a) Signal**
- 8-oxoG formed during shearing is read as **G>T on R1 and C>A on R2** (read orientation), in context CCG > others
  (Costello et al. 2013 NAR 41:e67).
- FFPE and heat deamination give **C>T on R1 and G>A on R2** ("pre-adapter" artifacts).
- Hybrid-capture "bait-bias" artifacts are reference-strand biased, not read-number biased.
- Under alnbase's directional walk, R2 is shown in R1's frame. **Both oxoG classes therefore collapse to walk `read T /
  refr G`, with the control being `A/C`.** One pattern plus a 3-mer context (`~T~ / NGN`) covers what Picard computes
  in several branches.

**(b) Comparators**
- Picard 3.5.0 (bioconda `picard`; Java; easy).
- GATK4 4.6.2.0 `CollectF1R2Counts` (feeds `LearnReadOrientationModel` for Mutect2).
- fgbio 4.1.1 `ErrorRateByReadPosition` (bioconda `fgbio`; Java; easy).

**(c) Output**
- **CollectOxoGMetrics** (Picard metrics file; `CollectOxoGMetrics.java:140-200`), one row per 3-mer context. Columns:
  `SAMPLE_ALIAS LIBRARY CONTEXT TOTAL_SITES TOTAL_BASES REF_NONOXO_BASES REF_OXO_BASES REF_TOTAL_BASES
  ALT_NONOXO_BASES ALT_OXO_BASES OXIDATION_ERROR_RATE OXIDATION_Q C_REF_REF_BASES G_REF_REF_BASES C_REF_ALT_BASES
  G_REF_ALT_BASES C_REF_OXO_ERROR_RATE C_REF_OXO_Q G_REF_OXO_ERROR_RATE G_REF_OXO_Q`.
  - `OXIDATION_ERROR_RATE = max(ALT_OXO − ALT_NONOXO, 1)/TOTAL_BASES`.
  - Contexts centred on C, with G sites reverse-complemented.
  - Oxo means ref-G in read orientation on R1, or ref-C on R2.
- **CollectSequencingArtifactMetrics** writes several files:
  - `.pre_adapter_detail_metrics`: `SAMPLE_ALIAS LIBRARY REF_BASE ALT_BASE CONTEXT PRO_REF_BASES PRO_ALT_BASES
    CON_REF_BASES CON_ALT_BASES ERROR_RATE QSCORE`, with ERROR_RATE =
    max(1e-10,(PRO_ALT−CON_ALT)/(all four))
  - `.pre_adapter_summary_metrics` (TOTAL_QSCORE, WORST_CXT, ..., ARTIFACT_NAME)
  - `.bait_bias_detail_metrics` (FWD_CXT_REF_BASES ... ERROR_RATE)
  - `.bait_bias_summary_metrics`
  - `.error_summary_metrics`
- **fgbio ErrorRateByReadPosition**: `.error_rate_by_read_position.txt` with columns `read_number(0/1/2) position
  (1-based cycle) bases_total errors error_rate a_to_c_error_rate a_to_g_error_rate a_to_t_error_rate
  c_to_a_error_rate c_to_g_error_rate c_to_t_error_rate [g_to_a ... t_to_g when not collapsed] collapsed`
  (`ErrorRateByReadPosition.scala:318-335`).

**(d) Defaults**
- Picard OxoG:
  - `MINIMUM_QUALITY_SCORE 20`, `MINIMUM_MAPPING_QUALITY 30`, `MINIMUM_INSERT_SIZE 60`, `MAXIMUM_INSERT_SIZE 600`
    (reads outside are dropped, including unpaired)
  - `INCLUDE_NON_PF_READS true`, `USE_OQ true`, `CONTEXT_SIZE 1`
  - Non-primary and duplicate reads filtered; optional `DB_SNP` mask
  - SamLocusIterator does **no mate-overlap dedup**
- Picard SequencingArtifact: the same, but `INCLUDE_NON_PF_READS false`, `INCLUDE_DUPLICATES false`,
  `INCLUDE_UNPAIRED false`, `TANDEM_READS false`.
- fgbio: `--min-mapping-quality 20`, `--min-base-quality 0`, duplicates excluded, secondary/supplementary and QC-fail
  excluded, `--collapse true`.
  - Cycle = `readLength − offset` for reverse reads, so it is **read-orientation SEQ offset including soft clips**.
    That matches alnbase off_5p+1.
  - Non-ACGT bases skipped, indels not counted, optional `--variants` mask.

**(e) Expressibility: Yes, today.**
- Queries `read=N refr=N` with a 3-mer context (`~N~ / NNN`, mark `.+.` with `^` on flanks, or 16 context queries)
  give refr_base, read_base, flag and off_5p.
- DuckDB reproduces PRO/CON using read number: under the directional walk, R1 is shown as read, and R2 needs
  complementing to get "baseAsRead".
- No syntax extension needed. A mate predicate would help only for tag output.

**(f) Test data**
- Any small PE WGS BAM: nf-core modules `test.paired_end.sorted.bam` (chr22 subset), or GIAB HG002 downsampled.
- Ground truth: inject G>T into R1 only and C>A into R2 only at a known rate on ART-simulated reads (a 30-line pysam
  script).
- Real oxoG examples are in TCGA (controlled access), which is not needed.

**(g) Issues and flip scenarios**
1. **fgbio's default `--collapse`** folds G>T into C>A per read number, which **hides the R1/R2 oxoG asymmetry**. An
   oxoG-damaged library looks like elevated "C>A", indistinguishable from a real C>A process.
2. **Picard `USE_OQ=true`** silently thresholds on pre-BQSR qualities when the OQ tag exists.
3. **Picard's insert-size 60-600 filter** drops short cfDNA/FFPE fragments. Those are the fragments most affected by
   damage, so artifact rates are underestimated.
4. **No overlap dedup** in SamLocusIterator. In short-insert libraries, R1 G>T and R2 C>A at the same molecule base are
   counted twice as independent evidence.
5. Somatic calling impact: Mutect2's orientation-bias filter learns from F1R2 counts. A mis-estimated artifact prior
   flips low-VAF G>T calls between PASS and filtered.

#### 5b. Duplex / error-corrected sequencing (fgbio CallDuplexConsensusReads, NanoSeq)
- Base-level read-vs-reference applies **after consensus**. Error profiles of consensus BAMs use the same queries as
  5a. Single-strand versus duplex agreement per molecule is an aggregation over the `MI` tag (`-f MI`).
- Consensus calling itself is cross-read and out of scope.
- NanoSeq (Abascal et al. 2021 Nature) variant calls apply many read-level filters: clipping, fragment-end
  exclusion, mismatch-per-read limits. All are expressible as DuckDB predicates over off_5p/off_3p and per-record
  mismatch counts.
- **Rank low**: the comparator is a pipeline, not a single base-level table.

#### 5c. Truncation-based lesion and nick maps (CPD-seq, Damage-seq, HydEn-seq/emRiboSeq, GLOE-seq, END-seq)
- **(a) Signal**
  - Polymerase or ligation stops at a lesion or nick, so the lesion sits at a **fixed offset 5' of the read start**.
  - CPD-seq: the dipyrimidine immediately upstream of the read 5' end (Mao et al. 2016 PNAS).
  - HydEn-seq/emRiboSeq: the rNMP is the base immediately 5' of the read start, on the opposite strand, depending on
    the protocol (Clausen et al. 2015 NSMB).
  - Validation uses the dinucleotide composition at -2..-1 (for example TT enrichment for CPDs).
- **(b) Comparators**
  - No dominant single tool: in-house scripts, bedtools `genomecov -5` / `bamtobed`.
  - Ribose-Map (Gombolay et al. 2019 NAR) for rNMP mapping.
  - Weak comparators.
- **(c) Output**: BED/bedGraph of 5' end counts (0-based start), and composition tables at relative positions.
- **(d) Defaults** vary. Ribose-Map uses MAPQ filters and UMI dedup. bedtools uses the aligned start, so soft clips are
  ignored.
- **(e) Expressibility**
  - **Partially.** For single-end R1, `read=__N / refr=YYN` anchored at column 2 (or on a pad column in parquet)
    fires only at the walk start, and the walk start is R1's 5' end.
  - For PE R2 under the directional walk, the left pad is R2's 3' end.
  - `_` also denotes a contig edge, and a skipped soft clip sits between the flank and the first aligned column. The
    truncation site is therefore the **aligned** start, which matches bedtools convention but not every protocol.
  - **Minimal syntax extension: end-specific pad codes** (e.g. `<` = off the read's 5' end as sequenced, `>` = off
    its 3' end, distinct from a contig-edge pad). With these, "base immediately 5' of the read start" is
    expressible for any mate and strand, and in BAM-tag mode.
- **(f) Test data**: CPD-seq and Damage-seq SRA runs are available but large. Simulate by placing read starts 1-2 nt 3'
  of chosen dipyrimidines.
- **(g) Issues**: soft-clipped starts (adapter or UMI remnants) shift the inferred lesion by the clip length. Tools
  that use `pos` versus unclipped start disagree.

---

### Ranking inputs

- **aDNA damage (mapDamage2 / DamageProfiler / PyDamage)**
  - Popularity: very high in paleogenomics; mandatory authentication plot.
  - Comparator cleanliness: good. The TSVs are simple, and mapDamage and DamageProfiler share a header.
  - Ground truth: easy (gargammel).
  - Discrepancy likelihood: **high**, with source-verified causes:
    - three position conventions under soft clips
    - DamageProfiler counting dups, secondary and supplementary alignments
    - DamageProfiler's window-shrink field bug
    - PyDamage ignoring 3' ends
  - Expressibility: **Y**. Optional alignment-column offset output column; optional `--library as-sequenced`.
  - Overall: **top pick** in this family.
- **8-oxoG / FFPE artifact profiles (Picard OxoG / SequencingArtifact, fgbio ErrorRateByReadPosition)**
  - Popularity: high in clinical and somatic pipelines.
  - Comparator cleanliness: very good (flat metrics tables).
  - Ground truth: easy to simulate.
  - Discrepancy likelihood: **high** (fgbio collapse hides asymmetry; USE_OQ; insert-size filter; overlap double
    counting).
  - Expressibility: **Y**. It also shows off the directional-orientation trick: R1 G>T and R2 C>A become one pattern.
  - Overall: **second pick**.
- **SNV pileups (bcftools / samtools / pysamstats / GATK)**
  - Popularity: universal.
  - Comparator cleanliness: good, but BAQ and overlap tweaks make exact parity require flags.
  - Ground truth: easy.
  - Discrepancy likelihood: medium to high (BQ 1/13/0/20; 250 cap; pysamstats overlap and orphan behaviour).
  - Expressibility: **Y**.
  - Overall: an essential baseline, though less novel.
- **CRISPR base editing (CRISPResso2)**
  - Popularity: high in genome editing.
  - Comparator cleanliness: medium. Its own aligner and headers made of repeated amplicon bases mean alignments must
    be pinned via `--bam_output`.
  - Ground truth: very easy (amplicon simulation).
  - Discrepancy likelihood: medium to high (uncovered positions as `-`, window off-by-one PR #651, gap-incentive
    placement, `min_bp_quality_or_N`).
  - Expressibility: **Y** for base editing; **N** for prime-editing classification (out of scope).
  - Overall: a good third or fourth pick.
- **Bisulfite SNP/ASM (BISCUIT; Bis-SNP; Revelio)**
  - Popularity: medium.
  - Comparator cleanliness: BISCUIT is OK (VCF plus vcf2bed). Bis-SNP is poor (GATK3).
  - Ground truth: moderate (Sherman plus SNP genome).
  - Discrepancy likelihood: high (SEQ-orientation end trimming, R2 drop, cross-sample redistribution, 5% cutoff).
  - Expressibility: **Y** via cross-strand join. The optional strand/mate predicate matters only for tag output.
  - Overall: a strong methylation-adjacent demo, and it overlaps other sections.
- **Duplex/NanoSeq error profiles**: popularity growing. Comparator is a pipeline, so cleanliness is poor. Simulation
  is moderate. Expressibility **Y** on consensus BAMs. Overall: low.
- **Truncation lesion maps (CPD-seq, HydEn-seq, ...)**: popularity niche. Comparator weak (scripts, bedtools).
  Simulation easy. Expressibility **partial**; needs **end-specific pad codes** for PE and tags. Overall: low as a
  validation target, but a good motivating example for the pad-code extension.

### Sources
- mapDamage2: https://github.com/ginolhac/mapDamage (tags 2.2.3, master); Jónsson et al. 2013 Bioinformatics 29:1682. Issues #27, #30, #38, #41.
- DamageProfiler: https://github.com/Integrative-Transcriptomics/DamageProfiler ; Neukamm et al. 2021 Bioinformatics 37:3652 (https://academic.oup.com/bioinformatics/article/37/20/3652/6247758); docs https://damageprofiler.readthedocs.io/en/latest/contents/output.html
- PyDamage: https://github.com/maxibor/pydamage ; Borry et al. 2021 PeerJ 9:e11845.
- Oliva et al. 2021 Brief Bioinf 22:bbab076, https://academic.oup.com/bib/article/22/5/bbab076/6217726
- nf-core/eager test data: https://github.com/nf-core/test-datasets/tree/eager
- gargammel: Renaud et al. 2017 Bioinformatics 33:577.
- CRISPResso2: https://github.com/pinellolab/CRISPResso2 (args.json, CRISPRessoCORE.py, plots/data_prep.py, CHANGELOG.md); Clement et al. 2019 Nat Biotech 37:224 (https://www.nature.com/articles/s41587-019-0032-3); demo data README https://github.com/eric-erki/CRISPResso2
- BE-Analyzer: Hwang et al. 2018 BMC Bioinformatics 19:542, https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2585-4
- Broad BE validation pipeline: https://broadinstitute.github.io/be-validation-pipeline/
- bcftools/samtools/htslib source (develop branch; mpileup.c, bam_plcmd.c, sam.c) and NEWS (1.13 BAQ revamp, PR https://github.com/samtools/bcftools/pull/1474); Danecek et al. 2021 GigaScience 10:giab008; Li 2011 Bioinformatics 27:1157 (BAQ).
- pysamstats: https://github.com/alimanfoo/pysamstats ; pysam libcalignmentfile.pyx (pileup defaults).
- GATK CollectAllelicCounts: https://github.com/broadinstitute/gatk (CollectAllelicCounts.java, LocusWalker.java).
- BISCUIT: https://github.com/huishenlab/biscuit (pileup.c, bisc_utils.h, vcf2bed.c); Zhou et al. 2024 NAR.
- Bis-SNP: Liu et al. 2012 Genome Biol 13:R61 (https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-7-r61); https://github.com/dnaase/Bis-tools
- Revelio: Nunn et al. 2022 BMC Genomics 23:477 (https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08691-6); https://github.com/bio15anu/revelio
- Picard CollectOxoGMetrics / CollectSequencingArtifactMetrics source: https://github.com/broadinstitute/picard ; Costello et al. 2013 NAR 41:e67.
- fgbio ErrorRateByReadPosition: https://github.com/fulcrumgenomics/fgbio
- CPD-seq: Mao et al. 2016 PNAS 113:9057 (protocol https://link.springer.com/protocol/10.1007/978-1-0716-0763-3_7); Damage-seq: Hu et al. 2017 PNAS (https://www.pnas.org/doi/full/10.1073/pnas.1706522114); review https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2022.1102593/full ; HydEn-seq: Clausen et al. 2015 NSMB; Ribose-Map: Gombolay et al. 2019 NAR.
- NanoSeq: Abascal et al. 2021 Nature 593:405; Duplex-seq: Schmitt et al. 2012 PNAS.


---

## Family 07 — RBP-fusion editing (TRIBE / HyperTRIBE / STAMP)

### RBP-fusion editing: TRIBE / HyperTRIBE (ADAR A>G) and STAMP (APOBEC1 C>U)

Added to the survey because it is a large, growing family (TRIBE, HyperTRIBE, STAMP, scSTAMP, TRIBE-STAMP, "expanded palette" base editors) whose entire signal is a read-vs-reference substitution, whose comparators are small script pipelines with hard-coded choices, and which ships its own test data.

#### (a) Signal

- **TRIBE / HyperTRIBE**: an RBP fused to the ADAR catalytic domain (HyperTRIBE uses the E488Q hyperactive mutant) deaminates adenosines near RBP binding sites: RNA A>I, read as **A>G in RNA-sense orientation**. HyperTRIBE's ADARcd E488Q has a weaker neighbour preference (5' U / 3' G for wild-type ADAR, i.e. UAG) and edits more sites. Background = wild-type RNA or genomic DNA from the same strain (Drosophila S2 / fly neurons in the original work; McMahon et al. 2016 Cell; Xu et al. 2018 RNA).
- **STAMP** (Brannan et al. 2021 Nat Methods, https://pubmed.ncbi.nlm.nih.gov/33963355/): RBP fused to rat APOBEC1 → **C>U (read C>T in RNA-sense orientation)**, APOBEC1 prefers A/U-rich flanks (weak ACA/UCA preference). Signal compared against an APOBEC1-only control.
- Library strandedness: TRIBE libraries were originally unstranded TruSeq-like, strand taken from gene annotation (A>G on + genes, T>C on − genes). STAMP/SAILOR assume **reverse-stranded (dUTP / TruSeq stranded)** by default.

#### (b) Comparators

| Tool | Version | Install |
|---|---|---|
| HyperTRIBE (rosbashlab) https://github.com/rosbashlab/HyperTRIBE | scripts, no releases ("1.0.0" docs, https://hypertribe.readthedocs.io) | Perl + **MySQL server** + Trimmomatic + STAR + Picard; not on bioconda. High complexity. |
| hyperTRIBER (Rennie et al.) https://github.com/sarah-ku/hyperTRIBER | R package + perl mpileup parser | samtools mpileup + R/DEXSeq; moderate. |
| SAILOR (Washburn et al. 2014; Deffit et al. 2017) https://github.com/YeoLab/sailor | 1.0.4 singularity (original A>I), updated SAILOR inside FLARE (last commit 2024-03-27) https://github.com/YeoLab/FLARE | Snakemake + Singularity containers (samtools 1.3.1, bcftools 1.2). Not on bioconda. Moderate-high. |
| FLARE (Kofman et al. 2023, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10544219/) | downstream cluster caller on SAILOR sites | as above |

#### (c) Output formats

- HyperTRIBE `sam_to_matrix.pl` writes `*.matrix.wig`: `exp  timepoint  chr  pos  A  T  C  G  N  total` with **pos = SAM POS + offset, i.e. 1-based**, counts over the forward reference strand, no strand split. Loaded into MySQL.
- HyperTRIBE `find_rnaeditsites.pl` output (tab, header): `Chr  Edit_coord  Name  Type  A_count  T_count  C_count  G_count  Total_count  A_count_gDNA/wtRNA  T_count_gDNA/wtRNA  C_count_gDNA/wtRNA  G_count_gDNA/wtRNA  Total_count_gDNA/wtRNA  Editbase_count  Total_count  Editbase_count_gDNA/wtRNA  Total_count_gDNA/wtRNA`. `Type` = EXON/INTRON from refFlat; strand from the gene; `Editbase` = G for + genes, C for − genes. Then `convert_editsites_to_bedgraph.py` → bedGraph (0-based start) with edit % and `filter_by_threshold_without_header.pl`.
- SAILOR final: `*.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed` — BED6: `chrom  start(0-based)  end  name  score  strand` where name = `{cov}|{ref}>{alt}|{edit_frac}` built from the CONF line and score = confidence `1 - betainc(G+α, A+β, edit_fraction)` (FLARE/workflow_sailor/scripts/rank_edits.py `as_bed`, `process`). Strand = '+' if ref/alt equals edit_type, '−' if complemented, '0' for multiallelic.

#### (d) Implicit filters / defaults

HyperTRIBE (CODE/trim_and_align.sh, sam_to_matrix.pl, find_rnaeditsites.pl):
- Trimmomatic `HEADCROP:6 LEADING:25 TRAILING:25 AVGQUAL:25 MINLEN:19` (first 6 bases always removed); STAR `--outFilterMismatchNoverLmax 0.07 --outFilterMultimapNmax 1` (unique only); `samtools view -q 10`; Picard `REMOVE_DUPLICATES=true`.
- **`sam_to_matrix.pl` only counts reads whose CIGAR matches `^\d+M$`, `M N M`, `M N M N M` or `M N M N M N M`**; any read with I, D, S, H, or >3 introns is silently skipped (the `else` branch is commented out). No per-base quality threshold at all after trimming.
- Site calling SQL: gDNA/wtRNA `totalcount > 9`; background edit base fraction `< 0.005`; background non-edit base fraction `>= 0.8`; RNA edit count `> 0`; restricted to gene bodies from refFlat and strand from annotation (sites in overlapping antisense genes are evaluated for both). Thresholds: `rnaedit_*.sh` uses edit ≥ 5 % (bedgraph) and read ≥ 20; `Threshold_editsites_20reads.py` uses ≥ 20 reads and ≥ 10 %.
- hyperTRIBER README: `samtools mpileup --max-depth 50000 -Q 30 --skip-indels` (other example `-Q 40`), then DEXSeq-style differential test vs control.

SAILOR as shipped in FLARE (workflow_sailor/Snakefile defaults, scripts/*.py):
- `split_strands.py` routes reads **by exact FLAG value only**: `[16, 83, 163]` vs `[0, 99, 147]` (swapped for reverse-stranded). Every other flag — improper pairs (97/145/81/161), mate-unmapped (73/137/89/153), secondary, supplementary, duplicates-marked — is dropped silently.
- `samtools rmdup -S` (single, default) / `rmdup` (paired).
- `filter_reads.py`: drop reads with any `I` or `D`; drop reads whose junction overhang < 10 (`junction_overhang`); drop reads whose first or last MD mismatch lies < 5 nt from the (soft-clip-trimmed) end (`edge_mutation`, applied to the whole read, not the mismatch); drop reads with > 1 non-target mismatch (`mm_tolerance`), computed from the MD tag; drop secondary.
- `samtools mpileup -R -d 100000000 -E -p -t DP,DV,DPR,INFO/DPR,DP4,SP -g -I` (samtools 1.3.1): BAQ recomputation (`-E`), default **BQ ≥ 13** for DP4, MAPQ 0 allowed, indels skipped; `bcftools view -v snps`; `bcftools call -c -A` (consensus caller — a site must be *called* as a variant by a diploid germline model).
- `filter_variants.py`: DP4 total ≥ 5 (`min_variant_coverage`), ref must be edit_type[0] on the sense split or its complement on the antisense split; `filter_known_snp.py` removes dbSNP BED3 positions; `rank_edits.py` scores with Beta(G+α, A+β), α=β=0, `edit_fraction` 0.01.

#### (e) Expressible in alnbase today?

- **Yes for the site signal.** For a reverse-stranded (dUTP) library with `--library directional`, every read is walked antisense to the RNA, so RNA A>G is `read = "G"`, `refr = "A"` (and STAMP C>U is `read = "A"`, `refr = "G"`); for a forward-stranded library write the sense form. APOBEC1/ADAR neighbour context can be added as flanking columns (written reverse-complemented for dUTP). `qual`, `mapq`, `flag` columns + DuckDB reproduce BQ/MAPQ/coverage/fraction thresholds, and annotation-based strand (TRIBE, unstranded) is a DuckDB join on `refr_pos` (hits on both walk orientations must then be queried: `read G / refr A` or `read C / refr T`).
- SAILOR's **read-level filters** (≤ 1 non-target mismatch per read, no mismatch within 5 nt of read ends, no indels, junction overhang ≥ 10) are per-read aggregations over hits of a "any mismatch" query (`read "/"` relational code) plus `off_5p/off_3p` — DuckDB, not a syntax gap. Junction overhang needs a "distance to nearest junction" which is expressible only crudely (pattern `,~~~~~~~~~~` at read edges); CIGAR can also be parsed in DuckDB from a `-f cigar` field.
- **Nice-to-have syntax**: a `--library reverse` (RNA dUTP) mode so queries are written in RNA-sense orientation; an unstranded mode that walks reads in reference orientation (for TRIBE-style annotation-stranded libraries).

#### (f) Test data

- HyperTRIBE repo `examples/`: `HyperTRIBE_rep1_chr2L.sort.sam.gz` (33 MB) and `S2_wtRNA_chr2L.sort.sam.gz` (30 MB), Drosophila dm6 chr2L, with `rnaedit_wtRNA_RNA.sh` driver — ideal: tiny, a treated/background pair, ships the exact pipeline.
- SAILOR repo: `CWL-SINGULARITY-pipeline-building-code/example/ce11_example_single_end.bam` (~10k reads, C. elegans chrI) + `ce11.chrI.fa` + `ce11_known_SNPs.bed` (paths per README; the current master has moved files, the frozen 1.0.4 singularity bundle copies them out on first run).
- STAMP: GEO GSE155729 (Brannan et al. 2021, RBFOX2-/TIA1-STAMP HEK293T, stranded); FLARE paper data.

#### (g) Correctness issues / contestable choices

1. **HyperTRIBE drops every read with an indel or soft clip** (`sam_to_matrix.pl` CIGAR regex). With STAR's default local soft-clipping, a large fraction of reads is discarded, and it is *not* random: reads ending near a mismatch cluster (hyper-edited Alu-like regions, 3'UTR edit clusters) are preferentially soft-clipped by STAR, so heavily edited regions are depleted → edit % underestimated exactly where binding is strongest. alnbase with no CIGAR restriction would show the difference directly.
2. **SAILOR PE sense bug candidate**: `filter_reads.py` sets `sense = read.is_reverse` (reverse-stranded) without checking read 1 vs read 2, while `split_strands.py` correctly groups 83+163. For paired-end data read 2 (flag 163, forward) of a + strand transcript is labelled antisense, so its genuine A>G edits count as "non-target" mismatches and reads with ≥ 2 edits are discarded (mm_tolerance 1). This suppresses highly-edited read-2s in PE data — worth confirming on a PE STAMP dataset (the pipeline defaults to single-end).
3. **SAILOR relies on `bcftools call -c` to emit a site**: low-fraction editing (1–10 %, the typical STAMP regime) at moderate depth is often genotyped hom-ref by the consensus caller and never reaches scoring, even though `edit_fraction` is 0.01. An alnbase count table has no such gate; the discrepancy is a sensitivity difference that changes which transcripts are "targets".
4. **Edge filter is read-level**: SAILOR drops the whole read if *any* mismatch is within 5 nt of an end, which also removes the read's valid internal edits.
5. **Soft-clip parsing bug (verified by running the function)**: `filter_reads.py::get_softclip` uses `[\w\d]{1}(\d+)S` to detect a right clip, which also matches a *left* clip of ≥ 10 bases (`15S85M` → `(left=0, right=15)`; `12S80M500N8M` → `(0, 12)`; `5S95M` → `(5, 0)` is correct). The read sequence is then trimmed at the wrong end, so MD mismatch offsets index the wrong read bases and the per-read non-target-mismatch count is computed on shifted sequence; reads with long 5' clips (common with STAR end-to-end off, adapters, template-switch oligos) are kept or discarded essentially at random. alnbase walks the CIGAR, so it cannot make this mistake.
6. **Exact-FLAG routing** in SAILOR silently discards improper pairs; with STAR these are a few % of reads and enriched at junctions / 3' ends.
7. HyperTRIBE ignores base quality after trimming (no BQ threshold), while SAILOR uses BQ ≥ 13 with BAQ; the two pipelines disagree at low-quality read tails.

#### Ranking inputs

- Popularity: moderate and growing (TRIBE/HyperTRIBE ~ several hundred citations; STAMP high-profile, used in scSTAMP / "expanded palette" 2023–2024).
- Comparator cleanliness: poor (MySQL, singularity, snakemake), but deterministic and small test data exist.
- Ground truth / simulation: easy (inject A>G into simulated reads; HyperTRIBE background vs treated).
- Discrepancy likelihood: high (CIGAR-restricted counting; PE sense handling; genotype-caller gate).
- Expressibility: Y (site signal); library-mode convenience for RNA orientation.
