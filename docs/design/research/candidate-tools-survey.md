# Candidate tools for the finding hunt (survey, 2026-09-18)

A literature-and-tooling survey for targets we have **not** yet read. Already catalogued
and out of scope here: Bismark, MethylDackel, BISCUIT, ALLCools/methylpy, ScaleMethyl,
meRanGs, BSBolt, BSMAP, HISAT-3N.

The bar this list is screening for: an issue that is **more than a nuisance**. Per
`finding-hunt-method`, that means large enough to matter and small enough to survive the
field's own rate-based QC (conversion efficiency, spike-ins, mapping rate). The strongest
entries below are the ones where the field's QC number is *computed by the same code path
as the result*, so an error moves both in lockstep and cannot be caught.

Citation counts are Europe PMC `citedByCount`, which runs 1.5–3× below Google Scholar —
treat them as floors. Stars as of 2026-09-18. Nothing below has been cloned, read or run.

## Priority candidates

| # | Tool | Assay | Lang | Usage | Mol? | The comparison | Where a bug could hide |
|---|---|---|---|---|---|---|---|
| 1 | **mapDamage2.0** (`ginolhac/mapDamage`) | ancient-DNA deamination | Py+R | 52★, **1225 cites**, nf-core/eager | no | C→T and G→A frequency **by distance from the read's 5′/3′ end**, per strand | Literally alnbase's `off_5p`/`off_3p`. Is the offset measured on the original read including hard clips, or only the aligned block? The 3′ G→A mirror needs a strand flip. **The misincorporation plot *is* the field's QC**, so an error is self-confirming |
| 2 | **Hyper-editing pipeline** (`hagitpt/Hyper-editing`) | A-to-I hyper-editing | Perl | 7★, **218 cites**, basis of the hyper-editing literature | per-read clusters | A→G transform of reads *and* genome, realign, back-map to original sequence | Structurally the same class as a bisulfite aligner: three-letter-space realignment plus coordinate/strand back-transform. 25 KB of Perl, untouched since 2018 |
| 3 | **SlamDunk + Alleyoop** (`t-neumann/slamdunk`) | SLAM-seq 4sU T→C | Py | 49★, 145 cites, bioconda + nf-core + Galaxy | no | per reference T, reads showing T>C, BQ≥27 on the converted base | Counts on the *reference* T while the chemistry acts on the *transcript* strand; sense-only filter; mate overlap undocumented. SLAM-seq's QC is the global conversion rate, computed by this same path |
| 4 | **WASP** (`bmvdgeijn/WASP`) | allele-specific mapping-bias correction | Py+C | 114★, **547 cites**, near-universal in QTL work | yes (rewrites reads) | finds each het SNP's base inside each read and swaps it to the other allele | Maps a *reference* coordinate into *query* coordinate through the CIGAR. An off-by-one swaps the wrong base, which then remaps fine and looks like a legitimate bias filter |
| 5 | **REDItools v1 / v2 / v3** | A-to-I RNA editing | Py2 / Py+MPI / Py3 | 74/60/15★, 312+134+50 cites, the field default | no | per-site base counts vs reference, after strand *inference* | Strand is guessed, not declared. `-a/-A` trims N bases per read end. v2 splits the genome across MPI intervals (boundary double-count/drop). **Three reimplementations of one logic by one group — a disagreement between them on one BAM is a result needing no source reading** |
| 6 | **GATK ASEReadCounter** | allele-specific expression | Java | the ASE format others imitate | no | reads supporting REF vs ALT at het sites | `COUNT_FRAGMENTS_REQUIRE_SAME_BASE` is a *third* distinct mate-overlap rule; deletion-spanning reads are filtered, changing the denominator |
| 7 | **RNAEditingIndexer** (`a2iEditing/`) | global ADAR activity, one number per sample | C+Java+Py2.7 | 48★, **195 cites**, Nat Methods 2019 | no | mismatches over Alu intervals / total reference A covered | `_RefSeqThenMMSites_` decides strand **from the mismatch pattern it is trying to measure** — circular. Output is one number, so a strand flip or overlap double-count leaves no per-site trace. Highest invisibility in the survey, lowest inspectability (ships a PyInstaller binary) |
| 8 | **bam-readcount** (`genome/`) | generic per-base metrics | C++ | 326★, ubiquitous filtering step | no | `reference_base` + per-base counts with plus/minus strand splits | `avg_pos_as_fraction` and `avg_distance_to_effective_3p_end` are computed "after clipping" and must also strand-flip; `-i` excludes insertion-containing reads from base counts |
| 9 | **Picard `CollectSequencingArtifactMetrics`** | OxoG / pre-adapter / bait-bias | Java | 1077★, in every GATK best-practices run | no | read base vs reference base by read number, orientation and trinucleotide context | Docs concede OxoG "reverse-complements all the contexts". `CollectSamErrorMetrics` hand-rolls deletion bookkeeping and shows **no fragment-level overlap correction** |
| 10 | **SAILOR** (`YeoLab/sailor`, archived) | A-to-I, basis of STAMP/FLARE | Py+CWL | 40★ | no | mismatches from the **MD tag**; A>G on `+`, T>C on `−`; drops variants within 5 nt of read edges | MD is reference-consuming and blind to soft clips and insertions — an edge distance computed in MD space is wrong for every clipped read |
| 11 | **snp-pileup** (`mskcc/facets`) | ref/alt counts for FACETS CNV | C++ | 528★, mandatory front end | no | per pileup base: matches REF, matches ALT, else silently "error" | The REF/ALT/error trichotomy **hides** strand and coordinate mistakes in the error bucket. **528 lines in one file — the most readable serious target found** |
| 12 | **Bullseye** (`mflamand/Bullseye`) | DART-seq / TRIBE C→U, incl. single-cell | Perl | 16★; DART-seq = 407 cites; maintained 2026 | per-barcode | per-position nucleotide counts vs genome FASTA or a control | Unstranded mode collapses both strands at a position; no described mate-overlap dedup. Two Perl scripts do the real work |
| 13 | **DSBS analyzer** (`tianguolangzi/DSBS`) | hairpin BS-seq hemimethylation | Py+awk | 5★, Brief Bioinform 2021 | per-fragment | compares read1 to revcomp'd read2 at the same reference position, then both to the reference | The read1↔read2 correspondence is itself a coordinate alignment; any indel or soft-clip asymmetry shifts one strand and mislabels hemimethylation. Three coordinate systems |
| 14 | **dynast** (`aristoteleo/dynast-release`) | scRNA metabolic labeling | Py | 19★ | **yes** (one row per conversion) | every conversion vs reference, then UMI consensus | Docs hand you the hazard: bases are "relative to the forward genomic strand… a read on a reverse-strand gene should be complemented". Double-complement or missed complement. Cleanest source in its domain |
| 15 | **bam2bakR / fastq2EZbakR** | TimeLapse-seq, SLAM-seq, TUC-seq | **bash/awk** | bakR on CRAN; EZbakR PLoS CB 2025 | yes (`cB` table) | tallies user-chosen `[ref base][read base]` pairs; `nT` = reference Ts covered | The only per-base comparator here written in awk. `nT`'s clipping and indel semantics are undefined in the docs |
| 16 | **GLORI-tools** (+ independent `jhfoxliu/GLORI_pipeline`) | m6A via A→G readout | Py | Nat Biotech 2023 | no | A-vs-G at each reference A after converting reads and genome | Single-end only, and the README requires *the user* to reverse-complement reads first — an orientation contract enforced by documentation, not code. **Two independent implementations exist, so a differential test is cheap**; a third-party fork already lists fixes |
| 17 | **asTair** / **rastair** | TAPS (mC→T) | Py / **Rust** | the standard TAPS toolkit; rastair new (Aug 2026) | asTair has per-read utilities | `--method mCtoT` vs `CtoT` **inverts which base means modified** | A strand error silently flips the meaning of every call rather than erroring. `--start_clip/--end_clip` is coordinate arithmetic; the authors' own `IDbias` module suggests they know indels are a problem |
| 18 | **BS-SNPer** (`hellbelly/BS-Snper`) | SNVs in bisulfite data | Perl | 42★, 75 cites | no | distinguishes a real C→T SNP from bisulfite conversion, per strand | The whole tool is one strand-dependent disambiguation; a strand error converts every C/T SNP call on one strand |

## A cross-cutting test worth doing once

JACUSA2, SAILOR, GRAND-SLAM and parts of REDItools all derive read-vs-reference mismatches
from the BAM **`MD` tag** rather than from the reference FASTA. MD is reference-consuming
and blind to soft clips and insertions, and MD/CIGAR disagreement around indels is a known
cross-aligner problem. One construction — a clipped, insertion-containing read — probably
breaks several at once. SAILOR's "within 5 nt of the read edge" filter measured in MD space
is the sharpest instance.

## Overlap rules are mutually incompatible across the field

At least four distinct mate-overlap rules appear in this survey, which is a figure in
itself: htslib's quality-sum (one mate nullified to BQ 0; retained quality = sum if the
mates agree, 0.8× the higher if they disagree), `perbase`'s MAPQ/BQ winner, GATK's
require-same-base fragment drop, and Picard's apparent absence of any.

## Bench (real candidates, lower priority)

Building blocks: `samtools/bcftools mpileup` (heavily scrutinised, novel finding unlikely),
pysam `pileup()` vs `count_coverage()` (issue #945: the two APIs disagree on covered
positions; defaults differ from the CLI), `perbase` (151★ Rust; its README already claims
three behavioural disagreements with competitors), sambamba `depth base` (613★ D; counts
REF_SKIPs toward depth), `pileup.js` (282★, own coordinate conventions, no htslib safety
net), `Rsamtools::pileup()` (has `left_bins` "from the 5′ end regardless of strand" *and*
strand-aware `query_bins` side by side), cancerit `alleleCount`, `allelecounter` (parses
mpileup *text* including `+N`/`-N` indel tokens).

Ancient DNA: PMDtools (346 cites), DamageProfiler (112 cites, nf-core/eager default),
metaDMG-cpp, epiPALEOMIX, DamMet — same `off_5p`/`off_3p` question as mapDamage.

RNA editing and CLIP: PARalyzer/PARpipe (**306 cites**, most-cited in PAR-CLIP; source
availability unverified), wavClusteR, SPRINT, L-GIREMI, RES-Scanner2, RED-ML (hardcoded to
GRCh37 plus a synthetic exon-junction contig that must be lifted back), JACUSA2, phASER.

Methylation callers still outside the catalogue: DNMTools/MethPipe `methcounts`, gemBS
`bs_call`, MOABS `mcall`, CGmapTools, BS-Seeker2. Novel chemistries: PRAISE (two near-
duplicate realignment scripts, ships a modified biopython), m5C-UBSseq, RNA-m5C, NT-seq
(bacterial 6mA+4mC+5mC, where A→G on one strand is T→C on the other), BACS, PP5mC,
HBS-tools, BSPAT.

## Excluded after checking

kb-python `--workflow=nac` (pseudoalignment, no per-base comparison); hts-nim-tools,
rust-bio-tools, alignmentSieve, mosdepth, megadepth (coverage, not base-vs-reference);
nanodisco, 6mASCOPE, SMAC, 6mA-Sniper (compare current/IPD to a model); RDDpred, DeepRed,
RNAEditor (classifiers over another tool's calls); ribosome-profiling P-site tools
(read-position-to-codon, not base comparison). CUT&RUN/ATAC: none found doing
read-base-vs-reference — reported as "none found", not "verified absent".

## Not verified

Google-Scholar-scale citation counts for the newest chemistry papers; ASEQ (Linux binaries
only, no repo); **BID-seq has no official repository** — only a third-party
reimplementation (`b-psid`) exists, which is arguably itself a finding; eTAM-seq (no repo
found); CeU-seq, RBS-seq, HPoxBS (likely email-the-authors software); the asTair Bitbucket
page would not render.

## Two notes beyond the bug hunt

**Prior art for the paper's related work.** `wgbs_tools bam2pat` (198★; underpins the
Loyfer *Nature* 2023 methylation atlas, 556 cites) turns each read into a string of C/T
calls indexed by a **global CpG index** rather than a genomic coordinate, and `jvarkit
sam2tsv` (525★) emits one row per aligned read base ↔ reference base. These are the two
closest existing things to alnbase's per-read primitive and should be read as prior art.

**A competitor, not just a target.** `rastair` (`bsbludwig/rastair`) is a Rust integrated
SNP + methylation caller for mC→T chemistries, created August 2026 with a bioRxiv preprint.

**Supporting evidence for the generality selling point.** Molecule-level output is rare in
every domain surveyed. Outside a handful of tools — MARINE, L-GIREMI, dynast, fastq2EZbakR,
NASC-seq2, phASER, wgbs_tools, sam2tsv — essentially everything here is pileup- or
gene-level.
