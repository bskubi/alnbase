# Methylation callers: hardcoded behaviours, documentation status, suspected bugs

Research date: 2026-09-16. This builds on `methylation-tools.md` (bulk extractors) and `singlecell-formats.md` (single-cell formats) in this folder. Those files record configurable defaults and output formats; this one covers only what **no option can change**.

The summary table, the cross-tool table and the ranked bug list are at the end of the file (they are written last). The per-tool sections come first.

**Rules used in this file:**
- **Hardcoded** means no CLI option, config file or function argument changes it. When an option controls only part of a behaviour, the uncontrollable part is hardcoded. Option defaults go in each tool's short "Configurable defaults" list.
- **Documented?** "Undocumented" means the behaviour is absent from the CLI help/usage text, the README, the online docs or manual, and the paper; the table says where we looked. Anything else is cited.
- **Class:**
  - **bug**: output contradicts the tool's own documented intent, or basic correctness.
  - **contestable**: defensible, but open to dispute.
  - **neutral**: a convention.
- **Scope:** directional libraries only. Items that matter only for non-directional or PBAT libraries are listed per tool under "Future: non-directional".
- **Mate overlap** is excluded from the planned comparisons, so each tool gets one line on it. Supplementary and chimeric handling is recorded in more detail, because it matters for Methyl-HiC.
- Citations are `repo/path:line` at the commit named in each section heading. "Verified by execution" is used only where the tool's own code was run; everything else comes from reading the code. No tool was compiled or run in Docker when the inventory was written; the MethylDackel bugs were later verified against the source build (`validation/demos/methyldackel-bugs`).
- Prose says "CG" rather than "CpG", except when quoting a tool.

---

## Bismark (Perl v0.25.1 scripts; Rust suite 3.1.0; commit `ec5cb58`)

Scope: the aligner `bismark` (the only place where context and call are decided, written into `XM/XR/XG`), `bismark_methylation_extractor` (parses XM only), and the hardcoded parts of `bismark2bedGraph` / `coverage2cytosine` that change which calls survive. Citations are `Bismark/<path>:<line>` at `ec5cb58`. Docs checked: `Bismark/docs/src/content/docs/**` (source of felixkrueger.github.io/Bismark), `CHANGELOG.md`, the `print_helpfile` text in each script, and the paper (Krueger & Andrews 2011, *Bioinformatics* 27:1571, PMC3102221; the full text was read, and it covers only the CG/CHG/CHH classification and the +/- state column). "Verified by execution" means the Perl subroutine or script was run on the synthetic input shown (in a local `bmtest/` directory). No aligner was run.

**Output granularity:** molecule-level. The aligner writes a per-read call string (`XM:Z:`) into the BAM. The extractor writes one line per call with the read ID (`<read_id> <+|-> <chr> <pos> <Z|z|X|x|H|h>`), so per-read and per-fragment (via the shared QNAME of the two mates) reconstruction is possible. Pileup outputs (`.bismark.cov`, bedGraph, cytosine report) come from `bismark2bedGraph` / `coverage2cytosine`.
**Single-cell:** not native. There is no barcode or cell handling. The docs recommend one run per cell (`docs/.../faq/single-cell-pbat.md`). `deduplicate_bismark --barcode` only adds a UMI from the read name to the dedup key.
**Reference matching:**
- *Aligner.* It reads the FASTA files in the genome folder: all `*.fa`, else `*.fa.gz`, else `*.fasta`, else `*.fasta.gz` (`bismark:5031-5046`). The contig name is the first whitespace-delimited token after `>` (`bismark:5149-5159`). Sequence is uppercased with `uc` and no other normalisation (`bismark:5103`). Duplicate names die (`:5080-5082`).
- The bowtie2 hit's RNAME `<name>_CT_converted` / `_GA_converted` is stripped (dies if the suffix is absent, `:2763-2769`) and looked up **by name** in that hash.
- There is **no length or MD5 check** between the FASTA and the index built from it. If the FASTA is edited after indexing but keeps the same names, calls are silently made against the new sequence.
- If the name is missing from the hash, Perl `substr` on undef yields an empty sequence. The read is then discarded with a per-read warning and counted as "genomic sequence could not be extracted" (`:3126-3131`). Rust returns an error instead (`rust/bismark/src/aligner/methylation.rs:108-113`).
- *Extractor.* It uses no reference; RNAME is copied verbatim.
- *coverage2cytosine.* It uses the same first-token-by-name rule (`coverage2cytosine:1648-1751`), again with no length check. A `.cov` contig absent from the genome is silently omitted (the Perl loop `while ($chromosomes{$chr} =~ /([CG])/g)` on undef only emits an "uninitialized" warning, `:237`). Covered positions that are not C/G in the supplied genome are silently dropped.

**Mate overlap (brief):** R1 wins. R2 is walked from its 5' end and the walk *returns* at the first position reaching R1's end, computed after R1's `--ignore*` trimming (`bismark_methylation_extractor:2396-2417`, drop test e.g. `:2905/2987`). The rule is quality-blind. Dovetailed R2 parts beyond R1 are also dropped. It is on by default and turned off with `--include_overlap`, but the rule itself is fixed.
**Supplementary/chimeric:**
- *Aligner.* bowtie2 PE is forced to `--no-mixed --no-discordant` (`bismark:8044-8045`), so no mixed or discordant pairs, and bowtie2 emits no supplementary records. In minimap2 mode (`--secondary=no`, no `-Y`), extra records with the same QNAME after the first one are read and discarded by the "discard inferior alignments until next read ID" loop (`bismark:2862-2878`). The SA tag is never consulted. Chimeric reads therefore keep only their first-reported (soft-clipped) segment. This was established by reading, not execution.
- *Extractor.* It never reads FLAG or SA. Hard clips are not handled. If the `H` lies at the read-orientation 5' side it **dies**. If it lies at the 3' side, the record is processed normally, because the XM index never reaches the H slots (verified by execution: `3H7M` fwd dies, `7M3H` fwd gives correct positions, `7M3H` rev dies).
- In PE mode, pairing is purely by line adjacency plus a QNAME-equality check. A supplementary or secondary record with the same QNAME is taken as "R2". The desync is normally caught one pair later by a QNAME mismatch (`die`), but it can pass silently if the extra records come in even numbers or sit at end of file.

### Hardcoded behaviours

| ID | Behaviour | Source | Documented? | Class |
|---|---|---|---|---|
| BM-1 | **Context comes from the read-projected reference.** CIGAR `D` bases are not appended, so a C next to a deletion takes its context from the first reference base after the deleted block (e.g. `CAATGG` with `AAT` deleted is called CG). | `bismark:4327-4385` (SE), `:4530-4702` (PE); `methylation_call` `:4833-4855` | Yes: `docs/.../faq/changing-context.md` ("context would effectively have changed from CHH to CG") | contestable |
| BM-2 | **An insertion or soft clip next to a C gives unknown context.** `I` and `S` bases are padded with `X`. An `X` at +1 or +2 (read orientation) gives `U/u`, which the extractor drops (BM-17). | `bismark:4337-4356`, `:4839-4843`, `:4847-4850` | Insertions: yes (`faq/changing-context.md`, "insertions are padded with `X` ... context `Unknown`"). Soft clips: not stated (checked the same FAQ, `usage/alignment.md` local-alignment note, help text, paper) | contestable |
| BM-3 | **N in the reference at +1 or +2 gives `U/u`** (CN, and CHN including CNG), so the call is dropped by the extractor. | `bismark:4839-4852`, GA mirror `:4913-4931` | Yes: XM legend "C in Unknown context (CN or CHN)" (`usage/alignment.md:106-107`; extractor help `bismark_methylation_extractor` help text) | neutral |
| BM-4 | **IUPAC ambiguity codes other than N in the reference context count as H.** Genome preparation converts `[^ATCGN]` to N for the index (`bismark_genome_preparation:463`), but the aligner's in-memory genome only applies `uc` (`bismark:5103`). So C-S-x gives CHH, C-R-G gives CHG, C-Y... gives CHH, although S or R may be G. | `bismark:5103`, `:4833-4865` | Undocumented (checked: `usage/genome-preparation.md`, `options/genome-preparation.md`, aligner help, `faq/*`, paper). Verified by execution: `methylation_call('CAAA','CSAAAA','CT')` gives `H...`; `('CAAA','CRGAAA')` gives `X...`; `('CAAA','CNGAAA')` gives `U...` | **bug** (low severity) |
| BM-5 | The reference is uppercased, so soft-masked (lowercase) bases are treated as normal. | `bismark:5103`; `bismark_genome_preparation:458`; `coverage2cytosine:1709,1724` | Undocumented (checked docs, help, paper) | neutral |
| BM-6 | The read is uppercased before conversion and calling. | `bismark:5597`, `:2444` | Undocumented | neutral |
| BM-7 | **Only read C/T at reference C (CT reads), or G/A at reference G (GA reads), produce a call.** A read `N` or any other mismatch gives `.`. A reference `N` or IUPAC code at the site itself is never called. | `bismark:4833-4899`, `:4907-4973` | Partly: "`.` - not a C or irrelevant position" (`usage/alignment.md:108`); N-in-read not stated | neutral |
| BM-8 | **Context past the read end comes from the reference.** 2 reference bases are appended on the read's 3' side: after the aligned end for OT, before POS for OB, then reverse-complemented. With a 3' soft clip, those 2 bases are separated from the last aligned base by the `X` padding, so they are never used. | `bismark:4314-4323`, `:4387-4397`, `:4523-4541`, `:4606-4702` | Code comment only (`:4295-4296`). Undocumented in docs, help and paper | neutral |
| BM-9 | **Chromosome-edge discard of the whole read or pair.** If the 2 padding bases cannot be extracted (OT: aligned end within the last 2 bp; OB: POS ≤ 2), the entire read or pair is dropped, with every call in it, not only the edge Cs. The run warns per read "Chromosomal sequence could not be extracted" and reports a count. | `bismark:4317-4323`, `:4390-4395`, `:3126-3131`, `:3858-3867`, report `:2036` | Only the report line ("Sequences which were discarded because genomic sequence could not be extracted"); the reason is undocumented (checked `docs/` all pages, help, paper). Rust ports it unchanged | contestable |
| BM-10 | **PE-only off-by-one in the R1 5' padding guard.** R1 uses `($pos_1-2) > 0` while R2 and SE use `>= 0`, so an OB (or CTOB) pair whose R1 has POS = 3 is discarded although 2 padding bases exist. Rust deliberately reproduces it (`rust/.../aligner/methylation.rs` `walk_mate`, `strict_5p` comment). | `bismark:4535` vs `:4622` and `:4317` | Undocumented | **bug** (verified by execution; see bug list) |
| BM-11 | CIGAR ops `H`, `P`, `=`, `X` die at extraction of the genomic sequence (aligner side). | `bismark:4380-4385`, `:4598-4603` | Error message only | neutral |
| BM-12 | **Uniqueness.** A read whose best AS ties at a different position or strand is discarded. When 2 strands hit the same `chr:pos` with the same AS, the first index processed (OT/OB order) is kept. | `bismark:2800-2830`, `:3045-3100`; PE `:3740-3800` | Yes: paper and `usage/alignment.md` ("unique best alignment") | neutral |
| BM-13 | **Strand of origin is fixed at alignment time** (XR/XG from the index that produced the unique best hit). The call is made in read orientation, and XM is stored in SEQ (reference) orientation. | `bismark:4399-4450`, SAM writer `:8489+` | Yes (XR/XG described in `usage/alignment.md`) | neutral |
| BM-14 | **Extractor ignores FLAG, MAPQ and base quality.** No secondary, supplementary, duplicate, QC-fail, unmapped or proper-pair filter exists. A record without XM is silently skipped (SE, `if ($meth_call)`). | `bismark_methylation_extractor:1576-1640`, `:1867-1945` (no flag or qual parsing) | Undocumented as a statement (checked `options/methylation-extraction.md`, `usage/methylation-extraction.md`, help, paper) | contestable (safe for Bismark BAMs) |
| BM-15 | **Extractor orientation comes from XR/XG, not FLAG 0x10.** SE: XM is reversed iff OB/CTOT. PE: R1's XR plus XG pick the index for both mates; R1's XM is reversed for OB/CTOT, otherwise R2's. R2's XR is never validated (see BM-22). | `:1611-1640`, `:1915-1943` | Undocumented | neutral (wrong for foreign BAMs whose FLAG disagrees) |
| BM-16 | **PE pairing by adjacency.** The first line of each 2-line block is "R1" whatever FLAG says. A per-pair QNAME check (stripping `/1`, `/2`) dies on mismatch. Only the first 100,000 pairs are pre-tested. | `:1887-1893`, `:2418-2432`, `:1371-1432` | Partly: runtime error text ("use an unsorted file ... `samtools sort -n`"); docs say input must be Bismark output | contestable |
| BM-17 | `U/u` calls are silently ignored. Any character outside `.ZzXxHhUu` dies. | `:2940-2943` (and the other context branches), `:2972` | Yes: CHANGELOG ("These methylation calls are simply ignored") | neutral |
| BM-18 | **Hard clips in the extractor** (see Supplementary above): H on the read-5' side dies; H on the 3' side works. | `check_cigar_string` `:4228-4346` | Undocumented | neutral (Bismark emits no H) |
| BM-19 | **`--ignore` / `--ignore_3prime` position arithmetic counts soft-clipped bases as reference-consuming**, which shifts all call coordinates of that read. | SE `:1650-1672` (`+`), `:1765-1790` (`-`); PE `:2016-2037`, `:2091-2118`, `:2213-2237`, `:2272-2300` | Undocumented. The docs admit it was never checked whether `--ignore` interacts with soft clipping (`usage/alignment.md:119`) | **bug** (verified by execution) |
| BM-20 | **Overlap rule (when on):** R1 wins wholesale; R2 is truncated from the first position reaching R1's (trimmed) end; quality-blind; applied at every XM index (including `.` and soft-clip slots). | `:2396-2417`, `:2905`, `:2987`, `:3576`, `:3657`, `:3745`, `:3826` | Yes: help and `options/methylation-extraction.md:34` ("only methylation calls of read 1 are used") | contestable |
| BM-21 | **M-bias position** is the index in the XM string after `--ignore` trimming (so positions shift by N) and includes soft-clipped slots. Overlap-dropped R2 bases are not counted. | `:2931` etc. (`$mbias_2{..}->{$index+1}`), return at `:2905` | Undocumented (checked `usage/methylation-extraction.md` M-bias section, help) | contestable |
| BM-22 | `$second_read_conversion = s/\r//;` substitutes on `$_` (the whole R2 line) and stores the match count; R2's XR is never used. There is no output effect (a PE run with a correct R2 record was verified). | `:1905` | Undocumented | neutral (cosmetic defect) |
| BM-23 | **Dead sorted-file check.** `test_positional_sorting` dies only for a header line starting with `@SO`, which never occurs in SAM (sort order is `@HD ... SO:coordinate`). Sorted PE files are caught only by the QNAME check. | `:1398-1402` | Undocumented. Verified by execution: `@HD VN:1.6 SO:coordinate` prints "...passed!" | neutral (defect without output effect) |
| BM-24 | **Call-file semantics.** The `+`/`-` column is methylation state, not strand. Position is 1-based and, for G-strand reads, is the G coordinate. | `:2910-2960` etc. | Yes (paper; `usage/methylation-extraction.md`) | neutral |
| BM-25 | **bismark2bedGraph pooling.** All calls at the same `chr:pos` are pooled regardless of read, strand file or context letter. Lines whose `+`/`-` disagrees with the letter case are warned about and skipped. The percentage is written as an unformatted Perl float. | `bismark2bedGraph:445-470`, `:553-588`, `:590-620` | Partly (format examples in `usage/methylation-extraction.md`) | neutral |
| BM-26 | **coverage2cytosine context comes from the reference only.** `^CG` gives CG, `^C.G$` gives CHG, `^C..` gives CHH, so **N or IUPAC in the H slots counts as H** (CNG gives CHG, CNN gives CHH). This is inconsistent with the aligner (BM-3). | `coverage2cytosine:349-362` | Partly: `faq/changing-context.md` says c2c is "purely assigned based on the reference"; the N rule is undocumented | contestable |
| BM-27 | **coverage2cytosine skips unextractable trinucleotides.** Top-strand Cs in the last 2 bp and bottom-strand Gs in the first 2 bp are skipped, and so is any G at the last base, even if covered. | `coverage2cytosine:291-345` | Partly: CHANGELOG issue #127 (edge fix) | contestable |
| BM-28 | **Genome file discovery.** Only the first matching extension class is loaded (`.fa` shadows `.fasta`). The same rule is used in genome prep, so the index and in-memory genome agree. A mismatched FASTA after indexing is not detected (see Reference matching). | `bismark:5031-5046`; `bismark_genome_preparation:610-622` | Partly: help says files with `.fa`/`.fasta` (±`.gz`) are expected; precedence undocumented | neutral |

28 hardcoded items; 3 true-bug candidates (BM-4, BM-10, BM-19) plus 2 defects without output effect (BM-22, BM-23).

### Configurable defaults (not inventory)
- Aligner:
  - `--directional` (default; `--non_directional`, `--pbat`).
  - bowtie2 end-to-end (no soft clips) unless `--local`.
  - `--dovetail` on for PE unless `--no_dovetail`.
  - `--score_min L,0,-0.2`, `-N 0 -L 20`.
  - HISAT2 `--no-softclip`.
  - minimap2 `-x map-ont|sr|map-pb`.
- Extractor:
  - `--no_overlap` ON for PE (`--include_overlap`).
  - `--ignore`, `--ignore_r2`, `--ignore_3prime`, `--ignore_3prime_r2` all 0 (nf-core methylseq sets `--ignore_r2 2`).
  - Strand-specific 12 files (`--comprehensive`, `--merge_non_CpG`).
  - bedGraph `--cutoff 1`, CG only unless `--CX`.
  - `--zero_based` off.
- coverage2cytosine: `--coverage_threshold 0` (all Cs reported), CG only unless `--CX`, `--merge_CpG` off.
- filter_non_conversion: `--threshold 3`.

### True-bug candidates (with reproduction plans)

**BM-19 (high confidence; verified by execution on the Perl extractor). `--ignore*` shifts coordinates of soft-clipped reads.**
- *Mechanism.* When `--ignore N` trims the 5' end of a `+`-orientation read, or `--ignore_3prime N` trims the left (SEQ-start) end of a `-`-orientation read, the code shifts `start` by `N + D + N_skip - I`. It treats soft-clip (`S`) slots as reference-consuming. The rebuilt CIGAR still begins with the remaining `S` slots, and `check_cigar_string` subtracts one per `S`, so every call in the read lands N bp too far right. The same arithmetic exists for PE R1 and R2 (`--ignore_r2` on an OB-pair R2, which is forward; `--ignore_3prime_r2` on an OT-pair R2, which is reverse).
- *Minimal repro, SE.* Header `@SQ SN:chr1 LN:100`. Record `r1 FLAG 0 chr1 POS 11 MAPQ 42 CIGAR 3S7M SEQ AAACAATAAA XM:Z:...Z..z... XR:Z:CT XG:Z:CT`.
  - `bismark_methylation_extractor -s --comprehensive`: `Z` at 11, `z` at 14 (correct).
  - `... --ignore 2`: **`Z` at 13, `z` at 16** (wrong).
  - Reverse analogue: `FLAG 16 ... XG:Z:GA` with `--ignore_3prime 2` gives 13/16 instead of 11/14.
- *PE repro (nf-core default `--ignore_r2 2`).* R1 `FLAG 83 chr1 POS 50 10M XR:Z:CT XG:Z:GA`; R2 `FLAG 163 chr1 POS 11 3S7M XM:Z:...Z..z... XR:Z:GA XG:Z:GA`. With `-p`, calls are at 11/14; with `-p --ignore_r2 2`, **13/16**.
- *Impact.* Any Bismark `--local` run (recommended for scBS/PBAT, `usage/alignment.md:114-119`; nf-core `local_alignment`) combined with any `--ignore*` option. Calls are written at wrong coordinates, possibly outside CGs. By reading, the Rust 3.x extractor (`rust/.../extractor/call.rs`, which uses `iter_aligned` reference positions) is **not** affected, so Perl and Rust diverge on this path despite the byte-identity claim. Not executed.

**BM-10 (high confidence; verified by execution of the Perl subroutine). PE R1 edge guard off by one.**
- `extract_corresponding_genomic_sequence_paired_end` needs `($pos_1-2) > 0` to prepend R1's 2 padding bases (index 1 = CTOB, 3 = OB). R2 and SE use `>= 0`.
- *Repro.* Contig `chr1` of 30 bp, OB pair (index 3): R1 `FLAG 83 POS 3 10M`, R2 `FLAG 163 POS 3 10M`, consistent CT (R1) / GA (R2) sequences.
  - Expected: both mates get a 12-base genomic window and are called.
  - Observed: R1 window length 0, so the pair is discarded ("could not be extracted").
  - The SE OB read at POS 3 is accepted (length 12). At POS 4 both work.
- End-to-end repro: bowtie2-align a synthetic OB pair to a contig so that the reverse-mapped R1 starts at base 3. `--dovetail` (default) permits R2 POS ≥ R1 POS. Check "Sequences which were discarded because genomic sequence could not be extracted: 1" and the absence of the pair in the BAM.
- *Impact.* Tiny: one position per contig start, only for OB/CTOB pairs. Relevant to small spike-in contigs (lambda, pUC19). Rust reproduces it on purpose.

**BM-4 (medium confidence as a "bug"; mechanism verified by execution). IUPAC codes in the reference context count as H.**
- Genome preparation declares ambiguity codes to be N (`bismark_genome_preparation:461-463`), and the XM legend defines U as "Unknown context (CN or CHN)". But the calling genome keeps the codes, so a C followed by `S` (C/G) or `K` (G/T) is called CHH instead of U.
- *Repro.* Reference `chr1: AAAACSAAAAAAAA...`; OT read `FLAG 0 POS 5 CIGAR 10M` with `C` at position 5.
  - Expected (by the tool's own N rule): `U` (dropped).
  - Observed: `XM` starts with `H`. With `CRG` it gives `X` (CHG), whereas `CNG` gives `U`.
- *Impact.* Only for references that contain IUPAC codes (some assemblies, consensus or strain references). The extractor then outputs CHH/CHG calls where the tool's intent is "unknown".

### Future: non-directional
- CTOT/CTOB index handling (`--non_directional`, `--pbat` index modifier `bismark:4308-4312`). BM-10 also affects CTOB pairs.
- Non-directional tie-breaking when OT and CTOB hit the same position with different AS (`bismark:2797-2800` "overwrite" comment).
- PBAT uses only the GA-converted read files (`bismark:535`).

---

## BISCUIT `pileup` and `epiread` (1.10.3-dev, commit `0a5ceae`, 2026-08-18)

Source: `biscuit/` at `0a5ceae`. Docs: the `gh-pages` branch of the same repo (`biscuit-gh-pages/`, commit `1c0e9ec`, 2026-08-18), which is the source of https://huishenlab.github.io/biscuit/. Paper: Zhou et al. 2024, NAR 52(6):e32, https://pmc.ncbi.nlm.nih.gov/articles/PMC11014253/ (full text read through WebFetch; it gives no filter-level detail). "Help" means the usage strings in `src/pileup.c:1012-1059`, `src/epiread.c:1174-1222`, `src/vcf2bed.c:307-330` and `src/mergecg.c:139-157`, which are mirrored in `docs/subcommands/*.md`. "Undocumented" means the behaviour is absent from the help, README.md, every `docs/**/*.md` file (grepped) and the paper.

**Output granularity:**
- `pileup` writes a per-site aggregate (VCF, one row per reference position, one sample column per input BAM).
- `epiread` writes **per-read** rows (epiBED: one row per alignment record). R1 and R2 are separate rows. Supplementary records get their own rows. Mates are never merged into fragments (biscuiteer `readEpibed` can collapse them later; `docs/epiread/epibed_format.md`).

**Single-cell:** Not native. `pileup`/`epiread` ignore the `CB`/`RX` tags that `biscuit align -9` writes (grep finds no `CB` in `src/pileup.c`, `src/epiread.c`). Per-cell output requires one BAM per cell; `pileup` accepts many BAMs and writes one VCF column each.

**Reference matching:**
- The contig name comes from the BAM header (`bam_hdr->target_name[w.tid]`, `src/pileup.c:743`; `src/epiread.c:545`) and is looked up in the FASTA by name through htslib faidx (`faidx_seq_len`/`faidx_fetch_seq`, `src/refcache.h:69-117`). faidx names are the first word of each FASTA header.
- No `LN` or `M5` check is made.
- If the contig is missing from the FASTA, the program exits with "cannot retrieve reference" (`refcache.h:88-95`). This happens for **every @SQ contig, even one with no reads**, because windows are dispatched over all header contigs (`pileup.c:1236-1246`, `epiread.c:1361-1372`).
- If the FASTA contig is shorter than @SQ LN, it dies with a fatal "cannot retrieve"/"outside range" error when a window or read goes past the FASTA end (`refcache.h:57-66, 158-167`).
- If the FASTA contig is longer, or has the same name but a different sequence, the wrong sequence is used **silently**.
- Multi-BAM `pileup` uses only the **first** BAM's header and queries every BAM by that tid (`pileup.c:716-719, 753, 1164`) (BP-7).
- `mergecg` looks up the BED chrom by name (`mergecg.c:190`).

**Mate overlap (brief):**
- Always on unless `-d`.
- Only records flagged READ2 (0x80) lose bases. The dropped bases are those at reference positions in `[max(pos,mpos), min(own_end, mate_end)]`.
- `mate_end` is taken from the `MC` tag, or set to **the record's own reference length when MC is absent** (`pileup.c:784-827`, `epiread.c:681-697, 756-762`).
- The rule is blind to quality and to whether R1's base survived filtering. It does not check that the mate is on the same contig.

**Supplementary/chimeric:**
- 0x800 is **never filtered**, and no option exists (`pileup.c:764-769`, `epiread.c:630-635`). The SA tag is not read.
- A supplementary record passes if it has MAPQ ≥ 40 and, when paired, 0x2 set. Whether 0x2 is set depends on the aligner. BISCUIT's own paired path emits supplementary records (ALT hits) without mate info, so the proper-pair filter removes them (`lib/aln/mem_alnreg_format.c:689-697`).
- A supplementary record takes part in overlap removal **only by its 0x80 flag, its own MPOS and its MC**:
  - An R2 supplementary lying inside the R1 primary's interval is dropped.
  - An **R1 supplementary that lies inside R2's primary interval is never dropped**, so Methyl-HiC-style chimeras are double-counted in one direction.
- Hard-clipped supplementary records (the default for BISCUIT align without `-Y`, `mem_alnreg_format.c:288-289`) are read with a **shifted SEQ index** in `pileup` (BC-2) and **abort** `epiread` (BC-3).

### Hardcoded behaviours

| ID | Behaviour | Source | Documented? | Class |
|---|---|---|---|---|
| BC-1 | Reference matched by @SQ name through faidx. No LN/M5 check. Missing contig is fatal even without reads. Shorter FASTA contig is fatal. Longer contig or different sequence is silent. | `refcache.h:69-117`, `pileup.c:743-746`, `epiread.c:545,579` | Partly. `docs/pileup.md:28`: "You should use the same reference for `biscuit pileup` as was used to create the alignment." Mismatch behaviour is undocumented. | contestable |
| BC-2 | CIGAR `H` advances the query index (`qpos += oplen`) in the pileup walker, `cnt_retention` and `infer_bsstrand`. After a leading H, every base and quality is read from `SEQ[i+H]`, and the last H positions read past SEQ/QUAL into adjacent record bytes. | `pileup.c:872-874`, `bisc_utils.c:112-114,196-198` | Undocumented (help, docs, paper) | **bug** (verified by execution; `validation/demos/biscuit-bugs`) |
| BC-3 | `epiread`'s CIGAR walker has no `H` case, so it hits `default: abort()` ("Unknown cigar 5") on any hard-clipped record that passes the read filters. | `epiread.c:706-1002` | Undocumented | **bug** (verified by execution; `validation/demos/biscuit-bugs`) |
| BC-4 | CIGAR `N` or `P` aborts all walkers ("Unknown cigar") | `pileup.c:875-877`, `epiread.c:999-1001`, `bisc_utils.c:115-118` | Undocumented | neutral |
| BC-5 | Supplementary (0x800) never filtered; SA ignored; no switch | `pileup.c:764-769`, `epiread.c:630-635` | Undocumented for pileup/epiread. Help lists only secondary/dup/improper switches. The paper's "mapped reads include all reads except those … secondary or supplementary" is about QC counting. | contestable |
| BC-6 | BS strand only from tags, in priority `YD` (f/r) > `ZS` (+/−) > `XG` (CT/GA) > inference. FLAG, R1/R2 and orientation are never used. | `bisc_utils.c:208-238` | Partly. `docs/alignment/understand_bam.md:33` describes YD. ZS/XG fallback and order are undocumented. | contestable |
| BC-7 | Inference with no tag: count high-quality (≥ `-b`) ref C→read T and ref G→read A over the whole read. `nC2T >= nG2A` gives BSW, so **ties (including 0/0) become BSW**. | `bisc_utils.c:163-206` | Undocumented. `docs/alignment/QC.md:114-120` documents a *different* rule (0/0 = `u`, ratio ≤ 0.5, else `c`) for `biscuit bsstrand`. | contestable |
| BC-8 | `cnt_retention` has the strand polarity **inverted**. For BSW reads (bsstrand 0) it counts ref G = read G; for BSC reads it counts ref C = read C. Those are non-informative matches, not retained cytosines. It drives `-t` in pileup and epiread and `Rret` in verbose VCF. | `bisc_utils.c:94-98` vs the correct polarity in `bsconv.c:88-94` and `pileup.c:842-858` | Intent documented: help "-t INT Maximum cytosine retention in a read"; `understand_bam.md:19` "retained C's for OT/CTOT reads or G's for OB/CTOB". Behaviour contradicts it. | **bug** (verified by execution; `validation/demos/biscuit-bugs`) |
| BC-9 | Retention count spans all contexts (including CG) over the whole read, with no base-quality or end-distance filter | `bisc_utils.c:76-122` | Undocumented (help says only "cytosine retention") | contestable |
| BC-10 | `-5`/`-3` distances are measured in **SEQ left→right** coordinates (`qpos` 1-based, including soft clips; `rlen = l_qseq`). They are not measured from the read's 5'/3' ends, so the ends swap for reverse-strand reads. | `pileup.c:432`, `epiread.c:739` | Help and `epiread_format.md:82` say "5' end"/"3' end" of the read, which contradicts the code. A TODO (`pileup.c:430-431`) mentions only the soft-clip issue. | **bug** (effective only with asymmetric `-5`/`-3`; verified by execution; `validation/demos/biscuit-bugs`) |
| BC-11 | With no `MC` tag, the overlap interval assumes the mate's reference length equals this record's own length. If R1 is shorter than R2, non-overlapping R2 bases are dropped. If R1 is longer, overlapping R2 bases are double-counted. Bismark BAMs never carry MC. | `pileup.c:784-796`, `epiread.c:681-697` | Undocumented in help/docs/paper (docs only say "avoided by default"). Code comment `pileup.c:812-816`. GitHub issue #25 (2022) reported this; the fix added MC support only. | **bug** (verified by execution; `validation/demos/biscuit-bugs`) |
| BC-12 | Overlap rule: R2 always loses, whatever its quality. The rule ignores whether R1's base (or the whole R1 record) survived MAPQ/BQ/end filters, so the position can lose its call entirely. | `pileup.c:822-827`, `epiread.c:756-762` | Partly. Help: "-d Double count cytosines in overlapping mate reads (avoided by default)". The R2-drop rule is only in code comments. | contestable |
| BC-13 | Overlap test compares `pos`/`mpos` without checking `mtid == tid`. With `-p`, an R2 whose mate sits at overlapping coordinates on another contig loses bases. | `pileup.c:781,822-827`, `epiread.c:676,756` | Undocumented | **bug** (low impact; verified by execution; `validation/demos/biscuit-bugs`) |
| BC-14 | The AS filter applies only if an `AS` tag exists (Bismark BAMs pass). The threshold is absolute, not length-scaled. | `pileup.c:774-775`, `epiread.c:640-641` | Partly: help "Minimum alignment score (from AS-tag)" | contestable |
| BC-15 | Reference uppercased everywhere (soft-masking ignored) | `refcache.h:153,162-166` | Undocumented | neutral |
| BC-16 | Informative base: BSW reads at ref C (read C = methylated, T = unmethylated). BSC reads at ref G (G/A). Other read bases (A/G/N) are not counted. | `pileup.c:841-858`, `epiread.c:784-899` | Documented: `epiread/epibed_format.md` note 4 | neutral |
| BC-17 | Context comes from the **reference only**. Read bases, deletions and insertions next to the C are irrelevant: a C whose next reference G is deleted in the read is still called CG. | `bisc_utils.c:33-72`, `epiread.c:784-899` | Partly: VCF header `N5` "5-nucleotide context, centered around target cytosine" (`pileup.c:947`). Deletion behaviour undocumented. | contestable |
| BC-18 | Windows are half-open `[beg,end)` with `end = contig length`, so **the last base of every contig is never piled up**. With `-g chr:a-b`, **position b is excluded**. epiread never emits reads whose POS is the contig's last base (whole-BAM mode). | `pileup.c:808,890,1224-1246`, `epiread.c:1343-1371` | Undocumented (help: "-g Region") | **bug** (verified by execution; `validation/demos/biscuit-bugs`) |
| BP-1 | `pileup` context is a reference 5-mer. `N` **anywhere** in it (including 2 bp upstream, and the N padding at contig positions 1, 2, L−1, L) makes `CX=CN`. Non-N IUPAC codes (R, Y, …) count as H. | `bisc_utils.c:33-72` | Undocumented. The VCF header lists "CG, CHH or CHG" only; no CN anywhere in docs. | contestable |
| BP-2 | A reference `N` at the site gives no VCF row | `pileup.c:470-471` | Undocumented | neutral |
| BP-3 | SNP veto: sample is methylation-callable at ref C only if opposite-strand (BSC) read T count is 0, or `T/C < 0.05`. C includes BSW retentions plus BSC C (and Y→C redistribution only when no T). Mirrored at ref G. Threshold hardcoded. `-r` cannot change the outcome (redistribution only happens when T = 0). A single BSC error read at a low-coverage or fully unmethylated C (C = 0) vetoes the site. | `pileup.c:517-534,388-419` | Partly. `docs/pileup.md:10-12` "may revert the methylation call if the SNP interferes". Threshold and rule undocumented for pileup; epiBED docs mention "allele frequency < 0.05" for epiread. | contestable |
| BP-4 | `BT` written as `%1.3f`. `vcf2bed -c` and `mergecg` rebuild `M = round(beta·cov)`, which is **not the true count when CV > 1000**. Verified: 1 of 2001 retained gives `M=0 U=2001` from both. | `pileup.c:666`, `vcf2bed.c:173-176`, `mergecg.c:116-118` | Undocumented (help: "-c Output Beta-M-U") | **bug** (low impact; verified by execution; `validation/demos/biscuit-bugs`) |
| BP-5 | Multi-BAM pileup assumes identical headers. Contig names and tids come from BAM 1; every BAM is queried by tid. | `pileup.c:716-719,753,1164` | Undocumented in help/docs (code comments only) | **bug** (verified by execution; `validation/demos/biscuit-bugs`) |
| BP-6 | `DP` counts every alignment covering the position, including those whose base failed filters, R2 overlap bases and hard-clip-shifted records. `CV` is post-filter. | `pileup.c:621`, `660-666` | Documented: `docs/pileup.md` "DP: raw number of reads covering that position" | neutral |
| BP-7 | `vcf2bed -t cg` keeps rows whose `CX` string is exactly `CG`, so CN-context CGs are dropped. `-t c` keeps any ref C/G. | `vcf2bed.c:157-161` | Partly (`methylextraction.md` lists types; CN undocumented) | contestable |
| BP-8 | `mergecg` pairs rows by reference adjacency (C at i, G at i+1), whatever CX says. A lone C or G in CG context is widened to 2 bp. A row at contig position 1 calls `refcache_getbase(rc, 0)` and aborts ("outside range"). | `mergecg.c:102-103,192-196,206-215` | Merging documented in `methylextraction.md`; widening and position-1 crash undocumented | **bug** (fatal error, exit 1, low impact; verified by execution; `validation/demos/biscuit-bugs`) + contestable |
| BP-9 | Per-observation fields are packed as `qpos`/`rlen` `uint16`, `qual` 7 bits and `cnt_ret` `uint8`, so reads over 65535 bp wrap the end filter | `pileup.h:67-77` | Undocumented | neutral (short-read tool) |
| BE-1 | `epiread` CG test uses only the **adjacent** reference base (C→next G, or G→previous C). N elsewhere and contig edges are irrelevant, so epiread calls M/U at CGs that pileup labels `CN` (dropped by `vcf2bed -t cg`). | `epiread.c:785-899` | Contradicts `epibed_format.md` "Methylation and SNPs match the output from `biscuit pileup` in version 2.0" and help ("should have the same results in the epiBED as in pileup") | contestable (divergence from documented equivalence) |
| BE-2 | epiread has no internal SNP veto. It uses `-B` BED only, and parses **only lines with exactly 9 columns** (`strcount_char(..,'\t')==8`). A multi-sample `vcf2bed -t snp` BED is **silently ignored**. It needs `vcf2bed -s ALL` or several named samples; the default `-s FIRST` writes 9 columns. | `epiread.c:1091` | Veto via BED documented (help and `epibed_format.md` note 1, "unfiltered SNP BED"). The 9-column requirement is undocumented, while `methylextraction.md` says columns 6-9 repeat per sample. | **bug** (silent no-op; verified by execution; `validation/demos/biscuit-bugs`) |
| BE-3 | BED callability is `ref C & (alt != T or AF1 < 0.05)`, where AF1 is `T/(C+T)` printed `%1.2f`. Pileup uses unrounded `T/C < 0.05`, so boundary sites disagree (T=1, C=21: pileup callable, epiread `x`). | `epiread.c:1128-1137`, `pileup.c:524,654` | Threshold documented (`epibed_format.md` "allele frequency is < 0.05"); the formula difference contradicts the documented match | **bug** (low impact; verified by execution; `validation/demos/biscuit-bugs`) |
| BE-4 | One epiBED row per alignment record; mates are not merged; supplementary records produce extra rows with the same name and read number | `epiread.c:1011-1026` | Documented (`epibed_format.md` columns; biscuiteer collapse) | neutral |
| BE-5 | epiBED `start = POS − leading soft clips`. `end = start + l_qseq + n_del − n_ins` (includes trailing soft clips). The RLE has one symbol per query base plus deletions (`P` soft clip, `i` insertion, `d` deletion), so decoded RLE length = `end−start+n_ins`. | `epiread.c:1008-1009`, `948-996` | Documented (`epibed_format.md` notes 2-3). **Doc typo**: it states `end - start = rle_length + n_insertions`, but the code gives `rle_length − n_insertions`. | neutral |
| BE-6 | Reads whose soft-clip-adjusted start is ≤ 0 are skipped with a warning | `epiread.c:236-239` | Undocumented | neutral |
| BE-7 | Base-quality, end-distance and R2-overlap exclusions all encode as `F` (indistinguishable). An N or other non-C/T read base at a CG encodes as `x`, the same as a non-CG base. | `epiread.c:729-762`, `894-896,936-942` | `F` and `x` documented (`epibed_format.md`); the N-at-CG → `x` case is undocumented | neutral |

Counts: 32 hardcoded items (BC 18, BP 9, BE 5 excluding neutral duplicates). **10 true-bug candidates**: BC-2, BC-3, BC-8, BC-10, BC-11, BC-13, BC-18, BP-4, BP-5, BE-2 (plus low-impact BE-3 and the BP-8 crash). All twelve (BC-2, BC-3, BC-8, BC-10, BC-11, BC-13, BC-18, BP-4, BP-5, BP-8, BE-2, BE-3) are verified by execution at `0a5ceae`. BC-19 (the aligner-side `YD:A:u` finding, below) is an eleventh candidate, found by source reading and not yet reproduced.

### Configurable defaults (not inventory)
Both subcommands (`bisc_utils.h:98-116`) and help:
- `-b 20` minimum base quality
- `-m 40` minimum MAPQ
- `-a 40` minimum AS (only when AS is present)
- `-t 999999` maximum retention
- `-l 10` minimum read length (SEQ length)
- `-5 3 -3 3` end distance
- Improper pairs, duplicates and QC-fail are dropped (`-p`, `-u`; qcfail has no switch in help)
- Secondary is dropped (`-c`, hidden in epiread)
- Overlap removal on (`-d`)
- `-@ 3`, `-s 100000`

pileup only:
- Y/R redistribution on (`-r`)
- Genotyping priors `-E 0.001 -M 0.001 -x 0.001 -C 0.01 -P/-Q 0.33333`
- vcf2bed `-k 1 -t cg`; mergecg `-k 0`

epiread only:
- `-L 302` maximum read length (longer reads are fatal)
- Empty-epiread filter on (`-E`)
- modBAM `-y 0.9`

Note that QC-fail filtering has no CLI switch, so it is effectively hardcoded (`filter_qcfail=1`, no getopt case).

### True-bug candidates (with reproduction plans)

All plans assume:
- A FASTA `ref.fa` with one contig `chrT`, faidx-indexed.
- A hand-written SAM with header `@SQ SN:chrT LN:<len>`, converted to BAM, coordinate-sorted and indexed.
- Qualities all `I`.
- Unless stated, pileup runs as `biscuit pileup -m 0 -a 0 -5 0 -3 0 -b 0 ref.fa in.bam` followed by `biscuit vcf2bed -t c`.

**1. BC-2 Hard-clip query-index shift (pileup). Verified by execution (`validation/demos/biscuit-bugs`): BT 0.000 at the C, plus spurious variant rows 5 bp upstream.**
- Ref `chrT` (30 bp): `TTTTTTTTTTCGTTTTTTTTTTTTTTTTTT` (C at 11).
- Read: FLAG `2048` (supplementary, single-end, so the proper-pair filter does not apply), POS 1, MAPQ 60, CIGAR `5H25M`, SEQ `TTTTTTTTTTCGTTTTTTTTTTTTT` (C retained), tag `YD:Z:f`.
- Expected: `chrT 11` `CV=1`, `BT=1.000`.
- Buggy: base read from SEQ index 15 (`T`), giving `BT=0.000`. The quality is also taken from index 15.
- Control: the same read with CIGAR `5S25M` and SEQ prefixed by 5 bases gives BT 1.000.

**2. BC-3 epiread aborts on hard clips. Verified by execution: "Unknown cigar 5", exit 134.**
- Same BAM as bug 1: `biscuit epiread -m 0 -a 0 -5 0 -3 0 -b 0 ref.fa in.bam`.
- Expected: one epiBED row `chrT 0 25 r1 1 + x10Mx14 . x25`.
- Buggy: stderr "Unknown cigar 5" and SIGABRT, with no output.

**3. BC-8 Retention polarity inverted. Verified by execution (the demo uses the clean split: with `-t 5` the fully converted read on a G-rich contig is dropped and the fully retained read on a G-free contig is kept). Impact only when `-t` is set.**
- Ref (40 bp): `AAAAAAAAAAGGGGGGGGGGCCCCCCCCCCAAAAAAAAAA`.
- Read A: FLAG 0, POS 1, 40M, `YD:Z:f`, SEQ = ref with Cs→T (fully converted, 0 retained C).
- Read B: a second record with SEQ = ref unchanged (10 retained C).
- Run pileup with `-t 5`.
- Expected: A kept (CV at positions 21-30 from A: BT 0), B dropped.
- Buggy: A has 10 G/G matches > 5, so it is dropped; B's count is also 10 (its G/G), so it is dropped too.
- For a clean split, use read B′ with ref `AAAAAAAAAAAAAAAAAAAACCCCCCCCCCAAAAAAAAAA` (no G) and SEQ unchanged. Expected: dropped (10 retained). Buggy: count 0, kept, `BT=1.000` at 21-30.

**4. BC-11 Overlap with no MC tag (Bismark-style PE). Verified by execution (both layouts, and the MC control).**
- Ref of 200 bp with a CG at 150-151 and no other C/G nearby.
- R1: FLAG 99, POS 101, 30M, `XG:Z:CT`, no MC.
- R2: FLAG 147, POS 111, 80M, `XG:Z:CT`, no MC, C at 150 retained.
- Expected (R1 covers 101-130): position 150 `CV=1` from R2.
- Buggy: mate end is assumed to be 101+80−1 = 180, so R2 bases 111-180 are dropped and position 150 has no CV/BT.
- Adding `MC:Z:30M` to R2 restores it.
- Reverse case: R1 POS 101 80M, R2 POS 151 30M, CG at 160. Expected CV=1. Buggy: the interval is empty (`rmend=130`), giving `CV=2` from one fragment.

**5. BC-18 Last base and region end dropped. Verified by execution (whole-BAM last base, and `-g chrE:1-20` vs `1-21`).**
- Ref (20 bp): `AAAAAAAAAAAAAAAAAAAC`. Read FLAG 0, POS 1, 20M, `YD:Z:f`, SEQ = ref.
- Expected: a VCF row at `chrT 20` (CX=CN, CV=1).
- Buggy: no row.
- Region variant: ref with C at 10. `-g chrT:1-10` gives no row at 10; `-g chrT:1-11` gives the row.

**6. BP-5 Multi-BAM header order. Verified by execution (the demo uses contigs `chrI`/`chrJ`: the read is reported on the other contig as variant rows).**
- `ref.fa` with `chrA` (all `A`, 100 bp) and `chrB` (a C at 50).
- BAM1 header `chrA, chrB` with no reads.
- BAM2 header `chrB, chrA` with a read on `chrB` (tid 0) covering 50, `YD:Z:f`, C retained.
- Run `biscuit pileup ref.fa bam1.bam bam2.bam`.
- Expected: sample 2 has `chrB 50 CV=1 BT=1`.
- Buggy: BAM2's tid-0 read is processed in the `chrA` window against `chrA` reference `A`. No call is made, the chrB row is missing, and SNP support columns at chrA:50 show a spurious C.

**7. BC-10 5'/3' frame. Verified by execution in `pileup` and `epiread` (the demo uses a `YD:Z:r` read and its CG Gs at 6 and 56).**
- Ref 60 bp with CGs at 5 and 55. Read FLAG 16, POS 1, 60M, `YD:Z:f`, both Cs retained.
- Run `-5 10 -3 0`.
- Expected (the 5' end of a reverse read is at the right): position 55 excluded, position 5 kept.
- Buggy: position 5 excluded, position 55 kept.

**8. BE-2 Multi-sample SNP BED silently ignored. Verified by execution (the BED comes from `vcf2bed -s ALL -t snp`; the default `-s FIRST` writes 9 columns).**
- Take a `vcf2bed -t snp` output from a 2-sample VCF (13 columns) with a C>T SNP at AF 0.5 at a CG, and run `epiread -B`.
- Expected: `x` at that CG.
- Buggy: M/U emitted, the same as with no BED.
- Cutting to the first 9 columns restores `x`.

**9. BP-4 Count reconstruction. Verified by execution (the demo uses 1 retained of 2001: `M=0` from both `vcf2bed -c` and `mergecg -c`; the planned 500 of 1001 is right by chance in `mergecg`, whose `rintf` rounds half to even).**
- One site with CV=1001 and M=500 (BT 0.4995 printed `0.500`).
- `vcf2bed -c` gives M=501, U=500 instead of 500/501.
- Build it synthetically with 1001 single-base reads (or check the arithmetic against `vcf2bed.c:173`).

**10. BC-13 Cross-contig overlap (with `-p`). Verified by execution.**
- R2 on `chrT` POS 100 50M, FLAG 177 (paired, reverse, mate reverse, R2, not proper), RNEXT `chrU`, PNEXT 100, `MC:Z:50M`.
- Expected: counted.
- Buggy: every base is dropped as "overlap".

**Lower priority:**
- BE-3 AF rounding: pileup callable but epiread `x` at T=1, C=21. Verified by execution.
- BP-8 `mergecg` fatal error on a `vcf2bed -t c` row at contig position 1. Verified by execution (exit 1).

### Future: non-directional
- BC-6/BC-7: strand never taken from FLAG. This is library-agnostic, but for tagless BAMs inference ties (0/0) default to BSW, and PBAT/CTOT/CTOB reads with little conversion evidence get mis-stranded.
- `YD:Z:u` is ignored in pileup/epiread (`allow_u=0`), so the code falls through to ZS/XG/inference (`bisc_utils.c:217`).
- NOMe-seq (`-N`) context rules (HCG/GCH, `mergecg -N`) are out of scope here.

---

### `biscuit align`: the strand the aligner knew and discarded (read 2026-09-18)

The `pileup`/`epiread` inventory above treats `YD` as an input. This section reads the
other side — where `biscuit align` writes it. The strand demo
(`validation/demos/biscuit-strand`) exercised the tag, not the aligner that produces it.

**The aligner knows the bisulfite strand exactly.** It is not inferred from observed
conversions at all: it falls out of which of the two converted indices the read aligned
to, and where in the concatenated forward/reverse pac it landed.

```c
#define mem_getbss(parent, bns, rb) ((rb>bns->l_pac)==(parent)?1:0)   /* memchain.c:282 */
...
reg->bss = mem_getbss(parent, bns, reg->rb); /* set bisulfite strand */  /* memchain.c:867 */
```

`reg->bss` is trusted everywhere inside the aligner: pairing refuses to mate two regions
whose `bss` differ (`mem_pair.c:86`) and packs `bss` into the sort key as the top bit
(`mem_pair.c:161`), and the scoring matrix is selected from `parent`
(`mem_alnreg_format.c:70`). Regions whose start and end disagree are dropped outright as
"cross boundary" (`memchain.c:871-874`).

**At output it is suppressed whenever the read shows no conversion.** `bis_bwa_gen_cigar2`
walks the alignment to build MD/NM and counts conversion events on the way; a read T at a
reference C (on a `parent` alignment) or a read A at a reference G (on a non-`parent`
one). Anything else, including a read A at a reference G *on a parent alignment*, is an
ordinary mismatch. Then:

```c
if (n_conv_ct == 0 && n_conv_ga == 0) *bss_u = 1;   /* bwa.c:415-416 */
else *bss_u = 0;
```

and the tag is written as

```c
if (p.bss_u) kputc('u', str);            /* mem_alnreg_format.c:436-437 */
else kputc("fr"[p.bss], str);
```

So `p.bss` — correct, and computed without reference to conversion evidence — is present
in the same struct on the same line, and is not printed. `u` does not mean "the aligner
could not tell"; it means "the read supplied no corroboration", which is a different
claim.

**The one rescue is mate-dependent.** Immediately before the tag is written:

```c
// if the mate has certain bss so should the target read
if (m0 && m0->bss_u == 0) p.bss_u = 0;   /* mem_alnreg_format.c:251-252 */
```

A mate with conversion evidence clears the flag, after which `"fr"[p.bss]` prints the
correct strand. (The comment's intent is to inherit the mate's strand; the code only
clears the flag. That happens to be harmless, because `p.bss` was already right and, for a
proper pair, equals the mate's.) The consequence is that whether a read's strand survives
to the BAM depends on its **mate's** sequence: single-end reads, and pairs where neither
mate converted anything, keep `u`.

**Downstream the `u` is not read as "unknown" — it is ignored and re-guessed.** Every
caller in the repository passes `allow_u = 0`:

| caller | call |
|---|---|
| `pileup` | `get_bsstrand(rs, b, conf->filt.min_base_qual, 0)` (`pileup.c:758`) |
| `epiread` | `get_bsstrand(rs, b, conf->filt.min_base_qual, 0)` (`epiread.c:645`) |
| `qc_coverage` | `get_bsstrand(..., 0)` (`qc_coverage.c:440`) |
| `cinread` | `get_bsstrand(d->rs, b, 0, 0)` (`cinread.c:67`) |
| `bsconv` | `get_bsstrand(d->rs, b, 0, conf->filter_u)` (`bsconv.c:56`) — the only one that can honour `u` |

With `allow_u = 0` the `YD` arm of `get_bsstrand` matches neither `f` nor `r` and falls
through (`bisc_utils.c:212-218`). `ZS` and `XG` are absent on BISCUIT's own BAM, so the
read reaches `infer_bsstrand`, which re-derives the strand from the record
(`bisc_utils.c:163-206`), ending in BC-7's tie rule:

```c
if (nC2T >= nG2A) return 0;   /* forward */
else return 1;                /* reverse */
```

**The re-derivation cannot succeed, by construction.** The condition that produced the `u`
is exactly "zero conversion events", so `nC2T` is 0 for these reads. The outcome is
therefore decided entirely by `nG2A` — read A at reference G — which for a `parent`
alignment is precisely what the aligner classified as *mismatches*:

- no A-at-G mismatch → `0 >= 0` → **BSW (forward)**, whatever the true strand;
- one or more A-at-G mismatches → **BSC (reverse)**.

A single G→A sequencing error, or one heterozygous A/G SNP, flips the whole read to the
opposite conversion strand. Nothing about this is evidence of the read's strand; it is
mismatch noise, and it overrides an answer the aligner had already computed correctly.

### BC-19 (true-bug candidate)

| ID | Behaviour | Source | Documented? | Class |
|---|---|---|---|---|
| BC-19 | `biscuit align` writes `YD:A:u` in place of the known `f`/`r` whenever an alignment has zero conversion events, discarding `reg->bss`. Whether a read keeps its strand depends on its mate's sequence. `pileup`, `epiread`, `qc_coverage` and `cinread` all pass `allow_u = 0`, so `u` is not recognised and the read falls through to `infer_bsstrand`, whose result is forced to a 0-vs-`nG2A` comparison: forward on a clean read, reverse as soon as one A-at-G mismatch appears. | `bwa.c:415-416`, `mem_alnreg_format.c:251-252, 436-437`, `memchain.c:282, 867`, `bisc_utils.c:212-218, 203-205`, `pileup.c:758`, `epiread.c:645` | Undocumented. `docs/alignment/understand_bam.md:33` describes YD as f/r; `u` is not described anywhere, nor is the fact that callers discard it. | **bug** (not yet reproduced) |

**Which reads get `u`, and which of those actually matter (narrowed 2026-09-18).** The
`u` class is reads with no conversion event on their own strand: the unconverted fraction
(bisulfite conversion failure), reads over stretches with no cytosine on the read's own
strand, and reads whose cytosines are all genuinely methylated. Three of those four
sub-classes are largely inert for a methylation caller, and the entry previously
overstated the reach by listing them as if each contributed:

- *Conversion failure* is excluded by the experiment's own QC in any study that reports a
  conversion rate, the same reason meRanGs was rejected.
- *All cytosines methylated* needs every C in the aligned span methylated. At full read
  length that is ~20 cytosines in a mammalian genome, most of them CH and therefore
  converted in ordinary tissue, so it is negligible except over short spans.
- *No cytosine on the read's own strand* is the one that survives — but not for the
  obvious reason. Such a read makes no legitimate calls, so mislabeling looks harmless;
  in fact the mislabel is what gives it a voice. A true BSC read over a window with no
  reference G gets `u`, is relabeled BSW, is then read at reference **C** positions where
  it does have bases, and those bases are unconverted C (bisulfite never touched the top
  strand of that molecule). It emits a **methylated call at every cytosine under it**,
  manufactured from a read that should have been silent.

So the residual condition is **strand-asymmetric cytosine composition within the aligned
span**: own-strand cytosines absent, opposite-strand cytosines present. Short and clipped
alignments are not a separate sub-class but a *multiplier* on this one — a full-length
read needs an unusual G-free window (poly-pyrimidine tracts, (CT)n, subtelomeric CCCTAA,
mostly low-MAPQ and filtered), whereas a 20 bp aligned span needs no G in only 20 bases,
which happens in ordinary sequence roughly one time in a hundred.

What keeps the finding interesting at this reduced size is its shape rather than its
volume: the damage is one-directional (the tie-break defaults to forward, so true-reverse
reads are systematically relabeled while true-forward reads are relabeled only when a
variant or sequencing error intervenes), and the fabricated calls are always "methylated",
always forward-strand, and spatially clustered in G-poor windows — systematic and local,
so it survives averaging in a way scattered noise would not. The plausible venues are
amplicon/targeted bisulfite (one region carries the conclusion), single-cell (one read
decides a site) and CH quantification (not protected by CG strand merging). This is a
sound catalogue entry, not a candidate to overturn a published number.

**Why it is invisible to rate-based QC.** A `u` read has no conversion events by
definition, so it never enters a conversion-efficiency numerator. Placed on the wrong
strand it becomes informative at the opposite reference base (BC-16), where its bases are
unconverted, so it contributes **methylated calls at correct coordinates on the wrong
strand** for positions it never informed about. Global conversion rate, mapping rate and
spike-in figures are all untouched. This is the screening shape we are after: wrong
strand, right coordinate, every rate intact.

CG merging (`mergecg`) hides part of it, since the two strands of a CG are usually
concordant. CH calls, strand-asymmetric analyses and allele-specific work do not get that
protection.

**Reproduction plan (not yet run).** Two single-end reads on a contig chosen so that one
strand carries cytosines and the other does not:

1. A true BSC (reverse-strand) read placed over a window with no reference G in its span,
   so `n_conv_ga == 0` and the aligner writes `u`. Confirm `YD:A:u` in the BAM, then
   `biscuit pileup` and check the calls land at reference C positions (forward strand)
   rather than reference G.
2. The same read with one base changed to create an A-at-G mismatch, showing the call
   flip to reverse — the strand decided by a single mismatch.
3. Controls: the same fragment as a *pair* whose mate carries one C→T conversion, which
   clears `bss_u` via `mem_alnreg_format.c:252` and restores the correct label. This is
   the sharpest form of the demonstration: identical read, correct answer or wrong answer
   depending on its mate.

Two small BISCUIT BAMs already at hand (`bsq/dir.bam`, `bsq/nd.bam`) carry only
`f` and `r`, so the reproduction needs purpose-built input.

## MethylDackel `extract` and `perRead` (0.6.1, commit `3c77bda`)

Docs checked:
- `README.md`
- usage text: `extract.c:570-720`, `perRead.c:227-270`
- GitHub issues #102, #157, #163, #168, #170
- MethylDackel has no paper.

Citations are `MethylDackel/<file>:<line>`. Configurable defaults and output formats are in `methylation-tools.md` §2.

**Output granularity:**
- `extract` is pileup only (per-site bedGraph, methylKit or cytosine report).
- `perRead` is **molecule-level**: one row per alignment record (`read name, chrom, 0-based POS, CG methylation %, n informative`). Mates are not merged. Only CG calls are counted.

**Single-cell:** not native. There is no barcode handling, so it needs one BAM per cell.

**Reference matching:**
- The FASTA sequence is fetched by the BAM `@SQ` name through htslib faidx (`faidx_fetch_seq(fai, hdr->target_name[tid], ...)`, `extract.c:381`, `perRead.c:181`, `MBias.c:147`). faidx takes the first word of each FASTA header.
- **No LN or M5 check.**
- `extract` with a missing contig prints "faidx_fetch_seq returned ... Note that the output will be truncated!" and skips that chunk, then continues (`extract.c:382-387`), so the exit status is 0.
- `perRead` with a missing contig does not check `seqlen`. `isCpG` sees a negative `seqlen` and returns 0, so every read on that contig is written with `0.0` and 0 informative bases. Only htslib's faidx warning appears.
- A FASTA contig shorter than LN gives silently no calls past its end. A longer contig, or a different sequence under the same name, is used silently.

**Mate overlap (one line):** `extract` only. Mates are paired by QNAME in a custom pileup constructor, then per-base quality arbitration runs (agree: winner +20%, loser 0, tie goes to the later mate; disagree: the qualities are subtracted and both can drop below `-p`). `perRead` has none.

**Supplementary/chimeric:**
- 0x800 is dropped by the configurable default `-F 0xF00` in `extract`, and **kept** by `perRead` (`--ignoreFlags` default 0).
- If supplementaries are kept in `extract`, the overlap constructor pairs **any two records with the same QNAME** that are paired and not unmapped. It does not check R1/R2 or 0x800 (`overlaps.c:121-139`). The first two records to enter the pileup get arbitrated and the hash entry is deleted, so a third record (e.g. the real mate after a supplementary) is double-counted.
- SA is never read.
- Hard clips are handled correctly (htslib pileup; `perRead` treats `H` as consuming neither sequence).

### Hardcoded behaviours

| ID | Behaviour | Source | Documented? | Class |
|---|---|---|---|---|
| MD-1 | Context comes from the **reference only**. Read bases, deletions and insertions next to the C have no effect. | `common.c:49-82`, `extract.c:407-418` | Partly: README "Methylation Context" classifies by reference; the indel case is undocumented (checked README, help) | contestable |
| MD-2 | **N in the H slots counts as H** (CNG gives CHG, CNC gives CHH). A C too close to a contig end for its context counts as CHH. | `common.c:49-82` | Yes: README "Methylation Context" | contestable |
| MD-3 | Other IUPAC codes also count as H, including ones that may be G: `CSG` gives CHG, `CS` gives CHH. | `common.c:49-82` | Undocumented (README defines H as "any nucleotide other than G" and mentions only N) | contestable |
| MD-4 | Lowercase (soft-masked) reference bases are accepted as C/G. A reference `N` at the site is never called. | `common.c:49-82`, `extract.c:427-436` | Undocumented (README, help) | neutral |
| MD-5 | **Strand of origin:** `XG` is used only if its value starts with C or G. Otherwise it comes from FLAG under a directional assumption (R1 fwd or R2 rev = OT; R1 rev or R2 fwd = OB). bwa-meth `YD` is ignored. | `common.c:84-116` | Undocumented (README, help; issue #102 discusses it). Code comment "Can't handle non-directional libraries!" | contestable (correct for directional) |
| MD-6 | A paired record with neither 0x40 nor 0x80 set gets strand 0, which triggers `assert` (abort) in `extract` and `mbias`. | `common.c:94`, `:120-124` | Undocumented | neutral |
| MD-7 | Only read C/T at reference C (OT) or G/A at reference G (OB) count. Read N and other bases are ignored. Reads of the other strand at a column feed only the (off-by-default) variant filter. | `common.c:118-134`, `extract.c:427-441` | Partly (README says "evidence for a methylated C") | neutral |
| MD-8 | **`--OT/--OB` left bound excludes position A.** `trimAlignment` masks SEQ indices `0..A-1` (1-based 1..A). The help says "--OT 5,0,0,0 would include all but the first 4 bases", i.e. that position 5 is kept. The right bound is inclusive as documented. mbias's suggested `lthresh = i+2` (`svg.c:281`) carries the same shift. | `common.c:137-172`, `common.c:11-43` (no −1 applied) | Contradicts help `extract.c:676-686` | **bug** (verified by execution; `validation/demos/methyldackel-bugs`) |
| MD-9 | Trimming (`--OT`, `--nOT`) indexes **SEQ as stored** (reference left to right, including soft clips), not read 5' to 3'. For a reverse-strand read, "position 1" is its 3' end. mbias bins the same way, so mbias suggestions agree with `extract`. | `common.c:137-208`, `MBias.c` (`plp->qpos`) | Undocumented: help says "1-based position on a read ... read #1". Open issues #163, #102 | contestable |
| MD-10 | Trimmed bases are set to N with quality 0 **before** overlap arbitration, so the untrimmed mate's base wins the overlap. | `common.c:446-459` | Code comment only | contestable (sensible) |
| MD-11 | Overlap arbitration details (see the one-liner above): a tie goes to the later mate, and a winner's quality is capped at 255. | `overlaps.c:54-119` | README "A note on overlapping reads" describes the rule, but says "-p defaults to 10" (code: 5) and says nothing of the tie rule | contestable (excluded from comparison) |
| MD-12 | Overlap pairing is by QNAME only (see Supplementary). There is no R1/R2 check and no limit of two records. | `overlaps.c:121-139` | Undocumented | contestable (Methyl-HiC relevant) |
| MD-13 | `is_del` and `is_refskip` pileup entries are skipped. Hard and soft clips never contribute (htslib). | `extract.c:423-424` | Undocumented | neutral |
| MD-14 | bedGraph column 4 is `(int)(100*M/(M+U))`, which **truncates**. The README says "The methylation percentage rounded to an integer". | `extract.c:50` | Contradicts `README.md:77` | **bug** (verified by execution; `validation/demos/methyldackel-bugs`) (cosmetic) |
| MD-15 | A missing FASTA contig skips the chunk with a warning and exit status 0. There is no LN/M5 check. | `extract.c:381-387` | Undocumented | contestable |
| PR-1 | `perRead` counts CG only (reference C followed by G, both from the reference), per alignment. It reports 0-based leftmost POS. Reads with no informative base are still written as `0.0 0`. | `perRead.c:16-37`, `:57-73` | Partly: usage lists columns; the 0-based convention and zero rows are undocumented | neutral |
| PR-2 | `perRead` applies **no** trimming, overlap handling, conversion filter or singleton/discordant filter. None has an option. | `perRead.c:185-193` | Undocumented (usage lists none of these) | contestable |
| PR-3 | **`perRead` never applies the NH filter.** Its usage says "By default, if an NH tag is present and its value is >1 then an entry is ignored as a multimapper", and offers `--ignoreNH`. But no NH check exists, and `--ignoreNH` is not in `lopts`, so passing it gives "Invalid option". | `perRead.c:262-264` (usage), `:296-303` (lopts), `:185-193` (no check) | Contradicts usage text | **bug** (verified by execution; `validation/demos/methyldackel-bugs`) |
| PR-4 | **`perRead` low-quality skip does not `continue`.** After skipping a base with Q < `-p`, the *next* base is evaluated in the same iteration, with no quality check and no CIGAR-boundary check. | `perRead.c:57-63` | Contradicts usage "-p Minimum Phred threshold to include a base" | **bug** (verified by execution; `validation/demos/methyldackel-bugs`) |
| PR-5 | **`perRead` `-l BED` and `--keepStrand` act only at chunk level.** A chunk (default 1 Mb) that overlaps any BED interval has all its reads written. `keepStrand` is parsed into the BED structure but no read-level check exists. | `perRead.c:158-172`, `:185-193` | Contradicts usage ("regions for inclusion"; "only metrics from the top strand will be output") | **bug** (verified by execution; `validation/demos/methyldackel-bugs`) |
| PR-6 | Only reads whose POS lies in the current chunk are processed, against reference fetched to chunk end + 10 kb. Alignments spanning more than 10 kb get wrong values. | `perRead.c:176-189` | Yes: usage "incorrect values for alignments spanning more than 10kb" | contestable |
| PR-7 | Missing contig: all reads written with 0 informative bases and no error (see Reference matching). | `perRead.c:181`, `common.c:49-51` | Undocumented | contestable |

Counts: 22 hardcoded items (MD 15, PR 7); **5 bugs**, all verified by execution at `3c77bda` (MD-8, MD-14, PR-3, PR-4, PR-5): `validation/demos/methyldackel-bugs/run.sh`.

### Configurable defaults (not inventory)
- `extract`:
  - `-q 10`, `-p 5`, `-F 0xF00`
  - NH>1 dropped (`--ignoreNH`)
  - singletons and discordant pairs dropped
  - CG only (`--CHG`/`--CHH`)
  - `-d 1`
  - SNP filter off (`--minOppositeDepth 0`)
  - `--minConversionEfficiency 0`
  - no trimming
- `perRead`: `-q 10`, `-p 5`, `-F 0`, `-R 0`, `--chunkSize 1000000`.

### True-bug candidates (reproduction plans)

**PR-4 (verified by execution: the `10M5S` record and the two-consecutive-low-quality variant each give `100.000000 1`).** Reference `chrT` (40 bp) = `AAAAAAAAAACGAAAAAAAAAAAAAAAAAAAAAAAAAAAA` (C at 1-based 11, G at 12).
- Record: FLAG 0, POS 1, MAPQ 60, CIGAR `10M5S`, SEQ `AAAAAAAAAACAAAA`, QUAL `IIIIIIIII#IIIII` (the 10th base has Q2).
- Code path: the 10th base is skipped, then the soft-clipped base 11 (`C`) is evaluated at reference position 11, a CG.
- Expected `perRead`: `r1 chrT 0 0.0 0`. Buggy: `r1 chrT 0 100.000000 1`.
- Variant with CIGAR `10M2D5M`: a read base after the deletion is scored against a deleted reference base.
- Variant with two consecutive low-quality bases: the second is counted despite `-p`.
- Variant with an even-length read whose last base is low quality: `bam_seqi` reads one nibble past SEQ (the high nibble of `qual[0]`; Phred 32-47 decodes as `C`).

**PR-3 (verified by execution).** Two records with the same QNAME: primary FLAG 0 and secondary FLAG 256, both with `NH:i:2` and covering a CG.
- Expected: no rows, as the usage says.
- Buggy: two rows. `MethylDackel perRead --ignoreNH ref.fa in.bam` fails with "Invalid option".

**PR-5 (`-l` verified by execution with reads at 0-based 100 and 900,000 on one contig; `--keepStrand` not run).** Reference with CGs at 100 and 900,000. BED `chrT 50 150`. Reads at POS 100 and POS 900,000.
- Expected `perRead -l bed`: one row. Buggy: two rows.
- `--keepStrand` with BED strand `+` and an OB read inside the interval: expected no row, buggy one row.

**MD-8 (verified by execution: `--OT 3,0,0,0` gives starts 4, 6, 8; `--OT 0,5,0,0` keeps 0, 2, 4, so the right bound is inclusive as documented).** Reference `chrT` = `CGCGCGCGCG` + 30 A. OT read FLAG 0, POS 1, `20M`, SEQ = reference (all C retained).
- Run `extract --OT 3,0,0,0` (help: "include positions from 3").
- Expected: rows at C positions 3, 5, 7, 9 (0-based 2, 4, 6, 8).
- Buggy: rows at 5, 7, 9 only. Position 3 is masked.

**MD-14 (cosmetic; verified by execution: `chrA 0 1 66 2 1`).** At one CG, 2 methylated and 1 unmethylated read. Expected column 4 `67` (README "rounded"). Buggy: `66`.

### Future: non-directional
- Without XG, PBAT and non-directional reads are mis-stranded (MD-5), and there is no YD support.
- The XG path gives 4 strands. The CTOT/CTOB trimming coordinates follow the same SEQ-orientation rule as MD-9.

---

## CGmapTools `bam2cgmap` (commit `afa6e40`), brief

- Pileup only; not single-cell.
- Strand comes from mapping orientation, which is wrong for Bismark PE R2. The usage text documents this.
- Context comes from the reference: `--` for N or within 2 bp of an edge.
- **Pileup entries next to an indel are not counted** (`CGmapFromBAM.c:555-760`, `pl->indel != 0`). This is undocumented and contestable.
- No MAPQ or base-quality filter.
- Contigs are matched by name via faidx.
- Mate overlap is off unless `-O` (first QNAME wins).
- Supplementary records are kept (samtools-0.1.18 default mask).

Details are in `methylation-tools.md` §5.

---

## ScaleMethyl `bin/met_extract.py` (ScaleBio, commit `2a4973c`)

The script was read in full (500 lines). Also checked:
- `modules/dedup_and_extract.nf`, `modules/alignment.nf`, `nextflow.config`
- `docs/*.md` (notably `analysisParameters.md`, `outputs.md`, `sc_dedup.md`)
- `README.md`
- bwa-meth `bwameth.py:461-509`, to establish the tag and orientation semantics
- BSBolt `ea4870e` `bsbolt/External/BWA/bwa.c:280-340`, where XB is generated

No ScaleMethyl paper describing `met_extract.py` was found; the method is described only in the docs above. Citations are `ScaleMethyl/<path>:<line>`.

**Pipeline context:**
- The default aligner is **bwa-meth** (`nextflow.config:20`), run as paired-end with the FASTQs swapped (`bwameth.py ... ${pairs[1]} ${pairs[0]}`, `modules/alignment.nf:44`).
- Upstream `sc_dedup` (a closed-source binary) removes duplicates by barcode plus leftmost position and applies `--min-mapq 10`.
- `met_extract.py` then runs once per contig that has reads (`:392-405`).
- The `--aligner bsbolt` path parses BSBolt's `XB` call string. The `bwa-meth`/`parabricks` path calls methylation itself against the FASTA.

**Output granularity:** pileup per cell. `<sample>.met_{CG,CH}.parquet` has columns `barcode, chr, pos, strand, methylated, unmethylated`, aggregated per (barcode, pos) (`:51-61`). Per-read calls exist only transiently in memory, so there is **no molecule-level output**.

**Single-cell:**
- **Native.**
- The barcode is the regex `:([ACGT]+\+[ACGT]+\+[ACGT]+)$` on QNAME (`:41`). A QNAME that does not match gives a null barcode, and those rows are aggregated under null.
- Outputs are the per-cell `cellInfo.txt` and optional per-cell ALLC, Bismark `.cov` and Amethyst files.

**Reference matching:**
- *bwa-meth path.* The reference is `pyfaidx.Fasta(ref)[chr][start:end]`, keyed by the BAM contig name (pyfaidx uses the first whitespace-delimited header token) (`:235`, `:275`). There is no LN or M5 check. A contig missing from the FASTA raises `KeyError` in a worker, and the Pool re-raises it, so the run fails. A FASTA contig that is shorter or different is used silently: pyfaidx truncates slices at the contig end.
- *BSBolt path.* No reference is used; context comes from the aligner's XB tag.

**Mate overlap (one line):** `unique(subset=["qname","pos"])` keeps **one arbitrary row** per QNAME and position (polars `keep="any"`; `:35-39`). Mates count once, but which mate's call survives when they disagree is not defined.

**Supplementary/chimeric:**
- `met_extract.py` applies **no FLAG filter at all**: secondary, supplementary, QC-fail and duplicate records all pass.
- bwa-meth's chimera heuristic marks chimeric pairs 0x200 but leaves them mapped (`bwameth.py:485-508`), so they are extracted unless `sc_dedup` removes them. `sc_dedup`'s FLAG handling is undocumented (`docs/sc_dedup.md`).
- Supplementary records of the same QNAME are deduplicated with the primary only where they cover the same position. Elsewhere their calls are added.
- Hard clips are handled correctly (pysam `get_aligned_pairs`).
- SA is never read.

### Hardcoded behaviours

| ID | Behaviour | Source | Documented? | Class |
|---|---|---|---|---|
| SM-1 | **bwa-meth path: strand is `read.is_reverse`.** Forward reads are called at reference C (C/T); reverse reads at reference G (G/A). For **paired-end directional** bwa-meth BAMs, the second-in-pair mate maps reverse on OT fragments (and forward on OB) but carries C/T at reference C. Every such mate is therefore called on the wrong strand: every reference G it covers reads as "methylated". `YD` and R1/R2 are ignored. | `:276`, `:188-227` | Undocumented (docs, README, script help) | **bug** (verified by execution; `validation/demos/scalemethyl-bugs`) |
| SM-2 | **bwa-meth path: the CH read filter tests the running ratio after every call** and discards the read at the first moment `mCH/CH > threshold`. A read whose first (leftmost) CH call is methylated is dropped even when its whole-read ratio is far below the threshold. The BSBolt path correctly uses whole-read counts (`:102-106`). | `:287-303` | Contradicts `docs/analysisParameters.md:70` ("Reads with greater than this percentage CH methylation will be discarded") and the `--threshold` help (`:477`) | **bug** (verified by execution; `validation/demos/scalemethyl-bugs`) |
| SM-3 | **Positions are 0-based** (`read.reference_start + y`; BSBolt path: pysam aligned-pair reference positions). They are written unchanged into the parquet, and by the exporters into ALLC, Bismark `.cov` and Amethyst files (`write_allc.py:24-31`, `write_bismark.py:25`, `write_amethyst.py:40`). The ALLC and Bismark `.cov` formats are 1-based. | `:118`, `:225` | Contradicts `docs/analysisParameters.md:50` (ALLC columns "as described here", linking the ALLCools spec, which says 1-based); the base is undocumented anywhere | **bug** (verified by execution; `validation/demos/scalemethyl-bugs`) |
| SM-4 | The ALLC context column is the literal `CG` or `CH`, not a reference k-mer. ALLCools' `generate-dataset` patterns such as `CGN` or `CHN` expand to 3-mers and do not match it. | `write_allc.py:26` | Contradicts the same doc line (the spec's context is e.g. `CGT`) | **bug** (format; verified by execution; `validation/demos/scalemethyl-bugs`) |
| SM-5 | bwa-meth path context comes from the **reference only**: CG if pos+1 is G, CHG if pos+2 is G, otherwise CHH. The reverse strand is mirrored. Output keeps only CG versus CH. N or IUPAC at pos+1/pos+2 counts as H. The reference is uppercased. | `:135-151`, `:192-211`, `:275`, `:49` | Undocumented | contestable |
| SM-6 | **A read whose `reference_start` is below 2 is skipped entirely** (the 2-bp upstream padding cannot be fetched). Forward-strand Cs whose pos+2 is past the contig end are skipped. | `:269-271`, `:197-198` | Undocumented | contestable |
| SM-7 | Deletions and insertions next to a C do not affect context (reference only). Only M/=/X positions are called (`matches_only=True`). | `:278` | Undocumented | contestable |
| SM-8 | Only read C/T (forward) or G/A (reverse) at the reference base counts. A read N is ignored. | `:212-218` | Undocumented | neutral |
| SM-9 | No FLAG filter (secondary, supplementary, QC-fail, duplicate); no base-quality filter. MAPQ and duplicates are left to upstream `sc_dedup`. | `:260-266`, `:93-96` | Undocumented (`docs/sc_dedup.md` covers only duplicate logic) | contestable |
| SM-10 | Overlap and supplementary dedup on (qname, pos) with arbitrary survivor (see the one-liner above). | `:35-39` | Code comment only | contestable (excluded) |
| SM-11 | The CH filter is applied **per mate**, not per fragment, with a strict `>`. `CH_high` counts dropped records per barcode. | `:104-106`, `:300-303` | Partly (`analysisParameters.md:70`) | neutral |
| SM-12 | BSBolt path: strand comes from the first character of `YS` (`W` = +, `C` = −; a missing YS counts as +). The context letter comes from XB (`x/X` = CG, anything else = CH). XB-to-reference mapping uses digit and letter cumulative sums plus the leading soft clip, which **is correct** for BSBolt's query-based XB counts (`bsbolt/External/BWA/bwa.c:285-316`). Records without XB, or with an all-numeric XB, are skipped. | `:95-122` | Undocumented | neutral |
| SM-13 | Output is aggregated by (barcode, pos) with `strand` and `context` taken from the first row. cellInfo `<ctx>_Cov` counts **calls**, not sites. | `:51-58`, `:438-449` | Partly (`outputs.md`) | neutral |
| SM-14 | The Bismark-format `.cov` export writes 5 columns (no end column) with a 0-1 fraction to 2 decimals where Bismark writes `chr start end percent M U`. MethSCAn's default `bismark` format, which `write_bismark.py` links to, fails on it with `IndexError`. | `write_bismark.py:23-29` | Contradicts `docs/analysisParameters.md:49` ("bismark .cov format", `percent_methylated`) | **bug** (format; verified by execution; `validation/demos/scalemethyl-bugs`) |

Counts: 14 hardcoded items; **5 bugs** (SM-1, SM-2, SM-3, SM-4, SM-14), all verified by execution at `2a4973c`: `validation/demos/scalemethyl-bugs/run.sh`. SM-1 was verified on reads aligned by bwa-meth 0.2.7, whose FLAG and `YD` values confirm the orientation assumption.

**Interaction of SM-1 and SM-2.** In a mis-stranded mate, every reference G in CH context reads as methylated CH. So the incremental filter (SM-2) usually discards that mate at its first CH G, and its spurious CG calls are dropped with it. The visible symptom is then roughly half of all mates counted in `CH_high` and CG coverage from one mate only, rather than inflated mCG. Mates that are mis-stranded but have no CH G positions (short reads in G-poor sequence) would leak spurious 100%-methylated calls at bottom-strand CGs.

### Configurable defaults (not inventory)
- `--aligner bwa-meth`
- `--chReadsThreshold 50` (percent, passed as 0.5)
- `--contexts CG,CH`
- `sc_dedup --duplicate-key Leftmost --min-mapq 10`
- `--subprocesses 4`

### True-bug candidates (reproduction plans)

**SM-1 (verified by execution on simulated OT and OB pairs aligned by bwa-meth: at `--threshold 1.0` CH reads 50% in fully converted cells; at the default 0.5 every second-in-pair mate is counted in `CH_high` and dropped).** Reference `chrT` (60 bp) with a CG at 0-based 30 and no other G within ±10 bp.
- Pair: R1 FLAG 99, POS 21, `20M`, `YD:Z:f`, C at 30 retained; R2 FLAG 147, POS 25, `20M`, `YD:Z:f`, SEQ in reference orientation with C at 30 retained. QNAME `x:AAA+CCC+GGG`.
- Run `met_extract.py --aligner bwa-meth --threshold 1.0` (threshold 1 disables SM-2).
- Expected: `pos 30` methylated=1 (the two mates dedup to one).
- Buggy: an extra row `pos 31, strand "-"`, methylated=1. R2 is tested at reference G 31 and its unconverted top-strand G reads as a methylated bottom-strand C. R2's own call at 30 is never made.

**SM-2 (verified by execution).** Forward read, `40M`, `YD:Z:f`, reference with 10 CHH Cs. Only the leftmost C is retained in the read; the other 9 are converted.
- Run with `--threshold 0.5`.
- Expected: kept (ratio 0.1). Buggy: dropped, `CH_high` = 1, no CG calls from the read.
- Control: the same read with the rightmost C retained instead is kept.

**SM-3 / SM-4 / SM-14 (verified by execution; ALLCools `extract-allc --mc_contexts CGN` writes 0 rows, a tabix query at the 1-based C returns nothing, MethSCAn `prepare` raises `IndexError`).**
- Any single forward read covering a CG at 1-based position 101 gives parquet `pos 100`, and ALLC row `chrT 100 + CG 1 1 1`.
- Expected per the ALLC spec: `chrT 101 + CGx 1 1 1`.

### Future: non-directional
- SM-1's orientation rule is also wrong for PBAT and non-directional reads. The BSBolt path relies on `YS` and handles them.

---

## premethyst `bam-extract` (Adey lab, commit `fa5a47b`); sciMETv2 `sciMET_BSBolt2cellCalls.pl` has the same logic

Checked:
- the usage text in `premethyst_commands/bam_extract.pm:21-66`
- `README.md`
- sciMETv2 `README.md`

The sciMETv2 paper (Nichols et al. 2022, *Nat Commun*) describes the pipeline only at workflow level. Citations are `premethyst/premethyst_commands/bam_extract.pm:<line>`.

**Output granularity:**
- Per-cell pileup: `<barcode>.CG.cov` and `.CH.cov` with columns `chr, pos (1-based), pct, t, c`, plus cellInfo.
- Calls are aggregated **per fragment** first (mates pooled in a hash), so a site covered by both mates counts once.
- The temporary `<barcode>.meth` (one line per fragment with the raw call strings) is deleted, so there is **no retained molecule-level output**.

**Single-cell:**
- **Native.**
- The barcode is the QNAME up to the first `:` (`:155`).
- Input must be name-sorted with each barcode's reads contiguous. The per-barcode `.meth` file is reopened with `>` truncation whenever the barcode changes (`:164-177`).

**Reference matching:** none. Context comes only from the aligner's call string: BSBolt `XB` by default, Bismark `XM` with `-B`. RNAME is copied verbatim. A reference mismatch at alignment time propagates silently.

**Mate overlap (one line):** mates are pooled per fragment. A coordinate called in both mates counts once, and it is methylated if **either** mate says methylated (`:257-268`, `:330-369`).

**Supplementary/chimeric:**
- There is **no FLAG filter**. Pairing is by adjacency of identical read IDs (QNAME up to `#`, `:156`, `:183-197`).
- A supplementary or secondary record adjacent to its primary is paired with it as if it were the mate. The real mate then becomes a "single". Chimeric segments' calls are all pooled into fragments.
- Hard clips: in BSBolt mode no issue arises, because XB covers only the aligned query. In `-B` mode see PM-2.

### Hardcoded behaviours

| ID | Behaviour | Source | Documented? | Class |
|---|---|---|---|---|
| PM-1 | **BSBolt mode ignores CIGAR.** `pos` starts at POS and moves by 1 per call letter and by the XB digit counts. BSBolt's XB digits count **query** bases between calls, including inserted bases and excluding deleted ones (`BSBolt bsbolt/External/BWA/bwa.c:285-316`: insertion `meth_pos += len`, deletion adds nothing). Every call after an insertion is shifted right by its length; every call after a deletion is shifted left by its length. | `:352-381` | Undocumented (usage, README, sciMETv2 README) | **bug** (verified by execution; `validation/demos/premethyst-bugs`) |
| PM-2 | **Bismark mode (`-B`) ignores CIGAR.** `pos++` for every XM character, including soft-clipped and inserted bases. Leading soft clips shift all calls right by S; insertions shift right; deletions shift left. | `:334-351` | Usage: "-B Methylation call field is in Bismark format (XM:Z:) Eg. Bismark or UA-Meth as aligner"; position handling undocumented | **bug** (verified by execution; `validation/demos/premethyst-bugs`) |
| PM-3 | **Fragments with zero CH calls are dropped entirely** (all their CG calls too) and counted in cellInfo column 10. Fragments above `-M` are dropped silently, without being counted. | `:255-273` | Usage documents only "-M Max allowed fraction mCH sites methylated"; the zero-CH drop is undocumented | **bug** (contradicts the documented filter's intent; biases against CG-dense, CH-poor fragments; verified by execution, `validation/demos/premethyst-bugs`) |
| PM-4 | The documented option `-m` ("Minimum chromosome size to retain (def = 10000000)") is **never used**: `$minSize` is set but never read, so no contig is excluded. | `:13`, `:36` | Contradicts usage | **bug** (verified by execution; `validation/demos/premethyst-bugs`) |
| PM-5 | Context classes: BSBolt `x/X` = CG; `y/Y/z/Z` = CH. Bismark `z/Z` = CG; `x/X/h/H` = CH; `u/U` ignored. CHG and CHH are merged into CH. | `:337-369` | Partly (the `.cov` naming implies it) | neutral |
| PM-6 | Pairing by adjacent identical QNAME prefix, no FLAG check (see Supplementary). An unmapped record (POS 0) contributes nothing, but still makes its fragment a "pair". | `:183-197`, `:332` | Usage says "rmdup & filtered bam file, name-sorted"; pairing rule undocumented | contestable |
| PM-7 | Overlap union rule: methylated wins on mate disagreement. | `:257-268` | Undocumented | contestable (excluded) |
| PM-8 | A barcode that reappears non-contiguously truncates and overwrites its earlier `.meth` (possibly after a worker has already started on it). | `:164-177` | Partly (name-sorted input required) | contestable |
| PM-9 | The mCH filter pools both mates (`<=` threshold keeps). Positions are 1-based; strands are not merged (a G-strand call sits at the G). `pct` is a percentage with 2 decimals. | `:255-256`, `:280-289` | Partly (`-M` in usage) | neutral |
| PM-10 | `%CALL_CONV` (BSBolt-to-Bismark letter map) is defined but unused. | `:103` | n/a | neutral |
| PM-11 | **The last record of a barcode is discarded when it has no mate.** A record is held until the next one shows whether it is a mate; on a barcode change (`:175`) and at the end of the input (`:204-206`) the held record is dropped without being written. With single-end input every read is lost; with paired-end input a cell loses its last fragment when one mate was filtered (for example by `bam-rmdup`'s per-record MAPQ filter). | `:164-206` | Contradicts cellInfo's `singles` column and `fastq-align`'s optional `-2` | **bug** (verified by execution; `validation/demos/premethyst-bugs`) |

Counts: 11 hardcoded items; **5 bugs** (PM-1, PM-2, PM-3, PM-4, PM-11), all verified by execution at `fa5a47b`. PM-11 was found while building the demo.

### BSBolt XB context (the aligner side, brief; `NuttyLogic/BSBolt` `ea4870e`)
- Context is computed at alignment from the **reference** window (`cseq`, with 2 padding bases on each side), so deletions do not change context (unlike Bismark).
- Bases in the read at reference C (or G) that are neither C/T (or G/A) give no call.
- Templates whose 2-bp padding crosses a contig boundary get no CIGAR (`bwa.c:266` comment).
- The XB digits are query-base counts (see PM-1).
- **FLAG 0x10 is the converted reference strand, not the read orientation** (`bwamem.c:849-869`, `:890`); SEQ is still reverse-complemented by orientation (`:937-950`). Single-end and first-in-pair records of a directional library are consistent; second-in-pair records are not (r2 of an original-top fragment: FLAG 129, SEQ reverse-complemented, `YS:Z:W_G2A`). Tools that take strand from the FLAG, including alnbase, mis-strand those records (verified; `validation/demos/premethyst-bugs`).

### Configurable defaults (not inventory)
- `-M 0.4`
- `-N 10000` (with `-C`), `-P 100`, `-p 0`
- BSBolt `XB` unless `-B`
- `-H` (skip CH output)
- Upstream `bam-rmdup`: MAPQ ≥ 10, dedup key (barcode, chr, POS).

### True-bug candidates (reproduction plans)
All inputs are SAM text piped through `samtools view` (premethyst calls samtools itself). The logic of `load_read` is pure Perl and can be exercised directly.

- **PM-1 (high, from reading plus BSBolt XB semantics; verified by execution, `validation/demos/premethyst-bugs`).**
  - Reference CGs at 1-based positions 105 and 115.
  - BSBolt read: POS 101, CIGAR `8M2I10M`, both Cs methylated.
  - XB = `4X` + `11` + `X` (4 query bases, call, 11 query bases including the 2 inserted, call).
  - Expected `.CG.cov` rows: 105 and 115. Buggy: 105 and **117**.
  - With `8M2D10M` and CGs at 105 and 117, the buggy output gives 105 and **115**.
- **PM-2 (high; verified by execution, `validation/demos/premethyst-bugs`).**
  - Bismark read: POS 101, CIGAR `3S20M`, XM `.......Z...............` (the Z at SEQ index 7, i.e. reference position 105).
  - Expected row 105. Buggy: 108.
- **PM-3 (high; verified by execution, `validation/demos/premethyst-bugs`).**
  - A single-end read covering only CG sites and C-free sequence otherwise, with XB `X4X4X` and no y/z.
  - A single-end read is itself lost to PM-11, so the demo uses a pair whose strand has Cs only in CG context.
  - Expected: 3 CG calls. Buggy: no rows; cellInfo column 10 = 1.
- **PM-4 (high; verified by execution, `validation/demos/premethyst-bugs`).**
  - A read on a 5-kb contig with default `-m`.
  - Expected: excluded (per the usage). Buggy: included.

### Future: non-directional
None beyond the above: strand is implicit in the aligner's call string.

---

## ALLCools `bam-to-allc` (lhqing `c9f7be2` = PyPI 1.1.1; DingWB fork `4a1f8a1`), with methylpy and YAP

What was checked:
- **Documentation:**
  - the docstring and `_doc.py` option help (`ALLCools/_doc.py`), which the CLI shows;
  - the ALLC format page `docs/allcools/start/input_files.md`;
  - the command-line notebooks in `docs/allcools/command_line/`.
- **Upstream behaviour:** the `samtools mpileup` manual (htslib.org) and samtools `bam_plcmd.c`.
- **Paper:** ALLCools has no standalone methods paper. It is described in Liu et al. 2021 (*Nature*), which was not checked for these details.

The two forks are identical in every item below, except the `convert_bam_strandness` default (see §4 of `methylation-tools.md`). Citations are `ALLCools-lhqing/ALLCools/<file>:<line>`.

**Output granularity:** pileup only (ALLC, one row per C with coverage). No molecule-level output.

**Single-cell:** by convention there is one BAM and one ALLC per cell. `bam-to-allc` does no barcode handling itself; `generate-dataset` builds cell-by-region matrices from many ALLCs.

**Reference matching:**
- By name. `samtools idxstats` contigs are checked against the `.fai` index, and **any BAM contig missing from the FASTA raises `IndexError`** ("Make sure you use the same genome FASTA file...") (`_bam_to_allc.py:416-425`).
- No LN or M5 check.
- mpileup (`-f`) and the Python context lookup both read the same FASTA by name. Python seeks to the `.fai` offset and trims `LINEWIDTH-LINEBASES` characters from every line (`:106-117`).

**Mate overlap (one line):** samtools mpileup's default overlap detection (htslib `tweak_overlap_quality`). ALLCools passes no `-x`, so this is hardcoded: for agreeing bases, a hash of QNAME picks the mate that keeps the summed quality.

**Supplementary/chimeric:**
- mpileup's default `--ff UNMAP,SECONDARY,QCFAIL,DUP` **keeps supplementary** (0x800) records. No option is exposed.
- htslib pairs overlap candidates by QNAME among proper-pair records (`htslib/sam.c:6022-6056`). A proper-pair supplementary record can take the mate's slot, leaving the real mate unpaired and double-counted.
- For snm3C-seq, YAP maps split fragments as separate records (see the YAP note below).

### Hardcoded behaviours

| ID | Behaviour | Source | Documented? | Class |
|---|---|---|---|---|
| AC-1 | **`^` plus MAPQ character not stripped from mpileup column 5.** mpileup writes `^` then `chr(MAPQ+33)` at every read start (`samtools bam_plcmd.c:167-170`). The parser strips only indels and then counts `.`/`T` (ref C) or `,`/`a` (ref G). So a read **starting** at a C column with MAPQ 13 (`.`) adds a spurious methylated count, and MAPQ 51 (`T`) a spurious unmethylated count. At a G column, MAPQ 11 (`,`) and 64 (`a`) do the same. MAPQ 11 and 13 pass the default `min_mapq 10`; bwa-meth, bwa and hisat-3n emit such values. The same parser is in methylpy. | `_bam_to_allc.py:207-250`, `:283-285`; DingWB `:303`; methylpy `call_mc_se.py:1548-1600` | Undocumented (docstring, docs, methylpy README) | **bug** (verified by execution; `validation/demos/allcools-bugs`) |
| AC-2 | **Strand = alignment orientation**: `.`/`T` at ref C gives `+`, `,`/`a` at ref G gives `−`. `--convert_bam_strandness` only rewrites FLAG from XG or YZ beforehand. It picks Bismark or hisat-3n mode from the **first read's** tags and raises if neither tag is present, so bwa-meth BAMs are unsupported. | `:241-285`, `:62-90` | Partly: `_doc.py:42-46` ("Set this parameter to True if you are doing PE mapping with bismark or hisat-3n") | contestable (correct for SE Bismark and hisat-3n-converted) |
| AC-3 | Context comes from the **reference** (uppercased): `seq[pos-up : pos+down+1]`. Read bases and indels are ignored. An N in the context is kept literally (`CNG`), so it is not counted as CG or CH by `generate-dataset` patterns like `CGN` unless N matches the pattern. | `:245`, `:275-280` | Partly: `input_files.md` (context column); N handling undocumented | contestable |
| AC-4 | **Reverse-strand IUPAC asymmetry.** The complement dict covers only `ACGTN`. At a ref G whose context contains another IUPAC code (R, Y, S, ...), `KeyError` is swallowed by `except: continue`, so the `−` row is silently dropped. The equivalent `+` row is written with the IUPAC code in its context. | `:178`, `:274-282` | Undocumented | **bug** (low impact; verified by execution; `validation/demos/allcools-bugs`) |
| AC-5 | Rows are skipped when the context is truncated (contig edges) or `cov == 0`. The `try/except` around slicing never fires; truncation is caught by `len(context) == context_len`. | `:251`, `:286` | Undocumented | neutral |
| AC-6 | mpileup is run with hardcoded `-B` (no BAQ) and otherwise default filters: `--ff UNMAP,SECONDARY,QCFAIL,DUP`, **anomalous (non-proper) pairs skipped** (no `-A`), overlap detection on (no `-x`), **max depth 8000 per file** (no `-d`). Deep sites (spike-ins, amplicons, merged pseudo-bulk BAMs) are silently capped. | `:149`, `:164` | Partly: docstring says "via samtools mpileup"; the flags themselves are not in ALLCools docs (samtools manual only) | contestable |
| AC-7 | A reference base N at the site is skipped. Only C/G sites are considered. | `:203-204` | Undocumented | neutral |
| AC-8 | Column 7 (`methylated`) is always `1`. | `:262`, `:297` | Yes (`input_files.md`: "1 if no test") | neutral |
| AC-9 | `--cpu > 1` raises `NotImplementedError`; the TODO notes that region-split ALLCs would overlap. | `:443-444` | Undocumented (option exists in help) | neutral (fails loudly) |
| AC-10 | **`extract-allc --strandness merge` never writes the final buffered row** (no flush after the loop). A `−` row that is last on a chromosome is written unconverted (still `−`, at its own position), while every other lone `−` row is moved to `+` at pos−1. | `_extract_allc.py:19-58` | Undocumented | **bug** (verified by execution; `validation/demos/allcools-bugs`) |
| AC-11 | FASTA reading assumes a trailing newline on every sequence line (`line[:tail]` with a negative tail). If the file has no final newline, the last base of the last contig is cut off, so contexts there are truncated and the rows skipped. | `:106-117` | Undocumented | neutral (edge) |

Counts: 11 hardcoded items; **3 bugs** (AC-1, AC-4, AC-10), all verified by execution at `c9f7be2` with samtools 1.24.

### methylpy `call-methylated-sites` (`5bb95e0`, 1.4.7), differences only
- The mpileup parser is the same as ALLCools', so **AC-1 applies**.
- **Soft-masked (lowercase) reference C/G are skipped.** mpileup prints the reference base as stored (`bam_plcmd.c:444`), and methylpy compares `fields[2] == "C"` without uppercasing (`call_mc_se.py:1545`). This is documented: README "methylpy only considers cytosines that are in uppercase in the genome fasta file (i.e. not masked)". Class: contestable.
- A contig missing from the FASTA: `get_chromosome_sequence` returns None and every pileup line on it is skipped silently (`:1382-1397`, `:1543`). There is no error. Class: contestable.
- The `chr` prefix is stripped by default (`remove_chr_prefix=True`, configurable).
- Column 7 comes from a binomial test against `unmethylated_control` when one is given (configurable).

### YAP / Luo-lab snmC pipeline (`lhqing/cemba_data` `788e83c`)
- The caller **is ALLCools `bam-to-allc`**, with no extra call-level filter.
  - Bismark path (`mapping/Snakefile_template/mc.Snakefile:55-93, 171-179`): R1 and R2 are mapped **single-end separately**, with R1 `--pbat`. Filters are `samtools view -q 10` and Picard dedup. `bam-to-allc` runs without strandness conversion, which is correct because Bismark SE sets FLAG 0/16 by information strand (`Bismark/bismark:8521-8538`).
  - hisat-3n path (`hisat3n/snakefile/mc.smk:154-251`): PE `--directional-mapping-reverse --unique-only`, `samtools view -q 10`, Picard MarkDuplicates, then `bam-to-allc --convert_bam_strandness`.
  - snm3C (`m3c.smk`): split-read chimeric fragments go through the same `bam-to-allc`.
- AC-1 and AC-6 therefore apply to YAP output. hisat-3n MAPQ values are 0/1/60, so AC-1 is largely avoided there. The Bismark SE path can emit MAPQ 11.

### Configurable defaults (not inventory)
- ALLCools: `--min_mapq 10`, `--min_base_quality 20`, `--num_upstr_bases 0`, `--num_downstr_bases 2`, `--convert_bam_strandness False` (lhqing) / `True` (DingWB).
- methylpy: `min_mapq 30`, `min_base_quality 1`, `remove_chr_prefix True`, `sig_cutoff 0.01`, `min_cov 2`.

### True-bug candidates (reproduction plans)
- **AC-1 (minor practical impact: needs a read starting on the C or G with MAPQ 11, 13, 51 or 64; verified by execution: MAPQ 13 gives `1 2` where MAPQ 12 gives `0 1`; 51 at the C and 11 and 64 at the G also add a count; `validation/demos/allcools-bugs`).**
  - Reference `chrT` = `AAAAACGAAAAA...` (C at 6).
  - SE read: FLAG 0, POS 6, **MAPQ 13**, `10M`, SEQ `TGAAAAAAAA` (C converted, i.e. unmethylated), BQ 40.
  - mpileup column 5 at pos 6 is `^.T`.
  - Expected ALLC row `chrT 6 + CGA 0 1 1`. Buggy: `chrT 6 + CGA 1 2 1`.
  - Control: MAPQ 12 (`^-T`) gives the correct row.

  - **Magnitude measured 2026-09-18, on the pipeline path that matters.** The affected MAPQ
    values are exactly 11, 13, 51 and 64, because samtools writes `^` plus `chr(MAPQ+33)` and
    ALLCools counts `.`/`T`/`,`/`a` straight off the raw string (`_bam_to_allc.py:248-249`,
    `:283-284`, no `^` stripping). Bismark cannot reach them often: it discards ambiguous reads
    rather than assigning low MAPQ, and on a repeat-graded reference emitted only
    {32,34,35,36,38,39,42}. bwa-meth can, and bwa-meth is the YAP/snmC path into ALLCools --
    on the same reference it emitted a near-continuous 0-60 spectrum including MAPQ 51 at
    0.087% of records, with 10, 12 and 15 all populated so 11 and 13 are plainly reachable.
    Running `bam-to-allc` at its defaults (`min_mapq=10`, so 11 and 13 pass the filter) over
    30,000 bwa-meth records and diffing its counts against a correct parser on the *same*
    mpileup: **5 of 27,717 cytosine positions wrong (0.018%)**, all from MAPQ 51.
  - **The sign depends on which value occurs.** 11 (`,`) and 13 (`.`) inflate the *methylated*
    count; 51 (`T`) and 64 (`a`) inflate the *unmethylated* count. In the run above only 51
    occurred, so methylated inflation was 0 and coverage inflation 5.
  - **Conclusion: genuinely minor, and not a conclusion-flipping candidate.** The spurious
    count is one per affected read start, so relative inflation of a cell- or region-level
    aggregate is P(affected MAPQ) / read length -- about 0.001% here. It perturbs individual
    sites at the 0.018% rate and cannot accumulate into a cell-level bias. The catalogue's
    original "minor practical impact" judgement is now a measurement rather than a guess.
- **AC-10 (verified by execution).** An ALLC whose last line is a CG `+` row with no `−` partner: `allcools extract-allc --strandness merge --mc_contexts CGN` loses that row. A lone `−` row that ends an earlier chromosome is written unconverted.
- **AC-4 (verified by execution; low impact).** Reference `...CRG...`. A reverse read covering the G gives no `−` row; a forward read covering the C gives `+ CRG`.

### Future: non-directional
- Orientation-based strand (AC-2) also mis-assigns non-directional CTOT/CTOB PE mates. The XG/YZ conversion handles them.

---

## `hisat-3n-table` (HISAT-3N, `DaehwanKimLab/hisat2` branch `hisat-3n`, read 2026-09-18)

Source read: `alignment_3n_table.h`, `position_3n_table.h`, `hisat_3n_table.cpp`. This is the
caller that consumes the `YZ` tag whose polarity `../../../validation/demos/hisat-3n-strand/`
checks; that demo explicitly did not run the caller, and this section closes that gap by
reading it. Nothing here has been executed yet.

The design is unusual and worth stating, because it is what makes the first finding reachable:
`hisat-3n-table` does not pile up by coordinate and then resolve mates. Every position keeps a
`vector<uniqueID>` of the reads that have already contributed to it, keyed by a hash of the read
name, so that the two mates of a pair cannot both be counted where they overlap
(`position_3n_table.h:131-165`). Per-read identity is therefore carried all the way into the
pileup, which is the right shape. The defect is in what it does with it.

### Hardcoded behaviours

- **Mate disagreement voids the position for that read, rather than preferring a mate.** When the
  second mate reports a different conversion status than the first at the same position, both are
  meant to be discarded and the read marked `removed`, so it can never contribute again.
- **Uniqueness is the literal string `"1"` in the MAPQ column** (`:122-128`). Not a numeric
  comparison and not a threshold.
- **Reads whose CIGAR spans more than 500,000 reference bases are dropped** (`:190`), described in
  the comment as an intron-length cutoff. Not in the help text.
- **`--CG-only` marks a position `'?'` when the reference base is a C or G that is not part of a
  CG, and such positions are skipped silently** (`position_3n_table.h:474-478`).

### True-bug candidates (reproduction plans)

- **H3-1: the mate-disagreement removal erases a base by quality-character value, not by read,
  so it frequently removes nothing and leaves the discordant base counted.** A position stores
  its evidence as two plain strings, `convertedQualities` and `unconvertedQualities`, one
  character per contributing read (`:66-67`); the output columns are those strings and their
  `.size()` (`:341-344`). When the second mate disagrees, the code sets `removed = true` and then
  hunts the stored string for a character *equal to the new mate's quality* and erases the first
  one it finds (`:141-158`):

  ```cpp
  if (uniqueIDs[index].isConverted != InBase.converted) {
      uniqueIDs[index].removed = true;
      if (uniqueIDs[index].isConverted) {
          for (int i = 0; i < convertedQualities.size(); i++) {
              if (convertedQualities[i] == InBase.qual) {
                  convertedQualities.erase(convertedQualities.begin()+i);
                  return false;
  ```

  The character that the *first* mate actually deposited is its own quality, which is a different
  observation and in general a different character. So there are three outcomes, and only one is
  the intended one. If the two mates happen to share a quality character, the right number of
  observations is removed. If some *other* read at the position carries the new mate's quality
  character, that read's character is erased instead — harmless for the count, since the string
  is consumed as a multiset, but it corrupts the emitted quality string. **If no read at the
  position carries that character, the loop finds nothing, erases nothing, and the first mate's
  base stays in the count** even though the pair has been judged discordant and the read marked
  `removed`. That third case is the bug: the discarding the code believes it performed did not
  happen.
  - Which case dominates depends on how many distinct quality values the instrument emits, which
    is the uncomfortable part: on a NovaSeq, with four binned quality values, a covered position
    almost always contains a match and the counts come out right; on an unbinned MiSeq or HiSeq
    run, or at low coverage anywhere, the match often fails and the discordant base survives.
    **The same data processed on two instruments gives different answers, and the direction of
    the error is set by the sequencer's quality binning rather than by anything biological.**
  - Reproduction: two hand-written overlapping mates disagreeing at one C, with distinct quality
    characters at that base and no other read covering it. Expected count 0/0; predicted output
    keeps the first mate's base. Then repeat with the qualities made equal to show the count
    correcting itself.
  - Rate-invisible in the sense that matters here: it moves individual site counts without
    touching any global conversion rate, since the retained base can be converted or unconverted
    with roughly equal likelihood.

- **H3-2: the overlap check is dead code.** `checkOverlap()` (`:161-170`) sets an `overlap` flag
  from `location + sequenceCoveredLength >= mateLocation`, and both the function and the flag have
  no callers or readers anywhere in the three files. The per-read `uniqueIDs` bookkeeping runs
  unconditionally instead, which is why H3-1 is on the live path for every read rather than only
  for pairs detected as overlapping.

- **H3-3: `unique` is `mapQ != "1"` as a string comparison** (`:122-128`), so MAPQ 0 is counted as
  a unique alignment and `--unique-only` keeps it. HISAT-3N's own MAPQ 1 for multi-mappers is
  caught, so this is latent rather than active on HISAT-3N output; it is the reason the table
  cannot be pointed at another aligner's BAM.

- **H3-4: the `NM` tag is parsed into the `NH` member** (`:140-141` and again `:151-152`).
  `NH` is never read, so nothing downstream is wrong today, but the SAM tag actually needed for a
  uniqueness decision is being consumed into a dead field under the name of the tag that would
  have provided one.

- **H3-5: latent out-of-range if `YZ` is the last tag on the line.** The field loop ends when
  `find("\t")` returns `npos`, and the trailing block re-runs the tag tests at the final field
  with `endPosition` still `npos`. The `MD` and `NM` arms pass `npos` as a `substr` length, which
  clamps harmlessly; the `YZ` arm evaluates `line->at(endPosition-1)` (`:154`), which throws.
  HISAT-3N's own tag order puts something after `YZ`, so this does not fire on unmodified output.

### Future: non-directional
- `--directional-mapping-reverse` is undocumented and unexercised; whether `YZ` and the table
  agree under it is unmeasured.

---

## BSMAP `methratio.py` (2.90), brief

Read 2026-09-18 to answer the question left open in
`../../../validation/demos/bsmap-strand/README.md`: whether BSMAP's own caller inherits the
aligner's copy-strand misassignment or re-derives a strand.

**It inherits it.** `get_alignment` asserts the presence of `ZS:Z:` and returns `strand[0]` — the
conversion-strand character — with no reference to the FLAG and no re-derivation
(`methratio.py:46-49`, `:89`). `BS_conversion` then selects `('C','T','G','A')` or
`('G','A','C','T')` from that one character alone. So every cytosine of a misassigned
copy-strand read is counted on the wrong strand, at the right coordinate, and the aligner's
`NM` inflation is the only trace. A `-n 1` run carries the aligner's error through to the
output table unchanged.

Otherwise the script is more careful than its age suggests, and nothing else in it rose to a
finding on this reading. Worth recording as checked rather than as suspect:
- The CIGAR fixup handles only `I` and `D` (`:56-67`), which would shift every base of a
  soft-clipped read, but BSMAP emits `%uM`, `%dM%dD%dM` or `%dM%dI%dM` and never a clip
  (`align.cpp:587-589`, `pairs.cpp:362-364`, `:412-414`), so the gap is unreachable from its own
  output. It is reachable from anything else that writes `ZS:Z:`.
- Depth and methylation counters are `array('H')`, so a position saturates at 65535 and the
  reported ratio silently freezes at whatever the first 65535 reads said (`:148-160`). Only
  amplicon-depth data reaches it, and the estimate is already precise there.
- The minus-strand context classifier indexes `refcr[i-1]` and `refcr[i-2]` and correctly guards
  `i == 0` and `i == 1` first (`:227-232`), so it never wraps to the end of the chromosome
  through Python's negative indexing.
- The CT-SNP correction divides by `d1` only inside `if m1 != d1`, which excludes `d1 == 0`
  (`:212-217`).

---
