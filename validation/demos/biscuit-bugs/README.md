# BISCUIT bugs, reproduced

Twelve BISCUIT behaviours that give wrong calls or counts, drop calls or crash, each shown
on a few hand-written reads. They are all the bug candidates found by reading the source
(`docs/design/research/methylation-callers-hardcoded.md`, IDs BC-2, BC-3, BC-8, BC-10,
BC-11, BC-13, BC-18, BP-4, BP-5, BP-8, BE-2, BE-3); this directory confirms each one by
running the tool.

Tool: BISCUIT 1.10.3-dev built from GitHub master `0a5ceae` (`../../envs/biscuit.yaml`,
`../../envs/build_biscuit.sh`), with its bundled htslib 1.18. Checked on 2026-09-17 with
alnbase 0.1.10.

```
BISCUIT_ENV=... PY=python-with-duckdb ALNBASE=... ./run.sh [workdir]
```

`expected.txt` is the output of that run. Every case prints what is promised, then what
BISCUIT does, then alnbase's hits for the same reads from the query file `c.toml` in this
directory, where alnbase computes the same thing. `pileup` runs with BISCUIT's default filters, except that the `-5`/`-3` end
filters are set to 0 so that every aligned base can be called.

| ID | Command | Expected | Observed | Cause |
|---|---|---|---|---|
| BC-2 | `pileup` | a hard clip is not in SEQ, so a `5H25M` read that retains the C at chrH 11 gives `CV=1 BT=1.000` and no variant rows, the same as the `5S25M` control | `BT=0.000` at 11, plus spurious `T>C` and `T>G` rows at 6 and 7 and an `AB=Y` row at 12: every base is read 5 positions further along SEQ | the CIGAR walker advances the SEQ index on `H` (`pileup.c:872-874`; also `cnt_retention` and `infer_bsstrand`, `bisc_utils.c:112-114, 196-198`) |
| BC-3 | `epiread` | one epiBED row for the same read | `Unknown cigar 5`, SIGABRT (exit 134), no output | the CIGAR walker has no `H` case and falls to `default: abort()` (`epiread.c:706-1002`) |
| BC-8 | `pileup -t 5` | help: "Maximum cytosine retention in a read"; docs: "retained C's for OT/CTOT reads". A fully converted read (0 retained) is kept, a fully retained read (10) is dropped | the converted read is dropped and the retained read is kept | `cnt_retention` has the strand test inverted: for a C-strand read it counts reference G read as G, for a G-strand read reference C read as C (`bisc_utils.c:94-98`; compare the correct test in `bsconv.c:88-94`) |
| BC-10 | `pileup`, `epiread` | help: "-5 INT Minimum distance to 5' end of a read", "-3 INT … 3' end". On a reverse-strand read, whose 5' end is on the right, `-5 10` drops the call 5 bases from that end (G 56) and keeps G 6 | `-5 10` keeps G 56 and drops G 6; `-3 10` drops G 56. `epiread` filters the same end (`F10x45Mx4`) | both filters count along SEQ from the left, 1-based and including soft clips, for either strand (`pileup.c:432`, `epiread.c:739`) |
| BC-11 | `pileup` | help: "-d Double count cytosines in overlapping mate reads (avoided by default)", so each fragment counts once at every covered C. Pair p (R1 30M at 101, R2 80M at 111): `CV=1` at 150 and 160, from R2. Pair q (R1 80M at 101, R2 30M at 151): `CV=1` at 160 | with no `MC` tag, pair p has no rows at all, and pair q has `CV=2` at 160. Adding `MC` to pair p restores both rows | without `MC`, the mate's reference length is assumed to equal the record's own, and R2 drops every base in `[max(pos, mpos), min(end, mpos + own length − 1)]` (`pileup.c:784-796, 822-827`; the same code in `epiread.c:681-697, 756-762`) |
| BC-13 | `pileup -p` | `-p` keeps improperly paired reads; a mate on another contig cannot overlap, so `CV=1` at chrI 120 | no row when the mate is at chrJ 101; the control with the mate at chrI 1001 gives the row | the overlap test compares `pos` and `mpos` without checking that `mtid == tid` (`pileup.c:781, 822-827`) |
| BC-18 | `pileup`, `pileup -g` | a row at every covered C, including the contig's last base (chrE 30); region `chrE:1-20` includes 20 | no row at chrE 30; `-g chrE:1-20` gives no row at 20, `-g chrE:1-21` does | windows are half-open with `end` set to the contig length or the region end (`pileup.c:808, 1224-1246`) |
| BP-5 | `pileup` with two BAMs | usage: `<in1.bam> [in2.bam …]`, one sample column per BAM; the second BAM's read at chrI 120 gives its sample `CV=1` there whatever the argument order | BAMs listing the contigs in opposite orders: the read is reported on **chrJ**, as `A>C` and `A>G` variant rows at 120-121 with no methylation call. With the BAMs swapped, the call is right | the first BAM's header names the contigs, and every BAM is queried by that header's contig index (`pileup.c:716-719, 753`; a code comment calls the equal-header assumption "probably okay") |
| BP-4 | `vcf2bed -c`, `mergecg -c` | help: "Output Beta-M-U". 1 of 2001 reads retains the C at chrS 21: `M=1 U=2000` | both write `M=0 U=2001`. The VCF's own `SP` field still says `C1Y2000` | the counts are rebuilt as `round(beta × CV)` from `BT` printed to 3 decimals (`pileup.c:666`, `vcf2bed.c:173`, `mergecg.c:116-117`); from CV ≈ 1000 the rounding error can reach a whole count |
| BP-8 | `mergecg` | help: the input is `vcf2bed` output; the `vcf2bed -t c` row for the C at chrQ 1 is merged or written | `Error retrieving base 0 outside range chrQ:1-40`, exit 1, no output | the base before the row is fetched at 1-based position 0; the guard `end - 1 < 0` is written for 0-based positions (`mergecg.c:192-196`) |
| BE-2 | `epiread -B` | docs: `vcf2bed -t snp` "columns 6-9 are repeated for each sample"; epiread help: `-B` takes that BED. A C with a T allele at frequency ≥ 0.05 is not callable (`epibed_format.md`), so the SNP is shown and the C is not called | with the two-sample BED from `vcf2bed -s ALL -t snp` (13 columns), the output equals the output with no BED: `M` at 21, no SNP. Cutting the BED to 9 columns gives the veto (`x`) and the SNP string | lines are parsed only if they have exactly 8 tabs; other lines are skipped without a message (`epiread.c:1091`) |
| BE-3 | `epiread -B` vs `pileup` | epiread help: with an unfiltered BISCUIT SNP BED, "should have the same results in the epiBED as in pileup". 21 C-strand reads retain the C at chrS 21, 1 G-strand read shows T | `pileup` calls it (`CV=21 BT=1.000`); `epiread` writes `x` at 21 for every C-strand read | `pileup` tests unrounded `T/C < 0.05` (1/21 = 0.048); `epiread` tests `AF1 < 0.05`, where `AF1 = T/(C+T)` = 1/22 printed as `0.05` (`pileup.c:524, 654`; `epiread.c:1128-1137`) |

## What this means for comparisons

- **BC-2** is silent and corrupts every call on a hard-clipped record: methylation calls,
  base qualities and SNP evidence all come from the wrong read position, and the last H
  positions are read past the end of SEQ. `pileup` never filters supplementary records,
  and bwa mem and BISCUIT's own `align` hard-clip them unless run with `-Y`
  (whether bwa-meth's `-M` output is affected is not checked). Paired supplementaries from `biscuit align` have no mate information and are
  removed by the proper-pair filter, so exposure is highest for single-end data, for
  `-p`, and for aligners that set the proper-pair flag on supplementaries. alnbase reads the
  hard-clipped read correctly: the retained C at chrH 11 (0-based 10) with `off_5p` 15,
  counting the 5 hard-clipped bases from the read's 5' end.
- **BC-3** makes `epiread` unusable on any BAM with a hard-clipped record that passes its
  filters. It fails loudly, so a comparison would notice.
- **BC-8** matters only when `-t` is set; the default (999999) turns the filter off. When it
  is set, it counts matched Gs on C-strand reads, so it drops reads by G content instead of
  by incomplete conversion: a typical 150 bp read carries dozens of Gs and fails any small
  threshold, while an unconverted read over a G-poor region passes. alnbase leaves per-read
  filters downstream; the per-read retained and converted C counts it prints are a
  `group by` over its hits.
- **BC-10** has no effect at the default `-5 3 -3 3`, which is symmetric. It matters when
  the trimming is asymmetric, the usual response to an M-bias plot: for every
  reverse-strand record the trim lands on the wrong end, so the biased bases stay in and
  unbiased bases are dropped. Soft clips also count towards the distance, as a code TODO
  notes. alnbase reports `off_5p` and `off_3p` from the read's 5' and 3' ends as sequenced,
  so the intended filter is `off_5p >= 10` on either strand.
- **BC-11** affects every Bismark paired-end BAM, because Bismark writes no `MC` tag while
  `biscuit pileup` reads Bismark's `XG` strand tag. The overlap is right only when the two
  mates have the same aligned reference length. After adapter and quality trimming, which is
  standard before Bismark, they often differ: a longer R2 loses non-overlapping bases, a
  shorter R2 is double-counted in the overlap. BAMs from aligners that write `MC` (bwa mem,
  `biscuit align`) are unaffected. alnbase resolves overlap in its separate `overlap` command,
  from the reads' sequences rather than mate coordinates; it is shown separately, not here.
- **BC-13** needs `-p` and a mate on another contig whose position happens to fall inside
  the record's own interval, so it is rare.
- **BP-5** is silent. It needs BAMs whose headers list contigs in different orders, as when
  samples come from different pipelines or references. Calls move to the wrong contig and
  show up as variant rows. BAMs with identical headers are unaffected.
- **BP-4** changes counts only at CV ≥ about 1000, by up to CV/2000 per site. That depth is rare
  in single-cell or typical WGBS data but routine in deep amplicon or targeted data, where
  low-level methylation (1 in 2001 here) is rounded to zero. `vcf2bed -c` is the documented
  way to get Bismark `.cov`-like counts (`methylextraction.md`), so a comparison of counts at
  deep sites should expect these differences.
- **BP-8** needs a `vcf2bed -t c` row at a contig's first base. `-t cg` cannot produce one,
  because a 5-mer at position 1 contains padding N and is labelled `CN`. It fails loudly.
- **BE-2** needs `vcf2bed -s ALL` (or several named samples); the default `-s FIRST` writes 9
  columns and works. A cohort-wide SNP BED is a natural choice for epiread across samples,
  and with it the SNP veto and SNP strings silently disappear.
- **BE-3** affects only sites where the T allele is near 5%, for example T=1 against C=20-21.
  There the two tools disagree on whether a methylation call is made, contrary to the
  documented match.
- **BC-18** loses one position per contig in whole-BAM mode, where contig ends are usually N,
  so the effect there is small. With `-g`, the last position of every region is lost, so
  a run sharded into adjacent regions (`chr1:1-1000000`, `chr1:1000001-2000000`, …)
  silently loses one position per shard. alnbase reports the Cs at 0-based 19 and 29.
