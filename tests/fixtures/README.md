# Hand-written fixtures

The validation suite (alnbase-validation) rests on two legs, and they share a blind spot.
Simulated truth (alnbase-validation `sims/`) says what the reads were made from, but a simulator only produces
the cases it models, and BSReadSim models none of the ones below. Mirroring says alnbase
and Bismark or MethylDackel agree call for call, which is worth a great deal — and says
nothing at all if the two are wrong the same way.

These fixtures are the third leg. Each is a few hand-typed SAM records over a reference of
tens of bases, chosen so that the right answer follows from the SAM specification and the
rule in [`context.toml`](context.toml) rather than from another tool's output. Every
`expected.tsv` here was written before alnbase was run on it; all fourteen matched on the
first run.

```
PY=/path/to/python-with-duckdb ALNBASE=/path/to/alnbase ./run.sh [case ...]
```

`run.sh` builds each case's BAM with samtools, runs one `alnbase query`, and diffs the
result. A case with an `overlap.args` file is queried twice, once before and once after
`alnbase overlap` runs on it. Exit status 0 means every case matched. Three of the steps
are checks in their own right: `samtools view` refuses a record whose CIGAR and SEQ lengths
disagree, `alnbase index` refuses a reference it cannot digest, and the scan line alnbase
prints is compared against `expected_scan.txt` where a case has one.

## The cases

| case | the hazard | what it pins |
|---|---|---|
| [`read_end_cg`](read_end_cg/) | a CG whose partner base is past the end of the alignment, on each strand | the context comes from the reference through pad columns, so the call survives |
| [`contig_end_cg`](contig_end_cg/) | a cytosine at a contig's first and last base | context off the contig is undetermined, not CHH and not dropped |
| [`n_context`](n_context/) | `CNG` and `CNN` beside `CAG`, `CAA` and `CG` | IUPAC-strict: N is not H, so the context is undetermined |
| [`indel_reference_rule`](indel_reference_rule/) | a deletion and an insertion inside a context window | the reference rule gives the same answer all three ways; Bismark's does not |
| [`flagged_records`](flagged_records/) | duplicate, secondary, supplementary, QC-fail, unmapped | the first four are walked; the unmapped record is skipped and counted |
| [`low_quality_tail`](low_quality_tail/) | the same CG at Phred 40 and Phred 2 | both called, each carrying its own quality; no built-in threshold |
| [`unconverted_read`](unconverted_read/) | a read that retains every cytosine | six calls, not a discarded read; the filter is a `group by` |
| [`hard_clip_offsets`](hard_clip_offsets/) | a read split in two, each half hard-clipping the other, on each strand | `off_5p` is the sequencing cycle, so the second half's C is cycle 29 and not SEQ index 9 |
| [`soft_clip_context`](soft_clip_context/) | a soft-clipped base standing where a CG's G is | the clip emits a flank column carrying the real reference base, so the call lands; offsets count the clip |
| [`cigar_edge_deletions`](cigar_edge_deletions/) | `10M2D` and `2D10M` | a first or last operation that consumes no read base changes no call |
| [`iupac_reference`](iupac_reference/) | `Y`, `R`, `S`, `W` in the reference beside `N` | subset matching: `CYG` is a CHG, `CRG` is undetermined |
| [`overlap_agreeing_mates`](overlap_agreeing_mates/) | two mates observing the same CG, with no `MC` tag | one molecule base is one row, at the summed quality; the table before resolution is checked too |
| [`overlap_mate_disagreement`](overlap_mate_disagreement/) | mates that read the same base as `C` and as `T` | the end that keeps the overlap and the base that wins are separate decisions; the survivor can be the loser, at quality 0 |
| [`overlap_cross_contig`](overlap_cross_contig/) | a ligation product whose mates' primaries sit on different contigs | the overlap lives on R1's supplementary, which is found without `MC` or `SA` and is the record clipped away |

Seven of these pin a deliberate choice rather than a correctness claim — which records to
walk, which IUPAC codes resolve a context, which of the two indel rules is in force,
whether a quality threshold exists, what happens to a base two mates disagree about, and
which end of a fragment gives up its copy.
That is the reason to write them down. A choice nobody recorded is indistinguishable from
an accident the next time a comparison disagrees, and most of them are exactly where
alnbase and another caller will differ on real data.

## Anatomy of a case

```
n_context/
  ref.fa         one contig of tens of bases, all of it visible at once
  reads.sam      SAM text with explicit FLAG, POS and CIGAR
  expected.tsv   one row per call; its header names the columns the case is about
  README.md      the hazard, and which tool behaviour it probes

overlap_agreeing_mates/          a case that runs `alnbase overlap` first
  overlap.args       extra options for it, possibly none, and the marker that it runs
  expected_raw.tsv   the calls before resolution -- the double count itself
  expected.tsv       the calls after it
```

`expected.tsv`'s header is the case's column list: `run.sh` selects exactly those columns
from alnbase's output, so a case about read coordinates asks for `off_5p` and `off_3p`
while one about contexts does not carry columns it has nothing to say about. Everything
committed is text, so a reviewer sees the whole case on one screen and can check the
expected table against the reference by eye. `refr_pos` is 0-based, as alnbase
reports it, while the READMEs describe the references in the 1-based coordinates SAM and
the FASTA use; the two differ by one on purpose, because that is the conversion the
fixtures exist to get right.

## Adding a case

Write `ref.fa` and `reads.sam`, derive `expected.tsv` **by hand from the reference and
`context.toml`**, and only then run `run.sh`. Deriving it from alnbase's output first
would make the fixture a regression test for whatever alnbase does today, which is not
what this directory is for. Add `expected_scan.txt` when the case is about records that are
skipped rather than called, since a skip leaves no row.

If a case needs a different call rule — Bismark's read-projected context, say — it belongs
in a separate directory with its own copy of `ref.fa` and `reads.sam` rather than a second
expected table beside the first, so that each directory still answers one question.
