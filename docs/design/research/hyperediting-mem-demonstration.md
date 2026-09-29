# Top demonstration candidate: hyper-editing detection with `bwa mem` + alnbase

Status: **planned, not started.** This is the leading candidate for the "meaningful
analysis alnbase makes easy that existing tools make hard" demonstration. Findings in the
pipeline itself are catalogued separately in `non-methylation-tool-findings.md` (HE-1..10);
this file is the experimental design only.

Target: the Porath/Levanon hyper-editing pipeline (`hagitpt/Hyper-editing`, `5eee279`),
the method behind the published human hyper-editing atlases.

## The claim

The pipeline aligns its transformed reads with `bwa aln -n 2 -o 0 -N` and back-maps them
with a flat `substr` that assumes read offset equals reference offset (HE-1). Switching to
`bwa mem` would recover reads the current aligner cannot place — but **cannot be done
without a correct CIGAR walk**, because soft clips, insertions and deletions immediately
break that back-mapping and shift every emitted coordinate with no diagnostic.

So the claim is not that alnbase makes an existing analysis tidier. It is that alnbase
**unlocks an analysis the pipeline currently cannot perform**, because the aligner upgrade
is gated on read-to-reference mapping that the pipeline does not have and alnbase supplies
off the shelf. That is a materially easier claim to defend than a convenience argument: the
reviewer does not have to accept that the analysis is *nicer*, only that it is *blocked*.

Note the framing constraint recorded in HE-1: we do **not** argue that the authors chose
`-o 0` in order to avoid CIGAR complexity. That is an inference about their reasoning and
is not supportable from the code. We test the consequence of the choice, not the motive.

## What `bwa mem` should recover

`aln -n 2 -o 0` is end-to-end with a two-difference budget and no gaps. Reads it cannot
place, in rough order of expected volume:

1. **Splice-junction reads.** This is RNA-seq aligned against a genome. An end-to-end
   aligner cannot place a junction-spanning read at all; `bwa mem` soft-clips and aligns
   the longer arm. A meaningful share of Alu editing sits in introns and UTRs, so this is
   not a corner case.
2. **Reads with more than two residual mismatches** — SNPs plus sequencing error plus
   non-canonical substitutions together exceeding the budget. `mem`'s affine scoring is far
   more permissive.
3. **Reads carrying a genuine indel**, which `-o 0` makes impossible by construction.

All three are exactly the reads a hyper-editing study most wants, because the pipeline's
input is *already* the set of reads that failed normal alignment.

## Measurements

Primary: the number of reads emitted as hyper-edited clusters, and the number of distinct
edited sites, for the shipped pipeline versus the `mem` + alnbase arm, on identical input.

### The pipeline ships its own negative control

This is what makes the result measurable without an external truth set. The 12-combo design
includes transform families that should carry little or no real editing signal — `G2C` and
`A2T` in particular. If `mem` + alnbase raises `A2G` yield substantially more than it raises
`G2C`/`A2T` yield, the gain is signal rather than noise. If all families rise together, it
is noise, and the demonstration has failed honestly and cheaply.

### Corroboration on real data

- **REDIportal overlap.** Newly recovered sites should hit the known-editing database at a
  rate comparable to the baseline's, not at background.
- **A-to-I motif signature.** The field's own signature is G depleted immediately upstream
  and enriched immediately downstream of the edited A. Newly recovered sites should carry
  it. The pipeline already computes these statistics (`detect_ue.pl:550-585`), so this comes
  free — modulo HE-11's strand-normalisation inconsistency, which must be fixed first or the
  pooled triplets mix orientations.
- **Junction reads specifically.** Sites recovered only from soft-clipped junction-spanning
  reads should fall near annotated splice sites. If they do, that is a clean, visual,
  single-figure result.

## The risk that could invert this

The alignment is against a **three-letter** reference. That alphabet is degenerate, and
global alignment with a tight mismatch budget is plausibly specificity engineering against
exactly that degeneracy. Local alignment with soft clipping on a 3-letter genome may
manufacture spurious hits faster than it recovers real ones — and a spurious local alignment
whose clipped portion is discarded is precisely how a false editing cluster would appear.

If that is what happens, the honest result is an argument *for* the authors' choice, and we
report it that way. The negative-control families tell us which way it goes early, before
any real data or whole-genome indices are paid for. **This question is settled first.**

Related: HE-2 means the pipeline's multi-mapping safeguard is already broken in the
false-positive direction, so the baseline's own specificity is not what it appears. Any
comparison must either fix HE-2 in both arms or report both arms with it present. Fixing it
in both arms is the fairer test and is what this plan assumes.

## Staging

**Stage 1 — simulated, one chromosome, ground truth (do this first).**
Simulate RNA-seq reads from one chromosome with known editing sites at known rates,
including junction-spanning reads, indel-carrying reads, and reads at the mismatch-budget
boundary. Run both arms. This gives recovery *and* false-positive rates against truth, and
settles the specificity question above. Everything after this is contingent on Stage 1
showing a real gain. Roughly a day's work.

**Stage 2 — real data, one chromosome.** A public RNA-seq dataset with heavy Alu editing
(a brain or a hyper-edited transcriptome sample). Six transformed indices of one chromosome
rather than the whole genome. Adds REDIportal and motif corroboration.

**Stage 3 — full reproduction**, only if Stages 1 and 2 hold. Six whole-genome transformed
indices; substantial disk and index time.

## What it needs

- The Perl scripts do not start on this machine: `analyse_mm.pl`, `detect_ue.pl`,
  `Get_orig_read.pl` and `sort_R_read.pl` all begin `#!/private/apps/bin/perl516 -w`, an
  author-machine absolute path, and are invoked directly. Patch the shebangs; record the
  patch as part of the reproduction.
- `filter_sam.pl`'s relative `unique_simple_repeats.txt` path (HE-4) must be resolved, or
  the repeat filter silently rejects every read or none depending on the Perl version.
- An alnbase query replacing `analyse_mm.pl`'s back-mapping, emitting the same hit-line
  format so `sort_R_read.pl` and `detect_ue.pl` run unchanged. Keeping the rest of the
  pipeline intact is the point: the swap must be surgical, or the comparison stops being
  about the back-mapping.
- A simulator for edited RNA-seq reads. Prefer an existing one per the validation
  conventions; check what `ecosystem-simulators-tooling.md` already surveyed before writing
  anything.
