# The hyper-editing pipeline with bwa mem, and what it takes to get there

The A-to-I hyper-editing pipeline of Porath, Carmi & Levanon (*Nat Commun* 2014,
`hagitpt/Hyper-editing` at GitHub `5eee279`) finds RNA editing in reads that no
aligner will map, by collapsing reads and genome to a three-letter alphabet,
realigning, and reading the edits back off. Its edit caller decides which
reference base each read base sits on by cutting a window out of the genome at
the alignment's start position and pairing position *i* with position *i*
(`analyse_mm.pl:230`). That is correct only for a CIGAR which is one run of `M`,
and the pipeline guarantees it by passing `-o 0` to `bwa aln` in twelve places
in `TransRun.sh` — three stages earlier, with nothing in between recording the
dependency.

The guarantee costs every read that is not perfectly flat on the genome:
junction-spanning, indel-bearing, clipped. This demo replaces the pairing with
alnbase's CIGAR walk, switches the aligner to `bwa mem` and then to a
splice-aware one, and measures what comes back against ground truth.

## Result

Simulation: 6,000 single-end 100 bp reads over a 300 kb contig with a four-exon
gene and an inverted pair of A-rich elements, 53,036 edited bases with known
genomic coordinates. Six arms: three aligners crossed with the published
flat-window caller and the alnbase one, sharing alignment stages where they can.

| arm | site recall | placement | misplaced | read recall |
|---|---|---|---|---|
| `aln-published` — the pipeline as published | 40.6% | 99.9% | 24 | 39.3% |
| `mem-naive` — swap the aligner, change nothing else | 44.7% | 99.1% | 227 | 43.2% |
| `mem-flat` — `bwa mem`, flat-window caller | 44.7% | 99.6% | 92 | 43.2% |
| `hisat2-flat` — splice-aware, flat-window caller | 42.8% | 99.7% | 58 | 41.5% |
| `mem-alnbase` — `bwa mem`, alnbase caller | 85.6% | 99.4% | 264 | 88.1% |
| **`hisat2-alnbase`** — splice-aware, alnbase caller | **88.8%** | 99.4% | 261 | 87.9% |

**2.19× more editing events recovered, with placement accuracy unchanged.**

Where it comes from, as recall by read class:

| arm | contiguous | junction | deletion | insertion | clipped |
|---|---|---|---|---|---|
| `aln-published` | 99.6% | **0.0%** | **0.0%** | **0.0%** | 0.5% |
| `mem-naive` | 99.7% | 0.0% | 0.7% | 0.0% | 28.0% |
| `mem-flat` | 99.7% | 0.0% | 0.6% | 0.0% | 28.0% |
| `hisat2-flat` | 99.2% | 0.0% | 0.6% | 0.0% | 16.2% |
| `mem-alnbase` | 98.3% | 30.5% | 96.0% | 97.5% | **97.9%** |
| `hisat2-alnbase` | 99.7% | **99.7%** | 90.9% | 88.2% | 43.3% |

The published pipeline does not merely lose these reads at some rate. It finds
**none** of the junction-spanning, deletion- and insertion-bearing reads, and
one in two hundred of the clipped ones. The missing reads were never a random
sample, which is why the aggregate moves so far.

### A better aligner, on its own, buys nothing

This is the clearest result in the table, and it runs the opposite way to the
obvious intuition. Read down the `flat` rows: `bwa mem` gains four points over
`bwa aln`, a splice-aware aligner gains two, and **both still find zero
junction-spanning reads.** HISAT2 places those reads correctly, with a proper
`N` in the CIGAR — and then the flat-window caller cuts a 100-base window at the
alignment's start position, compares a read that spans 29 kb of genome against
100 bases of the first exon, and throws the result away as too mismatched to be
editing. A splice-aware aligner is *worse* than `bwa mem` under the old caller.

`mem-naive` — the change someone would actually make after reading the paper —
gains four points and multiplies misplaced sites by nine, from 24 to 227: calls
that look ordinary and are on the wrong base.

**The aligner is not the improvement. The caller is. The aligner is what makes
the caller worth having.**

### The two aligners recover different things

`bwa mem` gets 97.9% of clipped reads and 30.5% of junction reads; HISAT2 gets
99.7% of junction reads and 43.3% of clipped ones. HISAT2 prefers to find a
splice or extend an alignment where mem prefers to clip, so each aligner is
strong exactly where the other is weak. Running both and taking the union would
be close to complete on every class; this demo does not do that, because the
point here is the caller and a union arm would confound it.

### Recall by how heavily the read was edited

| arm | 5-8 edits | 9-12 | 13-16 | 17+ |
|---|---|---|---|---|
| `aln-published` | 29.8% | 39.0% | 44.9% | 42.9% |
| `mem-alnbase` | 77.3% | 85.1% | 91.0% | 77.4% |
| `hisat2-alnbase` | 81.3% | 87.8% | 90.0% | **96.2%** |

The most heavily edited reads — the ones the method exists to find, and the ones
an ordinary editing screen is guaranteed to miss — are where the splice-aware
arm is strongest.

## What the port actually touches

The walkthrough that preceded this work said the swap was contained to one
stage. Running it showed that is not quite true, in two ways worth recording.

**The flat assumption is made twice.** `detect_ue.pl:538` turns an edit's offset
along the read into a genomic coordinate by adding it to the alignment's start
position — the same assumption as the caller, made again downstream. Fixing only
the caller leaves every recovered read placed at the wrong base, which is worse
than not recovering it: the first attempt scored 65% placement against the
original's 99.8%. So the hit line gains a CIGAR field and `detect_ue.pl` gains a
`ref_offset_of` subroutine. The conversion is the identity on a CIGAR of one
`M` run, so the `aln` arm is unaffected — checked.

**Fixing the walk reaches a latent crash.** `detect_ue.pl:218` divides by
`mm_count` without checking. Under the published pipeline `mm_count` is never
zero, because `bwa aln` aligns end to end and a restored read always carries its
own edits as differences. Once the caller stops comparing clipped bases against
whatever the genome happens to hold there, a read whose informative part was
clipped away has a clean core and nothing to count, and the script dies. This is
HE-10 in the findings catalogue, previously recorded as latent. The guard is
applied to every arm so the downstream code is identical across them.

One more, which is about `bwa mem` rather than about alnbase: `analyse_mm.pl`
counts a read's alignment positions from the `X0` and `X1` tags, which only
`bwa aln` writes. Under `mem` the count is zero, and `sort_R_read.pl` — which
reads exactly that many hit lines after each record — falls out of step with
the file and mis-parses everything after the first read. Every `mem` arm needs
that line whatever its caller.

### An interaction worth knowing about

`detect_ue.pl:228` chooses between a read's alternative alignments by ranking
them on `edit_count / mm_count` and requiring the winner to lead by 0.1. That
rule is partly powered by the flat window's invented mismatches: they differ
between the true locus and a decoy, which separates the ratios. With the
corrected caller, alignments of the same read at two loci both score near 1.0
and the rule refuses to choose.

This showed up as a large effect in a first simulation whose inverted repeat was
an exact reverse complement, making every read from it a perfect two-locus
mapper. Real Alu copies in a hairpin are roughly 12% divergent, and the
simulator now models that; with realistic divergence the rule discriminates
normally and the effect is small. It is recorded because it is real, and because
it means the rule's threshold was calibrated against a caller that no longer
behaves the same way.

## Layout

| file | what it does |
|---|---|
| `setup.sh` | fetch `hagitpt/Hyper-editing` at `5eee279`; portability patches only |
| `simulate.py` | the simulation and its ground truth |
| `run.sh` | run all four arms and score them |
| `score.py` | recall, placement, and three bias strata |
| `bin/TransRunAligner.sh` | the twelve transform alignments, aligner as a parameter |
| `prep_real.sh` | whole-genome indexes for the real-data run |
| `prep_chrom.sh` | the same for one chromosome, in minutes rather than a day |
| `run_real.sh` | the real-data run: whole-genome stage 1, one chromosome for stage 2 |
| `score_real.py` | the arms against the pipeline's own controls, with no truth set |
| `expected_real_chr21.tsv` | the chr21 real-data table, as committed |
| `bin/alnbase_mm.py` | the replacement edit caller |
| `bin/patch_analyse_mm.py` | CIGAR field on the hit line; the X0/X1 line |
| `bin/patch_detect_ue.py` | `ref_offset_of`, the three coordinate sites, the crash guard |
| `bin/probe_splice_motif.py` | how much losing GT..AG costs, against a counterfactual reference |
| `bin/probe_motif_protection.py` | whether exempting the motif from the collapse can work |
| `bin/probe_novel_pairing.py` | novel pairings of annotated sites, and what a pair list does to them |
| `queries/hyperedit.toml` | the twelve ordered mismatch pairs, as alnbase queries |

## Running it

```sh
./setup.sh                $WORK
python simulate.py        $WORK/sim
ALNBASE=/path/to/alnbase ./run.sh $WORK $WORK/sim
```

Needs `bwa` (0.7.18 here), `hisat2` (2.2.1), `samtools`, `perl`, `xa2multi.pl`
from bwa's own distribution, and a python with `pysam` and `pyarrow`.

### On a real library, one chromosome at a time

```sh
./prep_chrom.sh $WORK/real chr21 $WORK/real/gencode.v46.annotation.gtf.gz
./run_real.sh   $WORK/real chr21 $WORK/real/reads.fastq
```

`prep_real.sh` builds thirteen whole-genome indexes. `hisat2-build` on hg38 is
hours per transform and there are six, so the whole set is most of a day and
about 60 Gb. `prep_chrom.sh` does the same for chr21 in **four minutes and
1.3 Gb**, and every arm still sees exactly the same reference, which is what a
comparison between arms needs.

The two references are used for different things, and that split is the design
rather than a shortcut:

| stage | reference | why |
|---|---|---|
| 1 — what ordinary alignment cannot place | whole genome | run against one chromosome it would call most of the library unmappable and hand stage 2 a pool made of other chromosomes' reads |
| 2 — the twelve transform alignments | one chromosome | the expensive part, and the only part the subset touches |

Reads that survive stage 1 come from everywhere, so some of them will find a
spurious home on the subset. That is measured rather than waved away. Of the six
mismatch classes `detect_ue.pl` reports, only A2G can carry an A-to-I event;
A2C, A2T, C2A, G2A and G2C are noise channels by construction, and `score_real.py`
reports every A2G number against them. It also reads the ADAR signature straight
off the trinucleotide context the pipeline already writes — depletion of G 5′ of
the edited base, enrichment 3′ — which is a specificity check that needs no
external truth set and that misalignment does not reproduce.

The library is Illumina BodyMap 2.0 brain, run `ERR030890` — the dataset Porath
2014 used. 64.3 M single-end 75 bp reads; every third read is taken, so the
subsample spreads across the flowcell rather than starting at the first tiles.

### Splice-aware alignment needs an annotation here — but not for the obvious reason

None of the six three-letter transforms leaves a canonical splice motif intact:
`a2g` turns `AG` into `GG`, `t2c` turns `GT` into `GC`, `g2c` destroys both. A
transformed genome has no `GT..AG` anywhere.

The tempting conclusion — that junctions therefore cannot be found de novo — is
**wrong**, and measuring it says something more useful.
`bin/probe_splice_motif.py` sweeps anchor length, comparing the collapsed genome
against a counterfactual one with `GT..AG` written back in at the true intron
boundaries: four bases per intron, so alphabet, length and repeat content are
identical and the difference is the motif alone.

| anchor (bases on the short side) | collapsed, de novo | motif restored, de novo | junctions supplied |
|---|---|---|---|
| 4–11 | 0% | 0% | 84–100% |
| 12 | 0% | 68.3% | 100% |
| 13 | 15.7% | 82.2% | 100% |
| 14 | 69.5% | 84.2% | 100% |
| 15–18 | 81–86% | 81–86% | 100% |
| 20+ | 100% | 100% | 100% |

The motif effect is real and **three bases wide**. Below a 12-base anchor no
junction is found either way; above 15 the motif contributes nothing, because a
long anchor across a 29 kb intron is worth far more to the aligner than the
non-canonical penalty costs it. Losing `GT..AG` shifts the de novo floor from
about 12 bases to about 15.

Two things matter more than that:

**Supplying the junctions beats the motif everywhere.** It recovers reads
overlapping by as little as six bases, which no amount of motif ever would, and
about **28% of junction-spanning reads have an anchor under 15** — for a 100 bp
read the short side is roughly uniform on 1..50, so this is not an edge case.

**De novo splicing on a three-letter genome invents junctions.** Across the same
9,000 reads it placed **499 gaps at a wrong intron; with the junctions supplied,
zero.** For this pipeline that is the serious one. A misplaced junction does not
produce an obvious gap in the output — it produces a read aligned against the
wrong sequence, and therefore a burst of mismatches indistinguishable from dense
editing, which is the only signal the method has.

#### Exempting the motif from the collapse

The obvious alternative — transform the genome but leave the donor and acceptor
dinucleotides alone — comes in two forms that behave completely differently.
`bin/probe_motif_protection.py` measures the second.

**Targeted**, at annotated intron boundaries only: works, and costs about one
base per intron. It is the `motif-restored` reference in the table above, so it
buys back the 12–14 anchor band and nothing more. It also needs the annotation,
and anyone holding the annotation can hand the aligner the junction list
instead, which reaches 100% from six-base anchors. Dominated.

**Global**, at every occurrence in the genome: needs no annotation, which is the
whole appeal, and it destroys the method.

| reference | alignment rate | edited reads placed | junction gaps |
|---|---|---|---|
| plain a2g | 45.1% | 2,038 | 414 |
| every genomic `AG` protected | **4.0%** | **167** | 84 |

Leaving 6.1% of positions unconverted costs 92% of the edited reads. The reason
generalises past this pipeline: **the three-letter trick works because the
transform is context-free.** Each base maps to its image regardless of its
neighbours, and that is exactly what guarantees an edited A and an unedited A
collapse to the same letter. "Leave an A alone when the next base is G" is
context-dependent, and the context can itself be edited — a genomic `AG` whose
adenosine was deaminated reads `GG`, the genome keeps its `A`, and the position
becomes a hard mismatch at an edited base. The rule fires hardest against
exactly the reads the method exists to find. Under a2g only the acceptor needs
protecting, since the donor `GT` has no adenosine, so the measurement above is
the most favourable case for the idea.

#### Novel pairings of annotated sites, where the pair list turns harmful

A splice-site file is a list of donor-acceptor *pairs*. Exon skipping and
alternative site usage join two sites that are both annotated in a combination
that is not, so those junctions are absent from the list even though nothing
about them is novel except the pairing. Protecting the motif marks the sites
*individually*, so it should help exactly there.
`bin/probe_novel_pairing.py` measures it, on reads crossing a skipping junction
that joins the donor of one annotated intron to the acceptor of another.

| setup | anchor 10 | 12 | 14 | 16 | 20 | 30 | gaps at a wrong intron |
|---|---|---|---|---|---|---|---|
| collapsed, no list | 0% | 0% | 75% | 100% | 100% | 100% | 0 |
| collapsed, canonical list | 0% | 0% | 52% | 76% | 74% | 100% | **483** |
| motif protected, no list | 0% | **74%** | **100%** | 100% | 100% | 100% | 0 |
| motif protected, canonical list | 0% | 74% | 78% | 76% | 74% | 100% | 385 |
| skipping pair in the list | 100% | 100% | 100% | 100% | 100% | 100% | 0 |

Two results, and the second is the one that matters.

**Motif protection works, as intended and only where intended.** It moves the de
novo floor for the novel pairing from a 14-base anchor to 12, and takes anchor 14
from 75% to 100%. It cost three bases in this genome — under a2g only the
acceptor needs saving, because the donor `GT` has no adenosine — and all of them
are intronic, so a spliced read never contains one. On hg38 that is about 200 kb,
0.006%, against the 6.1% that made global protection collapse.

**An incomplete pair list is worse than no pair list, for reads crossing the
junction it is missing.** Recovery of the skipping junction falls from 100% to
about 75%, and 483 of 3,200 reads come back with a *confidently misplaced*
junction. The mechanism is verified rather than guessed: reads whose true intron
is 49,100 bases were given 19,100 — the annotated junction that shares their
donor. The aligner snaps the read onto the known pair, lands its second half on
the wrong exon, and reports a clean spliced alignment.

For this pipeline that is the dangerous failure. A read placed on the wrong exon
produces a burst of mismatches against sequence it never came from, which is
indistinguishable from the dense cluster the method calls hyper-editing. And it
does not wash out with longer anchors: at a 30-base anchor the misplacement rate
is still around 25%.

Motif protection does not rescue it — protected plus canonical list is 78/76/74%,
so the list's harm dominates. **The conclusion is that the junction list has to
be complete for the sample, which means deriving it from the sample.** A first
pass over the whole library against the untransformed genome yields the
junctions actually used, skipping events included. That was already the better
option on completeness grounds; this makes it the only safe one.

Caveats: one simulated genome with three introns and two skipping events, and
HISAT2-specific behaviour. A real annotation is far denser, so there are many
more shared-donor decoys to snap onto, but also far fewer genuinely missing
pairings. The direction is clear; the magnitude on real data is not.

The reason supplying coordinates is legitimate rather than a fudge: the
transform is a pure per-base substitution, so every coordinate is identical in
all seven genomes.

For real data, the annotation need not come from GENCODE. A first pass aligning
the whole library to the *untransformed* genome — which this stage nearly does
already — yields a sample-specific junction catalogue including novel junctions.
The hyper-edited reads fail that pass by construction, but they do not need to
pass it: the ordinary 95% of the library defines the catalogue.

A caveat this demo cannot settle: hyper-edited reads come mostly from inverted
Alus in introns and 3' UTRs, and intronic Alus sit in unspliced pre-mRNA where
there is no junction to cross. The junction-spanning hyper-edited population may
therefore be dominated by 3' UTR and exonized-Alu reads near annotated
junctions, which would make the annotation restriction cheaper than it looks.
The real-data run should be able to say.

#### Where the junction list should come from

If omission is what hurts, the question is which source omits least. The two
candidates fail in opposite directions: an annotation cannot contain a junction
nobody has catalogued, and the pipeline's own first pass cannot contain a
junction the sample expresses too weakly to be seen. `bin/probe_junction_list_source.py`
builds a gene with one instance of each -- two exon-skipping junctions absent
from the annotation, one annotated intron carried by almost no reads -- and
scores six ways of supplying the list.

```
junction found, by anchor      8    10    12    14    16    20    30   misplaced
no list                       0%    0%    0%   64%   88%  100%  100%          62
annotation only              60%   64%   58%   77%   88%   90%  100%         104
first pass only              83%   80%   80%  100%  100%  100%  100%           0
union                       100%  100%  100%  100%  100%  100%  100%           0
union + motif protection    100%  100%  100%  100%  100%  100%  100%           0
first pass + protection      83%   80%   92%  100%  100%  100%  100%           0
protection, no list           0%    0%   68%   91%   88%  100%  100%          62
oracle                      100%  100%  100%  100%  100%  100%  100%           0

junction found, by junction  intron1  intron2  intron3    skip1    skip2
no list                         43%      48%      57%      58%      47%
annotation only                100%     100%     100%      34%      47%
first pass only                100%     100%      57%     100%     100%
union                          100%     100%     100%     100%     100%
first pass + protection        100%     100%      66%     100%     100%
oracle                         100%     100%     100%     100%     100%
```

Three things come out of it.

The first pass is a far better source than its description suggests. It runs
against the untransformed genome, in four letters, with GT..AG intact, so it
discovers junctions at depths where de novo discovery on a collapsed genome is
hopeless: three spanning reads were enough to recover the lowly expressed
intron, and at two reads it was still found. It invented no junctions. That is
the asymmetry that makes the idea work -- a junction needs only one *ordinary*
read to be catalogued, and that read is not the hyper-edited one that failed to
map and sent the library down this path in the first place.

The union reaches the oracle, at every anchor and every junction, at every
first-pass depth tested. Neither source alone does. Since the cost is entirely
in omission, adding junctions can only help, and the two sources cover each
other exactly.

Motif protection is a weaker belt than claimed above. Against a list missing a
junction it recovers part of the gap and not the whole of it (intron3 57% ->
66%, anchor 12 80% -> 92%), because the missing junction's acceptor is shared
with a listed junction and the snapping still wins most of the contest. With no
list at all it helps below a 14-base anchor and does nothing above. It is worth
the three intronic bases it costs, but it is not a substitute for a complete
list and should not be described as one.

Caveats unchanged: one simulated genome, three introns, two skipping events,
HISAT2-specific snapping. A real annotation is far denser, which means both more
shared-site decoys to snap onto and fewer genuinely missing pairings; the
direction is clear and the magnitude is not.

## On real data

Illumina BodyMap 2.0 brain, run `ERR030890` — the library Porath 2014 used.
Every third read of 64.3 M, so 21.4 M single-end 75 bp reads. Stage 1 against
the whole genome leaves **1,675,400 reads that ordinary alignment could not
place** (7.8%), and that pool is what the three-letter stage sees. Stage 2 runs
against chr21.

### What the port buys

| arm | hyper-edited reads | control-channel reads | editing sites | control sites | signal:noise |
|---|---|---|---|---|---|
| `aln-published` | 467 | 2 | 3,314 | 9 | 368:1 |
| `hisat2-alnbase` | **1,441** | 54 | **9,803** | 301 | 33:1 |

**3.1× more hyper-edited reads**, and 3.0× more editing sites. There is no
ground truth here, so the number needs two independent checks, and both are
computed from what the pipeline already writes.

The **control channels** are the first. `detect_ue.pl` reports six mismatch
classes and only A2G can carry an A-to-I event; A2C, A2T, C2A, G2A and G2C are
noise by construction. They rise with the signal, from 2 reads to 54, which is
what a subset reference does — reads whose true locus is on another chromosome
have nowhere correct to go, and HISAT2 is readier to place them than `bwa aln`
is. Signal to noise falls from 368:1 to 33:1 and stays high.

The **ADAR signature** is the second, and it is the one that decides the
question. ADAR is depleted for G immediately 5′ of the edited adenosine and
enriched for G immediately 3′, against a genomic background near 21%:

| arm | G at −1 | G at +1 |
|---|---|---|
| `aln-published` | 17.0% | 33.0% |
| `hisat2-alnbase` | 18.1% | 32.2% |

If the extra 974 reads were misalignment they would carry no motif preference
and would drag both columns toward 21%. They barely move. The gain is editing.

What this comparison cannot separate is the aligner from the caller: the two
arms differ in both. The simulation says the caller is nearly all of it, and
that is not re-measured here.

### The junction list changes the alignments and not the output

This is the result that goes against the simulation, and it is the more useful
one. Four arms, identical in every respect except the splice-site file:

| list | spliced alignments reaching the caller | hyper-edited reads | editing sites |
|---|---|---|---|
| `aln-published`, for reference | **0** | 467 | 3,314 |
| none | 53,460 | 1,441 | 9,803 |
| first pass | 58,479 | 1,441 | 9,803 |
| GENCODE | 59,773 | 1,441 | 9,803 |
| **union** | **59,963** | 1,441 | 9,803 |

The published arm's zero is the premise of this whole demo, confirmed on real
data: `bwa aln -o 0` is gapless, so not one spliced alignment exists to be
examined.

The list ordering is exactly what `probe_junction_list_source.py` predicted —
**union > GENCODE > first pass > none** — so the list is doing on real data
precisely what it did in simulation. And the called sites are not merely equal
in number. They are the same sites: 9,022 distinct positions, identical across
all four arms, 100% overlap. Six and a half thousand extra spliced alignments
produce **zero** extra editing calls.

The reason is the thing the simulation could not model. Hyper-editing happens in
inverted Alu pairs, and those sit overwhelmingly in introns and 3′ UTRs — in
unspliced pre-mRNA, where there is no junction to cross. Only **0.15–0.21%** of
called sites fall within 10 bp of a junction. The simulation's gene had a sixth
of its reads spanning junctions, which is what made the effect large there.

So the honest summary of the splice work is: the junction list is free, it is
measurably correct, the union is the right way to build it — and on this library
it buys nothing. A previously unmeasured caveat is now measured, and it went
against the thing it was a caveat on.

### What the subset costs

The `bwa mem` arm is missing above, and the reason is worth recording. `mem`
soft-clips, so against a subset reference it places reads from outside the
subset rather than rejecting them; `filter_sam.pl` then grinds through the extra
records for **thirteen hours**, against twenty-five minutes for a HISAT2 arm on
the same pool. The aligner comparison is therefore subset-confounded in a way
the junction-list comparison is not, since those four arms are all HISAT2 and
differ only in the list. `mem` belongs to a whole-genome run.

Caveats: one library, one tissue, one chromosome for stage 2, and the junction
stratum on chr21 is about fifteen sites, which is too few to call on its own.
