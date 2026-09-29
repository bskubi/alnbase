# How aligners specify a read's original strand

Written 2026-09-17, to decide what alnbase's declarative strand handling has to express.
alnbase currently derives the strand from the SAM FLAG alone, under a directional-library
assumption (`src/tags.rs:140`). Reproducing premethyst's bugs turned up an aligner whose
FLAG does not mean what SAM says it means (`validation/demos/premethyst-bugs/README.md`),
so the assumption needs replacing with something the user declares. This note inventories
the cases first; no syntax is proposed here.

Citations are `path:line` in shallow clones of each tool, as in the other files here.

## Three different things get called "strand"

1. **Record orientation.** Whether SEQ, as stored in the BAM, is reverse-complemented
   relative to the reference. SAM says this is FLAG 0x10. alnbase needs it to walk the
   record and to count `off_5p`/`off_3p` from the right end of the read as sequenced.
2. **Conversion strand**, two-way. Which reference strand carried the conversion, i.e.
   whether the informative reference base is C (the read shows C or T) or G (G or A).
   This is all a caller needs to decide whether a read base is evidence, and it is what
   BISCUIT reduces every input to (`methylation-tools.md`, BISCUIT row).
3. **Strand of origin**, four-way: OT, CTOT, OB, CTOB. Adds which molecule of the
   converted duplex the read came from. Needed for Bismark-compatible `XR`/`XG` tags,
   for per-strand reporting (Bismark's 12 default output files), and for M-bias by strand.

In a directional library, 2 and 3 follow from 1 plus the first/last-in-template flag, which
is what alnbase assumes today. That assumption can fail independently at either end: the
library may not be directional, so read number no longer selects the strand; or the aligner
may not use FLAG 0x10 conformantly, so even orientation cannot be read off the flag.

alnbase already keeps 1 and 2 apart internally, and a replacement has to keep them apart.
The walk reports a record along its **conversion** strand — `is_ot = is_last_in_template()
== is_reverse()` (`src/alignment.rs:108`), which sends both mates of a fragment the same way
along the reference, complementing read and reference bases together on the bottom strand —
while `off_5p`/`off_3p` stay the sequencing cycle counted from the read's own 5' end,
whichever way the walk runs (`docs/reference/examples/03-walk/01-orientation`). So a
declared rule has to yield, for each record:

- **a** the conversion strand, which sets the walk's direction and which reference base is
  informative;
- **b** the read's own sequenced direction, for `off_5p`/`off_3p` — today FLAG 0x10, and the
  one thing BSBolt's output gets wrong;
- **c** a four-way label, for `[tag.XX.strand]` values and the `conv_strand` column
  (`src/tags.rs:140`, `src/record_columns.rs:277`, both hardcoded to directional today).

In a conformant directional BAM, **b** is FLAG 0x10 and **a** is `0x10 == 0x80`. Nothing
else in alnbase reads the FLAG for strand, so these three are the whole contract.

## What each aligner writes

Read from the sources listed in `README.md`'s clone directory; not run. Where a claim would
be a bug in the aligner, treat it as a hypothesis until a validator reproduces it.

| Aligner | Strand tag and values | Resolves | FLAG 0x10 | Library modes |
|---|---|---|---|---|
| BSMAP 2.6 / 2.90 / BSMAPz | `ZS:Z:` = `++`, `+-`, `-+`, `--`. First character is the conversion strand (Watson `+`, Crick `-`), second the read's orientation within it, so `++`=OT, `+-`=CTOT, `-+`=OB, `--`=CTOB (`param.cpp:234`, emitted `align.cpp:690`, `pairs.cpp:354,410,490`; README.txt:121-124) | all four, tag alone | conformant: set iff `ZS[0] != ZS[1]`, and SEQ is the precomputed reverse complement then (`align.cpp:564`, `align.h:252-254`) | `-n 1` non-directional. Even at the default `-n 0`, a PE run emits all four `ZS` values, because read 2 is searched against the complementary chain (`align.cpp:93-94`). No PBAT. Run 2026-09-17: see "What running BSMAP found" below |
| BS-Seeker2 | `XO:Z:` = `+FW`, `-FW`, `+RC`, `-RC` single-end (`bs_single_end.py:415-432`); `+FR`, `-FR`, `+RF`, `-RF` paired (`bs_pair_end.py:531-545`) | all four single-end; paired needs `XO` plus 0x40/0x80, because both mates carry one `XO` | SEQ is stored reference-forward and 0x10 follows the genome strand, not the record's orientation (`output.py:54-61`), so the two `RC` classes are inverted. Paired flags are literals `67/115` and `131/179` (`bs_pair_end.py:607-612`), so both mates always share one 0x10 state. Its own caller reads only `is_reverse` (`bs_seeker2-call_methylation.py:335-340`) | `-t Y` non-directional. No PBAT option; the README's recipe reverse-complements mate 1 beforehand and leaves no trace in the BAM |
| methylpy (bowtie2, 3-letter) | none: the converted read carries its C positions in the read name and bowtie2's own tags pass through (`call_mc_pe.py:851`, `call_mc_se.py:731-795`) | conversion strand only, and only after the flip below | the output BAM is conformant; the intermediate it pileups is not. `flip_read2_strand` flips 0x10 on every second-in-pair record without touching SEQ (`call_mc_pe.py:1238-1290`, called `:1343`), so that 0x10 alone gives the conversion strand for `samtools mpileup` | `--pbat` swaps the mate files (PE) or reverse-complements each read (SE) and marks nothing in the output. Non-directional unsupported |
| dnmtools (abismal) | `CV:A:` = `T` (T-rich read) or `A` (A-rich read), the read's own conversion as sequenced, written by abismal itself (`abismal.cpp:480,:660,:688`) and translated in from BSMAP `ZS[1]` or Bismark `XR` only inside `dnmtools format` (`bam_record_utils.cpp:833-913`, called only from `format-reads.cpp:339,348`) | conversion strand only, and **only when combined with 0x10** -- `CV` alone is the read's richness, which outside PBAT modes is the mate number (`abismal/README.md:112-117`). abismal takes the XOR itself (`abismal.cpp:1269-1270`) | non-conformant **before** `format`, not only after: abismal stores SEQ exactly as it appears in the input FASTQ on every record, reverse ones included (`abismal/README.md:108-109`), so 0x10 means what SAM says but SEQ does not. `format` then merges the two mates into one N-gapped record, drops 0x1 and writes `CV:A:T` on everything; dnmtools' own tools read SEQ backwards without complementing (`methcounts.cpp:299-325`) | `-P` PBAT and `-R` random PBAT; only under those does `CV` vary with the read. Run 2026-09-17: the rule as first written was wrong on every reverse record; see "What running abismal found" below |
| meRanTK (meRanGh, meRanGs) | `YG:Z:` and `YR:Z:`, each `C2T` or `G2A`: `YG` the converted genome index that won, `YR` the conversion applied to the read (`meRanGh.pl:1432-1433,3029-3060`) | **all four strands.** `YG` is the fragment's strand and `YR` the mate, so read 2 of a pair is a CTOT or CTOB record even though the *fragment* can only be OT or OB -- wrong-direction alignments are discarded before output (`meRanGh.pl:2030-2100`) | conformant, except meRanGs single-end with `-r` only, where the single-end registration takes no direction argument (`meRanGs.pl:2035-2064`) and SEQ comes out reversed on **every** record, `NM:i:0` and all | directional only, fixed by `-f`/`-r`. meRanT writes no conversion tag |
| asTair (TAPS) | none, in or out: it reads `MD`, `AS`, `XS`, `NM` and ignores `XG`, `XR`, `YD`, `ZS` (`filter.py:105`, `summary.py:156-163`) | conversion strand, from whole-FLAG equality | conformant, measured. The strand comes from a literal table of FLAG integers — 99/147 and 0 are OT, 83/163 and 16 are OB (`caller.py:235-267`) — and any other FLAG contributes nothing. It does not reach depth either: pysam's pileup drops a non-proper pair first, and the same fixture with and without one gives byte-identical `.mods` files | `--library reverse` flips the polarity globally, the nearest thing to PBAT. Non-directional is "under development", and TAPS gives it no help — the reference is unconverted, so copy-strand reads map normally and are counted on the wrong strand. `--ignore_orphans False` raises `TypeError` at `caller.py:247` and is swallowed per position by `except Exception: continue` (`caller.py:339`), so the run writes an empty `.mods` and exits 0 |

| Bismark 0.25.1 | `XR:Z:` = `CT`/`GA`, the conversion this record's read shows (`bismark:8611`, `:9059`); `XG:Z:` = `CT`/`GA`, the converted genome index that won, the same value on both mates (`:8615`, `:9064`). `YS:Z:` = `OT`/`CTOT`/`OB`/`CTOB`, but only under `--strandID` and only paired (`:8744-8761`) | single-end: all four, from `XR` and `XG` together. Paired by default: not recoverable, see below | single-end: 0x10 is set from the genome conversion, so it equals `XG:Z:GA` rather than SEQ's orientation (`bismark:8518-8546`) — OT 0, CTOB 16, OB 16, CTOT 0, with SEQ always reference-forward (`:8578-8582`). Paired by default: conformant (99/147, 83/163), but 0x40 and 0x80 are deliberately swapped on CTOT and CTOB pairs so that browsers do not call them discordant (`:8818-8867`). `--old_flag` restores honest mate bits and makes 0x10 wrong on one mate of every pair | `--non_directional` searches all four indexes, `--pbat` only CTOT and CTOB (`:7199-7240`). Neither changes how the FLAG is built, so a single-end `--pbat` run has a non-conformant 0x10 on every record |
| BSBolt 1.6.0 | `YS:Z:` = `W_C2T`, `W_G2A`, `C_C2T`, `C_G2A`: reference contig (Watson/Crick) then the conversion applied to the read (`bwamem.c:985-1008`). `XG:Z:` = `CT`/`GA` is the contig only, written for MethylDackel. Unmapped reads can carry a fifth value, `YS:Z:WC` (`bs_sorter.cpp:79`) | all four, tag alone | not conformant: `p->is_rev` is overwritten with "aligned to the Crick contig" (`bwamem.c:850,866`) while SEQ is reverse-complemented by the true orientation (`:938-963`), so the bit is inverted on every `*_G2A` record — read 2 of every directional pair. Also 0x20 is never set on a Watson record (`:854,868`) | `-UN` non-directional, with the conversion inferred per fragment from base composition and a hardcoded tie window (`bs_helpers.cpp:41-62`). No PBAT |
| bwa-meth 0.2.10 | `YD:Z:` = `f`/`r`, the converted contig hit (`bwameth.py:480`). `YC:Z:` = `CT`/`GA` is assigned from the mate number alone (`:205`) and carries nothing the FLAG does not | conversion strand only | conformant: bwa's flag is untouched and SEQ is restored from the original read, reverse-complemented iff 0x10 (`:500-503`) | directional only (`README.md:28`). Run 2026-09-17: holds on a directional library; see "What running bwa-meth found" below |
| BISCUIT 1.10.3 | `YD:A:` = `f`, `r`, or `u`, where `u` means no conversion event was seen anywhere in the alignment (`mem_alnreg_format.c:434-437`). `ZC:i:`/`ZR:i:` count converted and retained cytosines. It never writes `ZS`; it only reads it from BSMAP input | conversion strand alone; all four with 0x10, the truth table being in the source (`mem_alnreg.c:422-429`). Not recoverable for `YD:A:u` | conformant (`mem_alnreg_format.c:81-84`, `:338,357`) | `-b 0`, non-directional, is the default; `-b 1` directional. No PBAT mode: `scripts/flip_pbat_strands.sh` flips 0x10 with awk and leaves SEQ, POS and `YD` alone |
| HISAT-3N | `YZ:A:` = `+`/`-`, the 3N index hit (`alignment_3n.h:717-718`), with `Yf:i:` and `Zf:i:` counting conversions | conversion strand only | conformant (`aln_sink.h:2574-2575`, `:2689-2692`) | non-directional is the default; `--directional-mapping` and an undocumented `--directional-mapping-reverse` (PBAT) select two cycles each (`hisat2.cpp:1983-1994`). In non-directional mode `YZ` is inferred from conversion counts and falls back to the directional rule when the counts tie, so a fully methylated read is assigned a strand by assumption (`alignment_3n.h:299-308`). Run 2026-09-17: holds wherever one conversion is visible; see "What running HISAT-3N found" below |

Two sources of the same information are worth naming apart, because aligners split them
differently: a **genome tag** says which converted reference the read matched (Bismark `XG`,
BSMAP `ZS[0]`, BSBolt's `YS` prefix, meRanTK `YG`, bwa-meth `YD`, BISCUIT `YD`, HISAT-3N
`YZ`) and a **read tag** says which conversion the read itself shows (Bismark `XR`, BSMAP
`ZS[1]`, BSBolt's `YS` suffix, meRanTK `YR`, dnmtools `CV`). Either one alone gives the
conversion strand; the two together give the four-way strand without consulting the FLAG.

## The cases

**Every one of these aligners stores SEQ in the reference's forward orientation.** That is
what SAM requires and none of them break it, including the ones whose FLAG is wrong. So the
walk direction — output **a** — never has to come from 0x10 as such; it needs the conversion
strand, however that is expressed. What 0x10 is needed for is output **b**, which end of the
read the sequencer started at, and that is where the non-conformant aligners hurt.

1. **Conformant FLAG, no tag, directional library.** bwa-meth, BISCUIT with `-b 1`,
   HISAT-3N with `--directional-mapping`, methylpy's own output. The current rule —
   `0x10 == 0x80` for the strand, 0x10 for the orientation — is right, and it is the only
   thing available.
2. **Conformant FLAG, two-state tag.** bwa-meth `YD`, BISCUIT `YD`, HISAT-3N `YZ`. The tag gives the conversion strand without assuming the library is directional, and
   0x10 still gives the orientation. The four-way strand needs 0x10 as well (BISCUIT) or the
   mate flag, and for some inputs it is simply not there. dnmtools' `CV` looks like this case
   and is not: it names the read's own conversion, so the conversion strand is `CV` combined
   with 0x10 rather than read off the tag.
3. **Conformant FLAG, four-state tag.** BSMAP `ZS`, Bismark via `XR` plus `XG`, and meRanTK
   via `YG` plus `YR`. Everything is available from the tag; the FLAG only has to be trusted
   for **b**. The two-tag members of this case pair a genome tag with a read tag, and the
   read tag is the mate number rather than a measurement -- which costs nothing, because a
   record's strand of origin is exactly the fragment's strand combined with which mate it is.
4. **FLAG whose 0x10 means the conversion strand.** Bismark single-end (all of it; only
   `--pbat` and `--non_directional` runs actually produce records where this differs from
   the conformant answer), BSBolt paired-end read 2, BS-Seeker2's `RC` classes and all of
   its paired output. Here **a** can be read straight off 0x10 and **b** cannot be read off
   anything: it has to be reconstructed, which for a directional paired-end library means
   `0x10 xor 0x80`.
5. **FLAG deliberately rewritten downstream.** methylpy's `read2flipped.bam` and BISCUIT's
   `flip_pbat_strands.sh` both flip 0x10 without touching SEQ, producing case 4 on purpose;
   dnmtools `format` reverse-complements SEQ so that its own tools can read it backwards
   without complementing — and abismal, upstream of it, never complemented SEQ in the first
   place, so a raw abismal BAM is already in that state. These are intermediates, but they are
   also what a user may hand us.
6. **Mate bits rewritten.** Bismark's default paired output swaps 0x40 and 0x80 on CTOT and
   CTOB pairs, so on a non-directional library an OT record and a CTOT record are identical
   in FLAG, `XR` and `XG` together; only the order in the file, `--strandID`'s `YS`, or
   `--old_flag` tells them apart. A rule that reads the mate bit is reading a rewritten one.

Two consequences for whatever syntax follows. First, the rule has to be able to produce
**b** from something other than 0x10, because in case 4 the bit means something else.
Second, a rule must be allowed to say that the four-way strand is **unknown** for an input
that cannot support it: on bwa-meth, HISAT-3N or dnmtools output, only two of OT, CTOT, OB
and CTOB are ever determinable, so a `[tag.XX.strand]` table that names four values would be
writing a distinction its input does not carry.

## One thing this survey found in alnbase's own documentation

The Bismark example at `docs/reference/examples/04-outputs/01-bismark-tags/queries.toml`
(and the same fixtures in `src/bam_out.rs:637` and `src/golden.rs:109`) declares

```toml
[tag.XR.strand]             # read conversion, as the read was sequenced
CT = ["OT", "CTOB"]
GA = ["OB", "CTOT"]
```

which is right for the top-strand pair and inverted for the bottom-strand one. Real Bismark
output, checked on the BAMs under the earlier Bismark runs, pairs the flags and tags like
this — directional paired-end on the left, single-end on the right:

```
 99 XR:Z:CT XG:Z:CT   (OT)        0 XR:Z:CT XG:Z:CT   (OT)
147 XR:Z:GA XG:Z:CT   (CTOT)      0 XR:Z:GA XG:Z:CT   (CTOT)
 83 XR:Z:CT XG:Z:GA   (OB)       16 XR:Z:CT XG:Z:GA   (OB)
163 XR:Z:GA XG:Z:GA   (CTOB)     16 XR:Z:GA XG:Z:GA   (CTOB)
```

so `XR` is `CT` for OT and OB — the two reads that are converted originals, T-rich as
sequenced — and `GA` for CTOT and CTOB, their A-rich complements. That is what the example's
own comment says, and it is the sense dnmtools relies on when it turns `XR` into `CV`
(`bam_record_utils.cpp:873`). `XG` is correct in the example. The fix is to swap the two
strand lists under `[tag.XR.strand]`; the example's `expected.txt` has to be regenerated
with it.

## What running BSMAP found: a conversion strand it names wrongly

BSMAP 2.90, aligned and measured in `../../../validation/demos/bsmap-strand/`. At `-n 1`,
some CTOT and CTOB reads come back with `ZS` naming the opposite conversion strand — at the
correct position, with a unique MAPQ and only an inflated `NM` to show for it. In a sweep of
a 60 bp read from each strand of origin every 10 bases along a 1200 bp reference, 10 of the
228 copy-strand reads were misassigned and none of the 228 OT or OB reads was.

Every misassignment flips **both** characters of `ZS`: CTOT (`+-`) is called OB (`-+`), CTOB
(`--`) is called OT (`++`). That is the signature of the mechanism. BSMAP searches a read in
two chains, as sequenced and reverse-complemented (`align.cpp:92-93`), against two converted
reference chains, and the two indices become `ZS` (`align.cpp:269`). Flipping both places the
read at the same forward reference position, so the right and wrong candidates collide there,
and `AddHit` keeps whichever was found first (`align.h:232-249`):

```c
if(!hitset[_ghit.chr>>1].insert(_ghit.loc).second) return 0; //hit already exist
```

`>>1` drops the reference-chain bit, so the duplicate key is the forward position alone and
the two candidates are never compared on mismatches. The as-sequenced chain is searched first
(`align.cpp:211`), and for a copy-strand read that is the wrong chain, so whenever the wrong
chain fits within the mismatch budget — 5 for a 60 bp read at the default `-v 0.08`,
`align.cpp:501` — its hit is registered first and the perfect alignment is discarded as a
duplicate. Predicting a misassignment from that budget alone agrees with what BSMAP did on
226 of the 228 copy-strand reads. An OT or OB read is never affected, because the chain
searched first is its own.

This is a property of `-n 1`. A directional library at the default `-n 0` is safe: for a
single read only the as-sequenced chain is searched, so copy-strand reads are outside the
search space by design rather than by accident. Paired-end runs enable both chains for read 2
at `-n 0` too, but read 2 of a directional pair *is* a copy-strand read whose correct chain is
the second one — and the demo saw no misassignment in either paired run, where the mate
anchors the placement.

The consequence for a caller is the whole read: every cytosine in it is scored on the wrong
strand. No rule written over `ZS` can detect this, because `ZS` is the thing that is wrong.
What does detect it is counting read bases that differ from the reference after the walk —
alnbase's scan line goes from `0 of 240` to `5 of 480` between the two runs in that demo,
and those 5 are the misassigned read's own bisulfite conversions.

## What running bwa-meth found: a documented limitation that leaks

`validation/demos/bwameth-strand/` aligned the same four-strand fixture with bwa-meth 0.2.10.
On the directional library the table row above holds without qualification: all four records
carry the conversion strand and the sequenced direction their origin implies, `YD` is the same
on both mates of every pair, 0x10 agrees with how SEQ was really stored on every mapped record,
and SEQ is reference-forward throughout. `YC` was `CT` on every first-in-pair record and `GA`
on every second-in-pair one whatever the strand the record came from, which confirms by
measurement what the row asserts from the source: it names the conversion applied to the read,
not the strand the read came from.

The part worth recording is what "directional only" does in practice. bwa-meth converts read 1
C>T and read 2 G>A (`bwameth.py:205`) on the assumption that read 1 came from an original
strand, so a CTOT or CTOB read arriving as read 1 gets the wrong conversion applied. Sweeping a
60 bp read from each strand of origin every 10 bases along a 1200 bp reference, aligned
single-end:

| strand of origin | reads | mapped | conversion strand wrong | misplaced | clipped | median MAPQ |
|---|---|---|---|---|---|---|
| OT | 114 | 114 | 0 | 0 | 0 | 60 |
| CTOT | 114 | 15 | 15 | 0 | 0 | 56 |
| OB | 114 | 114 | 0 | 0 | 0 | 60 |
| CTOB | 114 | 17 | 17 | 0 | 0 | 52 |

93% of the copy-strand reads are dropped, which is the outcome "directional only" implies.
The other 7% is not a refusal: every copy-strand read that maps at all is on the opposite
conversion strand, at the right position, with a full-length match, MAPQ in the fifties and no
QC-fail flag. A read survives precisely because the wrong conversion made it look like a clean
hit on the other contig, so the survivors are selected for being indistinguishable from correct
ones. This is unlike the BSMAP case above, where the aligner found the right alignment and
threw it away; here the right alignment is never in the search space. The consequence for a
caller is the same, and so is the remedy: nothing downstream can tell these records apart, so
the library has to be filtered by design.

Two smaller facts from the same run. bwa-meth appends `YD` only to mapped records
(`bwameth.py:469-473` returns before the tag), so an unmapped record is `unknown` to any rule
over `YD` — three of the eight paired-end records, which alnbase skipped and counted. And
`directional.toml` over the same BAM produced the same hits at the same reference positions as
`bwameth.toml`, 0 disagreements over 29 and 37 hits, while additionally naming the strand of
origin. Both rules trace back to the contig bwa-meth chose, one through `YD` and one through
the 0x10 that same choice set, so they agree even on the record it placed on the wrong strand:
two rules agreeing establishes that the BAM is self-consistent, not that it is right.

## What running HISAT-3N found: a tag that is an inference, and what happens when it has nothing to infer from

HISAT-3N 2.2.1-3n-0.0.3, built from commit `f5dda37`, was run on the four-strand fixture in
`validation/demos/hisat-3n-strand/`. On ordinary partly converted reads the rule holds
exactly: all eight records called with the conversion strand and the sequenced direction their
origin implies, `YZ` the same on both mates of every pair, 0x10 agreeing with how SEQ was
really stored, SEQ reference-forward throughout.

What distinguishes `YZ` from the other two-state tags in this survey is that it is a
conclusion rather than a record. BISCUIT and bwa-meth write down which converted contig the
read was aligned to. HISAT-3N aligns to both three-letter indexes and then decides by counting
the alignment's C→T and G→A differences, so `YZ` is derived from the same evidence a caller is
about to read, and inherits that evidence's gaps. A sweep of 456 reads in four states measured
where the gap is:

| state | reads | conversion strand wrong | of those, carrying `Yf:i:0` |
|---|---|---|---|
| every CH converted | 456 | 0 | 0 |
| a single converted cytosine | 456 | 0 | 0 |
| nothing converted (fully methylated) | 456 | 228 | 228 |
| one conversion of each kind | 456 | 228 | 0 |

One converted cytosine is enough — the inference is not a vote that needs depth, so the tie is
its only failure mode. But when the counts do tie, HISAT-3N does not decline: it writes the
`YZ` a directional library would have implied, which is right for OT and OB and backwards for
every read from a PCR copy. That is half the sweep, at the right position, unclipped, MAPQ 60.
Aligning the fully methylated fixture with `--directional-mapping` produced the same `YZ` on
all eight records, which is the finding stated plainly: **in the absence of evidence the
default non-directional mode is the directional mode, applied without having been asked for.**

Two signals bear on it and neither covers both ties. `Yf` is the number of conversions of the
kind `YZ` names, so `Yf:i:0` says outright that the strand was assumed — it marks the fully
methylated case exactly, and says nothing about the equal-counts case, where `Yf:i:1` is
indistinguishable from a genuine single conversion. alnbase's count of read bases differing
from the reference, which is what surfaced BSMAP's wrong-strand records, reads `0 of 480` on
the fully methylated fixture while half its records are on the wrong strand: a read with no
conversion matches the reference whichever way it is walked. On the equal-counts fixture that
count is non-zero, but splitting the records into the ones the rule got right and the ones it
got wrong gives `4 of 240` for each half — it is counting the substitution, which both halves
carry, not the strand.

This is a third shape of problem, distinct from the two already in this document. BSMAP's `ZS`
is wrong about a read whose correct alignment it found and discarded; bwa-meth's limitation is
documented and its leak is a read it should not have been given. HISAT-3N's tag is wrong about
a read it aligned perfectly, in a mode it advertises for exactly that library, and the remedy
is neither a better rule nor a filter on the output but an argument at the aligner:
`--directional-mapping` when the library is directional.

One incidental finding: two comments in `alignment_3n.h` invert the rule the code implements.
The header of `makeYZ` (`:285-286`) and the declaration of `YZ` (`:95-96`) both say `+` is for
the smaller C→T count, while `:309` assigns `+` when it is the larger. The code is right, and
the run confirms it. A rule written from the comments rather than the code would have the tag
backwards on every record.

## What running abismal found: a rule that was wrong, and a SEQ convention no rule can express

`queries/strand/dnmtools.toml` was written from this table and read `CV` on its own: `CV:A:T`
as conversion strand `+`, `CV:A:A` as `-`. Running abismal 3.3.0 on the four-strand fixture
(`validation/demos/dnmtools-strand/`) showed that backwards on every reverse record. `CV` is
the conversion the **read** shows as sequenced, not the reference strand that carried it; the
conversion strand is `CV` combined with 0x10, which is how abismal computes it internally
(`abismal.cpp:1269-1270`). Under `-R`, where all four pairs map at `NM:i:0`, the old rule got
4 of 8 records wrong and the corrected one 0 of 8. The file now declares the combination, and
`queries/strand/records/dnmtools.sam` was corrected with it.

Two things about the tag were also wrong in the reading rather than the rule. The first is
who writes it: `standardize_format`, the function this table cited, is reachable only from
`dnmtools format`, so the BSMAP and Bismark translations never produce an *unformatted* BAM
carrying `CV`. Every such BAM comes from abismal. The second is how much the tag says: outside
PBAT modes abismal maps read 1 T-rich and read 2 A-rich and nothing else, so `CV` is the mate
number, carrying no more than bwa-meth's `YC`.

**The finding that outlives the rule is the SEQ convention.** abismal stores SEQ exactly as it
appears in the input FASTQ, on reverse records too, which its README states as a feature and
which its own caller compensates for by reading SEQ backwards without complementing. It is not
what SAM says, and it breaks the invariant every rule in `queries/strand/` depends on — before
`dnmtools format`, not only after, which is the opposite of what this table used to claim. A
strand rule cannot express it, because the rule declares which strand a record belongs to and
not how its bases are stored; the file has to be repaired first, which is a reverse complement
of SEQ and a reversal of QUAL on each 0x10 record.

alnbase's count of read bases differing from the reference separates all four combinations,
and by size rather than merely by presence: on the fixture the raw BAM reads 155 of 480 bases
differing (32%) whichever rule is used, the repaired BAM with the old rule reads 38 of 480
(8%), and the repaired BAM with the corrected rule reads 0. A third of all bases is the SEQ
convention; a twelfth is a wrong strand. That is the sharpest contrast yet with the HISAT-3N
round, where the same count could see nothing at all.

## What running meRanTK found: a comment that was wrong, and a bug in one of the two tools

The meRanTK row above had the tags right and their consequences wrong. Running meRanGh and
meRanGs on the four-strand fixture (`validation/demos/merantk-strand/`), three modes each,
called all twenty mapped records as the strand they were really sequenced from — the rule's
four-way table needed no change, unlike the dnmtools round before it.

What the row claimed and the run refuted is that "CTOT and CTOB never appear". They appear
constantly. `YR` is never a per-read measurement: paired-end it is the mate number, `C2T` on
read 1 and `G2A` on read 2 unconditionally, and single-end it is a run-level constant set from
whether `-r` was given. So half of every paired-end run's records and all of a `-r` run's are
CTOT or CTOB. This is the third tag in the survey to turn out to be the mate number, after
bwa-meth's `YC` and abismal's `CV`, and the first where it costs nothing: `YG` carries the
fragment's strand, `YR` carries the mate, and a *record's* strand of origin is exactly the two
combined. The tag pair is genuinely four-state; only the fragment is restricted to two.

**meRanGs single-end with `-r` alone writes SEQ the wrong way round on every record.** Its
single-end registration takes no direction argument (`meRanGs.pl:2035-2064`) and stores the
read as sequenced on a `C2T` hit and reverse-complemented on a `G2A` hit; the `-r`
wrong-direction filter keeps the opposite orientations, so both halves come out reversed.
meRanGh's equivalent does take the argument (`meRanGh.pl:1988-2018`) and is correct in all
three modes. Nothing inside the BAM contradicts itself — `NM` and `MD` are computed against
the converted genome, where the alignment really is perfect — so a reader that trusts `NM`
sees flawless records. alnbase's count of read bases differing from the reference reads 81 of
120 on that run and 0 of 120 on meRanGh's, and the two tools, which are meant to be
interchangeable, emit identical hit sets in the two modes where meRanGs is correct and 7
against 13 in the broken one. The seven are not a safe subset: they sit at the right
coordinates but take their base from the wrong end of the read.
