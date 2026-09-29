# meRanTK's `YG`/`YR`, checked against meRanGh and meRanGs — and a bug in one of them

[`queries/strand/merantk.toml`](../../../queries/strand/merantk.toml) is the ninth of the
twelve strand rules to be run against its aligner's real output. Its four-way table came
through unchanged — every one of the twenty records this demo produces is called the strand
it was really sequenced from — but the file's comment was wrong about how often three of the
four cells are used, and the run turned up a real bug in meRanGs.

```
MERANTK_SRC=<meRanTK checkout> MERANTK_ENV=<env> PY=<python with pyarrow> \
  ALNBASE=<alnbase> ./run.sh [workdir]
```

`MERANTK_SRC` is a clone of <https://github.com/icbi-lab/meRanTK>, run against commit
`c06d726` (meRanTK 1.3.0). The environment is [`../../envs/merantk.yaml`](../../envs/merantk.yaml).

meRanTK is the only RNA tool in this series: it maps bisulfite-treated RNA to find 5mC in
mRNA, so both of its genome tools are spliced aligners — meRanGh drives hisat2, meRanGs
drives STAR. The fixture has no introns, because what is under test is the tag.

## The table is right, in all three modes

`YG` names the converted genome index that won and `YR` the conversion applied to the read,
the same split as Bismark's `XG` and `XR`, so the two together name the strand of origin of
the record. Both tools were run three ways — paired-end, `-f` alone, `-r` alone — and every
record landed in the cell its origin implies:

```
       gh_pe   gh_f   gh_r   gs_pe   gs_f   gs_r
OT       1       1      -       1      1      -
CTOT     1       -      1       1      -      1
OB       1       1      -       1      1      -
CTOB     1       -      1       1      -      1
```

`strand of origin wrong: 0` on all six runs, twenty records in total. The `CTOT` and `CTOB`
input pairs do not map at all: meRanTK is directional, and its wrong-direction filter
(`meRanGh.pl:2030-2100`) discards them before output.

## What the comment got wrong: `YR` never varies with the read

The file said "in practice only OT and OB occur". That is true of one mode out of three, and
the reason is that `YR` is not a per-read measurement in any mode:

- **Paired-end**: `YR` is the mate number. Read 1 gets `C2T` and read 2 gets `G2A`,
  unconditionally, in both tools (`meRanGh.pl:3029-3060`, `meRanGs.pl:3056-3087`). So half of
  every paired-end run's records are `CTOT` or `CTOB` — 2 of the 4 above.
- **Single-end**: `YR` is a run-level constant, assigned once before the read loop from
  whether the run was given `-r` (`meRanGh.pl:1432-1433`, `meRanGs.pl:1483-1484`). A `-f` run
  is `OT`/`OB` throughout and a `-r` run is `CTOT`/`CTOB` throughout.

This is the third tag in the series to be the mate number rather than a measurement, after
bwa-meth's `YC` and abismal's `CV`. The difference here is that it does not matter: `YG`
carries the fragment's strand, `YR` carries the mate, and the strand of origin of a *record*
is exactly the combination of the two. The rule reads both and is right. What is misleading
is only the claim about which cells are visited.

There is a second consequence. Because the wrong-direction filter keeps `C2T`-index hits only
when 0x10 is clear and `G2A`-index hits only when it is set (and the reverse in a `-r` run),
`YG` and FLAG 0x10 are redundant with each other within any one run. In single-end the whole
rule therefore reduces to reading 0x10 — the same call `directional.toml` makes. The tags are
worth reading anyway, because they say which reduction applies without the caller having to
know how the run was invoked.

## The bug: meRanGs single-end `-r` stores SEQ the wrong way round

meRanGh passes the read direction into its single-end SEQ registration
(`meRanGh.pl:1988-2018`) and stores SEQ reference-forward in all three modes. meRanGs'
registration takes no such argument (`meRanGs.pl:2035-2064`): it stores the read as sequenced
on a `C2T` hit and reverse-complemented on a `G2A` hit, full stop. That is right for a `-f`
run, where the surviving `C2T` hits are forward and the `G2A` hits reverse, and backwards for
**both halves of a `-r` run**, where the filter keeps the opposite orientations.

alnbase's count of read bases differing from the reference is what finds it:

| run | | bases differing from the reference |
|---|---|---|
| `gh_pe`, `gh_f`, `gh_r` | meRanGh, all three modes | 0 of 240, 0 of 120, 0 of 120 |
| `gs_pe`, `gs_f` | meRanGs, paired-end and `-f` | 0 of 240, 0 of 120 |
| `gs_r` | meRanGs, `-r` | **81 of 120 (67.5%)** |

Nothing inside the BAM contradicts itself. Both broken records carry `NM:i:0 MD:Z:60`,
because those are computed against the *converted* genome, where the alignment really is
perfect; it is only the SEQ field written out afterwards that is reversed. A reader that
trusts `NM` sees a flawless record.

In hits, the two tools agree exactly in the two modes where meRanGs is correct, and the
broken mode is worse than it looks:

| mode | meRanGh | meRanGs | in common |
|---|---|---|---|
| paired-end | 29 | 29 | 29 |
| `-f` | 16 | 16 | 16 |
| `-r` | 13 | **7** | 7 |

The seven that survive are not a safe subset. They are at the same reference positions as
seven of meRanGh's thirteen, but the base each one reports comes from the wrong end of the
read, so they are methylation calls at the right coordinates from the wrong cytosines. The
six that vanish are the visible half of the damage; these seven are the invisible half.

## Caveats, and what is not exercised

- meRanT, meRanTK's transcriptome aligner, writes no conversion tag and no rule here covers
  it. It is not run.
- meRanCall, meRanTK's own caller, is not run. This demo is about the tag, and the meRanGs
  `-r` bug is measured in the BAM rather than traced through to a call file.
- No introns, even though both tools are spliced aligners, and no gaps, trimming, repeats or
  multimappers. `-samMM` is not used.
- MAPQ is not exercised: meRanGh writes 60 and meRanGs writes 255 on every record here.
- `-forceDir`, the internal base-composition directionality filter, is left off.
- The fixture's copy-strand pairs never map, so nothing here measures what meRanTK does with
  a non-directional library. It is not built for one.

## Scope of the bug, and whether it reaches the literature

Followed up 2026-09-18, because a silent orientation error in an m5C tool is the kind of
thing that could explain disputed sites in a field with a well-known false-positive problem.
It does not, and the reason is worth writing down.

**The bug is older and wider than this demo's version.** The same defect is in meRanGs in
every release, and was in meRanGh too until it was fixed. Verified here against the source
on disk: `meRanGh.pl`'s `registerC2Tsingle` reads `my $readDir = $_[4]` and branches on it,
while `meRanGs.pl`'s takes no such argument and assigns `$seq->{Seq}` unconditionally.
meRanGs *computes* `$readDirection` (`meRanGs.pl:1235-1239`) and threads it into
`filterSEalignments` (`:1802`, `:2076-2114`), which uses it correctly to decide which
alignments to keep -- it is simply never forwarded to the subs that write SEQ. The `CHANGES`
file records the fix twice for the sibling tools and never for meRanGs:

    1.2.1a (2019/01/29)  meRanGh:   bugfix for situations where only reverse reads are provided.
    1.2.1b (2019/04/08)  meRanCall: bugfix for situations where only reverse reads are provided.

The 1.2.1b entry lists meRanGs as well, for an unrelated dovetailing fix. So the same class
of bug was reported by an outside user, fixed in two of the three tools that had it, and the
third instance was left in place and is still there.

**But it cannot be hiding in published m5C data, and the test is cheap.** The corruption
applies to every record in the run, spike-ins included, and it drives apparent methylation to
roughly half at essentially every cytosine rather than producing a sparse false-positive set.
Any study that reports a global or spike-in C-to-U conversion rate above about 95% is
therefore excluded by its own QC number. That is a single published figure, and the m5C
literature reports it almost universally, precisely because incomplete conversion is the
field's acknowledged artefact. A survey of the open-access papers citing meRanTK found one
single-end meRanGs study, and it reports >99% spike-in conversion.

**What the demo is good for instead.** This is a clean example of a BAM that is internally
consistent and externally wrong: `NM` and `MD` are computed against the converted genome, so
they are correct; no structural check on the file can fail; and the only instrument that sees
it is a comparison of the read against the *reference*, which is what alnbase does by
construction. The argument the demo supports is about the class of error -- orientation
treated as an implicit property of an aligner's output path rather than an explicit
declaration -- not about any particular published result.
