# HISAT-3N's MD tag, probed (H3-1, H3-2, H3-3), reproduced

HISAT-3N aligns against a converted index but rebuilds the `MD` tag itself, against the
original reference, in `alignment_3n.h`. That rebuild is the only thing that tells a
downstream tool which reference base a read position sat on, so it is worth checking that
the string it writes describes the alignment it writes.

This also settles an open lead in the bam2bakR section: `mut_call.py` depends on `MD` and
raises `ValueError: MD tag not present` without it, and the HISAT-3N path of bam2bakR does
not request the tag. HISAT-3N writes it unconditionally, against the **original**
reference, so the conversions bam2bakR is counting are visible. See H3-1 below.

Full entries: `docs/design/research/non-methylation-tool-findings.md` (H3-1..H3-4).

Tool: `DaehwanKimLab/hisat2` branch `hisat-3n` at GitHub `f5dda37` (2022-10-13), built from
source with `make hisat2-align-s hisat2-build-s`. Checked 2026-09-18.

```
HISAT3N=<built hisat-3n checkout> PY=<python with pysam> ./run.sh [workdir]
```

`expected.txt` is the output of that run. `make_reads.py` writes a 20 kb random reference
and eleven 100 bp reads, each carrying a legal `T`→`C` conversion at read offset 49 followed
immediately by a 2 bp deletion. Sites where the reference has a `C` just downstream are
skipped, because HISAT-3N would otherwise shift the deletion to absorb the converted base
and align the read with no mismatch at all.

## 1. MD is built against the original reference — clean, and it closes the bam2bakR lead

`hisat-3n-build` writes the *converted* sequence into the FM index but deliberately
restores the original before writing the packed reference:

```c
if (threeN) {
    // save the unchanged reference in .3.ht2 and .4.ht2
    baseChange.restoreNormal();
    ...
    baseChange.restoreConversion();
}
```

`constructMD` reads that packed reference back and rebuilds `MD` from scratch, discarding
whatever HISAT2 computed in converted space. The run confirms it: every read comes out with
`MD:Z:49T...`, and `T` is the original reference base at the converted position.

## 2. A mismatch immediately before a deletion loses its separator — bug

Seven of the eleven reads:

| read | CIGAR | HISAT-3N | samtools calmd |
|---|---|---|---|
| s611_delGA | `50M2D50M` | `49T^GA50` | `49T0^GA50` |
| s3625_delGG | `50M2D50M` | `49T^GG50` | `49T0^GG50` |
| s8557_delTA | `50M2D50M` | `49T^TA50` | `49T0^TA50` |

`samtools calmd` reports each one:
`[bam_fillmd1] different MD for read 's611_delGA': '49T^GA50' -> '49T0^GA50'`.

The mismatch branch of `constructMD` emits a `0` when the previous character is a letter
(`alignment_3n.h:600-601`); the deletion branch (`:619-629`) does not. So a deletion that
follows a mismatch is appended straight onto the mismatch letter. `49T^GA50` is not a
production of the MD grammar `[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*` at all — there is no digit
between the mismatch letter and the `^`.

The other direction is handled correctly: `49^TA0T50` and `48^CT0T0T50` both carry the `0`
that a mismatch after a deletion needs, and agree with samtools.

Honest sizing: pysam reads all eleven correctly, because its MD parser is CIGAR-driven and
takes the deletion length from the `D` operator rather than from the MD string. bam2bakR and
SlamDunk go through exactly that path, so neither is affected. The exposure is a consumer
that parses MD on its own terms. What makes it worth recording anyway is the frequency: in
3N mode every conversion is a mismatch, so in bisulfite data essentially every deletion
preceded by a reference `C` produces a malformed tag, rather than the rare
mismatch-abuts-indel case an ordinary aligner would hit.

## 3. The repeat-index builder drops the match run before a deletion — bug

`constructRepeatMD` has the same missing separator and one more: its deletion branch
(`:504-508`) never flushes the pending match count before appending `^`. The count is not
reset either, so it carries into the run after the deletion. A `5M2D13M` alignment with no
mismatch at all comes out as `^CG18` instead of `5^CG13` — the leading `5` is gone and the
trailing `13` has absorbed it.

This one pysam does not survive:

```
CIGAR 5M2D13M  MD ^CG18
pysam reads the reference as  CGGTATACGTACGTACGT
the reference actually is     ACGTATACGTACGTACGT
```

The two deleted bases are handed to the first two aligned positions, so a caller reading
reference bases through `get_aligned_pairs(with_seq=True)` sees two mismatches that are not
there, at the start of every affected read.

Not reproduced against the aligner: the repeat path only fires for alignments to a repeat
index, `hisat-3n-build --repeat-index` is opt-in (`hisat2_build.cpp:92, :416`), and it
aborts on a random sequence with no repeats. `probe_repeat.py` therefore writes the MD
string that branch constructs and measures what the consumers do with it. When a repeat
index *is* present, it is used by default at align time (`hisat2.cpp:566`).

## 4. NM excludes legal conversions — by design, not a finding

`samtools calmd` disagrees about `NM` on all eleven reads (`2 -> 3`), and on the read with
two conversions (`s9790_delTT`, `2 -> 4`). That is not the same defect. HISAT-3N's `NM`
comes from the alignment in *converted* space, where a legal conversion is a match, and
`constructMD` adds back only conversions on the wrong strand (`:651-660`). So `NM` here is
"real mismatches plus indel bases", which for a 3N alignment is the more useful number and
is consistent with `XM:i:0` and `Yf:i:1` on the same records. Recorded so that the `NM`
disagreement in `expected.txt` is not mistaken for part of finding 2.
