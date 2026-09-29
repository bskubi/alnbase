# dnmtools' `CV`, checked against abismal — and a rule this run found wrong

[`queries/strand/dnmtools.toml`](../../../queries/strand/dnmtools.toml) is the eighth of the
twelve strand rules to be run against its aligner's real output, and the first the run found
**wrong**. It read `CV:A:T` as conversion strand `+` and `CV:A:A` as `-`; on a reverse record
that is backwards. The file now reads `CV` together with FLAG 0x10, and this demo is the
check on the corrected version as well as the record of what the old one did.

```
DNMTOOLS_ENV=<env> PY=<python with pyarrow> ALNBASE=<alnbase> ./run.sh [workdir]
```

The environment is [`../../envs/dnmtools.yaml`](../../envs/dnmtools.yaml): dnmtools 1.5.0,
which ships abismal 3.3.0 inside it.

## Who writes `CV`

The rule was written believing dnmtools translates `CV` in from BSMAP's `ZS` or Bismark's
`XR`. It does — but only inside `standardize_format`, which is called from exactly two places
(`format-reads.cpp:339,348`), both in `dnmtools format`. So a BSMAP or Bismark BAM that
carries `CV` has already been formatted, and every *unformatted* BAM carrying `CV` was written
by **abismal**, dnmtools' own aligner, which assigns the tag directly
(`abismal.cpp:480,:660,:688`). The rule's subject was never the tools it named.

## What `CV` says, and what the old rule made of it

`CV` is the conversion the **read** shows as it came off the sequencer: `T` for a T-rich read,
`A` for an A-rich one. abismal's README says so outright — "This tag is independent of the
strand the read was mapped to" (`README.md:112-117`) — and adds that outside PBAT modes the
first end is always T-rich and the second always A-rich. The run confirms both halves.

In the default directional mode `CV` is the mate number and nothing else:

```
pair  record from  flag  pos  CV  want   shipped   corrected
OT    read 1 OT    99    101  T   +      +         +
OT    read 2 CTOT  147   241  A   +      -         +
OB    read 1 OB    83    741  T   -      +         -
OB    read 2 CTOB  163   601  A   -      -         -
```

Under `-R` (random PBAT) all four pairs map at `NM:i:0` and `CV` varies with the read:

| | shipped rule (`CV` alone) | corrected rule (`CV` with 0x10) |
|---|---|---|
| default mode, 6 mapped records | 3 wrong | 2 wrong |
| `-R`, 8 mapped records | **4 wrong** | **0 wrong** |

The conversion strand is `CV` combined with the mapped strand: a T-rich read placed forward
carries C→T against the reference, so `+`; the same read placed reverse reads G→A against the
reference, so `-`. That is how abismal computes it itself (`abismal.cpp:1269-1270`, the XOR of
`read_is_a_rich` and `read_rc`), and it is what the file now declares.

The two records the corrected rule still misses in the default run are the `CTOT` pair, which
is copy-strand input to a directional run — the same leak
[`../bwameth-strand/`](../bwameth-strand/) measured. Here it is not silent: both records carry
`NM:i:5` and `NM:i:6` and one is soft-clipped, so the mismatch count marks them. The `CTOB`
pair does not map at all.

## The larger problem, which no rule can fix

**abismal stores SEQ exactly as it came off the sequencer, on reverse records too**
(`README.md:108-109`). That is not what SAM says, and it is the invariant every rule in
`queries/strand/` relies on. Measured on the fixture, the reverse records differ from the
reference in 42–49 of 60 bases as stored and in 5–13 once reverse-complemented. dnmtools' own
caller compensates by reading SEQ backwards without complementing
(`methcounts.cpp:299-325`); every other reader has to repair the file first. The demo does it
in five lines of `samtools view | python | samtools sort`.

alnbase's count of read bases differing from the reference is the diagnostic, and here —
unlike in [`../hisat-3n-strand/`](../hisat-3n-strand/), where it could see nothing — it
separates all four cases:

| BAM | rule | bases differing from the reference |
|---|---|---|
| raw | shipped | 155 of 480 (32.3%) |
| raw | corrected | 159 of 480 (33.1%) |
| repaired | shipped | 38 of 480 (7.9%) |
| repaired | **corrected** | **0 of 480** |

Only the repaired file read with the corrected rule comes back clean, and the two failures
are distinguishable by size: a third of all bases is the SEQ convention, a twelfth is the
strand.

## What it costs in hits

The same four runs, counting CG hits, with the repaired-and-corrected run as the answer:

| BAM | rule | hits | on the wrong conversion strand | of the 56 right hits, found |
|---|---|---|---|---|
| raw | shipped | 28 | 0 | 28 |
| raw | corrected | 35 | 0 | 35 |
| repaired | shipped | 57 | **29** | 28 |
| repaired | corrected | 56 | 0 | **56** |

The shipped rule on a repaired BAM is the worst of the four to receive downstream: it emits
*more* hits than there are, half of them anchored on the other base of the CG, at coordinates
one base away and on the opposite strand, with nothing to distinguish them. The raw BAM
simply loses hits, because a reverse record's stored SEQ does not match the reference at all.

## `dnmtools format`, for the record

Still out of scope, and more so than the old comment said. It does not merely normalise the
tag: it merges the two mates of a fragment into one record whose CIGAR bridges the gap with
an `N`, drops the `0x1` paired bit, and writes `CV:A:T` on everything.

```
CTOB  80   851  60M80N60M  NM:i:0  CV:A:T
CTOT  128  351  60M80N60M  NM:i:0  CV:A:T
OB    144  601  60M80N60M  NM:i:0  CV:A:T
OT    64   101  60M80N60M  NM:i:0  CV:A:T
```

A formatted BAM is recognisable by having `CV:A:T` everywhere and no `0x1`.

## Caveats, and what is not exercised

- Single-end abismal is not run. `CV` there comes from `-A`/`-P`/`-R` rather than the mate
  number, and `abismal.cpp:480` is the assignment; the rule is the same XOR but unmeasured.
- The `-P` (plain PBAT) mode is run only far enough to see that it maps the two copy-strand
  pairs and nothing else.
- The BSMAP and Bismark translations in `standardize_format` are read but not run, because
  reaching them means running `dnmtools format`, whose output no rule covers.
- No gaps, trimming, repeats or low MAPQ. abismal writes MAPQ 255 on every record
  (`abismal.cpp` passes it as a literal), so nothing here measures a mapping-quality filter.
- `dnmtools counts`, dnmtools' own caller, is not run. This demo is about the tag.
