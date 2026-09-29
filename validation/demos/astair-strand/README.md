# asTair's table of six FLAG integers, mirrored — and an option that cannot work

[`queries/strand/astair.toml`](../../../queries/strand/astair.toml) is the tenth of the twelve
strand rules to be run against its tool's real output, and the only one that transcribes a
*caller's* table rather than an *aligner's* tag. The table came through right: mirroring it in
alnbase reproduces asTair's own `.mods` file position for position. Three of the claims in the
file's comment did not, and the run found a documented asTair option that silently produces
nothing.

```
ASTAIR_ENV=<env> PY=<python with pyarrow> ALNBASE=<alnbase> ./run.sh [workdir]
```

`ASTAIR_ENV` is the environment from [`../../envs/astair.yaml`](../../envs/astair.yaml), with
asTair installed into it from a checkout of <https://bitbucket.org/bsblabludwig/astair> — run
against master `b7df279` (3.3.3). **Bitbucket, not GitHub**: there is no asTair repository on
GitHub, and the rule's comment said there was. That master's last commit deprecates asTair in
favour of [`rastair`](https://bitbucket.org/bsblabludwig/rastair), which has its own strand
handling and is covered by no rule here.

## The chemistry is inverted, and so is everything downstream of it

asTair is a TAPS tool. TAPS reads a **modified** cytosine as T and leaves an unmodified one as
C — the opposite way round from bisulfite — and it does not convert the whole library, so the
reference is not converted either. Two consequences run through this demo:

- The aligner is an ordinary `bwa mem` against an ordinary reference, not a three-letter
  aligner against a converted one.
- There is no converted index, so nothing at alignment time filters a pair whose orientation
  does not match a directional library. That only matters if you feed asTair one, which it
  does not support; see below.

The fixture is the same 1200-base reference the other strand demos use, with a CG planted every
17 bases, and **every other planted CG modified**. Modifying all of them would leave no exact
stretch longer than 16 bases — under `bwa mem`'s default `-k 19` seed — and nothing would map
at all; alternating also gives the fixture a mixture of modified and unmodified sites instead of
a uniform modification level of 1.0.

## Mirroring the table: 261 positions, none disagreeing

`taps.toml` asks two questions at one anchor — a CG cytosine read as `T`, and read as `C` — so
alnbase's two hit counts line up against asTair's `MOD` and `UNMOD` columns directly, rather
than the comparison being a check on strand labels alone.

|  | positions | mod | unmod |
|---|---|---|---|
| asTair | 261 | 710 | 1954 |
| alnbase | 262 | 711 | 1954 |
| in common | 261 | | |
| **disagreeing** | **0** | | |

The one extra position is alnbase's: reference position 1, the bottom-strand cytosine of the CG
planted at offset 0. asTair drops it because its `SPECIFIC_CONTEXT` is a trinucleotide (`CGA`,
`CGC`, `CGG`, `CGT`) and that trinucleotide runs off the start of the contig; alnbase's CG query
needs only the dinucleotide. The reference holds 266 CG cytosines and asTair's `TOTAL_POSITIONS`
is 265, which is that one position and no other; the remaining gap to its 261
`COVERED_POSITIONS` is four positions past the end of the fixture's read coverage, which
alnbase does not report either.

## Three things the rule's comment got wrong

**It cited GitHub.** See above.

**It said the rule reports an out-of-table record.** It did not — alnbase *aborted*:

```
no rule in [strand.conversion] matches read 'ff' (FLAG 65)
```

alnbase refuses to run a strand rule that leaves a record uncovered; the escape hatch is an
explicit `unknown` arm, and the shipped file had none, so it was never total. With the arm
added, the two FF records are skipped and counted, which is what the comment always described:

```
read 410 records, scanned 408, skipped 2 (given a null row)
  of those, 2 were declined by this run's strand rule: no strand, so no walk and no row.
```

**It said asTair counts such a record toward depth.** It does not. The demo aligns the same
fixture twice, with and without the FF pair, and diffs the two `.mods` files: **identical**. Not
one position's `MOD`, `UNMOD` or `TOTAL_DEPTH` changes. pysam's pileup drops a non-proper pair
before asTair's own FLAG table ever sees it, so it never reaches depth either.

## The bug: `--ignore_orphans False` silently produces nothing

`-io/--ignore_orphans` is a documented CLI option (`caller.py:61`). Turning it off is supposed
to admit the non-proper-pair FLAGs `97/145/81/161` alongside the six. The branch that does so
cannot run:

```python
# caller.py:247-248 (and 259-260, for reference == 'G')
if ignore_orphans == False:
    expected_flags.extend(97, 145)
    unexpected_flags.extend(81, 161)
```

`list.extend()` takes exactly one argument. Called directly, both halves of the option are
visible at once:

```
ignore_orphans=True  -> 8 desired, 8 undesired tuples
ignore_orphans=False -> TypeError: list.extend() takes exactly one argument (2 given)
```

What makes this worth recording is not the typo but what happens to it. `clean_pileup` wraps its
whole per-position body in `except Exception: continue` (`caller.py:339`), so the `TypeError` is
raised and swallowed once per position, for every position. The run then:

- logs `asTair modification finder finished running`,
- writes a `.mods` file containing its header and nothing else,
- writes a `.stats` file reporting `COVERED_POSITIONS 0` and `*` for every modification rate,
- and **exits 0**.

A user who passes the option gets a successful-looking run and an empty result, with no
indication that the option is the reason.

## Directional-only, tested out of spec

asTair is directional-only and says so: `--library directional|reverse` flips the polarity
globally, there is no third option, and non-directional is described upstream as under
development. The fixture includes two pairs a directional library does not produce -- pairs
whose **read 1** came off the copy strand -- and asTair calls both mates of each on the wrong
conversion strand, which is the expected failure:

```
pair  flag  sequenced from  informs about  rule says
CTOB  99    CTOB            -              +   <- wrong
CTOB  147   OB              -              +   <- wrong
CTOT  83    CTOT            +              -   <- wrong
CTOT  163   OT              +              -   <- wrong
OB    83    OB              -              -
OB    163   CTOB            -              -
OT    99    OT              +              +
OT    147   CTOT            +              +
conversion strand wrong: 4 of 8
```

The one thing worth understanding here is what the six integers actually claim, because it is
easy to misread them as being about individual reads. They are about **fragment orientation**:
`99` and `147` together mean "an FR pair whose read 1 is forward", which is a top-strand
fragment *and both of its mates* -- the OT read 1 and the CTOT read 2. So CTOT and CTOB records
are handled correctly and constantly; read 2 of every ordinary fragment is one. What the table
cannot survive is a pair whose orientation is inverted relative to its fragment's strand, at
which point both mates are counted on the wrong side.

## Caveats, and what is not exercised

- `--library reverse`, asTair's nearest thing to PBAT, is not run. It exchanges the two
  conditions globally; the rule file says so and the demo does not test it.
- Single-end is not run, so the `0`/`16` half of the table is exercised only by
  [`queries/strand/records/astair.sam`](../../../queries/strand/records/astair.sam), not by real
  output.
- `astair align` and `astair phred`, `astair filter`, `astair mbias` and `astair simulate` are
  not run. This demo is about the caller's FLAG table.
- No `-ks` known-SNP file, no indels, no soft or hard clipping, no repeats, no multimappers.
  MAPQ is 60 on every probe record and `-mq` is left at its default of 0.
- `-sc/--skip_clip_overlap` is left at its default. Mate overlap is a separate question and
  alnbase has its own `overlap` command for it.
- `rastair`, the successor named in asTair's own deprecation notice, is not covered by any rule
  in `queries/strand/` and is not run here.
