# Declaring how to recover a read's original strand

Status: **implemented** as of 0.1.18, except where a section says otherwise. A query
file's `[strand.*]` tables decide every record's strand, and they are the only thing that
does: `--library` was retired in 0.1.18, and a run that reads records without a rule is
refused. `queries/strand/directional.toml` is the retired flag's rule written out, and
reproduces it exactly.

The cases this has to cover are in `research/strand-determination.md`: eleven aligners, six
patterns of evidence, three of them with a FLAG whose 0x10 does not mean what SAM says.

## What a rule has to decide

Per record, alnbase needs three things:

- **a**, the **conversion strand**: which reference strand carried the conversion, which is
  also the direction the walk runs and which reference base is informative
  (`src/alignment.rs:108`, the `strand` column's `+`/`-`).
- **b**, the **sequenced direction**: whether the read as sequenced runs along the reference
  forward or backward, which is what `off_5p`/`off_3p` count from. Today FLAG 0x10.
- **c**, the **strand of origin**, four-way: OT, CTOT, OB, CTOB. Wanted for `conv_strand`
  and for `[tag.XX.strand]` values.

**b and a both follow from c, by biology rather than by convention.** A read of origin OT or
CTOB is a copy of the top-strand sequence, so as sequenced it reads along the reference
forward; OB and CTOT are bottom-strand sequences and read backward. OT and CTOT are the
converted top strand and its complement, so the informative reference base is C; OB and CTOB
give G. Since every aligner surveyed stores SEQ in the reference's forward orientation, c is
enough to place and orient a record without consulting the FLAG at all:

| c | a (walk) | b (as sequenced) | SEQ vs. the read |
|---|---|---|---|
| OT | + | forward | as sequenced |
| CTOT | + | reverse | reverse-complemented |
| OB | − | reverse | reverse-complemented |
| CTOB | − | forward | as sequenced |

This is the whole reason the design can be small: **a rule that yields c yields everything**,
and it does so without trusting 0x10 — so declaring the rule *is* the fix for Bismark
single-end, BSBolt read 2 and BS-Seeker2, whose 0x10 means the conversion strand.

What c cannot always be is *known*. On bwa-meth, HISAT-3N, dnmtools or BISCUIT output the
evidence distinguishes only two of the four. Then the rule has to give **a** and **b**
separately, and say c is unknown rather than invent it.

## The shape

Three tables, in the same TOML the queries live in, named by kind in the table key like
`[tag.XM.bases]`. Keys are outcomes; values are conditions on the record.

```toml
[strand.origin]             # four-way; a and b are derived from it
OT   = "not is_reverse and not is_last_in_template"
CTOT = "is_reverse and is_last_in_template"
OB   = "is_reverse and not is_last_in_template"
CTOB = "not is_reverse and is_last_in_template"
```

or, when only two states are distinguishable:

```toml
[strand.conversion]         # a
"+" = 'YD == "f"'
"-" = 'YD == "r"'

[strand.sequenced]          # b
forward = "not is_reverse"
reverse = "is_reverse"
```

`[strand.origin]` and the other two are alternatives: declaring `origin` together with
`conversion` or `sequenced` is an error, because a derived value is a theorem and a second
declaration of it can only disagree. Declaring `conversion` without `sequenced` is also an
error — the whole point is that "b is 0x10" is a claim about the input, not a default.

**Exactly one key must match each record.** No match, or two matches, is an error that names
the record and the rule, not a silent fallback. Two escapes, both explicit:

- an `unknown` key, whose records alnbase counts and skips (BISCUIT's `YD:A:u`, reads where
  no conversion event was observed, belong here);
- in `[strand.origin]` only, records may be claimed by `"+"`/`"-"` keys instead of a strand
  name, meaning "conversion strand known, origin not" — the same as writing a `conversion`
  table for those records. Then `b` still has to come from somewhere, so a `[strand.sequenced]`
  table is required alongside and applies to exactly those records.

Static checks at load: every key in a table is distinct, `origin` uses only strand names,
`sequenced` only `forward`/`reverse`, and a pair of conditions that are syntactically the
same is rejected. Conditions that overlap in general cannot be checked statically — that is
what the per-record exactly-one check is for.

## The condition language

Boolean expressions with `and`, `or`, `not` and parentheses, exactly as `where` combines
patterns today (`docs/reference/02-query-language.md`), over two kinds of term:

- **record fields that are already named** — `is_reverse`, `is_first_in_template`,
  `is_last_in_template`, `is_paired`, `is_proper_pair`, `is_secondary`, `is_supplementary`
  (`src/record_field.rs`). A bare name is the boolean field; `flags == 99` compares the whole
  FLAG for the aligners whose tables are written that way (asTair's is).
- **aux tag comparisons** — `XG == "CT"`, `ZS == "+-"`, `YD == "f"`. A bare two-character
  name is the tag read as a string; `XG:A` and the other type codes are accepted with the
  same spelling `FieldSpec` already parses. A tag that is absent from the record makes every
  comparison on it false, so the record falls through to no match and is reported.
- **the literals `true` and `false`**, which name every record and none. `true` is how a rule
  says that a property holds whatever the record is: an assay that converts nothing walks
  every record the same way, and `unconverted.toml` writes that as `"+" = "true"` rather than
  as a tautology (`is_reverse or not is_reverse`) over a field it does not otherwise care
  about. `false` is its complement and exists so the pair is symmetric; a key that can never
  match is a key that should be deleted, and the records test will fail a shipped rule that
  has one, since it requires every declared key to be exercised.

Nothing else: no arithmetic, no regular expressions, no reference to the read's bases. Every
rule below is expressible with equality alone, and a rule that needed more would be a rule
whose input should have been fixed upstream.

## One mechanism: files, not a menu

There is no built-in menu and nothing is compiled into the binary. alnbase ships a query
file per aligner in `queries/strand/` — `bismark.toml`, `bsbolt.toml`, `directional.toml`
and the rest below — and they are read from disk like any other query file. Supporting a
new aligner is writing a file, which anyone can do without a new release, and adapting one
is forking it.

That makes the file the whole interface: a user can read the rule alnbase is applying,
diff it against the aligner's documentation, and change it. It also means there is no rule
alnbase can express that a user cannot write, and no path through the code that a shipped
file does not take.

A strand rule is just tables in a query file, so it is passed with the same repeatable
option the queries use; two `[strand.*]` declarations in one run is an error, like two
declarations of the same tag. `--library directional` goes away, replaced by
`queries/strand/directional.toml`, which keeps today's behaviour byte for byte.

## The files to ship

These are now in the repository, one file per aligner under `queries/strand/`, with
`queries/strand/README.md` as the index; each names the aligner version its rule was read
from and what it would get wrong on another input. Every one is exactly what the survey
found, and the `verify` marks are claims no validator has run yet. Abridged here, without
those headers:

```toml
# queries/strand/directional.toml — bwa-meth, BISCUIT -b 1, HISAT-3N --directional-mapping,
# methylpy output, Bismark paired-end, BSMAP -n 0 single-end. The current rule.
# An unpaired read has is_last_in_template false, so it is read 1, as today.
[strand.origin]
OT   = "not is_reverse and not is_last_in_template"
CTOT = "is_reverse and is_last_in_template"
OB   = "is_reverse and not is_last_in_template"
CTOB = "not is_reverse and is_last_in_template"

# queries/strand/unconverted.toml — an assay that converts nothing. Every record walks the
# reference forward; 0x10 says only which way the read was sequenced.
[strand.conversion]
"+" = "true"

[strand.sequenced]
forward = "not is_reverse"
reverse = "is_reverse"

# queries/strand/bsmap.toml — ZS is the conversion strand then the read's orientation within it.
# Works for -n 1 (non-directional) unchanged.
[strand.origin]
OT   = 'ZS == "++"'
CTOT = 'ZS == "+-"'
OB   = 'ZS == "-+"'
CTOB = 'ZS == "--"'

# queries/strand/bsbolt.toml — YS is the contig (Watson/Crick) then the read's conversion. Four-way from
# the tag alone, which is what makes the FLAG's inverted 0x10 stop mattering: b comes
# from the derived table, not from the bit. Works for -UN unchanged.
[strand.origin]
OT   = 'YS == "W_C2T"'
CTOT = 'YS == "W_G2A"'
OB   = 'YS == "C_C2T"'
CTOB = 'YS == "C_G2A"'
unknown = 'YS == "WC"'          # unmapped records

# queries/strand/bismark.toml — XR is the read's own conversion, XG the genome index that won. Holds for
# single-end, where 0x10 means XG and must not be read, and for directional paired-end.
[strand.origin]
OT   = 'XR == "CT" and XG == "CT"'
CTOT = 'XR == "GA" and XG == "CT"'
OB   = 'XR == "CT" and XG == "GA"'
CTOB = 'XR == "GA" and XG == "GA"'

# queries/strand/bs-seeker2-se.toml — XO is the genome strand then the read's direction.  (verify)
[strand.origin]
OT   = 'XO == "+FW"'
CTOT = 'XO == "+RC"'
OB   = 'XO == "-FW"'
CTOB = 'XO == "-RC"'

# queries/strand/merantk.toml — YG is the genome index, YR the read's conversion, as Bismark
# splits XG and XR. Directional only: wrong-direction alignments are dropped, so only OT and
# OB occur, and the other two are declared so such a record is named rather than unmatched.
[strand.origin]
OT   = 'YG == "C2T" and YR == "C2T"'
CTOT = 'YG == "C2T" and YR == "G2A"'
OB   = 'YG == "G2A" and YR == "C2T"'
CTOB = 'YG == "G2A" and YR == "G2A"'

# queries/strand/{biscuit,bwameth}.toml — two-way tag, conformant FLAG, no four-way answer.
[strand.conversion]
"+" = 'YD == "f"'
"-" = 'YD == "r"'
unknown = 'YD == "u"'           # BISCUIT: no conversion event seen
[strand.sequenced]
forward = "not is_reverse"
reverse = "is_reverse"

# queries/strand/hisat-3n.toml — the same shape on a differently named tag: YZ is the 3N
# index hit. In the default non-directional mode YZ is inferred from conversion counts, so
# the file reports an assumption as fact on a read with none.
[strand.conversion]
"+" = 'YZ == "+"'
"-" = 'YZ == "-"'
[strand.sequenced]
forward = "not is_reverse"
reverse = "is_reverse"

# queries/strand/dnmtools.toml — CV is the read's own conversion, before `dnmtools format` normalises it.
[strand.conversion]
"+" = 'CV == "T"'
"-" = 'CV == "A"'
[strand.sequenced]
forward = "not is_reverse"
reverse = "is_reverse"

# queries/strand/astair.toml — no tag at all; the caller's own table is whole FLAG values.
[strand.conversion]
"+" = "flags == 99 or flags == 147 or flags == 0"
"-" = "flags == 83 or flags == 163 or flags == 16"
[strand.sequenced]
forward = "not is_reverse"
reverse = "is_reverse"
```

Two that do not fit, and should not be made to:

- **Bismark `--non_directional` paired-end.** The default output swaps 0x40 and 0x80 on CTOT
  and CTOB pairs, so a rule reading the mate bit is reading a rewritten one. `XR` and `XG`
  are per-record and appear to survive the swap, which would make `bismark.toml` above
  correct there too — this is the single claim in this document most worth checking against a
  real `--non_directional` BAM before that file claims to cover it. `--strandID`'s `YS:Z`
  gives the four-way directly and is the rule to prefer when it is present.
- **dnmtools after `format`.** It reverse-complements SEQ so its own tools can read reverse
  records backward without complementing, breaking the reference-forward invariant every
  other input holds. No strand rule can repair that; it is a different SEQ convention and
  would need its own declaration. Out of scope here, named so it is not mistaken for an
  oversight.

## What changes in alnbase

- `Library` (`src/tags.rs:125`) is replaced by a `StrandRule` holding the compiled tables;
  `Library::strand(record)` becomes `rule.call(record) -> StrandCall`, where a call carries
  a, b and an optional c.
- `src/alignment.rs:108` takes a from the call instead of `is_last_in_template() ==
  is_reverse()`, and the offsets take b instead of `is_reverse()`. These are the two places
  that read the FLAG for strand today; nothing else does.
- `conv_strand` (`src/record_columns.rs:277`) becomes null when c is unknown, and its
  description stops saying "for a directional library".
- A new `read_reverse` column: b, "the read as sequenced runs opposite the reference". When
  0x10 lies, this is the column that tells the user so, and it is the one `off_5p`/`off_3p`
  are counted with.
- `[tag.XX.strand]`'s requirement that all four strands be named stays — a value is still
  needed for whichever strands occur — but writing such a tag from an input whose c is
  unknown has to be an error at load time, not a wrong tag at run time. **Done** (0.1.19):
  `tags::strand_tags_need_origin` compares the run's declared tags against
  `StrandRule::names_origin()` while the query files are being read, so the pair is refused
  before the BAM is opened. `bam_out` keeps its per-record refusal as a backstop that should
  never speak.

No run-time checking of the declared rule against the FLAG. A counter for records whose
declared sequenced direction disagrees with 0x10 was proposed and rejected: there are too
many ways for an input to be wrong to chase each with its own detector, and reading the
aligner's source to find and document the bug — as was done for BSBolt in
`alnbase-validation/demos/premethyst-bugs/README.md` — is the productive form of the same work. The
shipped file per aligner is where that knowledge lands.

## Open questions

1. ~~The manifest should record which file the rule came from rather than `library`.~~
   **Done** (0.1.17): the manifest's `strand` key names the query file, whose whole text is
   already in `query_files` (0.1.18: always a file, since there is nothing else it can
   be). The underlying hazard stands — a
   strand rule is a property of the BAM, not of the queries, so a file carrying both can be
   pointed at the wrong input and will still run — but a hit table now says how its strands
   were called.
2. ~~Whether `unknown` records should be skipped or emitted with null coordinates.~~
   **Settled: skipped and counted** (0.1.17). The count is reported at the end of the run,
   so the loss is not silent, and it is what tells a user their rule does not cover their
   input — which they fix by widening the rule rather than by a flag alnbase provides. A
   tagged BAM still writes such a record through, untagged, as it does one the walk cannot
   walk; a hit table gives it no row, because a null row would claim the record was read
   and matched nothing.
3. How the shipped files stay honest as aligners change. Each names the aligner version and
   the source lines its rule was read from, so a file that was right for Bismark 0.25.1 says
   so. **Half done** (0.1.20): every file now ships `queries/strand/records/<name>.sam`,
   holding one record per strand the rule names, each record named for the call the rule
   must make, and a test runs all twelve on every build. That catches a rule which stops
   agreeing with its own documentation -- a swapped pair of conditions, a strand declared
   and never exercised -- and it is enforced both ways, so a rule cannot ship without
   records and records cannot outlive their rule. It cannot catch the aligner changing what
   it writes, because the records were written from the same reading of the source the rule
   was. That still needs the aligner's own BAM run in a validation demo, which is why
   `bs-seeker2-se.toml` stays marked `verify`.
4. Non-directional inputs with no tag — Bismark `--non_directional` under `--old_flag`,
   HISAT-3N's default mode when the conversion counts tie — have no evidence to declare over.
   The rule for them is honest only if it is `unknown`, which means the input is unusable;
   worth saying so in the docs rather than shipping a file that guesses.
