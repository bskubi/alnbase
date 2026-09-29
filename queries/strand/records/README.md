# Records for the strand rules

One SAM per rule in the directory above, holding a record for every strand that rule
names. `bismark.sam` is checked against `bismark.toml`, and so on.

**Each QNAME is the call the rule must make for that record.** That is the whole
convention: the file needs no separate table of expectations, and a record and its
expected answer cannot drift apart.

```
OT/pe-99      99  chr1  101  60  4M  =  301  200  ACGT  IIII  XR:Z:CT  XG:Z:CT
CTOT/pe-147  147  chr1  101  60  4M  =  301  -200 ACGT  IIII  XR:Z:GA  XG:Z:CT
```

A QNAME is one of:

| QNAME | What the rule must answer |
|---|---|
| `OT`, `CTOT`, `OB`, `CTOB` | that strand of origin, from which the conversion strand and the sequenced direction follow |
| `<conversion>:<sequenced>` — `+:forward`, `-:reverse`, … | that conversion strand and that sequenced direction, with the origin unknown. This is the shape for a rule built on `[strand.conversion]` |
| `unknown` | nothing: the rule declines the record, which alnbase skips and counts |

Anything after the first `/` is a label for the reader — the FLAG the record carries, the
tag value it was written for — and is ignored. It is what lets two records expect the same
call without sharing a name.

## What the test enforces

`strand_rule::tests::every_shipped_rule_calls_its_records_as_their_names_say` reads each
`*.toml` above with the ordinary query-file parser, runs every record of the matching
`*.sam` through it, and requires three things:

- every record is called exactly as its QNAME says;
- every key the rule declares is exercised by at least one record, so a file cannot claim
  a strand it is never asked about;
- every rule has a records file, and every records file has a rule
  (`no_records_file_outlives_its_rule`).

## What it does not enforce

These records are written from each aligner's source — the same reading the rule itself
came from — so they catch a rule that stops agreeing with its own documentation: an edit
that swaps two conditions, a condition that stops parsing the way it used to, a strand
declared and then forgotten. They cannot catch the aligner changing what it writes, because
they are not that aligner's output. Only running a real BAM through the rule can, which is
what the validation demos are for. `alnbase-validation/demos/bs-seeker2-strand/` is the first of
those, and it found something these records had wrong: they carried the FLAGs a conformant
aligner would write on the two `RC` classes rather than the inverted ones BS-Seeker2
writes. No call changed, because that rule reads only `XO` -- but a records file exists to
say what the aligner writes, so it was wrong.

## Adding a rule

Ship a `.sam` beside the `.toml` or the test fails. Copy the nearest existing file, change
the FLAGs and tags to what the aligner writes, and name each record for the answer you
expect. Keep the header as it is: `chr1` and the 4M alignment are there so the records are
valid SAM, not because anything reads them.
