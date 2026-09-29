# Changelog

## 0.1.27 — 2026-09-29

- **Licence and repository declared** (`Cargo.toml`). The crate is MIT-licensed, matching
  `LICENSE` at the repository root, and `repository` points at
  `https://github.com/bskubi/alnbase`. Both fields were commented out as undecided.

## 0.1.26 — 2026-09-21

Documentation aimed at rustdoc, starting with the page a reader lands on.

- **A crate-level doc comment** (`src/main.rs`). The crate root had no `//!` at all, so
  `cargo doc` produced a landing page with a module list and nothing to orient by. It now
  explains the three ideas the rest of the code assumes — a base is a set of possibilities
  rather than a letter, the column is the unit of work, a query is a boolean predicate over
  patterns rather than a pattern — then follows one record from `bam_io` through the strand
  rule, the walk, the automaton and the ring to an output, names the two scan pipelines and
  why there are two, gives the entry point of every subcommand, maps the modules by the job
  they do, and suggests a reading order. The audience is someone picking the codebase up
  cold, including its author after a year away.
- **`cargo docs`** (`.cargo/config.toml`, new). alnbase is a binary crate, so every module
  is private and a plain `cargo doc` renders almost nothing; the failure looks like broken
  documentation rather than a missing `--document-private-items`. The alias supplies the
  flag so the working command is the short one.
- **Broken and misleading intra-doc links fixed** (`src/bam_out.rs`, `src/parallel.rs`,
  `src/scanner.rs`). Three links did not resolve, and `cargo docs` is now warning-free.
  Two were pointing at `#[cfg(test)]` items that do not exist in a doc build: `bam_out`
  sent the reader to `tags::Library` for how a `strand` tag is chosen, which is the
  test-only oracle rather than the user's declared rule, so it now points at `strand_rule`;
  `parallel::scan_sharded_reader` described itself in terms of a test-only wrapper, and
  both functions now say which one production calls and why it takes an open reader.
  `Scanner::strand` carried a stray unresolvable link and now records what `Ok(None)` means.
- **Descriptions calibrated against the code they describe** (`src/main.rs`, `src/column.rs`).
  `PAD` was described as the reference past the end of a read, which is only its read-side
  meaning; it sits on whichever side has run out, and on the reference side it means a
  coordinate off the end of the contig. `Column` was described as a read base, a reference
  base and a coordinate, omitting the read offset and the base quality — and the offset it
  carries indexes `SEQ`, which is not the `off_5p` the outputs report, since `hits` adds the
  hard clip back on. `Column::refr` likewise documented only `GAP`, not the `PAD` and `SKIP`
  it can also hold. The walk does not emit one column per aligned position either: an
  inserted base is its own column and a long intron's interior collapses to one marker.
- **The same calibration applied along the query path** (`src/byg_search.rs`, `src/query.rs`,
  `src/dsl.rs`, `src/predicate.rs`, `src/column.rs`). `byg_search` claimed a step was "three
  operations", which omits the carry between words and the read-off of completed traces.
  `Column::read_off` was described as the offset "as sequenced, indexing SEQ", which are two
  different numbers — SEQ omits hard-clipped bases — and the field now says so and points at
  `hits`. `dsl` described a pattern's two sides as "equal-length strings" when the rule is
  equal *columns* after repeat counts expand, so `~8@GATCGATC` is well formed; it also said
  the anchor and captures come from the `mark` row without mentioning that the grammar still
  parses them inside a pattern in order to refuse them by name. `query` said the search and
  predicate layers meet at "a slice of `usize` and nothing else" (they also share the bitmap
  width), and did not mention that compiling is where a dead column or an always-satisfied
  query is refused. `predicate`'s popcount guards had no stated origin, and the count
  expression that produces them turns out to be reachable only from this module's tests,
  since no query-file syntax builds one — now said outright, so a reader does not go hunting
  for the syntax or delete the code as dead.
- **The driver and the query CLI** (`src/scanner.rs`, `src/cli_query.rs`). `scanner` stated
  an ordering rule the code does not have — that `ring.store` must precede `state.step` —
  when the two are independent and what matters is that both precede evaluation, the ring's
  because a query ending on a column reads its anchor back out of it. The header now says
  that, adds the span gate that stops a query firing before the ring holds its whole width,
  names the `--only-hits` flag rather than the internal field it sets, and records that a
  hitless row is also written for a record the walk never touches, which is what makes the
  file a denominator. `cli_query` claimed to carry "no tables of its own" while defining
  four enums; two mirror `parquet` crate types that cannot derive clap's `ValueEnum` from
  outside, one shadows an alnbase type to keep clap out of the writer, and one is a CLI
  convenience shadowing nothing — the header now distinguishes them and says why all the
  parquet-shaped ones resolve in one function. It also notes that query-file text is kept
  so a run can copy it into the BAM header.
- **A misattached doc comment** (`src/cli_query.rs`). `parquet_opts`'s three-sentence
  description had been glued onto `batch_rows_per_file`, whose own sentence trailed it,
  leaving `parquet_opts` undocumented and rendering the batch-size helper in rustdoc as a
  description of compression settings. Split, and `batch_rows_per_file` now states the
  convention it implements: `--batch-rows` is a total across output files, so raising the
  shard count does not quietly raise the run's memory.
- **Simplified Technical English throughout the doc comments.** The documentation written
  above had been written as literary prose: long sentences carrying several facts at once,
  metaphor, and asides set off with dashes. It is now rewritten in the controlled style used
  for technical writing — one idea per sentence, short sentences, active voice, literal
  wording — which reads more easily on a return visit and for a reader whose first language
  is not English. The facts and the level of detail are unchanged. Converted so far: every
  module header (`main`, `scanner`, `cli_query`, `bam_out`, `parallel`, `dsl`, `predicate`,
  `query`, `byg_search`, `column`), and the item-level docs in `main`, `scanner`, `query`,
  `byg_search` and `column`.
- **A second misattached doc comment** (`src/query.rs`). `QuerySet::compat_idx`'s
  description and its `#[inline]` attribute had come loose and were sitting on
  `QuerySet::shapes`, ahead of two different copies of `shapes`'s own description. So
  `compat_idx` was undocumented and not inlined, and rustdoc rendered `shapes` as a
  description of the compat table followed by the same paragraph twice. Each doc is now on
  its own item and `compat_idx` is `#[inline]` again.
- **A module header for the walk** (`src/alignment.rs`). The densest module in the program,
  and the one the crate root names as the place correctness is most fragile, had no `//!` at
  all — so it was the one blank row in the module table on the landing page. It now explains
  emitted order (the walk runs 5'→3' along the conversion strand, which is why a pattern
  written once works on both strands, and why a bottom-strand record walks the CIGAR
  backwards with both sides reverse complemented), gives a table of what each CIGAR operation
  produces, and covers the flank, the collapsed intron and its marker, the two cursors and
  the `read_off` invariant, and reading the reference through the contig map.
- **Every module's summary sentence rewritten** (37 of them, plus the crate root's). The
  module table on the landing page shows each module's first sentence and nothing else, and
  most of those were noun-phrase
  fragments: "Sharded parallel scan.", "Query files in TOML.", "The one place a record column
  is defined." They are now Simplified Technical English sentences that say what the module
  does — "Run the scan in parallel, and give each worker its own output files.", "Read a query
  file, which is written in TOML." Where the informative half of a summary had sat in a second
  sentence that the table never shows, it moved into the first.
- **The crate root's body converted too** (`src/main.rs`, `src/trace.rs`). The first pass
  converted the summary sentences and left the body of the landing page in the literary
  register it was drafted in, which the rendered page made obvious: "This is the spine of the
  program", a section headed "Three ideas to keep in mind", "The meaning of the pattern is the
  business of the user". Structural metaphor reads as ordinary technical writing and is still
  metaphor. The body is now Simplified Technical English throughout, including the section
  headings and the sentences that introduce a list.
- **A query that *fires*, defined rather than assumed.** The verb appears 285 times across the
  source, the help text and the outputs, so renaming it would be worse than keeping it — but on
  a landing page that never defines it, it reads as metaphor. The crate summary now says the
  program reports "every position where a query matches", and a sentence beside the first
  example defines firing as matching at a position.
- **The module pages below the crate root** (`src/seq.rs`, `src/alignment.rs`,
  `src/byg_search.rs`, `src/scanner.rs`). The same partial conversion that left the crate
  root's body untouched had left the module bodies untouched too: the summary sentence of
  each module was Simplified Technical English, and the page underneath it was not. `seq`
  and `alignment` are the two the crate's own reading order names first, and they carried
  the most of it — "two jobs pull in opposite directions", "a set that is nearly what was
  meant", "no symbol would say so honestly", "at the cost of breaking patterns that span
  one". Both are now converted throughout, header and rendered item docs, with every fact,
  link and table unchanged. `byg_search` said that interning a shared pattern "makes it
  free"; it now says the second use adds no cost, which is what the crate root already
  said. Each page was checked as rendered HTML rather than as source.
- **The origin of the matching algorithm, on the page where readers meet it**
  (`src/main.rs`). `byg_search` cited Baeza-Yates and Gonnet's "A New Approach to Text
  Searching" (CACM 1992) in its own header, but the crate root described the automaton
  without saying where it came from, and the crate root is the page a reader lands on.
- **"Interns" replaced by what actually happens** (`src/main.rs`, `src/byg_search.rs`,
  `src/query.rs`, `src/dsl.rs`, `src/trace.rs`, `src/cli_query.rs`). Interning is an
  analogy, and not one that clarifies: the reader has to know the analogy before the
  sentence says anything. All six sites now state the fact, which is that the run keeps
  one copy of a pattern two queries share.
- **A sentence that dropped a fact while getting shorter** (`src/main.rs`). Step 7 of the
  walkthrough read "`bam_out` writes tags onto the record, or `hits` builds a row" —
  builds a row where, and under what option. It now names `--parquet` as the path that
  reaches `hits`, and says the row is written into a parquet file.
- **Sections edited for what point they make** (`src/main.rs`). A pass at the level of the
  section rather than the sentence, cutting text that did not serve the point of the chunk
  it was in, even where the text was accurate. The opening carried a fenced code block for
  a two-word command and a paragraph justifying `--document-private-items`, where the
  point of the chunk is who the page is for; seventeen lines became ten. The flag table's
  explanation of `PAD` was then repeated in prose underneath it. The walkthrough's
  description of the automaton gave the instruction sequence, which `byg_search` already
  states, where the overview needs only the cost property. The note on `scanner` had grown
  a paragraph about why the list above it is ordered as it is, which is commentary on the
  page and not on the program.
- **`encoding` was filed under the walk** (`src/main.rs`). The module map listed it among
  `seq`, `column`, `alignment` and the rest of the walk. It describes how an output column
  is stored, so it belongs under Output.
- **Contractions, rhetorical questions and filler removed** (`src/trace.rs`, `src/dsl.rs`,
  `src/aref.rs`, `src/hits.rs`, `src/partition.rs`, `src/contig_map.rs`,
  `src/strand_rule.rs`, `src/predicate.rs`, `src/tags.rs`, `src/arrow_builder.rs`,
  `src/overlap.rs`, `src/overlap_apply.rs`, `src/record_columns.rs`, `src/explain.rs`,
  `src/cli_query.rs`, `src/batch.rs`, `src/hit_writer.rs`). Doc comments written as
  questions ("Is `c` usable as an alias?"), five uses of "simply" that changed no
  sentence's meaning, "worth knowing", "worth stating outright", "worth the trouble",
  "emphatically not", "blows up", "a nonsense overlap", "buys that back", "part company",
  "the cheapest way", "the main way this class of tool goes wrong", and `cli_query`'s "Two
  points on a curve" — the metaphor that prompted this whole style change, still in the
  source until now.
- **The section pass continued into the module pages** (`src/hits.rs`, `src/partition.rs`,
  `src/overlap.rs`, `src/encoding.rs`). `hits` opened by promising three properties and
  spent its third bullet forwarding to a section heading directly below it, which now
  carries itself. Its note on the two read offsets argued for four lines that deriving one
  from the other is error-prone; the point of the section is what the offsets mean. A
  closing aside about parquet dictionary-encoding empty lists belonged to `encoding` and
  did not follow from the sentence before it. `partition` and `overlap` lost their
  question-shaped headings.
- **Two sentences that said the wrong thing** (`src/overlap.rs`, `src/encoding.rs`).
  `overlap` introduces a table of *three* coordinate systems and then said a tool of this
  class fails by confusing "the two" — damage from an earlier rewrite in this same series,
  which replaced a correct "them". `encoding` said that an abandoned dictionary means "the
  bytes are saved twice over", where what happens is that the file pays for them twice,
  once in the dead dictionary and once in the plain values after it.
- **A wrong attribution caught by rustdoc.** A sentence drafted for the new `alignment` header
  said `cli_query` sizes the walk's flank. It is `scanner::WalkConfig::to_opts`; `cli_query`
  only carries the `--end-context` argument. Writing the claim as an intra-doc link rather
  than as plain prose is what surfaced it, since the link did not resolve.

## 0.1.25 — 2026-09-18

The condition language gains the literals `true` and `false`, and `queries/strand/` gains a
rule for assays that convert nothing and loses the one for PBAT libraries.

- **`true` and `false` in a strand rule's conditions** (`src/strand_rule.rs`). `true` names
  every record and `false` none. The need came from the file below: an assay with no
  conversion walks every record the same way, and saying so previously required a tautology
  (`is_reverse or not is_reverse`) over a field the rule does not otherwise care about, which
  is exactly the sort of line a user copies and then mangles. `true == 1` is rejected with a
  message saying the literal is a condition already. Two keys given `true` are caught by the
  existing identical-condition check at load.
- **`queries/strand/unconverted.toml`**, with its records file. alnbase is a tool for patterns
  over aligned read and reference columns, and a strand rule is mandatory for every run that
  reads records, so an assay that converts nothing — variants, motifs, RNA modifications read
  as mismatches, damage — still needs a rule and until now had none that fitted. This one
  walks every record along the reference forward, so a pattern written against the reference
  means the same thing on every record, reports the strand of origin as unknown because none
  exists, and keeps the read's own orientation in `read_reverse` so `off_5p` and `off_3p`
  still count from the sequencer's 5' end. Reaching for `directional.toml` instead makes the
  walk follow each fragment's orientation and reverse-complements roughly half the data.
- **`queries/strand/pbat.toml` is deleted**, along with its records file, on the evidence of
  `alnbase-validation/demos/pbat-strand/`. The rule was a correct reading of a conformant BAM from a
  PBAT library, but no such BAM is produced: both pipelines built for PBAT put the reads back
  into the directional arrangement first, Bismark `--pbat` by swapping 0x40 and 0x80 on every
  pair (under `--pbat` every pair is CTOT or CTOB) and methylpy `--pbat` by exchanging the
  mate files or reverse-complementing every read. On real Bismark `--pbat` output the rule
  named the opposite conversion strand on 402 of 402 paired-end records, putting 17.09% of
  read bases at odds with the reference against 0.00% under `directional.toml`. A file that is
  right in principle and wrong on every BAM anyone has is a trap. `bismark.toml` covers both
  `--pbat` layouts and its comment now says so; `directional.toml` covers the paired-end arm
  and its comment now says where it stops.

## 0.1.24 — 2026-09-17

`queries/strand/biscuit.toml` is the first rule that cannot name a read's strand of origin:
BISCUIT's `YD` distinguishes two states, so the file declares the conversion strand from the
tag and the sequenced direction from FLAG 0x10. That second half is a claim about the
aligner, not a default, and it is measured now.

- **`alnbase-validation/demos/biscuit-strand/`** builds four read pairs, one per strand of origin and
  each named for the strand its read 1 was sequenced from, and aligns them with BISCUIT 1.10.3
  twice: two pairs with `-b 1`, and all four with the default non-directional `-b 0`. Every
  record's conversion strand and sequenced direction is the one its strand of origin implies
  (0 of 4 wrong, then 0 of 8), and FLAG 0x10 agrees with how SEQ was really stored on all
  twelve records, measured by comparing SEQ against the fastq reads rather than by trusting
  the bit.
- The demo checks the walk rather than the label: a CG query anchors on a read cytosine, and
  the reference base under each of the 29 and 56 anchors is the informative one for the
  strand alnbase walked -- the C of a forward CG for `+`, the G of that CG for `-` -- with
  none wrong. `conv_strand` is null on every row and no record was given a strand of origin,
  so alnbase declines the split the input cannot make instead of reporting a guess as a
  measurement.
- **`queries/strand/biscuit.toml`**'s comment said the strand of origin is out of reach
  because "OT and CTOB both align to the f contig". That is the grouping by sequenced
  direction. The demo counted it: `YD:A:f` held OT and CTOT, `YD:A:r` held OB and CTOB -- a
  strand and its own PCR copy hit the same converted contig, which is the real reason a
  two-way tag cannot reach the four-way answer. Corrected, here and in
  `queries/strand/README.md`.
- The same slip in `alnbase query --help`, which offered `CT = ["OT", "CTOB"]` as a
  `[tag.XX.strand]` table "by conversion strand", is corrected to `CT = ["OT", "CTOT"]`,
  matching the reference docs and `queries/strand/bismark.toml`. That was the last surviving
  copy of the inversion `alnbase-validation/docs/research/strand-determination.md` found in the Bismark
  example; the example itself was fixed when the survey reported it.
- FLAG 0x2 is absent from every record in the demo, and the README says explicitly that this
  is not a BISCUIT behaviour: its bwa-derived aligner had too few pairs to estimate an insert
  size (`dir.log`). `queries/strand/records/biscuit.sam` keeps its conformant FLAGs, which
  are otherwise the ones BISCUIT wrote.

## 0.1.23 — 2026-09-17

`queries/strand/bsbolt.toml` exists because BSBolt's FLAG 0x10 is wrong. That had only ever
been read out of BSBolt's source; it is measured now.

- **`alnbase-validation/demos/bsbolt-strand/`** builds four read pairs, one per strand of origin and
  each named for the strand its read 1 was sequenced from, and aligns them with BSBolt 1.6.0
  twice: two pairs as a directional library, which already exercises all four `YS` values,
  and all four under `-UN`. Every record comes back as the strand it was sequenced from, and
  alnbase's scan line (0 of 240 and 0 of 480 compared bases differ from the reference) says
  they were placed and walked correctly rather than merely labelled consistently.
- The demo counts the inverted bit instead of asserting it. BSBolt sets 0x10 from the
  converted contig the read matched rather than from the read's orientation, and the run
  finds 0x10 contradicting how SEQ was really stored on every `*_G2A` record and no other —
  2 of 2 directional, 4 of 4 under `-UN` — while SEQ stays reference-forward throughout,
  which is what lets a `YS`-only rule place and orient every record without the bit.
- **`queries/strand/records/bsbolt.sam`** had shipped with conformant FLAGs (99/147/83/163).
  It now carries the ones BSBolt itself writes (65/129/113/177), where 0x10 and 0x20 both
  follow the Crick contig rather than orientation.
- **`alnbase-validation/demos/premethyst-bugs/`**' note that BSBolt paired-end BAMs need their
  read-2 FLAGs repaired before alnbase reads them is superseded: with `bsbolt.toml` they do
  not.

## 0.1.22 — 2026-09-17

`--non_directional` paired-end was the one layout `queries/strand/bismark.toml` declined to
claim. It is claimed now, on the strength of a Bismark BAM rather than a reading of Bismark.

- **`alnbase-validation/demos/bismark-nondirectional-pe/`** builds four read pairs, one per strand
  of origin and each named for the strand its read 1 was sequenced from, and aligns them
  with Bismark 0.25.1 under `--non_directional`. All eight records come back as the strand
  they were sequenced from, and alnbase's scan line (0 of 480 compared bases differ from
  the reference) says they were placed and walked correctly rather than merely labelled
  consistently.
- The demo also measures the hazard the rule exists to step around: on the CTOT and CTOB
  pairs, Bismark swaps 0x40 and 0x80 so that browsers do not call them discordant, so four
  of the eight records carry a mate bit naming the other mate. `XR` and `XG` stay with the
  record they describe, so the rule is unaffected; a rule that read 0x40 to decide which
  strand a mate came from would have those four backwards.
- **`queries/strand/bismark.toml`** no longer carries the "NOT yet claimed for
  `--non_directional` paired-end" caveat, and records what the run showed.

## 0.1.21 — 2026-09-17

The first shipped strand rule checked against its own aligner rather than against the
reading of that aligner's source it came from.

- **`alnbase-validation/demos/bs-seeker2-strand/`** aligns four reads — one per strand of origin,
  each named for the strand it was built from — with BS-Seeker2 2.1.8 and runs the BAM
  through `queries/strand/bs-seeker2-se.toml`. Every read comes back as the strand it came
  from, and alnbase's own scan line (0 of 240 compared bases differ from the reference)
  says the records were placed and walked in the right direction rather than merely
  labelled consistently. `../../envs/bs-seeker2.yaml` and `build_bsseeker2.sh` build the
  environment: BS-Seeker2 is Python 2 and writes its BAM through pysam's pre-0.9 API.
- **`bs-seeker2-se.toml` is no longer marked `verify`**, the last file that was. The run
  also confirms by measurement what its comment claims from the source: FLAG 0x10 follows
  the genome strand, so `+RC` carries 0 where a conformant aligner writes 16 and `-RC`
  carries 16 where it writes 0, while SEQ is reference-forward in all four cases.
- **`queries/strand/records/bs-seeker2-se.sam` corrected.** It had shipped with the FLAGs a
  conformant aligner writes on the two `RC` classes rather than the ones BS-Seeker2 writes.
  No call changed, since that rule reads only `XO` — but a records file exists to say what
  the aligner writes, and this is exactly the kind of error only the aligner can catch.

## 0.1.20 — 2026-09-17

Every shipped strand rule is now checked against records rather than only
described in prose, so a file that stops being true fails the build instead of
going on producing calls on the wrong strand.

- **`queries/strand/records/<name>.sam`, one per rule.** Each holds a record for
  every strand its rule names, and each record's QNAME *is* the call the rule
  must make — `OT`, `CTOT`, `OB`, `CTOB`, `unknown`, or `+:forward` and friends
  for a rule that reaches only the conversion strand. There is no separate table
  of expectations to drift out of step. `bismark.sam` is that file's own XR/XG
  table made executable, in both the paired-end and the single-end layout.
- **The test runs all twelve on every build**, and enforces three things: every
  record is called as its name says, every key a rule declares is exercised by
  some record, and rules and records exist one for one in both directions. So a
  rule cannot claim a strand nothing asks it about, cannot ship without records,
  and its records cannot outlive it.
- **`StrandRule::declared_keys()`** is what the coverage half reads: the keys of
  `[strand.origin]` or `[strand.conversion]`, then `[strand.sequenced]`'s, in
  declaration order.

These records are written from each aligner's source — the same reading the rules
themselves came from — so they catch a rule that stops agreeing with its own
documentation, not an aligner that changes what it writes. That still needs the
aligner's own BAM, which is why `bs-seeker2-se.toml` stays marked `verify`.

## 0.1.19 — 2026-09-17

A `strand` tag and a strand rule that cannot name a strand of origin are an
impossible pair, and now they are refused when the query files are read rather
than on the scan's first record.

- **`[tag.XX.strand]` on a two-way rule is a load-time error.** A `strand` tag
  writes one value per strand of origin, and a rule built on
  `[strand.conversion]` — or on a `[strand.origin]` table with a `"+"` or `"-"`
  key — reaches the conversion strand and stops, which does not say whether a
  `+` read is OT or CTOT. `tags::strand_tags_need_origin` compares the two
  before the BAM is opened, and the message names the tag, the rule's file and
  the table to drop. 0.1.17's per-record refusal in `bam_out` stays as a
  backstop that should never speak.
- **`StrandRule::names_origin()`** is the predicate: true when the rule's
  `[strand.origin]` table names only the four strands. A rule's `unknown` key
  does not make it false — a record the rule declines is skipped and counted,
  never written.

## 0.1.18 — 2026-09-17

`--library` is gone. A run gets its strand from a query file's `[strand.*]`
tables and from nowhere else, so there is no rule inside alnbase that can
disagree with the aligner that wrote the BAM — and no way to run without
saying which rule applies. Outputs are unchanged for anyone who was already
passing `queries/strand/directional.toml`, which is what the flag did.

- **A run with no strand rule is refused**, with an error naming a file to
  pass: `queries/strand/` ships one per surveyed aligner. The dry runs that
  read no records — `--explain`, `--list-codes`, a synthetic `--trace` — are
  exempt, since they call no strands. `--trace-records` reads a BAM, so it is
  not exempt.
- **`--library` and `LibraryArg` removed** from `query` and `extract`, and
  `StrandSource` with them: every path now takes one `Arc<StrandRule>`. The
  directional rule survives only as a test-only golden value, which is what
  the shipped `queries/strand/directional.toml` is checked against record by
  record.
- **`extract` recovers the rule from the BAM's header.** With no
  `--query-file`, the tag definition already came out of the stored query
  files; now the strand rule comes with it, so a re-extraction orients offsets
  the way the tagging run did instead of asking the user to remember.
- **`extract --query-file` with files that declare no tag** is an error naming
  the problem rather than a panic (`tag_extract.rs`).
- Every reference example now passes a strand file, and the examples' tagged
  BAMs store it: a run's strand rule travels in the header with its queries.
  Where a run stores more than one file, `dump-query --file N` picks one.

## 0.1.17 — 2026-09-17

The walk takes its strand from the rule a query file declares. A run handed
`queries/strand/bsbolt.toml` now calls strands the way BSBolt writes them; the
0.1.16 warning that a loaded rule was being ignored is gone because it no
longer is. A run that declares no rule still falls back to `--library
directional` and writes what 0.1.16 wrote.

- **One strand source for the whole run.** `StrandSource` is either the rule a
  query file declared or `--library`, resolved once and passed to the walk, the
  tag writer, the partition key and `trace`. A run's rows, its tags and its
  files cannot disagree about a record's strand because there is only one
  decision to disagree with. Retiring `--library` is now deleting a variant.
- **A record the rule declines is skipped and counted.** A rule's `unknown` key
  — BISCUIT's `YD:A:u`, BSBolt's unmapped `YS:Z:WC` — means the input does not
  say which strand a read came from, and there is no walk without one. The run
  reports the count at the end, which is what lets a user widen the rule's
  tables until it covers their input. A tagged BAM writes such a record
  through untagged, as it does one the walk cannot walk.
- **Strand columns are null when there is no call.** `strand`, `conv_strand`
  and `read_reverse` are null on the pass-through row of a record that was
  never walked. An unmapped read carries none of the evidence a rule reads, so
  it is not held to the rule at all — failing a run over records it was never
  going to scan would make every shipped rule unusable on a real BAM.
- **A strand tag from a two-way rule is an error, not a guess.** A rule naming
  only the conversion strand cannot say which of OT and CTOT a read came from,
  and `[tag.XX.strand]` names all four. This belongs at load time and is still
  open; the run-time refusal is the backstop until then.
- **The manifest records the rule, not the library.** `library: "directional"`
  became `strand`, naming the query file the rule came from (the file's whole
  text is already in `query_files`) or `--library directional`. The old field
  would have claimed `directional` for a run that used a rule, which is the
  silent wrongness this whole sequence of releases is about. Hit-table format
  version 9.

## 0.1.16 — 2026-09-17

The strand rule a query file declares now survives loading. Nothing consults it
yet — the walk still takes the strand from `--library` — so every output is
byte for byte what 0.1.15 wrote.

- **The run's rule is carried on `Resolved`**, named after the file it came
  from, so that a record no rule covers can say which rule failed to cover it.
- **Two `[strand.*]` declarations in one run is an error** naming both files.
  This is a check rather than a merge: two rules can only disagree about the
  same record, and there is no principled way to pick a winner.
- **A loaded rule warns that it is not yet applied.** The failure it guards
  against is silent: a run handed `queries/strand/bsbolt.toml` would otherwise
  call strands the directional way and write plausible output on the wrong
  strand. The warning goes when the walk takes the rule.

## 0.1.15 — 2026-09-17

The last three places that read the FLAG for strand now read the strand call
instead, and the call's third component gets a column of its own. Under
`--library directional` — still the only producer — every value is what it was.

- **`strand` and `conv_strand` are written from the call**, not re-derived. They
  were the last sites deriving a strand independently of the walk, so a rule
  that disagreed with the FLAG would have produced rows whose `strand` column
  contradicted their own `refr_pos`.
- **`conv_strand` is null when the rule names only the conversion strand.**
  Several aligners record two strands rather than four; reporting an invented
  `OT`/`CTOT` split would present a guess as a reading.
- **New `read_reverse` column** (in `all` and `bools`): whether the read as
  sequenced runs backwards along the reference, which is what `off_5p` and
  `off_3p` count from. `is_reverse` is the FLAG bit and a different question —
  the same question only for aligners that use the bit that way.
- **The partition key makes the same call.** It is hashed on the reader thread,
  before any worker has seen the record, so it cannot borrow the walk's call; it
  holds the rule and reaches the same conclusion, which is why partitioning by
  `strand` cannot split a file differently from how the column reads.

## 0.1.14 — 2026-09-17

The strand a record came from is now a value passed around the run, not
something three places each re-derive from the FLAG. Behaviour is unchanged:
the only producer is still `--library directional`, which now answers with a
`StrandCall` instead of a `Strand`.

- **`walk_alignment` takes the walk direction** rather than computing
  `is_last_in_template() == is_reverse()`. An aligner that writes the converted
  reference strand into 0x10 instead of the read's orientation would otherwise
  be walked backwards, silently, with every coordinate in the output wrong; a
  new test walks one flagless record both ways to show only the call decides.
- **`off_5p`/`off_3p` mirror on the call, not on read 2.** `HitBuilder` now
  holds `mirror` — whether the walk ran against the direction the read was
  sequenced in — which for a directional library is exactly read 2, and is
  asserted to be so for every FLAG the old rule distinguished.
- **`StrandCall::walk_reversed` and `mirrors_read`** are the two questions the
  rest of the program asks about a strand; everything else is presentation.
- **`Matcher` holds the library and calls the strand once per record**, handing
  the same value to the walk and to the row being written, so the walk direction
  and the written offsets cannot come apart.

## 0.1.13 — 2026-09-17

Query files can now declare a strand rule. Nothing applies it yet.

- **`[strand.origin]`, `[strand.conversion]` and `[strand.sequenced]` parse**,
  through the same typed, unknown-key-denying reader the other tables use, so a
  mistyped `[strand.orgin]` is an error rather than a table that does nothing.
  A rule that does not hold together is reported at `[strand]`, naming the key
  whose condition failed to parse.
- **A file may hold nothing but a strand rule**, which is what every file in
  `queries/strand/` is; the "no queries or tags in file" error became "no
  queries, tags or strand rule in file".
- **The tables are documented in `--help`** beside `[tag.XX.bases]`.
- All twelve shipped rules are now tested twice over: they compile in the
  engine's own tests, and they load through `parse_file` as the files they are.

## 0.1.12 — 2026-09-17

The strand-rule engine: `queries/strand/`'s files now compile and can answer for
a record. Nothing calls it yet, so no behaviour changes.

- **`src/strand_rule.rs`.** Parses the `[strand.origin]`, `[strand.conversion]`
  and `[strand.sequenced]` tables of `docs/design/strand-rules.md` into a
  `StrandRule`, and answers `rule.call(record)` with the conversion strand (the
  direction the walk runs), the sequenced direction (what `off_5p`/`off_3p`
  count from) and the strand of origin, which is `None` when the input
  distinguishes only two of the four. `unknown` is a call of its own.
- **Conditions are `and`/`or`/`not` and parentheses** over named record fields,
  `flags == 99`, and aux-tag equality such as `XG == "CT"`, with a bare
  two-character name read as a text tag. Precedence matches a `where` clause. A
  tag absent from the record makes every comparison on it false, so the record
  matches nothing and is reported by name rather than assigned a strand.
- **Exactly one key must match each record**; none or several is an error naming
  the read, its FLAG and the keys that matched. Load-time checks reject the two
  table forms together, `[strand.conversion]` without `[strand.sequenced]`, a
  four-way origin table that also declares the sequenced direction, and two keys
  carrying the same condition.
- **All twelve shipped rules are now tested**, as files: each compiles, and
  `queries/strand/directional.toml` is checked record by record against the
  hardcoded `Library::Directional` it will replace, including the walk direction
  and the sequenced direction, so the replacement changes nothing.

## 0.1.11 — 2026-09-17

Fixes the Bismark `XR` example, which named the wrong strands.

- **`[tag.XR.strand]`'s example mapping was inverted on the bottom-strand pair.**
  Every copy in the documentation and in the test fixtures read
  `CT = ["OT", "CTOB"]`, `GA = ["OB", "CTOT"]`. Bismark's `XR` is the conversion
  the read itself shows, so it is `CT` for the two converted originals, OT and
  OB, and `GA` for their complements, CTOT and CTOB — which is what the
  surrounding comment always said and what real Bismark BAMs contain
  (`alnbase-validation/docs/research/strand-determination.md`). The strand tables are user
  configuration, so no code changed and no behaviour did; what changed is the
  example everyone copies, the `bismark-tags` reference example's output, and
  the independent model it is checked against, which had the same inversion.
- **Strand rules, written but not yet implemented.** `queries/strand/` holds one
  file per surveyed aligner declaring how to recover a read's original strand
  from its flags and tags, in the syntax proposed in
  `docs/design/strand-rules.md`. Nothing reads them yet; `--library directional`
  is unchanged.

## 0.1.10 — 2026-09-17

Soft-clipped bases get their own flank symbol; ERR is gone and the indexer
refuses what it stood for.

- **New read symbol CLIP**, written `:` (`L` in a group or one-character
  output). Beside a soft-clipped read end, the flank columns nearest the aligned
  part stand for the clipped bases, one each, up to the `--end-context` width: the
  real reference base, the clipped base's own `off_5p`/`off_3p`, null quality.
  The rest of the flank is pads, past the read, with null offsets. A long clip
  therefore adds no columns. `~` matches a clip; `_` does not, so `_N` now fires
  on the end base only of an unclipped read end, `:N` on a clipped one, and an
  alias for `{_:}` either way. Queries that relied on `_` next to clipped reads
  change their hits. A `bases` tag never marks a clip (counted as unplaced), and
  the concordance check ignores clips.
- **ERR removed.** `alnbase index` refuses any sequence byte that is not an IUPAC
  base or `U`, in either case, with the line, contig and 1-based position, and
  writes no index. `-` and `.` (stored as a gap) and `0` (the empty set) are
  refused too. `%` is no longer a pattern code.
- **`.aref` format version 3**: rebuild indexes with `alnbase index`.
- **Hit table format version 8.**
- **Fix: capture offsets for a record with QUAL `*`.** `capture_off_5p` and
  `capture_off_3p` were null for every capture of such a record, because they
  followed the quality rather than the symbol. They are now null only on a
  deletion, intron column or pad.
- The unplaced-hits summary reads `N hits were anchored on a pad or a
  soft-clipped base, so had no aligned base to tag`.
- `--list-codes` lists the clip, and its junction note says k context bases, not
  k − 1 (stale since 0.1.8).
- The walk oracle derives clip columns independently and checks them.

## 0.1.9 — 2026-09-17

Pads carry no read offset; soft-clip record fields.

- **Pads have null `off_5p` and `off_3p`**, as the anchor of a hit row (they were
  already null in capture lists). A pad is not a base of the read: past a soft
  clip the offset it used to take was the clipped base's, beside a skipped
  leading insertion it was not, and on a supplementary alignment no offset
  describes it. `refr_pos` is unchanged.
- **Fix: a `bases` tag could mark a soft-clipped base.** A hit anchored on a pad
  next to a soft clip had the clipped base's offset, which lies inside SEQ, so the
  tag marked that base instead of counting the hit as unplaced. Queries anchored
  on a read base (all of `docs/bismark-xm.toml`) were not affected.
- **New record fields `soft_clip_5p` and `soft_clip_3p`**: soft-clipped bases at
  the read's sequenced 5' and 3' ends, already included in `off_5p`/`off_3p`, so
  `off_5p - soft_clip_5p - hard_clip_5p` counts from the first aligned base. In
  `all`, not in `core`.
- `--trace-records` prints `walk offset none` for a pad anchor.
- **Hit table format version 7.**
- The walk oracle (`tests/oracle/`) checks null pad offsets and the soft-clip
  fields.

## 0.1.8 — 2026-09-16

Requirements of the aggregation utility, and further reference safeguards.

- **Pads and intron context default to the widest query's span** (was span − 1).
  `--end-context` and `--splice-context` unset now show k columns past each read
  end and at each intron edge, where k is the widest query's span, so the widest
  query too can match windows made only of pads or intron context, like every
  narrower one. Runs with pad- or context-matching queries gain those hits (for
  example `mark = ".+"` over `read = "~_"` now counts 3 unplaced hits on one read,
  not 1); set both options to keep the old layout.
- **Run manifest.** Every file of `query --parquet` and `extract` carries plain
  parquet footer key-value entries (Arrow schema metadata for IPC):
  `format_version`, `coordinate_base` and `alnbase_manifest`, a JSON record of
  the run: alnbase version, `run_id`, command line, input and reference, library,
  output settings with the file names in slot order, the text of every query file
  (for `extract`, the tag and where its definition came from), the walk with its
  defaults resolved, and every `@SQ` contig with length, MD5 and whether the
  reference has it. `format_version` was previously visible only inside the
  serialized Arrow schema.
- **`{stem}.manifest.json` completeness marker**, written last beside OUT when a
  run succeeds: the manifest plus record counts, concordance totals and rows per
  file. A run removes any manifest already at that path before writing, and a
  failed run leaves none. Files beside OUT that follow its naming but are not the
  run's (left by an earlier run with more shards) are reported in a warning.
- **`coordinate_base = 0`** is stored in every hit table's metadata.
- **Hit table format version 6** (metadata added; no history of earlier versions
  is kept before the first release).
- **New record field `conv_strand`**: `OT`, `OB`, `CTOT` or `CTOB` under the
  directional library, by the same derivation as `strand` tags. In `all`, not in
  `core`.
- **New `query --require-m5`**: a record on a contig whose `@SQ` line has no `M5`
  stops the run, even with `--permissive`. For the tagged BAM, `--parquet` and
  `--trace-records`.
- **`--trace-records` checks reference concordance**: each traced record reports its
  rate, traced records comparing at least 1,000 bases above `--max-discordance`
  stop with the scan's error, and fewer bases above it print a warning.
- **`--explain` notes read-2 direction** under any query that places a pad:
  patterns run along the conversion strand, offsets from the ends as sequenced.
- Tests: the tag/extract round trip now covers a deletion, a hard-clipped
  supplementary alignment and `conv_strand`; new tests for `--require-m5`, the
  footer and stream metadata, and the explain note.

## 0.1.7 — 2026-09-16

- **Reverts the 0.1.3 firing rule.** A query again fires on any window of the
  columns the walk emits, including windows made only of pad or intron-context
  columns. Whether all-pad windows should count is the query's decision: exclude
  them with `and not` over an all-pad pattern (for example
  `[pat.all_pad] read = "__"  refr = "~~"` and `where = "cg and not all_pad"`).
- **New `--end-context N`**: the number of pad columns past each read end, which
  was always derived. Both it and `--splice-context N` default to (widest query
  span − 1). With the defaults, adding a wider query adds pad or intron columns,
  so a query that can match windows made only of those columns can gain hits;
  set both explicitly to make every query's hits independent of the other queries
  in the run. A randomised test checks that independence with fixed contexts, and
  that it fails with the derived defaults.

## 0.1.6 — 2026-09-16

- **The concordance check judges the whole run.** 0.1.4 counted disagreements per
  worker and stopped only when one worker alone had compared 100,000 bases, so a
  wrong reference could finish undetected when the input was small or spread over
  many threads (a one-base-shifted reference passed at 68.7% with `-@ 4`). The count
  is now shared by every worker; the run stops once the shared total reaches
  100,000 compared bases above `--max-discordance`, and a run that ends sooner is
  judged on its total (from 1,000 compared bases), removing its output on failure.
  A shifted reference now stops at about 100,000 bases with `-@ 1`, 4 and 16 alike.

## 0.1.5 — 2026-09-16

`overlap` now writes a valid, mate-consistent SAM. Outputs differ from 0.1.4 wherever a
template had an alignment clipped, dropped or promoted.

- **One primary line per read.** When a read's primary alignment has nothing left after
  clipping and one of its supplementary alignments survives, the survivor's alignment
  (RNAME, POS, MAPQ, strand, CIGAR) now moves onto the primary record, whose SEQ holds the
  whole read: the survivor's hard clips become soft clips (except on a side that
  `--clip-mode hard` itself clipped), NM/MD/AS come from the survivor, and the
  supplementary record is not written. 0.1.4 wrote the primary unmapped *and* cleared
  0x800 on the survivor, giving two primary lines for one read (flagstat: 11 primaries
  for 10 reads).
- **Unmapped reads follow SAM conventions.** A primary with nothing left and no surviving
  alignment is written with 0x4 set, 0x2 cleared, MAPQ 0, CIGAR `*`, SEQ/QUAL
  reverse-complemented back to sequencing orientation and 0x10 cleared, and RNAME/POS
  copied from its mapped mate (`*`/0 if the mate is unmapped too). 0.1.4 kept MAPQ, 0x2,
  0x10 and the old placement.
- **Mate fields are rebuilt.** In every template where an alignment was clipped, dropped or
  promoted, RNEXT, PNEXT, TLEN, flags 0x2/0x8/0x20 and `MC` (and `MQ` where already
  present) are recomputed from the final primary lines using `samtools fixmate`'s rules:
  TLEN from 5' end to 5' end, `MC:Z:*` for an unmapped mate, 0x2 cleared unless the
  primaries are mapped forward-reverse on one reference (never set). Supplementary and
  secondary lines get the same mate fields relative to the mate's primary. `samtools
  fixmate -m` on the output changes none of these fields. Templates that were not changed
  keep their mate fields as they came.
- **SA is rebuilt** on every surviving alignment of a read whose alignments changed (if
  the read carried SA): the other surviving alignments with their final position, strand,
  CIGAR (hard clips written as soft), MAPQ and NM, primary first; removed when no other
  alignment is left. Where `--stale-tags strip` removed a clipped record's NM, the SA
  entry's NM is recounted from the input MD over the bases that remain.
- **QUAL `*` stays `*`.** A base with no quality takes no part in the quality arithmetic
  (neither copy is quality-edited), and `--keep score` drops the quality term when either
  read has none. 0.1.4 wrote real values into part of the 255-filled string, producing an
  invalid QUAL.
- **`--compression-level` is honoured** (it was ignored) and must be 0 to 9.
- `--no-tag` help and cli-reference now say what it covers: the length tag and the prefixed
  overlap tags only. Tags describing a changed alignment (NM/MD/AS/XA per `--stale-tags`,
  SA, MC, MQ) are still kept correct, as before.
- The `overlap` module header no longer claims that two-copy molecules are refused; as the
  code and tests already did, the true fragment length is resolved and the decoy is
  passed over by the anchor test.

Docs: `docs/reference/05-overlap.md` gains an "Edge cases and how they are resolved"
section with diagrams (`docs/reference/figures/overlap-*.bob`) and a new example,
`examples/05-overlap/09-edge-cases-and-mate-fields`, which runs `samtools fixmate -m` and
`flagstat` on the output.

## 0.1.4 — 2026-09-16

Behaviour changes:

- **`off_5p` / `off_3p` are relative to the read, not the record.** Hard-clipped
  bases were sequenced, so they are now counted: a supplementary alignment with
  `60H20M` reports `off_5p = 60` for its first base, the same offset the primary
  alignment reports for that base. `off_5p + off_3p` is SEQ length plus hard clips
  minus 1.
- New record fields `hard_clip_5p` and `hard_clip_3p` (hard-clipped bases at the
  read's sequenced ends).
- **A missing quality (QUAL `*`) is null**: `qual` in hit rows (was 255),
  `qual_phred` (was a string of 0xff+33 characters), and `extract` rows.
- **Reference concordance check.** Every walked record counts read/reference
  disagreements over observed, unambiguous bases (a reference C read as T is not
  counted). The run summary reports the rate; a worker stops the run when the rate
  exceeds `--max-discordance` (default 0.25) after 100,000 compared bases. On
  correct data the rate is a few percent at most; a reference shifted by one base
  stops the run at about 63%. `--max-discordance 1` turns it off.
- Human-facing coordinates say their base: `--trace` hit lines give
  `refr_pos N (0-based)` and a "walk offset" rather than `read_off`, and read
  "matches ending at column"; `--trace-records` labels records with
  `(1-based POS)`; SEQ positions in error messages say `(0-based)`.

Added `tests/oracle/walk_oracle.py`, which checks every column's coordinates
against pysam and pyfaidx.

## 0.1.3 — 2026-09-16

- **Queries are independent.** A query fires only on a window that contains at
  least one column the read observed (a read base, an inserted base, or a
  deletion). Previously a query whose read row accepts pads (`~`, `_`) could fire
  on a window made only of flank columns; the flank grows with the widest query in
  the run, so adding an unrelated wider query added hits to a narrower one (a
  reference-context CG query: 1 hit alone, 3 beside a six-column query). The flank
  itself is unchanged. Intron-context columns and the `,@,` elision marker are not
  observations either: a query about the marker needs an explicit, small
  `--splice-context`. `--trace` and `--trace-records` apply the same rule. A
  randomised test checks that every query in a varied panel fires identically
  alone and in company, over all flag cases and CIGAR operations.

## 0.1.2 — 2026-09-16

Behaviour changes (outputs differ from 0.1.1):

- **`off_5p` / `off_3p` are measured from the ends of the read as sequenced for
  every record.** Previously read 2 (FLAG 0x80) was measured along the
  conversion strand, so its `off_5p = 0` was its last sequenced base. Pattern
  semantics are unchanged: patterns still run along the conversion strand, so
  for read 2 a pad written before a base (`_N`) is at its sequenced 3' end.
- **Soft-clipped bases are ignored entirely.** `--soft-clips` is removed; clipped
  bases never produce columns (they were previously emitted like insertions
  under `--soft-clips emit`). Offsets still count them, since they were sequenced.
- **`--batch-rows` is a total across all output files**, divided evenly over
  `--threads` x `--shards-per-worker` files, so it sizes memory directly.
  `--batch-rows` and `--row-group-rows` must be at least 1 (0 used to crash a
  worker).
- **Reference index format version 2.** Each contig's MD5 (the SAM `@SQ M5`
  digest) is stored; version 1 indexes are refused, so rebuild with
  `alnbase index`. `info` prints the digests.
- **`@SQ M5` is checked.** When the BAM header carries M5 for a contig, it must
  match the index, or the run stops: same name and length with different
  sequence is otherwise undetectable.
- `index` refuses a FASTA with duplicate contig names (the last one used to win
  silently) and writes to `OUT.partial`, renaming on success, so a failed run no
  longer truncates an existing index.

Fixes:

- A mapped record with SEQ `*` crashed a worker; it is now skipped and counted.
- Parquet/IPC runs reported every worker failure as
  `worker N stopped early (sending on a closed channel)`; the worker's own error
  (for example a missing output directory, or too many open files) is reported.
- A failed parquet/IPC run now removes every output file it may have written.
- When the OS refuses a thread (per-user thread limit), the run ends with an
  error naming `ulimit -u` instead of a panic, and a tagged-BAM run removes its
  partial output.
- `--parquet --trace-records N` traced nothing and ran a full scan; the trace now
  takes precedence like the other dry runs.
- `--list-codes` no longer advertises repeat counts (`C7`) and in-row markers,
  which query files reject; it describes the `mark` row instead.
- Parquet runs' summary said skipped records were "written untagged"; it says
  they were given a null row.

## 0.1.1 — 2026-09-16

Documentation and help text only; no change to what any command computes or writes.

- Output file naming documented as `{stem}_{worker}_{shard}.{ext}` everywhere
  (module docs said `hits_…`, user docs said `OUT_…`).
- Partitioning guarantee stated precisely: within a run, all records sharing a
  partition key are in exactly one file; across runs the same key lands in the
  same file only if the thread count and `--shards-per-worker` are unchanged.
- Hit column names in docs and help corrected to `off_5p` / `off_3p` (were
  `read_off` / `read_off_3p`); a nonexistent `expr` column removed from docs.
- `query --help` example uses `.aref`, not `.mm`.
- Error message for missing positional paths names the real dry-run flags
  (`--explain`, `--trace`, `--list-codes`) instead of a nonexistent `--emit`.
- Rustdoc link syntax no longer appears literally in `--soft-clips` /
  `--insertions` help.
- Stale references to removed modules (`visual`, the one-line query syntax,
  `cli.rs`) and a "writer does not exist yet" note removed from module docs.
- `docs/cli-reference.md`: `extract` reads standard input; `--stale-tags`
  default is `strip`; partial-output removal applies to `query` writing a BAM.

## 0.1.0

Initial prototype.
