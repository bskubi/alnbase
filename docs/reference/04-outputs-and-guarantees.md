# 04 — Outputs and guarantees

Scope: everything alnbase 0.1.5 **writes** and what a pipeline may rely on:
the tagged BAM (`bases` and `strand` tags, existing and damaged tags), the
query files stored in its header, `dump-query`, `extract`, the parquet / Arrow
IPC hit tables (every column, every record field, aux typing, file metadata
and the run manifest), rows and their
order, `record_id`, sharding and partitioning, run summaries, errors and exit
codes, and the reference concordance check. How columns are produced and which column each hit value comes from is
in [03 — walk and matching](03-walk-and-matching.md) (especially §3 offsets,
§6.2 hit values, §7 records not walked, §8 odd SEQ/QUAL); the query language is
in [02 — query language](02-query-language.md). This document cross-references
them rather than repeating them.

Conventions:

- `src/file.rs:N` cites the code. **(example NN)** means the claim was checked by
  running `examples/04-outputs/NN-*/run.sh` against the release binary
  (`alnbase` release build, samtools 1.24, duckdb,
  pyarrow); the output is that directory's `expected.txt` (temporary paths
  rewritten to `TMP/`). **(code only)** means read from source, not executed.
- Coordinates in hit tables are 0-based. SAM `POS` is 1-based.
- "Top"/"bottom" strand: OT and CTOT are top, OB and CTOB are bottom
  (`src/tags.rs:105-107`); the `strand` record column prints them as `+`/`-`.
- "CG" is used for the dinucleotide (not "CpG").

Examples (each has `ref.fa`, `reads.sam` whose `@CO` lines say what each record
tests, one or more query TOML files, `run.sh`, `expected.txt`). Run with
`PATH=<dir with samtools, python3+duckdb+pyarrow>:$PATH ./run.sh`; override the
binary with `ALNBASE=...`.

| NN | directory | shows |
|---|---|---|
| 01 | `01-bismark-tags` | `docs/bismark-xm.toml` + XR/XG tables vs an independent Bismark model on OT/OB/CTOT/CTOB; `@PG`/`@CO` storage; `dump-query` round trip; existing tags; re-tagging |
| 02 | `02-bases-tag-edges` | off-reference stop / `--permissive`; unplaced hits; unused-query warning; code clash; anchor-on-gap refusal; no tags; output name rules; misplaced options; exit codes |
| 03 | `03-parquet-schema` | flat and list schemas with nullability; `-f all`; every aux type; type mismatches; aliases; non-UTF-8; `--only-hits`; IPC and its schema metadata; row groups, codec, encodings, footer metadata and manifest keys; row order within a record |
| 04 | `04-sharding-ordering` | `--partition-by CB:Z -@ 4 --shards-per-worker 2`; file names; `shard` column; `record_id` density and order; empty files; thread-count invariance of rows and of the tagged BAM |
| 05 | `05-extract` | `extract` vs `query --parquet`; `refr_base` rules; hitless differences; `--tag`; strand-tag refusal; stdin; foreign tags; bad characters; no reference needed; unmapped records; aux text not checked |
| 06 | `06-failures-diagnostics` | QUAL `*` → null (255 before 0.1.4); failed parquet runs remove their files and report the worker's own error; BAM partial removal; degenerate option values (0 refused); `-` as a parquet prefix (and its manifest); `dump-query` errors and stdin |
| 07 | `07-idioms` | cell barcode as column and partition key; MAPQ/duplicate filtering downstream; context by query name vs by capture |
| 08 | `08-header-storage` | two runs, three stored files, escaping (TAB, backslash, CR, no final newline); `--pg`/`--file`; per-tag definition choice in `extract`; survival through `samtools sort`; a lost `@CO` line |

---

## 1. The outputs at a glance

| command | output | order | reference needed | pipeline |
|---|---|---|---|---|
| `query ... OUT` | one BAM (always BAM, BGZF), `-` = stdout | input order | yes | ordered (`src/ordered.rs`) |
| `query --parquet ... OUT` | `threads × shards-per-worker` parquet or IPC files beside `OUT`, then `{stem}.manifest.json` (§5.1) | see §5.5 | yes | sharded (`src/parallel.rs`) |
| `extract BAM OUT` | same file layout, schema family and manifest as `query --parquet` | see §5.5 | no | sharded |
| `dump-query BAM` | text on stdout | — | no | — |

`query` without `--parquet` writes **only** the tags declared by `[tag.*]`
tables; with `--parquet` the tag tables are ignored and every query writes rows
(`src/main.rs:378-380`, `src/cli_query.rs:320-325`). Options belonging to the
other mode are refused, not ignored: parquet options without `--parquet`, and
`--overwrite-tags`, `--writer-threads` with it
(`src/cli_query.rs:614-625`, `src/main.rs:524-529`) (example 02).
The `strand`, `conv_strand` and `read_reverse` columns are written from the
run's strand call — the same call the walk and the tag writer were given, so
a row and a tag from one record cannot disagree (`src/record_columns.rs:263-284`).

---

## 2. The tagged BAM

### 2.1 What reaches the file

- **Every input record, exactly once, in input order.** One reader numbers
  batches of 512 (`src/batch.rs:103`), any tagger tags a batch in place, the
  writer holds early batches in a `BTreeMap` and writes strictly by batch number
  (`src/ordered.rs:174-190`). Verified with 1500 records, `-@ 1` and `-@ 4`:
  `samtools view` output identical, and fields 1-11 identical to the input in
  order (example 04).
- A record is **walkable** when it is mapped (FLAG 0x4 clear), its RNAME is a
  contig the reference index has, and it has a SEQ (not `*`)
  (`src/scanner.rs:558-570`, called at `src/ordered.rs:262-263`). Walkable
  records get every declared tag. Unmapped records, and mapped records with
  SEQ `*`, are written **unchanged** (no tags at all) and counted (examples 01,
  02; SEQ `*` code only).
- A mapped record on a contig the index lacks **stops the run** with
  `record NAME is on contig X, which the reference does not have ...`
  (`src/batch.rs:37-43`, `69-73`) unless `--permissive`, in which case it is
  written unchanged and counted (`src/batch.rs:74`) (example 02). A header that
  merely *lists* contigs the index lacks is fine; the warning
  `1 of 2 BAM contigs are not in the reference; records on them will be skipped`
  is printed either way (example 02).
- With `--require-m5`, a record on a contig the index has whose `@SQ` line
  carries no `M5` **stops the run**, with or without `--permissive`:
  `record NAME is on contig X, whose @SQ line has no M5 checksum, and
  --require-m5 was given. ...` (`src/contig_map.rs:149-170`,
  `src/scanner.rs:563-564`, `src/batch.rs:45-53`, `75-80`). Contigs no record is
  on are not checked, and a contig the index lacks is handled as above (checked
  while writing this: a header whose only record-bearing contig lacks `M5` stops
  the run; with `M5` on it, an `M5`-less contig without records and a
  `--permissive` off-reference record pass). The partial BAM is removed like
  any failed run's.
- htslib itself rewrites a SAM record with FLAG 0 and RNAME `*` to FLAG 4 on
  conversion, so "mapped with no contig" does not reach alnbase from SAM text
  (observed with samtools 1.24 while building example 02).
- The header is the input header plus one `@PG` line and the stored query files
  (§3) (`src/bam_out.rs:485-556`).
- The file is opened with `hts_open(path, "wb")` (`src/bam_out.rs:545`): always
  BAM, never SAM/CRAM, whatever the name. Compression threads:
  `--writer-threads`, 0 = min(4, `-@`) (`src/main.rs:405-414`).
- The writer closes the file explicitly and reports a failed final flush
  (`src/bam_out.rs:572-578`) (code only).
- **Output path** (`src/bam_out.rs:226-238`): `-` (stdout), a name ending
  `.bam` (any case), or a name with no extension (FIFO, `/dev/stdout`).
  Anything else is refused before reading (`out.sam` → error) (example 02).
- An **empty input** writes a valid empty BAM (`src/ordered.rs:486-491`,
  a unit test).
- **No tags declared** → error `no tags to write` (`src/bam_out.rs:127-132`)
  (example 02). A query no tag uses still runs and prints
  `warning: query 'X' is used by no tag, so it writes nothing to the BAM output`
  (`src/bam_out.rs:157-166`) (example 02).

### 2.2 `bases` tags

Declared as `[tag.XX.bases]` (syntax and parse-time rules: 02). Resolved
behaviour:

| aspect | behaviour | source |
|---|---|---|
| tag name | two chars, letter then letter/digit | `src/tags.rs:42-52` |
| value type | `Z` string | `src/bam_out.rs:413`, `:455-460` (example 01) |
| length | exactly `seq_len()` (SEQ length; hard clips excluded) | `src/bam_out.rs:286-296` |
| fill | one printable non-space ASCII char, required; may not equal a code | `src/tags.rs:271-292`, `:327-339` |
| codes | one-char keys, printable non-space ASCII, each naming exactly one query; a query may be under only one code of one tag | `src/tags.rs:295-324` |
| placement | at the SEQ index of the query's **anchor** column | `src/bam_out.rs:301-328` |
| strand | on bottom strands (OB, CTOB) the walk offset is flipped: `pos = len-1-off`, so characters are in **SEQ order** (reference-forward, as stored) on every strand | `src/bam_out.rs:287-288`, `:320` (example 01) |
| orientation source | the run's strand rule, through `StrandCall::walk_reversed()` (§2.3) | `src/strand_rule.rs`, `src/bam_out.rs:286-288` |

**Which base gets the character.** The anchor's SEQ index. For bottom-strand
reads a CG call therefore sits on the **G** of the CG in SEQ (read `G`
protected, `A` converted), exactly as Bismark's `XM` (example 01, `ob`/`ctob`).
Inserted bases have SEQ indices too and can be marked when the walk emits them
(`--insertions emit`; 03 §4.4). Soft-clipped bases are never aligned columns:
they appear only as clip columns in the flank, and a hit anchored on one is not
marked (below; 03 §5.2).

**Anchors with no read base are refused at parse time.** A query used by a
bases tag whose anchor column's read side admits a gap or junction is rejected:
`code 'z': query 'q' can fire with its anchor on a deletion or intron ...`
(`src/tags.rs:232-245`) (example 02). The tagger re-checks at run time and
errors rather than marking a wrong base (`src/bam_out.rs:306-309`,
`:331-337`, code only).

**Anchors on a flank column** (a pad past either end of the read, or a clip
standing for a soft-clipped base, 03 §5) have no aligned base: the hit is **not
written and is counted** (`src/bam_out.rs:314-316`).
The run summary prints
`  N hits were anchored on a pad or a soft-clipped base, so had no aligned base to tag`
(it covers both ends) (`src/main.rs:699-704`) (example 02: query `end`,
`read = "N_"`, `mark = ".+"`, one per read). A query may also fire on windows
made only of pad columns (`_@_` fires on every pad), and with the default
`--end-context` their number follows the widest query in the run, so such counts
can change when a wider query is added ([`02-query-language.md` §8.6](02-query-language.md#86-pads-flank-contig-ends-and-record-boundaries)).

**Conflicts.** If two queries under **different codes of the same tag** mark
the same SEQ position of one record, the run **fails** at that record:
`record fwd: queries 'CG' and 'anyC' both mark SEQ position 8 (0-based) of tag
XM, as 'Z' and 'C'. ...` (`src/bam_out.rs:318-327`, `:339-355`) (example 02). There is no
precedence rule. Different tags are independent: the same base may be marked in
`XM` and in `XE`. The same query firing at one anchor once cannot conflict with
itself. Nothing checks for overlap statically; a query set that never collides
on the data at hand runs fine.

**No hits** → the tag is still written, all `fill` (example 02, `XE`).

### 2.3 `strand` tags

`[tag.XX.strand]`: each key is a value, each list the strands that write it;
all four strands must be assigned exactly once; strand names exact uppercase
`OT CTOT OB CTOB`; values non-empty printable ASCII, spaces allowed
(`src/tags.rs:346-390`). Written once per walkable record as `Z`
(`src/bam_out.rs:402`, `:413`) regardless of hits. Unmapped and (permissive)
off-reference records get none.

Every run declares a strand rule and a run that reads records without one is
refused (0.1.18); there is no default, since a default is a rule built into
alnbase that can disagree with the aligner that wrote the BAM. The table below
is `queries/strand/directional.toml` (FLAG only, nothing else consulted; a read
without 0x80 is read 1), which is what to pass for an aligner that writes no
strand tag:

| FLAG 0x80 (read 2) | FLAG 0x10 (reverse) | strand | top/bottom | typical FLAGs |
|---|---|---|---|---|
| 0 | 0 | OT | top | 0, 99, 65 |
| 0 | 1 | OB | bottom | 16, 83, 81 |
| 1 | 1 | CTOT | top | 147, 145 |
| 1 | 0 | CTOB | bottom | 163, 129 |

Secondary/supplementary/duplicate bits are ignored: such records are tagged like
any other (README "No filters").

### 2.4 Tags the input already has

Detection is by tag name only, any type (`src/aux.rs:164-183`,
`src/bam_out.rs:358-389`). Without `--overwrite-tags` the first such record
**stops the run**: `record ot already carries tag XM; pass --overwrite-tags to
replace it`, exit 1, output removed (example 01). With it, the old field is
removed and the new one **appended at the end of the aux block**, so tag order
changes (example 08: `XE` now precedes `XM`). Tags this run does not declare
are left alone (example 08: run 2 declared only `XM`; run 1's `XE` survives).

### 2.5 Damaged aux blocks

A block is damaged when a field partway through cannot be sized (unknown type
byte, truncated value, `B` subtype other than `cCsSiIf`) (`src/aux.rs:72-117`);
nothing after it is reachable. Default: the run stops
(`record X: its aux block is damaged ... Pass --permissive ...`,
`src/bam_out.rs:362-370`). With `--permissive` the record is still tagged, all
new tags written **at the front of the aux block** in one rebuild
(`src/bam_out.rs:408-410`, `:439-465`), counted as `damaged_aux` and reported
`N records had a damaged aux block; ...`. An existing same-named tag *after*
the damage cannot be removed, so the block may hold two; ours is first and is
what readers return (`src/bam_out.rs:430-438`). (code only: a damaged block
cannot be produced from SAM text; unit test `src/bam_out.rs:963-988`.)

### 2.6 Reproducing Bismark's `XM`, `XR`, `XG`

`docs/bismark-xm.toml` plus two strand tables yields Bismark's three tags
exactly. The full file used (example 01 `queries.toml`):

```toml
# Bismark's XM, including u/U for cytosines whose context the reference
# does not determine. Three columns: the cytosine and the two bases after it.

# --- the three determinable contexts, as patterns the u/U queries exclude ---
[pat.cg]
read = "~~~"
refr = "CG~"

[pat.chg]
read = "~~~"
refr = "CHG"

[pat.chh]
read = "~~~"
refr = "CHH"

# --- CpG ---
[query.Z]                   # protected
read = "C~~"
refr = "CG~"

[query.z]                   # converted
read = "T~~"
refr = "CG~"

# --- CHG ---
[query.X]
read = "C~~"
refr = "CHG"

[query.x]
read = "T~~"
refr = "CHG"

# --- CHH ---
[query.H]
read = "C~~"
refr = "CHH"

[query.h]
read = "T~~"
refr = "CHH"

# --- unknown context: a reference C whose next two bases settle nothing ---
[query.U]
read  = "C~~"
refr  = "C~~"
where = "U and not cg and not chg and not chh"

[query.u]
read  = "T~~"
refr  = "C~~"
where = "u and not cg and not chg and not chh"

[tag.XM.bases]
fill = "."
z = "z"
Z = "Z"
x = "x"
X = "X"
h = "h"
H = "H"
u = "u"
U = "U"

# --- Bismark's conversion tags (added for this example; not in docs/bismark-xm.toml) ---
[tag.XR.strand]             # read conversion, as the read was sequenced
CT = ["OT", "OB"]
GA = ["CTOT", "CTOB"]

[tag.XG.strand]             # genome conversion: which reference strand
CT = ["OT", "CTOT"]
GA = ["OB", "CTOB"]
```

Expected values were derived independently of alnbase by
`examples/04-outputs/01-bismark-tags/bismark_model.py`, a direct encoding of
Bismark's documented semantics: `XM` has one character per SEQ base in SEQ
order; top-strand reads (OT, CTOT) are called at reference C with context from
the next two reference bases; bottom-strand reads (OB, CTOB) at reference G with
context from the two preceding bases, complemented; read `C`/`G` uppercase
(methylated), `T`/`A` lowercase, any other read base `.`; an `N` in the context
gives `u`/`U`; `XR` is the read's conversion (OT CT, OB GA, CTOT GA, CTOB CT),
`XG` the genome's (OT CT, CTOT CT, OB GA, CTOB GA). Reads are 36M over
`chr1:3-38`, including a methylated CHG, a methylated CHH, an unknown-context C
next to a reference `N`, a mismatch at a reference C (no call), and a mismatch
at a non-C.

| record | FLAG | strand | `XM` | `XR` | `XG` |
|---|---|---|---|---|---|
| `ot` | 0 | OT | `.Z...X..h.H...U..xZ...x.........x...` | CT | CT |
| `ob` | 16 | OB | `..z....X...........Z....xH.....h..xh` | GA | GA |
| `ctot` | 147 | CTOT | `.Z...X..h.H...U..xZ...x.........x...` | GA | CT |
| `ctob` | 163 | CTOB | `..z....X...........Z....xH.....h..xh` | CT | GA |
| `unmapped` | 4 | — | no tags | — | — |

alnbase's output equals the model's output byte for byte (`IDENTICAL`), all
three tags are type `Z` (example 01). `ctot` carries the same bases as `ot` and
gets the same `XM` (calls at C, SEQ order) with `XR:GA`; `ctob` mirrors `ob`.
Scope of the verification: ungapped, unclipped alignments; Bismark's own
handling of indels inside a context window was not modelled (see 03 for how
alnbase's columns treat deletions and insertions).

---

## 3. Query files stored in the header

### 3.1 Format (version 1)

Every `--query-file` a tagged-BAM run reads is stored; a `--parquet` run
writes no BAM header, and stores the query files' text in its hit files'
manifest instead (§5.2). `QuerySource.path` is the path **as given on
the command line** (`src/cli_query.rs:596`); the text is the file read as UTF-8
(a non-UTF-8 file is an error: `std::fs::read_to_string`,
`src/cli_query.rs:582`).

```
@PG  ID:alnbase  PN:alnbase  [PP:<prev>]  VN:<alnbase version, e.g. 0.1.2>  CL:<argv joined by spaces>
@CO  alnbase:v1:file  PG:<pg id>  file:<k>  lines:<n>  path:<escaped path>
@CO  alnbase:v1:line  PG:<pg id>  file:<k>  line:<i>   <escaped line text>
```

(single tabs between fields; example 08 prints it with `cat -A`.)

- `@PG` is added by `sam_hdr_add_pg`, which picks a unique ID (`alnbase`,
  `alnbase.1`, ...) and chains `PP` to the previous program
  (`src/bam_out.rs:502-515`). Tabs, CR, LF and NUL in `VN`/`CL` become spaces
  (`src/bam_out.rs:168`). `CL` is `std::env::args()` joined with spaces
  (`src/main.rs:393`), so it contains the binary path as invoked and the
  temporary paths used (example 01).
- The stored lines' `PG:` is the first `@PG` ID absent from the input header
  (`src/bam_out.rs:517-541`).
- One `file` record per query file, numbered from 0 in command-line order, then
  one `line` record per line: the text is split on `\n` only, so a file ending in
  a newline has a final empty line (`lines:74` for a 73-line file, example 01;
  `lines:13` for a file without a final newline, example 08)
  (`src/header_query.rs:62-77`).
- Escaping (`src/header_query.rs:260-275`): `\` → `\\`, TAB → `\t`, CR → `\r`,
  LF → `\n`, other ASCII controls → `\xHH` (uppercase hex); non-ASCII UTF-8 as
  is. Unescaping rejects anything else (`src/header_query.rs:277-305`).
- Reading (`src/header_query.rs:116-211`): only `@CO` lines whose first field is
  exactly `alnbase:v1:file` or `alnbase:v1:line` are considered; others
  (including a future `alnbase:v2:*`) are ignored. Lines may be in any order.
  A missing tab before an empty line's text reads as an empty line. Runs are
  ordered by the position of their `@PG` line; runs whose `@PG` is gone sort last.
- Detected corruption is an **error, not a shorter file**: a lost line
  (`@PG alnbase.1, run2_calls.toml: 8 of its 9 lines are in the header`),
  a file record described twice, a line stored twice differently, lines with no
  file record, a missing file index (example 08; `src/header_query.rs:154-204`).
  Because every stored run is parsed, corruption in **any** run breaks
  `dump-query` without `--pg` *and* `extract` of every tag, including tags
  defined only by an intact older run (example 08, section 6).
- Standard samtools operations that keep `@CO` lines keep the definitions:
  verified through `samtools sort` (example 08).

### 3.2 Multiple runs

Re-tagging a tagged BAM (with `--overwrite-tags` if tags repeat) appends a
second `@PG` (`alnbase.1`) and a second set of `@CO` lines; the first run's
lines stay (examples 01, 08).

### 3.3 `dump-query`

`alnbase dump-query [--pg ID] [--file N] BAM` (`BAM` may be `-`)
(`src/cli_extract.rs:56-70`, `src/main.rs:808-813`, `src/header_query.rs:219-258`):

- Picks the run by `--pg` (default: the last run in `@PG` order) and the file by
  `--file` (required only when that run stored more than one, which every run
  does since 0.1.18: the strand rule is a query file too).
- Prints exactly `# <the @PG line>\n# file K of N: <path>\n` followed by the
  stored text. **The text after the two comment lines is byte-identical to the
  file the run read** (verified with `cmp`: a 74-line file, a file with a TAB
  and a backslash and no final newline, a file with a CRLF line) (examples 01,
  08). The whole output is itself a valid query file (the two lines are TOML
  comments). If the `@PG` line has been removed from the header the first line
  reads `# @PG\tID:<id> (no @PG record with this ID remains in the header)`.
- Errors (exit 1): `the BAM has no stored alnbase query file`;
  `@PG alnbase stored 2 files; choose one with --file (0: a.toml, 1: b.toml)`;
  `@PG alnbase stored 2 file(s); there is no file 2`;
  `no query file is stored under @PG nope; stored under: alnbase, alnbase.1`
  (examples 06, 08).

### 3.4 How `extract` chooses a definition

(`src/tag_extract.rs:92-183`)

- **With `--query-file`**: the header is not consulted at all. All given files'
  queries and merged tags form the definition (a tag declared in two given
  files is an error, `src/tags.rs:218-226`).
- **Without**: runs are walked newest first; **each tag** takes its definition
  from the newest run that declares it (tag by tag, not run by run). The
  definition's queries are all queries of all files of that run.
  Example 08: run 1 declared `XM` and `XE`, run 2 redeclared `XM`;
  `--tag XM` uses `@PG alnbase.1 (run2_calls.toml)`, `--tag XE` uses
  `@PG alnbase (run1_calls.toml, run1_end.toml)`.
- **Tag choice**: `--tag` if given; else the single `bases` tag defined; else
  error `several bases tags are defined (XM, XE); choose one with --tag` or
  `... defines no bases tag`. Naming a strand tag is an error
  (`XR is a strand tag in ...: it holds one value per read, not calls to
  extract`); naming an undefined tag lists the defined ones (examples 05, 08).
- The choice is printed to stderr:
  `extracting tag XM with the query file stored under @PG alnbase (calls.toml)`
  or `... with the given query files`.

---

## 4. `extract`

`alnbase extract [OPTIONS] BAM OUT` (`src/cli_extract.rs:11-54`,
`src/main.rs:741-805`). `BAM` may be `-`: the one reader reads the header,
chooses the definition, and scans on (example 05). Options: `--tag`,
`--query-file` (repeatable), `--permissive`, `-@`,
`--reader-threads`, and all *Parquet output* options (§5.8, §5.6). No
reference argument.

### 4.1 What a row is

For each record that carries the tag, one row per tag character that is not
`fill`, emitted in walk order (along the conversion strand: SEQ ascending on top
strands, descending on bottom) (`src/tag_extract.rs:324-357`). The row order
is the walk's; the offsets are not (below). Columns
(`src/hits.rs:261-296`):

| column | value |
|---|---|
| `name` | the query the character stands for |
| `off_5p` | the base's offset from the first base **as sequenced**: `k`, the index along the strand, mirrored for read 2 (FLAG 0x80), plus the hard-clipped bases at the sequenced 5' end, exactly as a direct run reports it (`src/hits.rs:280`, 03 §3.2) |
| `off_3p` | `seq_len + hard clips - 1 - off_5p` |
| `refr_pos` | reference coordinate of that SEQ base; insertions and soft-clipped bases take the nearest 5' coordinate along the strand (`src/tag_extract.rs:392-422`) |
| `qual` | the base's Phred quality; null when QUAL is `*` (**Fixed in 0.1.4**; it was 255) (example 06) |
| `read_base` | the SEQ base, **complemented on bottom strands** (as a walk reports it) |
| `refr_base` | the anchor's reference symbol if the query determines it (below), else **null** |

`refr_base` (`src/dsl.rs:707-753`, `src/tag_extract.rs:214-219`, `:342-346`):
every assignment of the query's operands that makes its `where` true is tried;
if every such assignment pins the anchor's reference side to the same single
symbol, that symbol (in walk orientation, e.g. `C` for a bottom-strand CG call);
if every one has `=` at the anchor, the read base; otherwise null (also null
beyond 20 operands). Example 05: `CG`/`TG` (`refr = "CG"`) → `C`; `eqG`
(`refr = "="`) → the read base `G`; `mm` (`refr = "/"`) → null, where the
direct run shows the walked base (`C`, `G`, `T`, `A`). The Bismark `U`/`u`
queries (§2.6) give `C`.

### 4.2 Capture layout is always flat

`extract` writes `read_base`/`refr_base`, never `capture_*` lists, even when the
defining queries mark `^` columns (`src/tag_extract.rs:266`,
`src/hits.rs:290-292`).

### 4.3 Differences from `query --parquet` on the same reads

Verified on a read set with a soft clip, an insertion, OB and CTOT reads, a read
with no CG and an unmapped read (example 05):

- **Hit rows agree** in `record_id, name, off_5p, off_3p, refr_pos, qual,
  read_base, refr_base, ref_name, strand` for the queries the tag's characters
  stand for, when those queries determine `refr_base` (`True` in example 05;
  also `src/tag_extract.rs:619-699`, which additionally compares `conv_strand`
  on all four strands and covers a deletion and a hard clip).
- **Only queries with a character** produce rows; the direct run reports every
  query (`spare`, `eqG`, `mm` rows appear only there) (example 05).
- **Hitless rows differ**: `extract` writes a hitless row for a record whose tag
  has no non-fill character, or that lacks the tag; the direct run only for a
  record where *no query at all* fired. Example 05: hitless `record_id`s
  `[3, 4]` from `extract`, `[4]` from `query --parquet`.
- **Row order within a record** differs: `extract` is strictly by strand index;
  the direct run emits a row when the walk reaches the query's *last* column
  (§5.5).
- **No reference**: a record on any contig is decoded (example 05, `chrZ`); no
  off-reference check exists, so `--permissive` has no effect on contigs.
- **Unmapped records are never decoded**: a tag on a FLAG-4 record gives a
  hitless row (example 05, section 8).
- **The tag is trusted, not validated**: a `Z` at a read `T` is written as a `CG`
  row with `read_base = T` (example 05, section 7).
- **Aux text is not checked**: a selected `Z` column with non-UTF-8 bytes is
  null with no error and no count; the direct run stops (or counts under
  `--permissive`) (example 05 section 9 vs example 03 section 1). A damaged aux
  block is not detected either; the tag lookup failing reads as "no tag"
  (`src/tag_extract.rs:304-308`, `:379-382`) (code only for damage).
- Errors that stop the run (exit 1): a character not in the definition
  (`record odd: tag XM has 'u' at SEQ position 2 (0-based), which its definition
  does not include`), a tag length different from SEQ length
  (`record short: tag XM has 3 characters, but the read has 4 bases`), a
  non-string tag (`tag XM is not a string`, code only) (example 05).
- The summary verb is `extracted from` (`src/main.rs:803`).
- A successful run writes `{stem}.manifest.json` like `query --parquet`, with
  `reference` and `walk` null, `extract` naming the tag and where its
  definition came from, the definition's query files (from the header when
  none were given), and `concordance` null in `result` (§5.1, §5.2; checked
  while writing this).

---

## 5. Hit tables: `query --parquet` and `extract`

### 5.1 Files

- **Count**: `max(1, -@) × --shards-per-worker` (`-@ 0` behaves as 1;
  `--shards-per-worker 0` is an error) (`src/main.rs:667`,
  `src/partition.rs:219-231`) (examples 04, 06).
- **Names** (`src/parallel.rs:211-243`): with an extension, `{stem}_{worker}_{shard}.{ext}`
  (`cb.parquet` → `cb_0_0.parquet` ... `cb_3_1.parquet`); without,
  `{name}_{worker}_{shard}` (`one` → `one_0_0`); in `OUT`'s directory, which must
  exist (`worker 0 opening TMP/nodir/x_0_0.parquet ... No such file or
  directory`). Both indices are always present, even for one file.
  The extension is not checked or changed: `--output-format ipc` into
  `s.parquet` writes IPC streams named `s_0_0.parquet`. `OUT = -` is **not**
  stdout: it writes a file literally named `-_0_0` (examples 03, 04, 06).
- **Every file is created and closed**, including files no record was routed to
  (0 rows, valid parquet with the full schema) (`src/parallel.rs:355-377`)
  (example 04: 1 record into 8 files).
- **A failed run removes every file it was set to write**, whether or not the
  file was complete (`src/parallel.rs:513-533`) (examples 03 §1,
  05 §6, 06 §2). This removes a same-named file left by an earlier run too.
  Errors found before any file is created (header, contig and M5 checks, a
  bad `--partition-by`) leave nothing to remove; an earlier run's same-named
  files then survive, but a header or contig error has already removed its
  manifest (below; checked while writing this with a length-mismatch BAM).
- **The run manifest** (new in 0.1.8; `src/manifest.rs:123-151`, `:209-246`,
  wired in `src/main.rs:540-609`, `741-805`, `642-657`). On success, after every data
  file is closed, the run writes `{stem}.manifest.json` beside `OUT`
  (`calls.parquet` → `calls.manifest.json`; `-` → `-.manifest.json`, example 06
  §4), first as `{stem}.manifest.json.partial` and then renamed, so a reader
  never sees half of one. Before the scan starts, any manifest already at that
  path is removed, so **a failed run leaves no manifest**, and one left by an
  earlier run cannot vouch for a partial output (checked while writing this).
  Its presence therefore marks a complete output. Its content is §5.2's
  per-file manifest plus a `result` object: `records_read`, `records_scanned`,
  `records_skipped`, `records_off_reference`, `records_without_seq`,
  `damaged_aux`, `non_utf8_aux`, `concordance` (`{compared, discordant}`, null
  for `extract`) and `rows`, the rows written to each file indexed by slot (the
  `shard` value), in the order of `output.files`. Errors found before the old
  manifest is removed (unreadable reference or BAM, a zero file count, bad
  parquet options) leave an earlier run's manifest and files as they were
  (`src/main.rs:540-604`, code only).
- Existing files of the same name are truncated (`File::create`,
  `src/hit_writer.rs:206`, `:221`); files from an earlier run with a larger file
  count are **not** removed by a later successful run, so a glob can pick up
  stale shards. Since 0.1.8 a successful run names them on stderr, after the
  summary, and the manifest's `output.files` lists this run's files, whose
  footers carry this run's `run_id` (§5.2) (`src/manifest.rs:248-293`, `src/main.rs:642-657`; checked
  while writing this, after a `-@ 2` run and then a `-@ 1` run into the same
  prefix):
  `warning: 1 file(s) beside the output follow its naming but are not part of this run (left by an earlier run that wrote more files?): o_1_0.parquet. pq/o.manifest.json lists this run's files.`
  Only names of the form `{stem}_{digits}_{digits}.{ext}` count; up to five are
  listed, then `...`.
- **Formats** (`src/hit_writer.rs:169-231`): parquet (default) or Arrow IPC
  stream (`--output-format ipc`), same schema and metadata; IPC ignores every
  parquet option's effect. `--row-group-rows 0` and `--batch-rows 0` are
  refused on the command line for either format (`src/cli_query.rs:408-413`,
  `699-705`) (example 06).

### 5.2 Schema

Column order (`src/hits.rs:368-398`): `shard`, `record_id`, the record fields
in selection order (§5.3), `name`, `off_5p`, `off_3p`, `refr_pos`, `qual`, then
the capture columns.

| column | Arrow type | nullable | null when |
|---|---|---|---|
| `shard` | UInt32 | **no** | never |
| `record_id` | UInt64 | **no** | never |
| *record fields* | §5.3 | yes (all) | §5.3 |
| `name` | Utf8 | yes | hitless row |
| `off_5p` | Int64 | yes | hitless row; anchor on a pad (03 §5.2) |
| `off_3p` | Int64 | yes | as `off_5p` (a record with SEQ `*` is not walked, so it only ever has a hitless row) |
| `refr_pos` | Int64 | yes | hitless row; (extract) never otherwise |
| `qual` | UInt8 | yes | hitless row; anchor with no read base (deletion, junction, pad); record with QUAL `*` (since 0.1.4) |
| **flat layout** | | | |
| `read_base` | Utf8 | yes | hitless row |
| `refr_base` | Utf8 | yes | hitless row; (extract) query does not determine it |
| **list layout** | | | |
| `capture_col` | List<Int32> (item nullable) | yes | hitless row (null list, not empty) |
| `capture_read` | List<Utf8> | yes | hitless row |
| `capture_refr` | List<Utf8> | yes | hitless row |
| `capture_qual` | List<UInt8> | yes | hitless row; item null where no read base |
| `capture_off_5p` | List<Int64> | yes | hitless row; item null where no read base |
| `capture_off_3p` | List<Int64> | yes | as `capture_off_5p` |
| `capture_refr_pos` | List<Int64> | yes | hitless row |

(verified with pyarrow for both layouts, example 03.)

- **Layout choice**: flat when every query in the run records only its anchor;
  list as soon as **any** query marks a `^` column, for every row and every file
  of the run (`src/hits.rs:141-162`, `QuerySet::flat_captures`) (example 03).
  `extract` is always flat. Values and meaning of each capture item: 03 §6.2.
- **Metadata** (since 0.1.8; `src/hits.rs:74`, `:448-453`,
  `src/hit_writer.rs:179-195`, `:216-220`, `src/manifest.rs:196-207`). Every
  file carries three plain key-value entries, written before any row:
  `format_version = "7"`, `coordinate_base = "0"` (hit-table coordinates count
  from 0) and `alnbase_manifest`, the run manifest as JSON. In a parquet file
  they are footer key-value metadata beside `ARROW:schema` (example 03 section 8:
  `kv keys ['ARROW:schema', 'alnbase_manifest', 'coordinate_base',
  'format_version']`), so any parquet reader sees them, e.g. DuckDB's
  `parquet_kv_metadata`; the Arrow schema serialized in `ARROW:schema` also
  carries `format_version` and `coordinate_base`. An IPC stream, which has no
  footer, carries all three as Arrow schema metadata (example 03 section 7).
  The manifest's top-level keys, in order (example 03 section 7):
  `alnbase_version, format_version, coordinate_base, run_id, command,
  command_line, input, reference, strand, output, query_files, extract, walk,
  contigs`. `run_id` (time and process id, hex) is the same in every file of a
  run and differs between runs; `command` is `query` or `extract`;
  `command_line` is the argv list; `input` and `reference` are the paths as
  given (`reference` null for `extract`); `strand` names the query file whose
  `[strand.*]` tables the run used, whose full text is in `query_files`; `output`
  holds `format` (`parquet`/`ipc`), `prefix` (OUT as given), `n_workers`,
  `shards_per_worker`, `partition_by`, `fields` (record columns in file order),
  `only_hits`, `capture_layout` (`flat`/`lists`) and `files` (file names beside
  OUT, indexed by slot); `query_files` holds each query file's `path` and full
  `text`; `extract` is `{tag, definition}` for `extract`, else null; `walk` is
  null for `extract` and for `query` holds the resolved `insertions`,
  `end_context`, `end_context_given`, `splice_context`, `splice_context_given`,
  `widest_span`, `max_discordance` and `require_m5`; `contigs` lists every `@SQ`
  in header order with `name`, `length`, `md5` (the index's MD5 for a contig it
  has, else the header's `M5`, else null) and `in_reference` (null for
  `extract`) (`src/manifest.rs:38-121`, `:153-184`). The per-file copy
  describes the run's settings only; the outcome (`result`) is only in
  `{stem}.manifest.json` (§5.1). Before 0.1.8 the only metadata was
  `format_version`, visible in parquet only inside `ARROW:schema`, and no query
  text, command line or alnbase version was stored in hit files.
- **Schema identity across files**: a function of the resolved field list and the
  layout only, so all files of a run share it (example 04: 8 schemas identical).

Values of the hit columns (anchor-based `off_5p`/`off_3p`/`refr_pos`/`qual`,
null offsets on pads, symbols `_`, `.`, `,`): see 03 §3, §5.2 and §6.2.
`off_5p` / `off_3p` are measured from the ends of the read **as sequenced** for
every FLAG, read 2 included (03 §3.2). Facts that belong to the output and bite
pipelines:

- `read_base` and `refr_base` on **bottom-strand** reads are **complemented**
  relative to SEQ and the reference: a protected CG on an OB read is reported
  `read_base = C`, `refr_base = C`, not `G` (example 03, `rev`).
- `refr_pos` of rows within a bottom-strand record **decreases** as `off_5p`
  increases (example 03). For read 2 it is the other way round: a
  forward-strand R2 (FLAG 129, bottom walk) has `off_5p` rising with
  `refr_pos` (03 example 01).
- An anchor on a flank column gives `off_5p`, `off_3p` and `qual` null,
  `read_base = "_"`, and its real `refr_pos` (03 example 06 `ends_on_C`).
- **QUAL `*`**: htslib stores 0xFF per base; since 0.1.4 alnbase treats that as
  no quality, so `qual` is null on every row of such a record in both
  `query --parquet` and `extract`, and the `qual_phred` record field is null
  (example 06). A `qual IS NULL` test therefore no longer means "no read base".
  **Fixed in 0.1.4**: up to 0.1.3 `qual` was 255 and `qual_phred` was `seq_len`
  copies of U+00FF `ÿ`.
- **Hard clips**: `off_5p` / `off_3p` count hard-clipped bases (since 0.1.4), so
  a supplementary `60H20M` reports `off_5p = 60` for its first base, as its
  primary does; `off_5p + off_3p = read_len + hard_clip_5p + hard_clip_3p - 1`
  (03 §3.2).

### 5.3 Record fields (`-f/--field`, `-F/--add-field`)

Selection (`src/record_field.rs:559-653`): comma-separated and/or repeated.
`-f` **replaces** the default `core`; `-F` **appends** to whichever base set
applies. Order is as given; groups expand in place; duplicates are dropped **by
column name**, first spelling wins, silently (`-F NM:i,NM:Z` keeps `NM:i`)
(example 03). Names are case- and hyphen-insensitive; aux tags and type codes
are case-sensitive.

Named columns (`src/record_field.rs:151-185`, types `:311-348`, extraction
`src/record_columns.rs:217-325`):

| column | type | value; null rule |
|---|---|---|
| `qname` | Utf8 | QNAME; invalid UTF-8 replaced lossily (U+FFFD) |
| `flags` | UInt16 | FLAG |
| `tid` | Int32 | header index; null when < 0 |
| `ref_name` | Utf8 | RNAME; null when tid < 0 |
| `pos` | Int64 | 0-based POS; **-1 (not null)** for an unplaced record (example 03 probe) |
| `end_pos` | Int64 | 0-based exclusive alignment end |
| `unclipped_start` | Int64 | `pos` minus leading S+H |
| `unclipped_end` | Int64 | `end_pos` plus trailing S+H |
| `mapq` | UInt8 | MAPQ |
| `cigar` | Utf8 | CIGAR text (`*` records give an empty string, code only) |
| `read_len` | Int64 | SEQ length (hard clips excluded) |
| `hard_clip_5p` | Int64 | hard-clipped bases at the read's sequenced 5' end: the CIGAR's leading `H` for a forward record, its trailing `H` for a reverse one; 0 when there are none (also for unmapped records); already included in `off_5p` (added in 0.1.4) |
| `hard_clip_3p` | Int64 | the same at the sequenced 3' end; already included in `off_3p` (added in 0.1.4) |
| `soft_clip_5p` | Int64 | soft-clipped bases at the read's sequenced 5' end: the CIGAR's first operation other than `H` if it is `S`, on a forward record; its last on a reverse one; 0 when there are none. Already included in `off_5p`, so `off_5p - soft_clip_5p - hard_clip_5p` counts from the first aligned base (added in 0.1.9) |
| `soft_clip_3p` | Int64 | the same at the sequenced 3' end; already included in `off_3p` (added in 0.1.9) |
| `is_paired` … `is_qc_fail` | Boolean | one FLAG bit each: `is_paired, is_proper_pair, is_unmapped, is_mate_unmapped, is_reverse, is_mate_reverse, is_first_in_template, is_last_in_template, is_secondary, is_supplementary, is_duplicate, is_qc_fail` |
| `strand` | Utf8 | the conversion strand the run's strand rule named, `+` or `-`; **not** FLAG 0x10; under `queries/strand/directional.toml` that is `+` when `is_last_in_template == is_reverse`. Also computed for unmapped records (`+` for FLAG 4) (example 03) |
| `conv_strand` | Utf8 | strand of origin `OT`, `OB`, `CTOT` or `CTOB`, as the run's strand rule named it (`src/record_columns.rs:275-283`): OB and CTOB are exactly the records whose `strand` is `-`, so it splits each `strand` value by mate. **Null** when the rule names only the conversion strand, which several aligners' output cannot refine — an invented OT/CTOT split would read as a measurement. Otherwise, like `strand`, filled for every record, unmapped ones included (code only). In `all`, not in `core` (added in 0.1.8; null case in 0.1.15; example 03 `all` schema; unit tests `src/record_columns.rs:900-915`, `src/tag_extract.rs:694-698`) |
| `read_reverse` | Boolean | whether the read as sequenced runs backwards along the reference, as the run's strand rule concluded, which is the direction `off_5p` and `off_3p` count from. **Not** FLAG 0x10 — `is_reverse` is that bit, and the two part company for any aligner that writes something else in it (BSBolt writes the converted reference strand). Under `queries/strand/directional.toml` the two agree on every record. In `all` and `bools`, not in `core` (added in 0.1.15; unit tests `src/record_columns.rs:918-955`) |
| `insert_size` | Int64 | TLEN |
| `mate_tid` | Int32 | null when RNEXT unset |
| `mate_ref_name` | Utf8 | null when RNEXT unset |
| `mate_pos` | Int64 | null when RNEXT unset |
| `seq_ascii` | Utf8 | SEQ letters |
| `qual_phred` | Utf8 | QUAL as phred+33; null when QUAL is `*` (since 0.1.4, see above) |

Groups (`src/record_field.rs:204-249`, `:508-514`): `core` = `ref_name,strand`
(the default); `all` = every named column above in that order; `bools` = the 12
flag-bit columns; `mate` = `insert_size,mate_tid,mate_ref_name,mate_pos,is_mate_unmapped,is_mate_reverse`;
`coords` = `tid,ref_name,pos,end_pos,unclipped_start,unclipped_end`.

Aliases (`src/record_field.rs:517-537`): `rname`, `refname` → `ref_name`;
`rnext`, `mrnm` → `mate_ref_name`; `pnext`, `mpos` → `mate_pos`; `tlen`,
`isize` → `insert_size`; `seq` → `seq_ascii`; **`qual` → `qual_phred`** (the
record column, not the hit column `qual`); `mtid` → `mate_tid`; `is_qcfail` →
`is_qc_fail`; `is_ot`, `refr_strand` → `strand`; `is_read1` →
`is_first_in_template`; `is_read2` → `is_last_in_template` (example 03,
section 5).

**Aux tags**: `TAG:type`, TAG two ASCII alphanumerics (a leading digit is
accepted here), type one of SAM's codes (`src/record_field.rs:569-584`,
`:63-70`). The column is named by the tag alone (`NM`, not `NM:i`).

| requested type | Arrow | stored `c C s S i I` | stored `f` / `d` | stored `A` | stored `Z` / `H` | stored `B` | absent |
|---|---|---|---|---|---|---|---|
| `i c C s S I` | Int64 | value | null | null | null | null | null |
| `f d` | Float64 | widened to float | value (f32 widened) | null | null | null | null |
| `Z A H B` | Utf8 | null | null | the char | text (null if not UTF-8) | elements comma-joined, `Display` formatting (`-3,0,7`, `1.5,-0.25`) | null |

(`src/record_columns.rs:327-386`; verified for `NM:i`, `XA:A`, `XC:i:-5`
stored as `c`, `XI:i:3000000000` stored as `I`, `XF:f`, `XH:H`, `XB:B:c`,
`BF:B:f`, and the mismatches `NM:Z`, `XA:i`, `XC:f`→-5.0, `XI:f`→3000000000.0,
`XF:i`, `XB:i`, `XH:i`, `XZ:i` — example 03.) A repeated tag reads its first
occurrence (`src/record_columns.rs:595-597`). An empty `B` array gives `""`.

**Trouble in aux blocks** (`query --parquet` only; `src/scanner.rs:352-365`):
aux columns are filled by one walk over the block (`src/record_columns.rs:578-612`).
If a selected aux tag's `Z`/`H` text is not UTF-8, or the block is damaged, the
run **stops** at that record
(`record badtext: text in an aux tag that is not UTF-8. Pass --permissive ...`)
unless `--permissive`, where the value is null and the record counted
(`N records had text in an aux tag that is not UTF-8, read as null`)
(example 03). Unselected tags are never checked. Unwalkable records are never
checked (`src/scanner.rs:382-398`). `extract`: never checked (§4.3).

**Encodings** (defaults; `src/record_field.rs:379-393`, `src/hits.rs:413-446`):
`qname` delta-text; `seq_ascii`, `qual_phred`, float aux tags plain; `pos`,
`end_pos`, `unclipped_*`, `mate_pos`, `insert_size`, `record_id`, `off_5p`,
`off_3p`, `refr_pos`, `capture_off_*`, `capture_refr_pos` delta; everything else
(including `shard`, `name`, `qual`, `read_base`, `refr_base`, text/int aux,
flags, bools, `strand`, `conv_strand`, `cigar`, `mapq`, `tid`) dictionary. Verified on disk
(example 03, section 8). Override with `--column-encoding COL=dict|plain|delta|delta-text`
(repeatable, comma-separable); the column name must exist in this run's schema
and `delta` needs an integer, `delta-text` a text column, else the run fails
before any file is written (`src/hit_writer.rs:106-144`). A list column's
override applies to its items.

### 5.4 Rows

- **One row per query firing** (direct) or per non-fill tag character
  (extract). Every row repeats the record's selected fields.
- **Hitless rows** (`src/scanner.rs:369-371`, `:382-398`,
  `src/tag_extract.rs:359-361`, `:367-376`): with the default (no
  `--only-hits`), a record that produced no row gets exactly one row: `shard`
  and `record_id` set, record fields filled (e.g. `ref_name` present for a
  mapped hitless record), every hit column null. This covers:
  1. walkable records where nothing fired (direct) / no tag character (extract);
  2. unmapped records (both);
  3. off-reference records under `--permissive` (direct; without it the run
     stops) (`src/parallel.rs:424-440`);
  4. mapped records with SEQ `*` (direct; 03 §7).
  Verified in examples 03, 04 (all 116 unmapped records of 1500 have one null
  row), 05.
- `--only-hits` drops all of them, including unmapped ones
  (`src/parallel.rs:746-766`) (example 03). `record_id` values keep their input
  ordinals, so they are then not dense.
- **`record_id`** (`src/parallel.rs:441-444`): assigned by the single reader,
  before routing, to every record read (walkable or not), starting at 0, in
  input order. With the default hitless rows the set of `record_id`s across all
  files is exactly `0..records_read` and equals the input ordinal (example 04:
  `qname rNNNN` ↔ `record_id NNNN` for 1500 records across 8 files).
- **Not deduplicated or filtered**: secondary, supplementary, duplicate, QC-fail
  and low-MAPQ records produce rows like any other.

### 5.5 Order

- **Within a file**: rows of one record are contiguous; records appear in
  increasing `record_id` (the reader sends each file's records in input order
  over one channel per worker, and a worker processes its channel in order,
  `src/parallel.rs:378-389`, `:441-462`) (example 04: non-decreasing in all 8
  files). A record's rows are never split across two batches or row groups
  (`src/scanner.rs:373-376`).
- **Within a record** (direct): in the column order of the walk (along the
  conversion strand, 03 §2), a row being emitted **when the walk stores the
  query's last column**; at one column, queries fire in declaration order
  (`src/scanner.rs:263-288`). Hence rows are **not** necessarily sorted by
  `off_5p`: example 03 section 9 gives `[('narrow', 4, 4), ('wide', 3, 3)]` for a
  1-column query at column 4 and a 3-column query anchored at column 3.
  Within `extract`, strictly by strand index.
- **Across files**: no order; read the dataset and sort by
  `record_id` (and `off_5p`) if order matters.
- A coordinate-sorted input gives files whose records are coordinate-sorted by
  `pos`; `refr_pos` of rows is not monotonic (bottom strands, flank columns,
  the rule above).

### 5.6 Sharding and partitioning

Routing (`src/partition.rs:213-267`, `src/parallel.rs:448`):

```
h      = FNV-1a over the key's extracted cells (in key order)
slot   = h % (workers * shards_per_worker)
worker = slot / shards_per_worker
shard  = slot % shards_per_worker        → file {stem}_{worker}_{shard}.{ext}
```

- The digest per field is: a type byte (null, bool, u8, u16, i32, i64, f64,
  str distinguished), then the value little-endian, text length-prefixed
  (`src/record_columns.rs:68-79`, `:428-446`). So an absent aux tag and an empty
  one route differently, and all records lacking the tag route together.
- **`--partition-by`** accepts exactly the `--field` vocabulary: any named
  column, alias, group, or `TAG:type` aux tag (`src/cli_query.rs:461-474`,
  `:503-509`). It **replaces** the default `qname` (not added to it), and is
  independent of the selected columns: partitioning by `CB:Z` does not put `CB`
  in the file (add `-F CB:Z`). The value hashed is the value the column of that
  name would contain (same extractor), so `CB:i` on a `Z` tag hashes null for
  every record and funnels everything into one file (code only, follows from
  the table in §5.3).
- **Unwalkable records** are hashed like any other and land with their key
  (`src/parallel.rs:445-448`).
- **Warning** (`src/partition.rs:179-204`, printed once, `src/main.rs:672-675`):
  any key other than exactly `qname` prints
  `warning: partitioning by CB rather than qname: a fragment's alignments can
  land in different files, so mate-overlap resolution cannot see both sides of
  one fragment within a file`; a key containing `qname` plus more prints
  `... stay in one file only if X is constant across the fragment`
  (examples 04, 07). `extract` prints it too.
- **`shard` column** = `worker * shards_per_worker + shard` = the file's slot,
  constant within a file (`src/hits.rs:104-110`) (example 04: all 8 files agree
  with their names).
- **Stability across runs**: FNV-1a is deterministic across machines and runs;
  a key keeps its `shard` value when `workers × shards_per_worker` is unchanged,
  and its file name only when both numbers are unchanged
  (`src/partition.rs:25-31`) (code only).
- **Memory** is sized by `--batch-rows`, which is a **total across all output
  files**: each file buffers up to its share, `--batch-rows` divided by
  `max(1, -@) × --shards-per-worker`, rounded up and at least 1
  (`src/cli_query.rs:523-527`, `src/main.rs:598`, `794`). Adding files does
  not multiply the buffered rows (it did before 0.1.2). The module comment at
  `src/parallel.rs:52-55` and the `--shards-per-worker` help
  (`src/cli_query.rs:485-487`) still describe the old per-file behaviour (D19).

### 5.7 Determinism across thread and file counts

- `query --parquet` with the default key, `-@ 1` vs `-@ 3` (2 files per worker):
  the multiset of rows, all columns except `shard`, is identical (2147 rows)
  (example 04). This holds for any key, because each record's rows depend only
  on the record (code only beyond the verified case). The files' metadata is not
  identical across runs: every run has its own `run_id`, and the manifest
  records the command line, `n_workers` and the file list (§5.2).
- Tagged BAM, `-@ 1` vs `-@ 4`: records identical, header identical except the
  `@PG CL` (which records the `-@` argument), so the files are **not**
  byte-identical (example 04).

### 5.8 Writer options

| option | default | effect (source) |
|---|---|---|
| `--batch-rows N` | 500000 | rows held in memory **in total across all output files** before they are written. Each file's builder is flushed once it holds ≥ its share, `ceil(N / (max(1, -@) × --shards-per-worker))` and at least 1, checked after each record (`src/cli_query.rs:523-527`, `src/hit_writer.rs:252-258`). Must be at least 1: `0` is a command-line error, exit 2 (`invalid value '0' for '--batch-rows <BATCH_ROWS>': must be at least 1`) (example 06). Not a row-group size: parquet accumulates flushed batches up to `--row-group-rows` |
| `--row-group-rows N` | 50000 | max rows per parquet row group (`src/cli_query.rs:411-413`); `2` gave groups `[2, 2, 2, 2]` (example 03). Must be at least 1: `0` is a command-line error, exit 2 (example 06); before 0.1.2 it panicked a worker |
| `--output-profile fast\|small` | fast | fast = LZ4_RAW + chunk statistics; small = ZSTD(`--compression-level`) + page statistics (`src/cli_query.rs:529-559`) (example 03: `ZSTD`) |
| `--compression none\|lz4\|snappy\|zstd` | profile | overrides the codec only |
| `--compression-level N` | 3 | validated (1-22) only when zstd is in force |
| `--statistics none\|chunk\|page` | profile | overrides statistics only |
| `--column-encoding COL=ENC` | per column | §5.3 |
| `--output-format parquet\|ipc` | parquet | IPC stream with end-of-stream marker written on close (`src/hit_writer.rs:294-308`) |

---

## 6. Run summaries, diagnostics, exit codes

**Summary** (stderr, after a successful scan only; `src/main.rs:687-732`):

```
read R records, <verb> S, skipped K[ (written untagged) | (given a null row)]
  U hits were anchored on a pad or a soft-clipped base, so had no aligned base to tag
  N records had text in an aux tag that is not UTF-8, read as null
  D records had a damaged aux block; tags after the damage could not be read
  of those, O were on contigs the reference does not have
  X of C compared read bases differ from the reference (P%; a reference C read as T is not counted)
  of those, Q were mapped but had no SEQ ('*')
```

- verb: `tagged` (BAM), `scanned` (`query --parquet`), `extracted from`
  (`extract`). `S` counts records walked (BAM, direct) or records seen by a sink's
  `scan` (extract, i.e. mapped records); `K` counts unmapped, off-reference
  (permissive) and, for `query`, mapped SEQ `*` records (`src/batch.rs:60-86`).
- The first line always prints (with `skipped 0`). When `K > 0` the suffix is
  `(written untagged)` for a tagged BAM and `(given a null row)` for a parquet
  or IPC output that keeps hitless rows; an `--only-hits` run, where skipped
  records get nothing, prints no suffix (`src/main.rs:692-697`,
  `src/parallel.rs:414`) (examples 03-06). Before 0.1.2 every output printed
  `(written untagged)`.
- The indented lines print only when non-zero. `of those, O ...` and
  `of those, Q ...` refer to `K`, though they print after the aux lines, and
  `of those, Q ...` after the concordance line too (`src/main.rs:717-731`;
  03 example 05; the order of the last two is code only, no example has both).
- `unplaced` is only ever non-zero for the BAM output.
- The concordance line (since 0.1.4, `query` only) reports the reference
  concordance check below. It prints when the run compared at least one base;
  `extract` reads no reference and prints no such line (examples 01-08; a run
  that scanned 0 records prints none, 03 example 05).
- After the summary, a successful `query --parquet` or `extract` run may warn
  about files from an earlier run that follow the output's naming (§5.1).

**Reference concordance check** (`--max-discordance FRACTION`, default `0.25`,
since 0.1.4). Every walked record counts, over its observed columns where the
read and the reference each show one unambiguous base, how many disagree; a
reference C read as T is not counted, so bisulfite conversion does not raise the
rate (`src/scanner.rs:103-155`). The tagged BAM and `query --parquet` stop on it,
and since 0.1.8 `--trace-records` does too (below); `extract` and `--trace`
compare nothing. Correct alignments disagree at a few percent at most (example
04: 4.58% over 16,608 bases); a reference shifted by one base stops the run at
about 63% (CHANGELOG), 68.2-68.6% on the random 20 kb reference used to check
this section. The value must be between 0 and 1 (`1.5` is refused by the
argument parser, exit 2); `--max-discordance 1` turns the check off, and the
summary line still reports the rate.

- Since 0.1.6 one count is shared by every worker (tagging thread, or
  parquet/IPC worker) (`SharedConcordance`, `src/scanner.rs:158-194`). After
  each record is added, the run stops with exit 1 once the shared total has
  compared **100,000 bases** at a rate above the threshold
  (`Concordance::MIN_COMPARED`), for any `-@`. Checked while writing this, on
  250,000 bases against a shifted reference: parquet `-@ 1`
  `Error: 68.3% of the 100000 read bases compared so far disagree with the reference (a reference C read as T is not counted), above --max-discordance 0.25. A wrong reference, a BAM aligned to a different assembly, or shifted coordinates cause this; pass --max-discordance 1 to turn the check off.`,
  parquet `-@ 4` `68.2% of the 100200 read bases compared so far ...` (the
  workers add whole records concurrently, so the count can pass 100,000), tagged
  BAM `-@ 4` `68.4% of the 100000 ...`.
- A run that ends before that is judged on its total, from **1,000 compared
  bases** (`Concordance::FINAL_MIN_COMPARED`; `src/parallel.rs:503-509`,
  `src/ordered.rs:322-326`): 10 reads of the same input (1,000 compared bases)
  stop with `Error: 68.6% of the 1000 read bases compared in total disagree ...`. A run
  comparing fewer than 1,000 bases is not judged, whatever its rate; read its
  summary line. (Before 0.1.6 each worker counted alone, so a small or thinly
  spread input finished with exit 0: N13.)
- The failure is handled like any other: a tagged BAM prints
  `removed the incomplete output PATH` and is removed, and a parquet/IPC run
  removes every file it was set to write and leaves no manifest (below; checked
  while writing this).
- **`--trace-records N`** (since 0.1.8; `src/main.rs:423-511`) prints, after
  each traced record, `reference concordance: X of Y compared read bases differ
  (P%; a reference C read as T is not counted)` (01 example 06). After the last
  record it judges the traced total like a finished run: at 1,000 or more
  compared bases above the threshold it exits 1 with `... of the 1200 read bases
  compared in the traced records disagree ...` (checked while writing this, 12
  records of 100 bases); on fewer bases above the threshold it prints
  `warning: 37.5% of the 8 read bases traced disagree with the reference, above --max-discordance 0.25; too few to judge (a scan judges from 1000 bases), but a wrong reference, a BAM aligned to another assembly, or shifted coordinates look like this.`
  and exits 0 (checked while writing this, on 01 example 06's `no_m5.bam`
  against `other.aref`). The traces are printed before either.

**`--require-m5`** (since 0.1.8; `query` only, for the tagged BAM,
`--parquet` and `--trace-records`): the first record on a contig the index has
whose `@SQ` line carries no `M5` stops the run with exit 1, whether or not
`--permissive` is given:
``Error: record on_chr10 is on contig chr10, whose @SQ line has no M5 checksum, and --require-m5 was given. Add M5 to the header's @SQ lines (`samtools dict REF.fa` prints them for the FASTA the BAM was aligned against; `samtools reheader` applies an edited header), or leave out --require-m5 to check names and lengths only.``
Contigs no record is on are not checked; records on contigs the index lacks
are handled by the off-reference rule, not by this one (`src/contig_map.rs:149-170`,
`src/scanner.rs:563-564`, `src/batch.rs:45-53`, `:75-80`, `src/main.rs:436-438`,
`:451-456`). A failed run removes its output as any other does (checked while
writing this, on 01 example 06's inputs). Whether it was given is recorded in the
manifest's `walk.require_m5` (§5.2).

**Other stderr lines**: parse warnings (`warning: in query ...`,
`warning: FILE: ...`), unused-query warning (BAM), contig warning
(`src/ordered.rs:152-154`, `src/parallel.rs:321-323`), fragment warning,
`extracting tag ...`, `removed the incomplete output PATH` (tagged BAM only;
a failed parquet run removes its files without a message), the stale-file
warning `warning: N file(s) beside the output follow its naming but are not part
of this run ...` after a successful parquet/IPC/`extract` run (§5.1,
`src/main.rs:642-657`), and the `--trace-records` concordance warning (above).

**Exit codes**: 0 success; **1** for every alnbase error (`main` returns
`anyhow::Result`; message printed as `Error: <msg>` with a `Caused by:` chain
when present), including misplaced options, thread-start failures and worker
panics; **2** for command-line syntax errors rejected by clap (unknown flag
such as the removed `--soft-clips`, bad value, `--batch-rows 0`,
`--row-group-rows 0`) (`src/main.rs:517-518`) (example 02, section 8; example 06
§4).

**Thread limits**: when the operating system refuses a thread (usually the
per-user process/thread limit), the run ends with exit 1 rather than a panic
(`src/batch.rs:93-99`; threads started at `src/ordered.rs:171-201`, `213-245`,
`src/parallel.rs:346-403`):

```
Error: could not start the worker 3 thread: <OS error>. The system refused a new thread; the per-user process/thread limit (`ulimit -u`) is the usual cause, so lower --threads or raise the limit.
```

`worker N` is `writer` or `tagger N` for a tagged BAM. A tagged-BAM run removes
its partial output and a parquet run removes its files, as for any other
failure (code only; `alnbase-validation/demos/resource-limits` describes the limits).

**Failure and partial outputs**:

- Tagged BAM: on any error the output **file is removed**
  (`removed the incomplete output ...`) unless it is `-` or not a regular file
  (FIFO, device) (`src/ordered.rs:328-332`, `:341-353`) (examples 01, 02, 06).
  The reported error is the first that is not a knock-on "stopped" error
  (`src/ordered.rs:330`).
- Parquet/IPC: on any error after the output files are set up, **every output
  file the run was configured to write is removed** (`src/parallel.rs:513-533`),
  so a glob over `OUT_*` finds nothing from the failed run (examples 03 §1,
  05 §6, 06 §2). No message names the removed files. No manifest is written,
  and one an earlier run left at `{stem}.manifest.json` was already removed
  when the scan started (§5.1). Before 0.1.2 nothing was
  removed and a failed run left unreadable files next to valid partial shards.
- Parquet/IPC error reporting: every worker is joined before the run returns,
  and a worker's own error wins over the reader's (`src/parallel.rs:474-501`).
  Example 06 §2: 80000 records with a bad record at index 10 report
  `record r10: tag XM has 'u' at SEQ position 2 (0-based), ...`, the same as 20 records.
  Before 0.1.2 the larger input reported
  `worker 0 stopped early (sending on a closed channel)`. A missing output
  directory is reported as `worker 0 opening PATH` caused by
  `No such file or directory` (example 06 §4).

---

## 7. Guarantees

Each item states the scope in which it holds.

1. **G1 — Tagged BAM is the whole input, in order.** Every record read is
   written exactly once, in input order, for any `-@`; fields 1-11 are unchanged
   (only aux fields are added/replaced). (examples 01, 02, 04; `src/ordered.rs:174-190`)
2. **G2 — Unwalkable records are untouched.** Unmapped records, and under
   `--permissive` off-reference records, are passed to the writer without any
   modification and get no tags (verified as unchanged SAM fields and no aux
   added; byte-level BAM encoding is htslib's). (examples 01, 02, 04)
3. **G3 — Every walkable record gets every declared tag**, `bases` tags of
   length `seq_len` in SEQ order and `strand` tags, all type `Z`; or the run
   fails. (examples 01, 02)
4. **G4 — Bismark equivalence.** `docs/bismark-xm.toml` plus the XR/XG tables of
   §2.6 produce Bismark's `XM`, `XR`, `XG` for OT, OB, CTOT and CTOB reads under
   the documented semantics, for ungapped alignments. (example 01)
5. **G5 — No silent tag loss or overlap.** A conflicting pair of codes, an
   anchor without a read base, an existing tag (without `--overwrite-tags`), a
   damaged aux block (without `--permissive`), an off-reference record
   (without `--permissive`) or, under `--require-m5`, a record on a contig whose
   `@SQ` has no `M5` stops the run rather than writing a partial answer;
   so does a reference concordance rate above `--max-discordance` once the run
   has compared 100,000 bases, or at its end from 1,000 compared bases (§6);
   hits with no SEQ position are counted and reported. (example 02)
6. **G6 — Failed runs leave no output file.** A tagged BAM at a regular-file
   path is removed; a parquet/IPC run removes every file it was set to write,
   including files from an earlier run with the same names, and leaves no
   `{stem}.manifest.json`. (examples 01, 02, 03, 05, 06; the manifest part
   checked while writing this, `src/main.rs:604`, `:801`)
7. **G7 — Query files travel with the BAM.** Every query file of a tagged-BAM run
   is stored under that run's `@PG` ID, earlier runs' files are kept, and
   `dump-query` returns each file's text byte-identically after two comment
   lines; loss or corruption of stored lines is an error. (examples 01, 08)
8. **G8 — `extract` uses each tag's newest definition** (tag by tag), or exactly
   the given `--query-file`s. (example 08)
9. **G9 — `extract` agrees with a direct run** on the rows of the tag's queries,
   given the same strand rule and walk defaults, for queries that determine
   `refr_base`; otherwise `refr_base` is null. (example 05)
10. **G10 — One schema per run.** All files of a run have the identical schema
    (a function of the field list and capture layout), including empty files;
    `shard` and `record_id` are never null. (examples 03, 04)
11. **G11 — Partition key locality.** Within one run, all records with the same
    `--partition-by` value (null included) are in exactly one file, whatever
    `-@` and `--shards-per-worker`. (example 04: 7 CB values plus null over 8
    files)
12. **G12 — `shard` is the file's slot** `worker × shards_per_worker + shard`.
    (example 04)
13. **G13 — `record_id` is the input ordinal**, assigned before routing, so it is
    dense over `0..records_read` across the dataset without `--only-hits`, and
    non-decreasing within every file. (example 04)
14. **G14 — Inventory.** Without `--only-hits`, every record read (mapped or not)
    has at least one row. (examples 03, 04, 05)
15. **G15 — Row content is independent of thread and file count** (only `shard`
    and file placement change). (example 04)
16. **G16 — All files exist after success**, each a complete parquet file (or
    IPC stream with end marker), including files that received no records.
    (example 04)
17. **G17 — Rows of one record are contiguous in one file** and never straddle a
    batch or row group. (`src/scanner.rs:373-376`, code only for the row-group
    part; contiguity by construction)
18. **G18 — A hit table describes itself** (since 0.1.8). Every file of a
    `query --parquet` or `extract` run carries `format_version`,
    `coordinate_base` and the run manifest (same `run_id` in every file) as
    parquet footer key-value metadata, or Arrow schema metadata for IPC; a
    successful run then writes `{stem}.manifest.json`, whose `output.files`
    lists exactly the files the run wrote and whose `result.rows` gives each
    one's row count. (example 03 sections 7-8; example 06 §4;
    `src/manifest.rs`, `src/hit_writer.rs:179-195`, `:216-220`)

## 8. Non-guarantees (do not rely on these)

1. **N1** Tagged BAM files are not byte-identical across `-@` values (the `@PG CL`
   differs, and BGZF block layout is not promised). (example 04)
2. **N2** No order of rows across files; within a record, rows are not sorted by
   `off_5p` (direct run). (example 03)
3. **N3** A key keeps its file across runs only if `-@` and
   `--shards-per-worker` are both unchanged. (`src/partition.rs:25-31`)
4. **N4** Stale shards from an earlier run that wrote more files (a larger
   `-@` or `--shards-per-worker`) are not removed, by a successful or a failed
   run, so a bare `OUT_*` glob can still mix runs. (Fixed in 0.1.2: a failed run
   no longer leaves its own partial files. **Partly fixed in 0.1.8**: a
   successful run warns about them, and they are told apart by `run_id` and by
   the manifest's `output.files`, §5.1.)
5. **N5** ~~The error message of a failed parquet run may be a knock-on
   `worker N stopped early` instead of the cause.~~ Fixed in 0.1.2: the
   worker's own error is reported. (example 06)
6. **N6** ~~`qual` is not null for records with QUAL `*` (255), and `qual_phred` is
   not empty.~~ Fixed in 0.1.4: both are null. (example 06)
7. **N7** `extract` does not validate tags against the reads or the reference,
   does not check aux text or damage, and does not decode unmapped records.
   (example 05)
8. **N8** ~~The parquet footer carries no query definitions, command line or tool
   version; only `format_version` inside the Arrow schema. Keep the tagged BAM
   or the query file alongside.~~ **Fixed in 0.1.8:** every file carries the
   run manifest (query file text, command line, alnbase version, walk settings,
   contigs) beside `format_version` and `coordinate_base` (§5.2). (example 03)
9. **N9** `pos` is not null for unplaced records (-1), and `strand` is not null
   for unmapped records. (example 03)
10. **N10** No filtering of secondary, supplementary, duplicate, QC-fail or
    low-MAPQ records in any output. (README; `src/scanner.rs:244-291`)
11. **N11** `--partition-by` other than `qname` does not keep fragments together.
    (`src/partition.rs:179-204`)
12. **N12** A corrupt stored run makes `extract` fail for every tag, even tags
    defined by an intact run. (example 08)
13. ~~**N13** The concordance check does not catch a wrong reference on a small
    input or a thinly spread one: a worker that compares fewer than 100,000
    bases never stops the run.~~ **Fixed in 0.1.6:** the count is shared across
    workers and the finished run is judged on its total (from 1,000 compared
    bases). Runs comparing fewer than 1,000 bases are not judged. (§6)

---

## 9. Idioms

**Carry a cell barcode into rows and partition by it** (example 07):

```bash
alnbase query --parquet -@ 8 --shards-per-worker 4 \
  -F CB:Z,mapq,is_duplicate --partition-by CB:Z --only-hits \
  --query-file calls.toml in.bam ref.aref calls.parquet
```

`-F` keeps `ref_name,strand` and adds the columns; `--partition-by` needs the
column added separately if you want it in the rows. Each barcode's rows are in
one file (per-cell processing can open one file), but a file holds many
barcodes; records without `CB` all share one file. Expect the fragment warning.
If mates will be overlap-resolved per file, keep the default `qname` key (or run
`alnbase overlap` upstream) instead.

**Carry MAPQ/flags and filter downstream** (example 07):

```sql
select CB, split_part(name, '_', 1) ctx, refr_pos,
       count(*) filter (where name in ('CG_meth', 'CHG_meth')) meth, count(*) calls
from read_parquet('calls_*.parquet')
where mapq >= 20 and not is_duplicate
group by all
```

Alternatively filter the BAM first (`samtools view -F 0xD04 -q 20`) and select
no extra columns. For `extract`, the same `-F mapq,flags` options apply.
(Beware SQL `LIKE '%_meth'`: `_` is a wildcard and also matches `_unmeth`; the
first draft of example 07 got this wrong.)

**Context by query name vs by captures** (example 07):

- *By name*: one query per context and state (`CG_meth`, `CHG_unmeth`, …, or
  Bismark's `z Z x X h H u U`). Flat layout, one low-cardinality `name` column,
  and the same file works for a `bases` tag and `extract`. Contexts must not
  overlap if they share a tag (§2.2).
- *By capture*: one query `read = "Y~~"`, `refr = "C~~"`, `mark = "+^^"`;
  every file switches to list columns, and the context is computed downstream
  from `capture_refr` (`[C, G, T]` → CG). Needed when a context does not fit
  one character per base; not usable through a tag (a tag records the anchor
  only).

**Read the definition back from a table** (since 0.1.8): every file's
`alnbase_manifest` footer entry, and `{stem}.manifest.json`, hold each query
file's text, so the queries behind a `name` travel with the rows:

```sql
select decode(value) from parquet_kv_metadata('calls_0_0.parquet')
where decode(key) = 'alnbase_manifest'
```

(`json_extract` on that value gives `query_files`, `command_line`, `walk`.)
Open the files the manifest lists (`output.files`) rather than a bare glob when
the prefix may hold an earlier run's files, and treat a missing
`{stem}.manifest.json` as an incomplete output. A tagged BAM still stores the
query files in its header (§3).

---

## 10. Five facts pipeline authors most often get wrong

1. Bottom-strand rows report **complemented** bases (`read_base = C` for a
   protected CG on an OB read), while `bases` tags are in plain **SEQ order** with
   the call on the G. (examples 01, 03)
2. A failed `--parquet`/`extract` run removes its files, but files from an
   **earlier** run with more shards are left alone, so a glob can still mix
   runs. Since 0.1.8 a successful run warns about them and lists its own files
   in `{stem}.manifest.json`, which only a successful run leaves. Check the exit
   status and the manifest, read the files it lists, or write each run to a
   fresh prefix or directory. (examples 04, 06; §5.1)
3. `extract` and `query --parquet` differ in **hitless rows** (tag character vs
   any query), in **which queries** appear, and in `refr_base` (null unless the
   query pins it). (example 05)
4. `--partition-by` **replaces** `qname`, and does not add the key as a column;
   `-f qual` means `qual_phred`. (examples 03, 04)
5. Rows within a record are not sorted by `off_5p`, `qual` is null on real
   bases for QUAL `*` (since 0.1.4; it was 255), and `record_id` is dense only
   without `--only-hits`. (examples 03,
   06)

---

## Discrepancies

Code vs comments vs docs (`README.md`, `docs/alnbase.md`,
`docs/cli-reference.md`) vs `--help`, with evidence.

| # | claim | where | reality | evidence | severity |
|---|---|---|---|---|---|
| D1 | ~~"Any error exits 1; there are no distinct codes"~~ | ~~`docs/cli-reference.md` Exit status~~ | **Fixed in 0.1.2** (documentation): `docs/cli-reference.md` now says command-line syntax errors exit 2 | example 02 section 8 | ~~low~~ |
| D2 | "A zero count prints nothing, so a quiet run means none of these happened" (records skipped listed among them) | `docs/cli-reference.md` Exit status | the summary line always prints `skipped K`, including `skipped 0` | examples 03-06 | low |
| D3 | ~~summary suffix "(written untagged)" vs "(given a null row)"~~ | ~~`src/main.rs:562-566`~~ | **Fixed in 0.1.2.** Parquet/IPC runs print "(given a null row)", the tagged BAM "(written untagged)", and `--only-hits` runs no suffix (`src/main.rs:692-697`) | examples 03-06 | ~~low~~ |
| D4 | "non-UTF-8 aux text: values read as null rather than lossily converted" | `docs/cli-reference.md` diagnostics table; `src/batch.rs:139-141` | `query --parquet` **stops** on it unless `--permissive`; `extract` nulls it silently and never counts it | examples 03 s1, 05 s9; `src/scanner.rs:352-365`, `src/tag_extract.rs:379-382` | medium |
| D5 | `extract --permissive`: "as for `query`"; "Warn and count instead of stopping on data the run cannot handle" | `docs/cli-reference.md` extract; `extract --help` | has no effect in `extract`: no reference means no off-reference check, and `TagExtractShard` ignores aux trouble | `src/main.rs:741-805`, `src/tag_extract.rs:298-367`; example 05 s7, s9 | medium |
| D6 | ~~`qual_phred`: "empty when QUAL is '*'"; walk: "A record with '*' for QUAL has no array ... fall back to -1"~~ | ~~`src/record_field.rs:473`, `--help`; `src/alignment.rs:125-129`~~ | **Fixed in 0.1.4.** `qual` and `qual_phred` are null for QUAL `*` (`--help` now says "null when QUAL is '*'"); up to 0.1.3 htslib's 0xFF gave `qual = 255` and `qual_phred` = `ÿ` × length | example 06 s1 | ~~medium~~ |
| D7 | ~~Only the BAM output's partial-file behaviour is documented~~ | ~~`docs/cli-reference.md`~~ | **Fixed in 0.1.2.** A failed parquet/IPC run removes its output files (`src/parallel.rs:513-533`), and `docs/cli-reference.md` says so | example 06 s2 | ~~medium~~ |
| D8 | ~~the sharded pipeline has no equivalent of the ordered pipeline's "first error that is not [Stopped]"~~ | ~~`src/ordered.rs:42-45` vs `src/parallel.rs`~~ | **Fixed in 0.1.2.** Every worker is joined and a worker's own error wins over the reader's closed-channel error (`src/parallel.rs:474-501`) | example 06 s2 | ~~medium~~ |
| D9 | ~~`FORMAT_VERSION` 5 "added `read_off_3p` (and `capture_read_off_3p`)"~~ | ~~`src/hits.rs:78-80`~~ | **Moot in 0.1.8.** The version history comment was removed: no history of format versions is kept before the first release (`src/hits.rs:70-74`); `FORMAT_VERSION` is 6 | example 03 metadata | ~~low~~ |
| D10 | "Each file stays sorted by position" | `src/parallel.rs:64-67` | true of record `pos` for a sorted input; hit `refr_pos` is not monotonic within a file (bottom strands descend; rows emit at a query's last column) | example 03 s2, s9 | low |
| D11 | `--permissive` help: "Without it, a record on a contig the reference does not have stops the run ... Nothing about what is written changes either way" | `query --help` (`src/cli_query.rs:332-337`) | it also governs damaged aux blocks (BAM) and non-UTF-8/damaged aux columns (parquet), which change what is written | `src/bam_out.rs:362-373`, `src/scanner.rs:352-365`; example 03 s1 | low |
| D12 | `strand` column: "reference strand the pattern matched" | `src/record_field.rs:491`, `--help` | also filled (`+`/`-` from the FLAG) on hitless and unmapped rows, where nothing matched | example 03 s2 | low |
| D13 | `docs/alnbase.md` lists `--output-format ipc` among `query` options without noting it | `docs/alnbase.md` query options table | refused without `--parquet` (`only used with --parquet`) | `src/cli_query.rs:614-625`, `src/main.rs:527`; example 02 s7 (same mechanism) | low |
| D14 | "Every record read gets a row, so ... the table lines up one-to-one with the BAM" (for `extract` too) | `docs/alnbase.md` extract | holds only without `--only-hits`; and a record's "matched nothing" means *no tag character* in `extract` but *no query fired* in `query --parquet`, so hitless rows differ between the two | example 05 s2 | low |
| D15 | ~~`--row-group-rows` "Number of rows per row group" with no range~~ | ~~`--help`, `docs/cli-reference.md`~~ | **Fixed in 0.1.2.** `0` is refused by the argument parser (`must be at least 1`, exit 2), for `--batch-rows` too (`src/cli_query.rs:699-705`) | example 06 s4 | ~~low~~ |
| D16 | Parquet `OUT` is "a path prefix" | `docs/cli-reference.md`, `--help` | `-` is not special: it writes a file named `-_0_0`; not stdout | example 06 s4 | low |
| D17 | ~~`format_version` "Written into the parquet file metadata"~~ | ~~`src/hits.rs:70-71`~~ | **Fixed in 0.1.8.** `format_version` (with `coordinate_base` and the manifest) is now a plain footer key-value entry as well as Arrow schema metadata (`src/hit_writer.rs:179-195`), so the comment (`src/hits.rs:70-72`) is true | example 03 s8 | ~~low~~ |
| D18 | ~~`off_5p` "counts from the 5' end of the read as sequenced"~~ | ~~`--help` for `--field`, `src/hits.rs:16-34`~~ | **Fixed in 0.1.2.** Offsets are as sequenced for every FLAG, read 2 included (`src/hits.rs:456-481`), so the help text is now true | 03 §3.2, 03 example 01 | ~~see 03~~ |
| D19 | "each fills to `batch_rows` rows before it drains — so the rows held in memory scale with the file count ... Many files wants a smaller `--batch-rows`"; "memory scales with the file count: budget roughly `--batch-rows` rows per file" | `src/parallel.rs:52-55` (module doc); `--shards-per-worker` help (`src/cli_query.rs:485-487`) | `--batch-rows` is now a total divided over the files (`src/cli_query.rs:523-527`, `src/main.rs:598`), so buffered rows do not grow with the file count | code only | low |
