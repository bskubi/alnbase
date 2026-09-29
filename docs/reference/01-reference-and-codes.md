# 01 — Sequence model, pattern codes, the reference index, and BAM↔reference matching

Reference for alnbase 0.1.2. It covers:

- the 8-bit value that every column holds,
- the codes a pattern can use and **exactly** how an observed value is tested against them,
- how `alnbase index` reads a FASTA, and the `.aref` file it writes,
- `alnbase info`,
- coordinate conventions,
- how a BAM's contigs are matched to the index, and what that check misses.

**How to read this document.**

- Statements are normative. Each one cites the source that implements it as `src/file.rs:line`.
- **(example NN)** means the claim was checked by running the release binary on the example in `docs/reference/examples/01-reference/NN-*/`. That directory's `expected.txt` holds the captured output.
- **(code only)** means the claim was read from the source and not run.
- "Observed value" means the symbol a column actually holds when a record is walked. "Code" or "pattern set" means what a query row says.

**Running the examples.** Each example directory has a `run.sh`. It copies its inputs to a temporary directory, converts SAM to BAM with samtools, runs `alnbase index` and the command being tested, and prints the result. It needs `samtools`, and `python` with `duckdb`, on `PATH`. The alnbase binary is taken from `$ALNBASE`, which defaults to `alnbase` on `PATH`.

| Example | Shows |
|---|---|
| `01-subset-matching` | Which codes match reference `N`, `R` and `S`, lowercase reference, a read `N`, `=` and `/`, and a reverse-strand read |
| `02-fasta-symbols` | Which FASTA bytes are accepted (`U u N n R Y` and the rest of IUPAC), how they are stored, printed and matched, and the error for `X * 7 0 . - E %` |
| `03-fasta-parsing` | Headers, comments, blank lines, duplicate names refused, gzip/BGZF/CRLF, CR-only files, error cases, no partial output |
| `04-aref-format` | The version 3 file layout decoded byte by byte, including the per-contig MD5 (compared with `samtools dict`); version, magic and truncation errors; an empty FASTA |
| `05-info-and-coordinates` | `info`, `info --seq`, 0-based hit coordinates, PAD past a contig's ends, an alignment that runs off the contig |
| `06-contig-matching` | Contig order, BAM-only and index-only contigs, length mismatch, no shared names, unmapped records, a wrong sequence with the same name and length (caught by `@SQ M5`, missed without it), `--permissive`, `overlap`, `--trace-records` |
| `07-pattern-codes` | Parsing of pattern rows: spellings, groups, aliases, lowercase, digits, relational codes |
| `08-idioms` | Selecting ambiguous and exactly-`N` reference symbols; capturing reference context |

---

## 1. The value: `Seq`, a set in one byte

Every column has two observed values, one for the read and one for the reference. Each is a `Seq(u8)` (`src/seq.rs:65`). A `Seq` is a **set**:

| bit | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| member | A | C | G | T | GAP | PAD | CLIP | SKIP |

The bits are defined at `src/seq.rs:74-82`. All 256 values are legal (`src/seq.rs:112-129`). The two masks are `BASES = 0b0000_1111` and `FLAGS = 0b1111_0000` (`src/seq.rs:108-110`).

### 1.1 What each member means

| Member | Meaning | Produced by |
|---|---|---|
| A C G T | A nucleotide | Read SEQ; the reference FASTA |
| GAP | No base here on this side: a deletion (read side), or an inserted read base under `--insertions emit` (reference side). | The walk (`src/alignment.rs`, `D` and `I`) |
| PAD | Outside what exists on this side. **Read side:** a flank column past either end of the read, beyond any soft clip. **Reference side:** a position `< 0` or `>= contig length`. A query fires on read-side PAD columns like any other, windows made only of them included; how many exist is `--end-context`, by default the widest query's span (span − 1 before 0.1.8), so see the edge case in [`02-query-language.md` §8.6](02-query-language.md#86-pads-flank-contig-ends-and-record-boundaries). | `emit_flank` puts PAD on the read side (`src/alignment.rs`). `Aref::fetch` puts PAD on the reference side off the contig (`src/aref.rs:321-345`). |
| SKIP | A reference base skipped by `CIGAR N` (read side), or the single collapsed intron marker, which is SKIP on both sides. A window of only SKIP columns fires if the query matches it (example 07: `J,@~~` over `,,@AA`). | `src/alignment.rs:207-256` |
| CLIP | **Read side only.** A flank column standing for a soft-clipped base: the read continues here, but the aligner did not align it. The column has the real reference base, the clipped base's offsets and no quality. One per clipped base at that end, nearest the aligned part, up to the flank width (`--end-context`); the rest of the flank is PAD. | `emit_flank` (`src/alignment.rs`); 03 §5.2 |

Until 0.1.10 bit 6 was ERR, which the indexer stored for any FASTA byte it did not recognise. The indexer now refuses such a FASTA (§5.1), and the bit is CLIP.

### 1.2 Named constants

| Constant | Bits | Canonical name |
|---|---|---|
| `EMPTY` | `0` | `0` |
| IUPAC `R Y S W K M B D H V N` | unions of base bits (`src/seq.rs:85-95`) | the letter |
| `GAP_PAD` | GAP\|PAD | `{._}` |
| `N_GAP_PAD` | ACGT\|GAP\|PAD = `0b0011_1111` (`src/seq.rs:103`) | `{N._}` |
| `ANY` | all eight bits = `0b1111_1111` (`src/seq.rs`) | `~` |

**`~` is every symbol**, clips included. It still cannot match the empty value `0` (a read `=` in SAM SEQ), because nothing can (§4). Verified in examples 01 and 07.

### 1.3 How a value is written: canonical names and one-byte forms

There are two renderings.

**Canonical name** — `Seq::name`, defined at `src/seq.rs:137-189`. This is the form used in parquet `read_base`, `refr_base`, `capture_read` and `capture_refr` (`src/hits.rs:221-228`), and in `--explain` and warnings.

- `0` is the empty set. `~` is ANY.
- The base nibble is written as its IUPAC letter (table at `src/seq.rs:68-70`).
- Flags are written as `.` (GAP), `_` (PAD), `:` (CLIP), `,` (SKIP).
- A value that mixes two or more categories is wrapped in braces, for example `{C.}` or `{N._}`. Inside braces, SKIP is written `J` rather than `,` (`src/seq.rs:181-185`).
- Every name parses back to the same value (`src/seq.rs:391-397`, test).

**One-byte form** — `Seq::to_ascii`, defined at `src/seq.rs:234-268`. This is the form used by `info --seq` and trace grids.

- GAP→`Z`, PAD→`X`, CLIP→`L`, SKIP→`J`, empty→`0`.
- A value that mixes bases with flags becomes `*`, so this form is lossy.

An index holds only base sets, so `info --seq` prints only IUPAC letters (example 02).

---

## 2. Where observed values come from

### 2.1 Read side

The walk copies htslib's 4-bit SEQ code straight into `Seq`, with no decode step (`src/alignment.rs:120-123`). htslib's order `=ACMGRSVTWYHKDBN` lines up bit for bit with the base nibble. Two consequences:

- **An observed read value is always a base set.** It is either one of the 15 IUPAC sets or the empty set. It never carries a flag bit from SEQ.
GAP, PAD, CLIP and SKIP on the read side come only from the walk: deletions, flanks (clips and pads) and introns.

When samtools or htslib converts **SAM text** to BAM, SEQ characters are mapped as follows (example 01, read `odd`: SEQ `A=GT.XU*acgt` is stored as `A=GTNNTNACGT`):

| SAM SEQ character | Stored as | Observed `Seq` |
|---|---|---|
| `A C G T` and IUPAC letters, either case | that code (lowercase is folded) | that set |
| `U` / `u` | T | T |
| `=` | code 0 | **empty set `0`**, which no pattern matches, not even `~` |
| `.`, `X`, `*`, any other character | 15 (N) | N |

In example 01, read `odd` has no row at `refr_pos` 1 (the `=`), even for `~@~`. Its `.`, `X` and `*` are observed as `N`, `U` as `T`, and `acgt` as `ACGT`.

### 2.2 Reference side

`alnbase index` converts each sequence byte with `Seq::FROM_FASTA` (`src/aref.rs`, table in `src/seq.rs`):

| FASTA byte | Stored `Seq` |
|---|---|
| `A C G T R Y S W K M B D H V N`, upper or lower case | that set. Case is folded; soft-masking is lost. |
| `U`, `u` | T |
| ASCII whitespace inside a sequence line | dropped |
| **anything else**, including `X`, `*`, `-`, `.`, digits, `E`, `J`, `Z`, `_`, `%`, `~` | **refused**: `index` stops with `line 4: 'X' at c2:3 (1-based) is not an IUPAC base; …` and writes no index (§5.1) |

All rows verified in example 02, and lowercase also in example 01. Up to 0.1.9 the indexer stored `-` and `.` as GAP, `0` as the empty set and every other byte as ERR, without a warning, so a pattern silently failed to fire there.

Past either end of a contig, the reference value is PAD, whatever coordinate is asked for (`src/aref.rs:321-345`). This applies to flank columns and also to **aligned** read bases, when a record's alignment runs past the contig end. Neither samtools nor alnbase rejects such a record (example 05, read `overhang`).

The intron marker column is SKIP on the reference side too (`src/alignment.rs:247-253`).

### 2.3 Orientation: the complement applies to both sides

A read is walked along its conversion strand. For a read the directional library puts on OB, the walk runs backwards along the reference and complements **both** sides (`src/alignment.rs:101-102`, `167-168`, `404-413`). The rule `is_ot = is_last_in_template == is_reverse` means that an unpaired or R1 reverse read is OB, and so is an **R2 forward** read (code only for R2).

`Seq::complement` swaps A↔T and C↔G and leaves the flag bits alone (`src/seq.rs:292-300`). On IUPAC codes this means `R`↔`Y`, `K`↔`M`, `B`↔`V` and `D`↔`H`, while `S`, `W` and `N` are unchanged.

So a reference `R` under a reverse read is **observed as `Y`**. A pattern `~@R` does not fire there (example 01, read `rev`, `refr_pos` 5). `refr_pos` still reports the forward-strand coordinate.

---

## 3. Pattern codes: what a query row may contain

In a query file, a pattern is two rows, `read` and `refr`. The rows are joined as `read@refr` and parsed by `dsl::parse_pattern_text` (`src/lower.rs:201`). Each side is scanned by `parse_side_with` (`src/dsl.rs:299-445`).

### 3.1 Single-character codes

A single character is parsed with `Seq::from_name` (`src/dsl.rs:422-427`, `src/seq.rs:192-224`).

| Code | Set | Notes |
|---|---|---|
| `A C G T` | one base | |
| `R Y S W K M B D H V N` | IUPAC union | `src/seq.rs:85-95` |
| `U` | T | example 07 |
| `.` `-` `Z` | GAP | example 07 |
| `_` `X` | PAD | example 07 |
| `,` `J` | SKIP | example 07 |
| `:` `L` | CLIP | example 07 |
| `~` | ANY = `{N._:J}` | every symbol |
| `=` `/` | relational (§3.4) | |
| lowercase `a`–`z` | a declared alias (§3.3). Otherwise an error: "'c' is lowercase; base codes are uppercase…" (`src/dsl.rs:411-421`) | example 07 |
| `E`, `O`, `P`, `%` and other non-codes | error "unknown code" | example 07 |
| `0`–`9` | **rejected in query-file rows**: "digit in the read row" (`src/lower.rs:175-186`). So `0` (the empty set) and repeat counts like `C7` cannot be written. | example 07 |
| whitespace | error (`src/dsl.rs:310-312`) | |
| `^` `+` | error: markers go in the `mark` row (`src/dsl.rs:480-491`) | |

The set of parseable letters does not depend on case: `from_name` uppercases its input (`src/seq.rs:198`). Outside braces, the side parser rejects lowercase before `from_name` is ever called.

### 3.2 Groups `{...}`

`{body}` is the union of the codes in `body`, parsed by `Seq::from_name` (`src/dsl.rs:314-346`).

- Codes can appear in any order. Braces inside the body are ignored by `from_name`.
- Nested `{` is an error, as are an empty group, a group that names only the empty set (`{0}`), and a comma inside a group (`src/dsl.rs:321-343`). Example 07.
- **Lowercase is accepted inside a group.** `{c.}` is `{C.}` (example 07).
- **An alias letter inside a group is not the alias.** It is read as whatever the letter means as a built-in code. With `j = "A"` declared, `{j}` means `J`, which is SKIP (example 07). By the same code path, `{x}`=PAD, `{z}`=GAP, `{u}`=T and `{a}`=A (code only).
- **In a query-file row, the row-width check counts characters, not columns** (`src/lower.rs:188-194`). `read = "{C.}~"` with `refr = "CG"` is rejected as "5 columns vs 2". A group only gets through when the other row happens to have the same number of *characters*: `read = "{C.}"`, `refr = "{NJ}"` works, and so does `refr = "{~~}"` (example 07). **In practice, use aliases for multi-member sets.**

### 3.3 Aliases

`[alias]` maps one lowercase letter to a set (`src/dsl.rs:41-130`, `src/lower.rs:145-159`).

- The name must be a single lowercase ASCII letter. An uppercase name is an error ("uppercase letters are base codes"). A non-letter is an error (`src/dsl.rs:82-99`). Example 07.
- The value is parsed with `Seq::from_name`, so braces are optional. `"c."`, `"{C.}"` and `"C."` are the same set, and lowercase is accepted (example 07).
- The value must be non-empty. It can be any set, for example `{_:}`, a pad or a clip (03 example 07).
- `=` and `/` cannot be alias values: "'=' is not a base set" (example 07).
- Declaring the same letter twice in one table is an error. Across several query files, the same letter may be declared again with the **same** set, but a different set is an error (`src/dsl.rs:104-119`).
- Aliases can be declared after the patterns that use them (`src/lower.rs:119-124`).

### 3.4 Relational codes `=` and `/`

- A relational code goes on **one** side. The column's set on that side becomes `N_GAP_PAD`, and the relation does the constraining (`src/dsl.rs:348-354`). Putting a relational code on both sides is an error (`src/dsl.rs:543-548`, example 07).
- `Rel::holds(read, refr)` requires that **each** observed value has no flag bit and exactly one base bit. Then `=` needs the two bases equal and `/` needs them different (`src/dsl.rs:210-224`).
- The other side's code still applies as a subset test. `Y@=` fires on C over C, but not on T over C (example 07).
- Neither code fires on GAP, PAD, CLIP, SKIP, empty, or **any IUPAC ambiguity on either side**. A read `A` over reference `N`, and a read `N` over reference `A`, satisfy neither `~@=` nor `~@/` (`src/query.rs:433-442`, test; example 01 at `refr_pos` 4, and read `allN`). They are **not complements**.
- If a relational column can never be satisfied — for example `.@=` — the query fails to compile with `DeadColumn` (`src/query.rs:290-307`, example 07).
- **Compile warnings on the idiomatic `~@=`.** When a declared member can never survive the relation, a warning is printed ("'{._:J}' on the read side is unreachable; the column reduces to 'N'", `src/query.rs:309-326`). So the idiomatic `~@=` and `~@/` print two warnings each (example 01). These warnings are harmless.

### 3.5 Where `^`, `+` and counts go

Digits are rejected in query-file rows (§3.1), and `^` and `+` belong in the `mark` row (`src/dsl.rs:480-491`). Since 0.1.2 `--list-codes` says so: its "MARK ROW" section lists `.`, `^` and `+` for the `mark` key and states that rows take no repeat counts (`src/codes.rs:78-86`). See D5 (fixed).

---

## 4. The matching rule (normative)

A column of a pattern matches an observed column `(read_obs, refr_obs)` when **all** of the following hold:

1. `read_obs` is a **non-empty subset** of the column's read set.
2. `refr_obs` is a **non-empty subset** of the column's reference set.
3. If the column is relational, the relation holds (§3.4).

This is how the compatibility table is built. For each column, alnbase enumerates every non-empty subset of the read set against every non-empty subset of the reference set, filters by the relation, and marks those pairs as compatible (`src/query.rs:275-287`, `Seq::bit_combinations` drops the empty set at `src/seq.rs:303-316`). At scan time, the observed pair is looked up by index (`src/query.rs:160-162`, `src/scanner.rs:268`). The documented contract is at `src/dsl.rs:229-230`.

**It is a subset test, not an intersection test.** An observed ambiguity matches only a code that contains *all* of its members.

### 4.1 Consequences

All rows below verified in example 01, except where noted.

| Observed | Codes that match it | Codes that do **not** |
|---|---|---|
| reference `A` | `A R W M D H V N ~` and any group or alias containing A | `C G T Y S K B` |
| reference `N` (FASTA `N`/`n`) | `N`, `~`, and supersets such as `{N.}` | **every** base code and **every** partial IUPAC code: `G`, `H`, `R`, `S`… and `=`, `/` |
| reference `R` (A/G), forward read | `R V D N ~` | `A`, `G`, `S`, `H`, `=`, `/` |
| reference `R` under a reverse (OB) read | observed as `Y`, so `Y B H N ~` | `R` |
| reference `S`, either strand | `S V B N ~` | `C`, `G`, `R` |
| reference lowercase `c` | same as `C` | — |
| read `N` | `N`, `~` | `C`, `Y`, `=`, `/` |
| read `=` in SAM SEQ (empty set) | **nothing**, not even `~` (read `odd`) | everything |
| read CLIP (flank beside a soft clip) | `:` (`L`), `~` | bases, `_` (03 examples 02 and 07) |
| reference PAD (off contig) | `_` (`X`), `~` | bases (examples 05 and 08) |

**Practical meaning.**

- A pattern such as `C~@CG` never fires where the reference has `N`, `R`, `S`, etc. in either column.
- A read `N` never counts as a C or a T.
- A `~` column in a query stops that query from firing across a read `=`, and only there.

---

## 5. `alnbase index FASTA INDEX`

The implementation is `aref::write_from_fasta` (`src/aref.rs:85-188`), called from `run_index` (`src/main.rs:296-297`). It takes no options.

### 5.1 Input

- **Compression.** Plain, gzip (not BGZF) and BGZF inputs are all accepted, detected by htslib's `bgzf` reader (`src/aref.rs:100-103`). Indexes built from `ref.fa`, `gzip -c ref.fa` and `bgzip -c ref.fa` are byte-identical (example 03).
- **Line endings.** Each line has trailing `\n` and `\r` stripped (`src/aref.rs:69-79`), so CRLF gives an identical index (example 03). A **CR-only** file (classic Mac) is one single "line":
  - If the file starts with `>`, it becomes **one contig of length 0**, and the rest of the file is treated as a header description. No warning is given (example 03).
  - If the file starts with `;`, it becomes **zero contigs** (example 03).
- **Blank lines** are skipped (`src/aref.rs:120`, example 03).
- **Comment lines.** Any line whose first byte is `;` is skipped, wherever it appears — before a header, or between sequence lines (`src/aref.rs:120`, example 03).
- **Header lines.** A line starting with `>` begins a contig. The **name** is the bytes after `>` up to the first ASCII whitespace (space, tab, etc.). The rest is discarded (`src/aref.rs:130-133`).
  - `>chr1 desc` gives `chr1`, and `>chr2\tdesc` gives `chr2`.
  - Punctuation is kept: `>chr1|x:y` gives `chr1|x:y` (example 03).
  - `>` alone, or `> c1` (whitespace straight after `>`), fails with "FASTA record with an empty name" (example 03).
- **Sequence before the first header** fails with "sequence data before the first '>' header" (`src/aref.rs:156-158`, example 03).
- **Within a sequence line**, whitespace is removed. Every other byte must be an IUPAC letter or `U`, in either case (§2.2). **Any other byte stops indexing** with the line number, the contig and the 1-based position within it, and no index is written (example 02). A protein FASTA, an alignment FASTA with `-` gaps or a damaged download is therefore caught here rather than silently indexed.
- **Line width** is irrelevant. Lines of any length can be mixed (example 03).
- **Duplicate contig names** are refused: "the FASTA has more than one contig named dup", and no index is written (`src/aref.rs:137-144`, example 03). Before 0.1.2 both records were stored and lookup by name found the last one.
- **Names that are not UTF-8** are written without complaint, but every later open of the index fails with "contig name is not valid UTF-8" (`src/aref.rs:268-269`, example 03).
- A missing input fails with "file not found: missing.fa" (example 03).
- An empty FASTA gives a valid index with 0 contigs (example 04).

### 5.2 Output behaviour

- The index is written to `INDEX.partial` beside the destination and renamed to `INDEX` only when it is complete (`src/aref.rs:85-97`). If any error happens, the `.partial` file is removed and nothing is left behind. An existing `INDEX` is untouched by a failed run (example 03). A file already named `INDEX.partial` is overwritten.
- The whole conversion is a single streaming pass. Memory is one line plus the contig table (`src/aref.rs:81-84`).
- **Each contig's MD5 is computed in the same pass** (`src/aref.rs:60-67`, `159`, `128`, `172`). It is the SAM `@SQ M5` digest: MD5 over the contig's sequence lines with every byte outside `!`..`~` (ASCII 33–126) removed and lowercase letters uppercased. It is computed over the **FASTA bytes**, not the stored `Seq` values, so it matches what `samtools dict` writes for the same FASTA (example 04), and it distinguishes bytes that alnbase stores identically (`U` and `T`). Comment and blank lines are not part of it. A contig of length 0 gets the MD5 of the empty string, `d41d8cd98f00b204e9800998ecf8427e` (example 03).
- Contig order is FASTA order. It is irrelevant to BAM matching (§9).

---

## 6. The `.aref` file format

Layout is documented at `src/aref.rs:1-24`. All integers are little-endian `u64`, and offsets are absolute from the start of the file:

```
[ data body: name0 seq0 name1 seq1 ... ]   names: raw bytes, no NUL; seqs: 1 byte per base (Seq bits; bases only, no flag bit)
[ contig table: 48 B per contig        ]   seq_size, seq_offset, name_size, name_offset (4 x u64), then md5 (16 raw bytes)
[ trailer: 24 B                        ]   n_contigs, version, magic
```

- **Magic** is `b"REFRMM\x1a\n"`, read as a little-endian u64 (`src/aref.rs:39`).
- **Version** is `3` (`src/aref.rs`). Each contig table entry is 48 bytes (`src/aref.rs:44`): the four u64s at entry offsets 0, 8, 16 and 24, then the contig's 16 MD5 bytes at offsets 32–47, raw (not hex) (`src/aref.rs:176-181`, `257-262`). The table therefore starts at `file_len - 24 - 48 * n_contigs`.
- Example 04 decodes a two-contig file by hand: the names, the sequence bytes (`N` = `0b00001111`, `R` = `0b00000101`), each entry's MD5, and the trailer. The MD5s equal `samtools dict`'s `M5` values.

### 6.1 Opening an index

`Aref::open`, at `src/aref.rs:220-274`, maps the file into memory and checks it in this order. The checks run on every `info`, `query` (tagging, `--parquet`, `--trace-records`) and `overlap --refr`.

1. The file must be at least 24 bytes: "file too short to hold a trailer" (example 04).
2. The last 8 bytes must equal the magic: "bad magic". A truncated file normally fails here (example 04).
3. **Version.** A version greater than 3 fails with "index format version 4, but this alnbase understands 3; rebuild the index or use a newer alnbase". Any other version other than 3 fails with "index format version 2, but this alnbase understands 3; rebuild it with `alnbase index`" (example 04). **An index built by alnbase 0.1.9 or earlier must be rebuilt** with `alnbase index`. There is no migration and no compatibility window.
4. The contig table must fit in the file: "contig table overruns the file" (`src/aref.rs:248-251`).
5. Each entry's name and sequence must lie inside the data body: "contig entry points outside the data body" (`src/aref.rs:264-267`).
6. Each name must be UTF-8 (`src/aref.rs:268-269`).

Two further points:

- An I/O error is reported without the path, e.g. "Error: No such file or directory (os error 2)" (example 04). In `query`, that message does not say whether the BAM or the index was missing.
- The index is memory-mapped. Modifying or truncating it while a run has it open is undefined behaviour (`src/aref.rs:215-219`, code only).

**What the index stores and does not store:** it stores one MD5 per contig (the `@SQ M5` digest of the FASTA sequence). It does not store a checksum of the index file itself (the MD5s are never re-verified against the stored `Seq` bytes), the source path, soft-masking, the FASTA description text, or the original line width.

---

## 7. `alnbase info [--all] [--seq CHR:FROM-TO] INDEX`

Implemented at `src/main.rs:130-217`.

### 7.1 Summary output

The summary prints the path, `format version`, `contigs`, `total bases` (the sum of stored lengths), and `file size` (`src/main.rs:177-190`). `format version` is `3`.

It then prints `contigs`, one `name  length  md5` line each, in index (FASTA) order: the name left-aligned, the length right-aligned in 12 columns, and the contig's MD5 as 32 lowercase hex digits, the value `samtools dict` writes as `@SQ M5` (`src/main.rs:192-216`, examples 02–05):

- **11 or fewer** contigs are all listed.
- **More than 11** are shown as the first 5, then `... N more, --all to list them`, then the last 5. `--all` lists every contig (examples 03 and 05).
- An index with 0 contigs prints the summary only (example 04).

Duplicate names cannot occur: `index` refuses them (§5.1).

### 7.2 `--seq CHR:FROM-TO`

Implemented by `print_region` (`src/main.rs:130-159`).

- **Coordinates are 1-based and inclusive**, as in `samtools faidx`. The output header `>c1:1-5` and its sequence were identical to `samtools faidx` in example 05.
- Commas in the numbers are removed: `c1:1,0-1,2` gives `c1:10-12` (example 05).
- `TO` past the end is clipped, and the header shows the clipped value: `c1:18-40` on a 20 bp contig prints `>c1:18-20` (example 05).
- Errors (example 05):
  - `FROM < 1` or `TO < FROM`: "FROM must be at least 1 and TO at least FROM".
  - `FROM > length`: "c1 is 20 bases, so 21 is past its end".
  - A missing `-TO`: "expected CHR:FROM-TO".
  - An unknown contig: "the index has no contig chrX".
- The sequence is printed with `to_ascii`, 60 characters per line (`src/main.rs:153-157`). Letters are uppercase, and a FASTA `U` prints as `T`, the only difference from `samtools faidx` beyond case (example 02).
- Only the name up to the first `:` is taken as the contig (`src/main.rs:131-133`). A contig whose name contains `:` (e.g. HLA alleles) cannot be addressed (code only).

---

## 8. Coordinates

| Where | Convention |
|---|---|
| `refr_pos` and `capture_refr_pos` in hit rows | **0-based**. A SAM `POS 1` read's first aligned base has `refr_pos 0` (example 05). |
| Off-contig reference positions | Emitted as-is: `-1`, `-2`… before the start, `LN`, `LN+1`… after the end (examples 05 and 08). Their reference value is PAD. |
| Flank columns (read side PAD) | `off_5p`/`off_3p` are null, as the anchor and in capture lists (example 05; since 0.1.9 for the anchor). |
| `Aref::fetch` / `slice` (internal) | 0-based, half-open `[start, end)`. `fetch` pads out of range with PAD; `slice` returns `None` (`src/aref.rs:311-345`). |
| `info --seq` | **1-based inclusive** (§7.2) |
| `--trace-records` label `contig:pos` | **1-based** (`rec.pos() + 1`, example 06); since 0.1.4 the label says so: `chr1:1 (1-based POS)`. The per-column `refr_pos` in the same trace is 0-based and printed as `refr_pos 0 (0-based)`. |
| `pos` / `end_pos` record fields (`-f`) | 0-based leftmost / 0-based exclusive end (from `--help` for `--field`) |

On a reverse (OB) read, `refr_pos` decreases as the column index increases, but each value is still the forward-strand 0-based coordinate (`src/alignment.rs:76-82`).

---

## 9. BAM ↔ reference matching

### 9.1 The algorithm

`ContigMap::build(header, aref)` is implemented at `src/contig_map.rs:59-147`. It is called by `query` (BAM output, `src/ordered.rs:147`; `--parquet`, `src/parallel.rs:316`; `--trace-records`, `src/main.rs:435`) and by `overlap --refr` (`src/overlap_apply.rs:50`).

1. For every `@SQ` entry (BAM tid) in header order:
   - The name must be valid UTF-8, or the run fails (`src/contig_map.rs:69-76`).
   - The name is looked up **by exact, case-sensitive byte match** in the index (`src/contig_map.rs:79`, `src/aref.rs:293-295`). No aliasing is done: `chr1`≠`1`, `chrM`≠`MT`, `chr1_KI270706v1_random`≠`KI270706.1`.
2. **Found:** if the `@SQ LN` differs from the index length, the mismatch is recorded. Up to 10 are listed (`src/contig_map.rs:81-90`).
3. **Found, and the `@SQ` line carries `M5`:** the header's value is compared with the index's MD5 for that contig, ignoring ASCII case. A disagreement is recorded; up to 10 are listed (`src/contig_map.rs:91-100`). The header's `M5` values are read from the header text, keyed by `SN` (`header_m5`, `src/contig_map.rs:225-236`).
4. **Not found:** the tid maps to `None`, and the name and length are kept for a warning (up to 1000 stored, 10 printed; `src/contig_map.rs:104-109`, `204-221`). An `M5` on such a contig is not looked at.
5. **Any length mismatch is fatal**: "the BAM and the reference disagree about contig lengths, so they are different assemblies" (`src/contig_map.rs:113-120`, example 06). This is checked before the M5 comparison, so a contig that differs in both reports the length.
6. **Any M5 mismatch is fatal**: "the BAM header's @SQ M5 checksums disagree with the reference, so a contig with the same name and length holds different sequence:", followed by `NAME: M5 <header value> in the BAM header, <index value> in the reference` lines and "Re-index the FASTA the BAM was aligned against." (`src/contig_map.rs:122-129`, example 06 §2).
7. **No shared names at all is fatal**, provided the header has at least one `@SQ` — even if the BAM holds no mapped records: "no contig in the BAM header is present in the reference by name…" (`src/contig_map.rs:131-144`, example 06). A header with **zero** `@SQ` lines is accepted (`src/contig_map.rs:363-369`, test).
8. If some but not all names are found, a warning is printed to stderr: "N of M BAM contigs are not in the reference; records on them will be skipped:" (example 06).
9. **`query --require-m5`** (0.1.8): `ContigMap::require_m5` marks every contig the index has whose `@SQ` line carries no `M5` (`src/contig_map.rs:155-170`; called at `src/ordered.rs:148-150`, `src/parallel.rs:317-319`, `src/main.rs:436-438`). `walkability` classifies a record on such a contig as `MissingM5`, after the unmapped and off-reference tests (`src/scanner.rs:558-570`), and the first one stops the run with "record NAME is on contig CHR, whose @SQ line has no M5 checksum, and --require-m5 was given. Add M5 to the header's @SQ lines (`samtools dict REF.fa` prints them …; `samtools reheader` applies an edited header), or leave out --require-m5 to check names and lengths only." (`src/batch.rs:46-53`, `76-80`), **whatever `--permissive` says**. Like an off-reference contig it is decided per record: a header contig no record is on is not checked, and a contig the index lacks is off-reference, not missing an `M5`. (Checked by running 0.1.8 on example 06's inputs, not in `expected.txt`: `reordered.bam`, whose `chr1` and `chr10` lack `M5`, stops at `on_chr10`, the first record on either; `chr1` holds no record. With `--permissive` it stops the same way.)

After that:

- Per record, `refr_tid(bam_tid)` indexes a `Vec` (`src/contig_map.rs:177-182`). A tid of `-1` or one past the end of the header gives `None`.
- The BAM's tid order and the index's contig order are unrelated, and matching is correct across any reordering (`src/contig_map.rs:272-285`, test; example 06, scenario 1).
- Index contigs that are not in the BAM are ignored without comment (example 06).

#### What the `M5` check catches

The comparison is between the header's `@SQ M5` string and the MD5 that `index` computed from the FASTA (§5.2). Exactly:

- **Caught:** a contig present in both, by name, whose `@SQ` line carries `M5`, where the FASTA that was indexed differs from the one the `M5` was computed from in any byte that survives the `M5` normalisation (printable ASCII, uppercased). That includes a same-name, same-length contig with different bases (example 06 §2), and also differences alnbase itself would not see, because the digest is over FASTA bytes: `U` vs `T`, `X` vs `*`, `.` vs `-`. Such a run stops even though the stored sequences are identical (code only).
- **Accepted:** hex case differences in the header value (`src/contig_map.rs:95`, test at `src/contig_map.rs:326-340`); FASTA differences only in case (soft-masking), line width, whitespace, comment lines or header description text, none of which enter the digest.
- **Not checkable:** a contig whose `@SQ` line has no `M5`. **Nothing can tell a same-name, same-length, different-sequence contig apart there, and the run passes silently** (example 06 §2, second run), unless `--require-m5` turns such a contig into an error for any record on it (step 9). The reference concordance check (`--max-discordance`, see `04-outputs-and-guarantees.md` §6) judges a finished run only from 1,000 compared bases; example 06 §2 compares 8 (37.50% differ) and exits 0. Many aligners do not write `M5`; `samtools dict` output and headers built from it do. `UR`, `AS` and `SP` are never read.
- **Not looked at:** `M5` on a contig the index lacks, and any `M5` in `extract` (no index) or in `overlap` without `--refr`.
- A malformed `M5` value (not 32 hex digits) can never be equal, so it stops the run (code only).

### 9.2 Per-record handling

| Record | `query` (BAM out) | `query --parquet` | `query --trace-records` | `overlap --refr` |
|---|---|---|---|---|
| unmapped (FLAG 4) | skipped, counted, written through untagged, **never an error** (`src/ordered.rs:262-263`, `src/batch.rs:60-86`) | skipped, counted, gets a null row (`src/parallel.rs:432-440`) | skipped silently (`src/main.rs:451-457`) | not a reference question |
| mapped, contig not in index, default | **run stops** at the first such record (`src/batch.rs:37-43`, `69-73`). The partial BAM is removed. | **run stops** (`src/batch.rs:69-73`), and every output file the run was set to write is removed (`src/parallel.rs:516-518`, `522-533`). | **skipped silently, no error** (`src/main.rs:451-457`) | NM/MD removed, no error (`src/overlap_apply.rs:370-374`) |
| mapped, contig not in index, `--permissive` | skipped, counted ("of those, N were on contigs the reference does not have"), written untagged | skipped, counted, null row | same as above | n/a (no `--permissive`) |
| mapped, contig in index, `@SQ` without `M5`, `--require-m5` (0.1.8) | **run stops** (`src/batch.rs:76-80`), with or without `--permissive`; the partial BAM is removed | **run stops**, files removed as above | **trace stops** with the same error (`src/main.rs:453-456`) | n/a (no `--require-m5`) |

All table cells verified in example 06, except `overlap`'s NM/MD removal (code only) and the `--require-m5` row (checked by running 0.1.8 on example 06's inputs for all three `query` modes; the BAM mode printed `removed the incomplete output`).

- A header error (length mismatch, M5 mismatch, or no shared names) is raised before any output is created (example 06 §5 for parquet, `src/parallel.rs:314-326`; code only for BAM, `src/ordered.rs:147-156`). `--require-m5` is not a header error: it stops at the first record on an unchecked contig, after the output exists, so the output is removed as for an off-reference record.
- `extract` and `query --trace` never open an index, so they do no matching at all (`src/main.rs:354-364`, `741-805`).
- `overlap` builds the same map, including the M5 check, only when `--refr` is given (`src/overlap_apply.rs:40-57`, example 06).
- Since 0.1.8 `--trace-records` also checks each traced record against the reference as a scan does: after each record's trace it prints `reference concordance: X of Y compared read bases differ (P%; a reference C read as T is not counted)` (example 06 §7: `0 of 3` for each record), and at the end the traced records' total stops the trace with exit 1 when at least 1,000 compared bases exceed `--max-discordance`, or prints a warning when fewer do (`src/main.rs:477-488`, `498-510`). On example 06's `no_m5.bam` against `other_assembly.fa` it printed `warning: 37.5% of the 8 read bases traced disagree with the reference, above --max-discordance 0.25; too few to judge (a scan judges from 1000 bases), …` and exited 0 (run by hand, not in `expected.txt`).
- In parquet mode, the end-of-run message says "(given a null row)"; "(written untagged)" is used only for a tagged BAM (`src/main.rs:687-698`). Before 0.1.2 parquet runs said "(written untagged)" too (D8, fixed).

### 9.3 Fragility analysis

**Caught (the run fails):**

| Mismatch | Evidence |
|---|---|
| Same name, different `LN` | example 06 §5 |
| Same name and `LN`, different sequence, **when the `@SQ` line carries `M5`** | example 06 §2 |
| Duplicate contig names in the FASTA (refused by `index`) | example 03 |
| No `@SQ` name exists in the index (e.g. a UCSC-style index with an Ensembl-style BAM, where every name differs) | example 06 §6 |
| A mapped record on a contig the index lacks (default mode) | example 06 §4 |
| A mapped record on a contig whose `@SQ` line has no `M5`, **with `--require-m5`** (any mode, `--permissive` included) | §9.1 step 9 (run on example 06's inputs) |
| Index file damaged, truncated, wrong version (including every version 1 index from 0.1.1), or with non-UTF-8 names | examples 03 and 04 |
| `@SQ` name not UTF-8 | code only, `src/contig_map.rs:69-76` |

**Skipped and counted (only with `--permissive`), or warned only:**

| Situation | Behaviour |
|---|---|
| BAM-only contigs **with no records** | warning only, even without `--permissive` (example 06 §3) |
| Records on BAM-only contigs under `--permissive` | counted as `records_off_reference` |
| **`--trace-records`** on records from BAM-only contigs | silently omitted from the trace (example 06 §7) |
| **`overlap --refr --stale-tags recompute`** on records from BAM-only contigs | silently loses NM/MD (code only) |

**Pass silently (no error, no warning, wrong or partial output):**

1. **Same name and same length, different sequence, with no `M5` in the header.** Examples: another patch release of the same assembly, a masked variant with substitutions, a FASTA with IUPAC codes instead of `N`, a lifted or edited contig. Every column pairs the read with the wrong base (example 06 §2, second run: a read `G` is reported over a reference `T`). The only sequence identity check is `@SQ M5` (§9.1), so such a header cannot be checked. Adding `M5` to the header (for example from `samtools dict` with `samtools reheader`) makes it checkable, and `--require-m5` (0.1.8) refuses to run on a record whose contig lacks it rather than passing silently.
2. **Same sequence, different case/masking** — harmless, because case is folded (and the `M5` digest ignores case too).
3. **Partial naming overlap.** If *any* name matches, the check passes. Two cases:
   - A BAM with `chr1…chr22, chrX, chrY, chrM` against an index with `chr1…chr22, chrX, chrY, MT`. Only `chrM` is unresolved, which gives a startup warning and then a fatal error at the first `chrM` record. That is caught in default mode but **silently dropped** under `--permissive`, apart from a count.
   - A decoy or alt set present in the BAM but absent from the index behaves the same way.
4. **A record aligned past the contig end** is walked normally, with PAD on the reference side of aligned bases (example 05).

Fixed and removed from this list: reference bytes that are not IUPAC, which were stored as ERR so queries silently failed to fire there (refused by `index` since 0.1.10, §5.1); duplicate FASTA names (now refused, §5.1); re-indexing over an existing index destroying it when the conversion fails (now written to `.partial` and renamed, §5.2); 0-row parquet shards left by a failed run (now removed, §9.2); a bogus `@SQ M5` being accepted (now checked, §9.1).

### 9.4 Mitigation implemented in 0.1.2 and 0.1.8

The mitigation this section used to propose is implemented:

1. **At `index` time** each contig's `M5` digest is computed over the FASTA bytes and stored in a 48-byte table entry, and a byte that is not an IUPAC base is refused (§5.1, §5.2, §6).
2. **In `ContigMap::build`** a resolved contig whose `@SQ` carries `M5` must agree with the stored value; a disagreement is fatal. Where `M5` is absent nothing is checked (§9.1).
3. **`info` prints each listed contig's MD5**, always (there is no `--md5` flag), so an index can be compared with a `.dict` file without a BAM (§7.1, §10).
4. **`query --require-m5`** (0.1.8) makes a missing `M5` fatal for any record on that contig, so a header that cannot be checked no longer passes silently when the flag is given (§9.1 step 9). It is off by default.

Still not implemented:

- Warn when a contig has length 0.
- Make `ContigMap::build` warn when *any* contig is unresolved **and** the unresolved names look like a known renaming (`chrM`/`MT`, a `chr` prefix difference).
- Name the offending path in `Aref::open` I/O errors.

---

## 10. Idioms

These are capabilities a reader might assume are missing.

- **Select positions where the reference is ambiguous.** Because matching is by subset, `~@N` alone also fires over plain bases. Exclude the single bases:

  `where = "n and not a and not c and not g and not t"`, with `n = ~@N`, `a = ~@A`, and so on.

  This fires on `N`, `R`, `S`, … on both strands (example 08).
- **Select positions where the reference is exactly `N`.** Exclude the three-base codes:

  `where = "n and not b and not d and not h and not v"` (example 08).
- **Capture reference context as data** instead of writing one query per context. Put `^` on neighbouring columns of a free (`~`) window. `capture_refr` then holds the observed reference symbols, including `_` and IUPAC codes, in read orientation (examples 05 and 08).
- **Tell a contig end from a read end.** On the reference side, `_` is true only off the contig. On the read side, `_` is true off the read. `read="N_" refr="~N"` fires at a read's last base only when the contig continues. `read="N_" refr="~_"` fires only when the read ends exactly at the contig end, or runs past it. **Put a read base (`N`), not `~`, in the first column:** `~` also admits a pad, and with the default of two pads per end for these span-2 queries, example 05's `after` (`~_@~_`) also fires on the all-pad windows off the contig (`start` at `refr_pos -2`, `end` at 8, `overhang` at 10) and `after_on_contig` (`~_@~N`) on all-pad windows over the contig (`end` at 3, `overhang` at 4, `start` at 3). With `N_` only the read-end hits remain: `end` at 7 and `overhang` at 9 for `N_@~_`, `start` at 2 for `N_@~N` (checked by running 0.1.8 on example 05's reads, not in `expected.txt`).
- **Match an ambiguity code on both strands.** On an OB read, the reference side is complemented (`R`→`Y`). Use complement-closed sets (`S`, `W`, `N`), or write both (`R` or `Y`), when the query must fire on both strands (examples 01 and 08).
- **`=`/`/` as a "both sides are real, unambiguous bases" filter.** `~@= or ~@/` excludes read `N`, reference `N`/IUPAC, gaps, pads, clips and junctions in one column (example 01). The compile warnings it prints are expected.
- **Multi-member sets in rows:** use aliases, not groups (§3.2). Remember that an alias letter inside `{}` is *not* the alias.
- **Compare an index with a FASTA:**

  ```
  diff <(alnbase info --seq R idx | tr a-z A-Z) <(samtools faidx fa R | tr a-z A-Z)
  ```

  Expect `T` where the FASTA has `U` (§7.2).
- **Check an index against a sequence dictionary** without running a query. `info --all` prints the digests `samtools dict` writes as `M5` (example 04 shows them agreeing):

  ```
  diff <(alnbase info --all idx | awk 'NR>7 {print $1, $3}' | sort) \
       <(samtools dict fa | awk -F'\t' '$1=="@SQ" {sub("SN:","",$2); sub("M5:","",$4); print $2, $4}' | sort)
  ```

---

## 11. Discrepancies

Differences between the code and its comments, the docs (`docs/alnbase.md`, `docs/cli-reference.md`, `README.md`) or `--help`. Rows marked **Fixed in 0.1.2** are kept, struck through, for reference. Severity: **High** = could lead to wrong conclusions or silent wrong output; **Medium** = misleading; **Low** = cosmetic or stale.

| # | Where | Says | Code / behaviour | Evidence | Severity |
|---|---|---|---|---|---|
| D1 | `src/contig_map.rs:13-18` (module doc) | A contig the reference lacks "is not an error… records on it cannot be scanned and are skipped and counted" | Skipping and counting happen only under `--permissive`. By default the first such record stops the run (`src/batch.rs:69-73`, called from `src/ordered.rs:263` and `src/parallel.rs:437`). `ContigMap::warning` also says "records on them will be skipped" (`src/contig_map.rs:210`) immediately before the run fails. | example 06 §4 | Medium |
| D2 | `docs/cli-reference.md` (`index`) | "Accepts plain or bgzip-compressed FASTA" | Plain gzip is accepted too (`src/aref.rs:81`, `100`). The code comment is right; the doc is incomplete. | example 03 | Low |
| D3 | ~~`src/seq.rs:20`~~ | ~~ERR — "no reference is available to say"~~ | **Fixed in 0.1.10.** ERR is gone: unrecognised FASTA bytes are refused by `index`, and bit 6 is CLIP. | example 02 | ~~Medium~~ |
| D4 | `--list-codes` (`src/codes.rs:1-7`: "Every character a pattern can contain"), `docs/alnbase.md` "Reference: codes", `query --help` code list | Lists `. - Z _ X , J ~ {} = /` and the IUPAC letters | `U` (=T) is also accepted (`src/seq.rs`). **Fixed in 0.1.10:** `%` (ERR) is gone, `~` matches every symbol, and `--list-codes` lists `: L` (clip). `~` still cannot match a read `=`, which `docs/alnbase.md` says. | examples 02 and 07 | ~~High~~ Low |
| D5 | ~~`--list-codes` "COUNTS AND MARKERS" (`src/codes.rs:78-85`)~~ | ~~`C7` repeats a column; `^` and `+` "attach to the column before them"~~ | **Fixed in 0.1.2.** The section is now "MARK ROW": `.`, `^`, `+` in the `mark` key, and "rows take no repeat counts" (`src/codes.rs:78-86`). | `--list-codes` output | ~~Medium~~ |
| D6 | `src/dsl.rs:53-54` | "only 19 of the 127 non-empty `Seq` values have a single-character form" | 255 non-empty values, 20 single-character names. `--list-codes` correctly says "20 of the 255". | `--list-codes` output | Low |
| D7 | ~~`src/aref.rs:161-163`, `236`, `241`~~ | ~~`version` / `n_contigs` "not yet built… TODO(4.1): `alnbase info`", marked `dead_code`~~ | **Fixed in 0.1.2.** The TODO comments and `dead_code` attributes are gone; the field doc says "reported by `alnbase info`" (`src/aref.rs:204`). | example 04 | ~~Low~~ |
| D8 | ~~`src/main.rs:558` end-of-run summary~~ | ~~"skipped N (written untagged)" in `--parquet` mode~~ | **Fixed in 0.1.2.** Parquet and IPC runs say "(given a null row)"; only a tagged BAM says "(written untagged)" (`src/main.rs:693-697`). | examples 06 §1, 04-outputs/03 | ~~Low~~ |
| D9 | `docs/cli-reference.md` and `docs/alnbase.md` (`info --seq`), `info --help` | "The output is uppercase… compare case-insensitively" with `samtools faidx` | True except that `U` prints as `T`. Since 0.1.10 an index holds nothing else that could differ. | example 02 | Low |
| D10 | `docs/alnbase.md` "The model", `--list-codes` | Pad is "past the end of the read" / "off the end of the read, or off the contig" | On the **reference** side PAD means only "off the contig". It can also appear under *aligned* read bases when an alignment overruns the contig end, which the docs do not mention. | example 05 | Low |
| D11 | `docs/cli-reference.md` "Data you did not expect", `docs/alnbase.md` `query` | Off-reference records stop the run; `--permissive` counts them | Not true of `query --trace-records`, which silently skips them (`src/main.rs:451-457`), or of `overlap --refr`, which silently strips NM/MD. | example 06 §7 | Low |
| D12 | ~~`README.md`, `docs/cli-reference.md` ("Use the FASTA the reads were aligned against"), `src/contig_map.rs:15-18`~~ | ~~Length disagreement is the assembly check~~ | **Fixed in 0.1.2.** `@SQ M5` is checked against a per-contig MD5 stored in the index, when the header carries it (§9.1); the user docs now say a header without `M5` cannot be checked. The module doc at `src/contig_map.rs:13-18` still mentions only lengths (see D16). | example 06 §2 | ~~Medium~~ |
| D13 | `query --help` positional `[REFR]` | (empty help text) | It is the `.aref` index, not a FASTA; passing a FASTA fails with "bad magic". | `query --help` | Low |
| D14 | ~~`docs/cli-reference.md` "Exit status"~~ | ~~"A failed `query` run writing a BAM removes its own partial output", silent about `--parquet`~~ | **Fixed in 0.1.2.** A failed parquet/IPC run removes every output file it may have written (`src/parallel.rs:513-533`), and `docs/cli-reference.md` says so. | example 06 §4; 04-outputs/06 §2 | ~~Medium~~ |
| D15 | `src/aref.rs:9-24` (module doc, on-disk layout) | Table entries are "32 B × n_contigs" with four u64 fields; the table is at `len - 24 - n_contigs*32` | Version 2 entries are 48 bytes: the four u64s, then 16 MD5 bytes (`src/aref.rs:44`, `176-181`, `257-262`). The comment was not updated for version 2. | example 04 | Low |
| D16 | `src/contig_map.rs:13-18` (module doc), `src/contig_map.rs:55-58` (`build` doc) | `build` "fails if the two share no contig at all, or if a shared name has two different lengths" | It also fails when a resolved contig's `@SQ M5` disagrees with the index (`src/contig_map.rs:91-100`, `122-129`). | example 06 §2 | Low |
