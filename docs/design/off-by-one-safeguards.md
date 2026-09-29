# Off-by-one and coordinate safeguards

Status: **resolved** (0.1.8). Items marked *done* shipped in 0.1.2. Outcomes of the
open risks below:

1. Per-query reach was rejected; queries fire on all-pad windows, with
   `--end-context`/`--splice-context` settable and `and not` over an all-pad pattern
   (0.1.7). Both default to the widest query's span (0.1.8).
2. A: concordance check, shared across workers (0.1.4, 0.1.6), also applied to
   `--trace-records` (0.1.8). B: `--require-m5` (0.1.8).
3. Offsets include hard clips; `hard_clip_5p`/`hard_clip_3p` fields (0.1.4). Pads carry no
   offsets, and `soft_clip_5p`/`soft_clip_3p` fields give offsets from the first aligned base
   (0.1.9).
   Soft-clipped bases are clip columns (`:`) in the flank, with their own offsets, so a pattern
   can tell a clipped end from a true read end (0.1.10).
4. Null quality for QUAL `*` (0.1.4).
5. Coordinate base labelled in human-facing output (0.1.4); `coordinate_base = 0` in
   every hit table's metadata, with a run manifest (0.1.8).
6. `--explain` note under queries that place a pad (0.1.8).
7. Walk oracle (0.1.4, `tests/oracle/`; pads, intron context and the marker added
   after 0.1.8, and run in CI); tag/extract round trip is
   `extraction_matches_a_direct_parquet_run`, extended to deletions and hard clips
   (0.1.8); Bismark tags are `bismark_xm_xr_xg_on_all_four_strands`; writer
   conventions belong to the aggregation utility's validation.

The proposals are kept below as written. Sources: the code-verified reference in `docs/reference/` (its Discrepancies
sections) and probes run against 0.1.2.

The failure mode these guard against is the worst kind: output that looks like data.
A shifted coordinate is a real coordinate, a wrong reference base is a real base, and
nothing crashes. So the aim is (1) make the conventions impossible to misread,
(2) detect mismatches at run time, cheaply and loudly, and (3) pin every convention with
a test that fails if it moves.

## Found and fixed in 0.1.2

| Risk | What happened | Fix |
|---|---|---|
| Read 2 offsets | `off_5p = 0` was read 2's **last** sequenced base (offsets followed the conversion strand). M-bias trimming or read-end logic on read 2 was mirrored. | *done*: offsets are as sequenced for every FLAG; golden test pins all four flag cases. |
| Same name, same length, different sequence | Undetectable; a masked analysis set or patched assembly silently gave wrong reference bases. | *done*: `.aref` v2 stores MD5 per contig and checks `@SQ M5` when present. |
| Soft clips as insertions | `--soft-clips emit` made a clipped base look exactly like an insertion, with the nearest 5' coordinate. | *done*: clipped bases are never aligned columns; they are clip columns in the flank (0.1.10). |

## Open risks, with proposed safeguards

### 1. A query's hits depend on the *other* queries in the run — decision needed

The flank (pad columns past each read end) is `widest query − 1`, shared by every query.
A narrow query can fire entirely inside the flank, so adding an unrelated wider query
changes its hits. Verified on 0.1.2: `read="_" refr="A"` fires **0** times alone and
**8** times on the same BAM once a 6-column query is added. The query file stored in
the BAM header does not record which other files were in the run's width calculation
in an obvious way, and a validation run that adds a query changes every other result.

**Proposal:** each query sees only its own reach. A window may start at most
`span − 1` columns before the first aligned column and end at most `span − 1` after the
last, whatever the run's widest query. Equivalently: a query fires only if its window
contains at least one aligned (non-pad) read column. Either makes every query's hits a
function of that query alone. (`--splice-context` has the same shape and would get the
same rule.)

### 2. Reference mismatch when the header has no M5 — flag proposed

Most aligners do not write `@SQ M5`, so the 0.1.2 check usually has nothing to compare.

**Proposal A (default on): a concordance check.** While walking, count aligned columns
where read and reference are both unambiguous bases and differ, *excluding* reference
C / read T (the conversion mismatch; the walk's orientation makes that the only one on
every strand). On correct data this is sequencing error plus variants, typically under
2%. A wrong reference, a wrong contig, or a coordinate shift gives about 75%. Report
the rate in the run summary; stop with an error when it exceeds a threshold after the
first 100,000 aligned bases. `--max-discordance F` (default 0.25) adjusts it, and `1`
disables it. Cost: one comparison per column, already in hand.

**Proposal B: `--require-m5`.** Refuse a BAM whose header lacks M5 for any contig that
has records, with a message pointing at `samtools dict` + `samtools reheader`. Off by
default.

### 3. Hard clips shift offsets on supplementary alignments — decision needed

Offsets count SEQ, which excludes hard-clipped bases. A supplementary alignment with
`60H20M` reports `off_5p = 0` at what is really the read's 61st sequenced base. For
Methyl-HiC, where supplementary alignments carry many calls, M-bias computed from
`off_5p` mixes true read ends with chimera junctions.

**Proposal:** add record fields `hard_clip_5p` / `hard_clip_3p` (hard-clipped bases at the
sequenced 5' and 3' ends) so `off_5p + hard_clip_5p` is the offset in the original read,
and say so in the `off_5p` help. Alternatively an `--offsets original-read` flag that
adds them in. The field is the smaller change and keeps one meaning per column.

### 4. Missing quality reads as 255 — fix proposed

A record with QUAL `*` reports `qual = 255` rather than null, so `qual >= 20` keeps it;
`qual_phred` is a string of `ÿ` characters. Deletion and pad columns already use null.

**Proposal:** null for QUAL `*`, matching the "no read base, no quality" convention.

### 5. Mixed coordinate bases in human-facing output — fix proposed

`refr_pos` and all parquet coordinates are 0-based; `info --seq` and the
`--trace-records` record label are 1-based (to match `samtools faidx` and SAM POS).
Each is defensible alone; together they invite an off-by-one when copying a position
from a trace into a query of the parquet.

**Proposal:** label the base wherever a coordinate is printed for people
(`pos 1234 (1-based)`), and store `coordinate_base = 0` in the parquet footer next to
`format_version`, so downstream code can assert it.

### 6. Pattern direction vs offset direction for read 2 — documentation

After 0.1.2, offsets are as sequenced but patterns still run along the conversion strand:
for read 2, `_N` matches at the sequenced **3'** end. This is intentional (contexts must
read along the strand that was converted) but is the next most likely confusion.

**Proposal:** state it in `--explain` output whenever a query contains a pad, and in the
reference (in progress).

### 7. Tests that pin conventions — to add

Unit tests pin individual decisions; these would pin the conventions end to end:

- **Walk oracle.** For randomly generated records (all four flag cases, random CIGARs
  with I/D/N/S/H), check every emitted column against an independent derivation using
  pysam's `get_aligned_pairs(with_seq=False)` and `pyfaidx`: `refr_pos`, the reference
  base (complemented on the bottom strand), the read base, `off_5p` against the
  sequenced read, and `off_5p + off_3p == read_len − 1`. Independent code is the point:
  a test that reuses alnbase's own arithmetic cannot catch it being wrong.
- **Tag / extract round trip.** For a bases-tag query set, the rows of
  `query --parquet` restricted to the tagged queries equal the rows of `extract` on the
  tagged BAM, for every flag case.
- **Bismark tag equivalence** (already an example in `docs/reference/examples/04-outputs`):
  promote to a test.
- **Writer conventions.** Every aggregation writer that emits 1-based or half-open
  coordinates is compared against the real tool's file on the same data, in the
  validation suite, not against our reading of the tool's documentation.

## Recommended order

1 and 2A first (they change results or catch silent corruption), then 4, 3, 5, the
tests, and 6 alongside the documentation.
