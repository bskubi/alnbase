//! Build the hits file, which holds one denormalized row for each query that
//! fired.
//!
//! Three properties of the schema explain the rest of it:
//!
//! - A row identifies its query by `name`. What that name means is in the
//!   query file the run stored in the BAM header, which `alnbase dump-query`
//!   prints -- so a row stays interpretable without the command line, without
//!   spelling the whole query out on every row.
//! - `read_5p` / `read_3p` / `refr_pos` / `qual` describe the **anchor**
//!   column, not the pattern's first column. With the default anchor those are
//!   the same column, so single-pattern output reads as you would expect.
//!
//! # The two read offsets
//!
//! `read_5p` counts from the 5' end of the read as sequenced, and `read_3p`
//! counts from the 3' end. The last sequenced base is therefore `read_3p == 0`.
//!
//! Both are written because they localise different artefacts. End-repair and
//! adapter read-through sit a fixed distance from the 3' end. The 5' bias from
//! random priming sits a fixed distance from the other end. A consumer can
//! derive either offset from the other and `read_len`, but only if it selected
//! that column and knows that the derivation runs over SEQ. One subtraction per
//! row here removes that step.
//!
//! Both count from the ends of the read *as sequenced*, whatever order the walk
//! emitted the columns in. The walk runs along the conversion strand, which is
//! the sequencing order only for a read sequenced in that direction. For the
//! rest, the walk offsets are mirrored when written; see `as_sequenced`. The
//! run's strand call says which reads those are, and for a directional library
//! it is read 2. `read_5p == 0` is therefore always the first base the sequencer
//! read, whatever the aligner put in the FLAG.
//!
//! Soft-clipped bases that the walk skipped are still counted, because the
//! sequencer read them. Both offsets are measured against the length of SEQ,
//! which excludes hard-clipped bases. For a hard-clipped record neither offset
//! refers to the full original read, and `read_5p + read_3p == read_len - 1`
//! still holds.
//!
//! A pad is not a base of the read, so a row anchored on one has both offsets
//! null. The same is true past a clip, and on a supplementary alignment.
//!
//! # Flat versus list captures
//!
//! The anchor is always recorded, so every query captures at least one column.
//! For most queries that is all they capture, because most do not mark an extra
//! `^`. The read offset, reference position and quality of that capture are the
//! anchor's, and the row already holds all three as scalars. Only the observed
//! bases are new, so the file gets two extra scalar columns and no lists.
//!
//! When any query in the run captures more, the whole file uses LIST columns.
//! The capture count varies by query, and one shard has one schema. Offsets and
//! positions are then written for each capture, and not derived from the anchor.
//! Across an indel the read offset and the reference position do not advance
//! together, so the distance from the anchor is not a constant.

use std::sync::Arc;

use arrow::array::{
    ArrayRef, Int32Builder, Int64Builder, ListBuilder, StringBuilder, UInt32Builder,
    UInt64Builder, UInt8Builder,
};
use arrow::datatypes::{DataType, Field, Schema, SchemaRef};

use crate::column::Column;
use crate::arrow_builder::ArrowBuilder;
use crate::encoding::ColumnEncoding;
use crate::column_ring::ColumnRing;
use crate::query::Query;
use crate::record_columns::RecordColumns;
use crate::record_field::RecordField;
use crate::seq::Seq;

/// Bumped whenever the column set, the meaning of a column, or the metadata
/// changes. Written into every file's metadata (the parquet footer and the
/// Arrow schema) beside `coordinate_base`, so a reader can tell what it has.
/// No history of earlier versions is kept before the first release.
pub const FORMAT_VERSION: &str = "10";

/// Capture columns, in whichever shape the query set allows.
enum Captures {
    /// Every query records only its anchor: two scalars, no lists.
    Flat { read: StringBuilder, refr: StringBuilder },
    /// At least one query records more: parallel lists.
    Lists {
        col: ListBuilder<Int32Builder>,
        read: ListBuilder<StringBuilder>,
        refr: ListBuilder<StringBuilder>,
        qual: ListBuilder<UInt8Builder>,
        read_off: ListBuilder<Int64Builder>,
        read_read_3p: ListBuilder<Int64Builder>,
        refr_pos: ListBuilder<Int64Builder>,
    },
}

pub struct HitBuilder {
    record_id: UInt64Builder,
    name: StringBuilder,
    read_off: Int64Builder,
    /// The same offset from the other end; see the module comment.
    read_read_3p: Int64Builder,
    refr_pos: Int64Builder,
    qual: UInt8Builder,

    captures: Captures,

    record_cols: RecordColumns,
    /// Which of the run's output files the row is in, counted across every
    /// worker: `worker * shards_per_worker + shard`, the slot `partition::Router`
    /// routed to. Constant within a file, but carried per row so the data says
    /// where it came from rather than the filename: a consumer reading several
    /// files off one stream needs to demultiplex them, and a file's rows are
    /// only ordered relative to each other.
    shard: UInt32Builder,
    shard_id: u32,
    rid: u64,
    /// SEQ's length for the record being emitted, which is what `read_off` is
    /// an offset into. Set once per record by `begin_record`, for the same
    /// reason `RecordColumns` extracts once per record: every row this record
    /// emits needs it, and reading it back off the record per row would be the
    /// same value fetched a dozen times.
    read_len: i64,
    /// Whether the walk ran opposite the direction the read was sequenced, in
    /// which case its walk offsets are mirrored before they are written. Taken
    /// from the strand call, not from the FLAG: for a directional library it is
    /// read 2, and for an aligner that says otherwise it is what the aligner
    /// says. Set by `begin_record` alongside `read_len`.
    mirror: bool,
    /// Hard-clipped bases at the read's sequenced 5' and 3' ends, added to the
    /// SEQ-based offsets so they are relative to the read rather than the record.
    hard_5p: i64,
    hard_3p: i64,
    rows: usize,
}

/// `ListBuilder` names its child field `item` and makes it nullable; the
/// declared schema has to say the same thing or `RecordBatch::try_new` rejects
/// the batch.
fn list_field(name: &str, item: DataType) -> Field {
    Field::new(name, DataType::List(Arc::new(Field::new("item", item, true))), true)
}

impl HitBuilder {
    /// `fields` is already resolved -- see `RecordField::resolve`.
    ///
    /// `flat` comes from `QuerySet::flat_captures` and must be the same for
    /// every worker, or the shards would not share a schema.
    pub fn new(fields: &[RecordField], cap: usize, flat: bool, shard_id: u32) -> Self {
        Self {
            record_id: UInt64Builder::with_capacity(cap),
            name: StringBuilder::with_capacity(cap, cap * 8),
            read_off: Int64Builder::with_capacity(cap),
            read_read_3p: Int64Builder::with_capacity(cap),
            refr_pos: Int64Builder::with_capacity(cap),
            qual: UInt8Builder::with_capacity(cap),

            captures: if flat {
                Captures::Flat {
                    read: StringBuilder::with_capacity(cap, cap * 2),
                    refr: StringBuilder::with_capacity(cap, cap * 2),
                }
            } else {
                Captures::Lists {
                    col: ListBuilder::new(Int32Builder::new()),
                    read: ListBuilder::new(StringBuilder::new()),
                    refr: ListBuilder::new(StringBuilder::new()),
                    qual: ListBuilder::new(UInt8Builder::new()),
                    read_off: ListBuilder::new(Int64Builder::new()),
                    read_read_3p: ListBuilder::new(Int64Builder::new()),
                    refr_pos: ListBuilder::new(Int64Builder::new()),
                }
            },

            record_cols: RecordColumns::new(fields, cap),
            shard: UInt32Builder::with_capacity(cap),
            shard_id,
            rid: 0,
            read_len: 0,
            mirror: false,
            hard_5p: 0,
            hard_3p: 0,
            rows: 0,
        }
    }

    /// Read every selected column off the record, once. Each row this record
    /// emits then copies from that one extraction.
    pub fn begin_record(
        &mut self,
        rid: u64,
        record: &rust_htslib::bam::Record,
        header: &rust_htslib::bam::HeaderView,
        strand: Option<crate::strand_rule::StrandCall>,
    ) {
        self.rid = rid;
        // `walk_alignment` measures `read_off` against this same length -- it
        // decodes SEQ and counts every base including skipped clips -- so the
        // two ends of the read are measured against one number.
        self.read_len = record.seq_len() as i64;
        // `None` is a record the run's rule declined to call. It reaches here
        // only as a row that carries no walk -- an unmapped record under
        // `--keep-unscanned` -- so there are no offsets to mirror, and the
        // strand columns are written null rather than guessed.
        self.mirror = strand.is_some_and(|s| s.mirrors_read());
        (self.hard_5p, self.hard_3p) = crate::alignment::hard_clips_as_sequenced(record);
        self.record_cols.extract(record, header, strand);
    }

    /// One row for a query that fired.
    ///
    /// `ring` is read through the back-indices `Query` resolved at compile time;
    /// no arithmetic happens here.
    pub fn push_hit(&mut self, q: &Query, ring: &ColumnRing) {
        let anchor = ring.get(q.anchor_back);

        self.shard.append_value(self.shard_id);
        self.record_id.append_value(self.rid);
        self.name.append_value(&q.name);
        // A pad is past the read: no offset into the read describes it, so it
        // gets none. A clip is a sequenced base and keeps its offset.
        if anchor.read == Seq::PAD {
            self.read_off.append_null();
            self.read_read_3p.append_null();
        } else {
            let (five, three) = as_sequenced(self.read_len, self.mirror, anchor.read_off, self.hard_5p, self.hard_3p);
            self.read_off.append_value(five);
            self.read_read_3p.append_option(three);
        }
        self.refr_pos.append_value(anchor.refr_pos);
        append_qual(&mut self.qual, anchor);

        // Copied out before the match so the borrow of `self.captures` below is
        // the only one taken.
        let (read_len, mirror, hard_5p, hard_3p) = (self.read_len, self.mirror, self.hard_5p, self.hard_3p);
        match &mut self.captures {
            // The anchor is the only capture, and we are already holding it.
            Captures::Flat { read, refr } => {
                read.append_value(anchor.read.name());
                refr.append_value(anchor.refr.name());
            }
            Captures::Lists { col, read, refr, qual, read_off, read_read_3p, refr_pos } => {
                for (&back, &col_idx) in q.capture_back.iter().zip(q.captures.iter()) {
                    let c = ring.get(back);
                    read.values().append_value(c.read.name());
                    refr.values().append_value(c.refr.name());
                    // Null, not -1, for a column with no quality: -1 would
                    // survive a `qual >= 30` filter.
                    append_qual(qual.values(), c);
                    // A column with no sequenced base has no offset here: a
                    // deletion, an intron column, a pad. A clip is a sequenced
                    // base. Decided by the symbol, not the quality, which a
                    // record with QUAL `*` lacks everywhere.
                    if c.read.0 & (Seq::GAP.0 | Seq::SKIP.0 | Seq::PAD.0) == 0 {
                        let (five, three) = as_sequenced(read_len, mirror, c.read_off, hard_5p, hard_3p);
                        read_off.values().append_value(five);
                        read_read_3p.values().append_option(three);
                    } else {
                        read_off.values().append_null();
                        read_read_3p.values().append_null();
                    }
                    refr_pos.values().append_value(c.refr_pos);
                    col.values().append_value(col_idx as i32);
                }
                col.append(true);
                read.append(true);
                refr.append(true);
                qual.append(true);
                read_off.append(true);
                read_read_3p.append(true);
                refr_pos.append(true);
            }
        }

        self.record_cols.append_row();
        self.rows += 1;
    }

    /// One row for a record that matched nothing, so the file stays a complete
    /// inventory of what was scanned.
    /// A hit read back out of a per-base tag rather than found by a walk: the
    /// anchor's position, quality and read base are known, and the reference
    /// base is known when the query determines it. Only for the flat capture
    /// layout -- a tag records the anchor and nothing else a query captured.
    ///
    /// `read_off` is in walk order, like a walked hit's, and is converted to
    /// the as-sequenced `read_5p` / `read_3p` the same way.
    pub fn push_decoded(
        &mut self,
        name: &str,
        read_off: i64,
        refr_pos: Option<i64>,
        qual: Option<u8>,
        read: &str,
        refr: Option<&str>,
    ) {
        self.shard.append_value(self.shard_id);
        self.record_id.append_value(self.rid);
        self.name.append_value(name);
        let (five, three) = as_sequenced(self.read_len, self.mirror, read_off, self.hard_5p, self.hard_3p);
        self.read_off.append_value(five);
        self.read_read_3p.append_option(three);
        self.refr_pos.append_option(refr_pos);
        self.qual.append_option(qual);
        match &mut self.captures {
            Captures::Flat { read: r, refr: f } => {
                r.append_value(read);
                f.append_option(refr);
            }
            Captures::Lists { .. } => {
                unreachable!("decoded hits are only written with flat captures")
            }
        }
        self.record_cols.append_row();
        self.rows += 1;
    }

    /// What the last [`begin_record`](Self::begin_record) found wrong with
    /// the record's aux block.
    pub fn aux_trouble(&self) -> crate::record_columns::AuxTrouble {
        self.record_cols.aux_trouble()
    }

    pub fn push_hitless(&mut self) {
        self.shard.append_value(self.shard_id);
        self.record_id.append_value(self.rid);
        self.name.append_null();
        self.read_off.append_null();
        self.read_read_3p.append_null();
        self.refr_pos.append_null();
        self.qual.append_null();
        match &mut self.captures {
            Captures::Flat { read, refr } => {
                read.append_null();
                refr.append_null();
            }
            // A null list, not an empty one: nothing was examined here, which
            // is different from a query that captured nothing.
            Captures::Lists { col, read, refr, qual, read_off, read_read_3p, refr_pos } => {
                col.append(false);
                read.append(false);
                refr.append(false);
                qual.append(false);
                read_off.append(false);
                read_read_3p.append(false);
                refr_pos.append(false);
            }
        }
        self.record_cols.append_row();
        self.rows += 1;
    }

    /// Order must match `arrow_fields`.
    pub fn finish(&mut self) -> Vec<ArrayRef> {
        let mut cols: Vec<ArrayRef> = Vec::with_capacity(13 + self.record_cols.width());
        cols.push(Arc::new(self.shard.finish()));
        cols.push(Arc::new(self.record_id.finish()));
        cols.extend(self.record_cols.finish());
        cols.push(Arc::new(self.name.finish()));
        cols.push(Arc::new(self.read_off.finish()));
        cols.push(Arc::new(self.read_read_3p.finish()));
        cols.push(Arc::new(self.refr_pos.finish()));
        cols.push(Arc::new(self.qual.finish()));
        match &mut self.captures {
            Captures::Flat { read, refr } => {
                cols.push(Arc::new(read.finish()));
                cols.push(Arc::new(refr.finish()));
            }
            Captures::Lists { col, read, refr, qual, read_off, read_read_3p, refr_pos } => {
                cols.push(Arc::new(col.finish()));
                cols.push(Arc::new(read.finish()));
                cols.push(Arc::new(refr.finish()));
                cols.push(Arc::new(qual.finish()));
                cols.push(Arc::new(read_off.finish()));
                cols.push(Arc::new(read_read_3p.finish()));
                cols.push(Arc::new(refr_pos.finish()));
            }
        }
        self.rows = 0;
        cols
    }

    pub fn len(&self) -> usize {
        self.rows
    }

    /// Order must match `finish`.
    pub fn arrow_fields(&self) -> Vec<Field> {
        let mut fields = vec![
            Field::new("shard", DataType::UInt32, false),
            Field::new("record_id", DataType::UInt64, false),
        ];
        fields.extend(self.record_cols.arrow_fields());
        fields.extend([
            Field::new("name", DataType::Utf8, true),
            Field::new("read_5p", DataType::Int64, true),
            Field::new("read_3p", DataType::Int64, true),
            Field::new("refr_pos", DataType::Int64, true),
            Field::new("qual", DataType::UInt8, true),
        ]);
        match &self.captures {
            // The anchor's own coordinates and quality are the scalars above,
            // so only the observed bases are added.
            Captures::Flat { .. } => fields.extend([
                Field::new("read_base", DataType::Utf8, true),
                Field::new("refr_base", DataType::Utf8, true),
            ]),
            Captures::Lists { .. } => fields.extend([
                list_field("capture_col", DataType::Int32),
                list_field("capture_read", DataType::Utf8),
                list_field("capture_refr", DataType::Utf8),
                list_field("capture_qual", DataType::UInt8),
                list_field("capture_read_5p", DataType::Int64),
                list_field("capture_read_3p", DataType::Int64),
                list_field("capture_refr_pos", DataType::Int64),
            ]),
        }
        fields
    }

    /// Each column's storage policy, in schema order and by schema name.
    ///
    /// Stated here rather than in the writer for the same reason the schema is:
    /// this is where the columns are, and a table living anywhere else would be
    /// a second list to keep in step. Columns absent from this list get the
    /// writer's default, so adding a column without an opinion is not a bug.
    ///
    /// The reasoning is the shape of the file. `shard` is constant within one,
    /// `name` has one distinct value per query, `qual` has ninety
    /// odd and the base columns five -- all dictionaries. `record_id` only ever
    /// increases, and a record's rows are contiguous and ordered 5'->3', so
    /// `read_5p`, `read_3p` and `refr_pos` each move in small steps within
    /// a record and jump only between records: delta.
    pub fn column_encodings(&self) -> Vec<(String, ColumnEncoding)> {
        use ColumnEncoding as E;
        let mut out: Vec<(String, ColumnEncoding)> = vec![
            ("shard".into(), E::Dictionary),
            ("record_id".into(), E::Delta),
        ];
        out.extend(self.record_cols.column_encodings());
        out.extend([
            ("name".into(), E::Dictionary),
            ("read_5p".into(), E::Delta),
            ("read_3p".into(), E::Delta),
            ("refr_pos".into(), E::Delta),
            ("qual".into(), E::Dictionary),
        ]);
        match &self.captures {
            Captures::Flat { .. } => out.extend([
                ("read_base".into(), E::Dictionary),
                ("refr_base".into(), E::Dictionary),
            ]),
            // The list columns' policies apply to the leaf inside the list, not
            // to the list itself; the writer derives that path from the schema.
            Captures::Lists { .. } => out.extend([
                ("capture_col".into(), E::Dictionary),
                ("capture_read".into(), E::Dictionary),
                ("capture_refr".into(), E::Dictionary),
                ("capture_qual".into(), E::Dictionary),
                ("capture_read_5p".into(), E::Delta),
                ("capture_read_3p".into(), E::Delta),
                ("capture_refr_pos".into(), E::Delta),
            ]),
        }
        out
    }

    pub fn arrow_schema(&self) -> SchemaRef {
        let meta = std::collections::HashMap::from([
            ("format_version".to_string(), FORMAT_VERSION.to_string()),
            ("coordinate_base".to_string(), crate::manifest::COORDINATE_BASE.to_string()),
        ]);
        Arc::new(Schema::new(self.arrow_fields()).with_metadata(meta))
    }
}

/// The two written offsets for a column at walk offset `walk_off`.
///
/// `read_5p` counts from the first base the sequencer read and `read_3p` from the
/// last, over the whole read: SEQ plus any hard-clipped bases, which are not in
/// the record but were sequenced (`hard_5p` at the 5' end, `hard_3p` at the 3'
/// end). So a supplementary alignment reports the same offsets for a base as the
/// primary alignment of the same read, and `read_5p + read_3p` is the read's length
/// minus 1. The walk emits columns along the conversion strand (see
/// `walk_alignment`), which is the sequencing order only when the read was
/// sequenced in that direction; `mirror` says it was not, and the walk offset
/// is counted from the other end. See `StrandCall::mirrors_read`.
///
/// Flank columns lie outside the read, and the arithmetic carries that through
/// without a special case: a column before the first sequenced base has a
/// negative `read_5p` and an `read_3p` past the far end.
///
/// A record with no SEQ (`*`, so `seq_len` is 0) has no ends to measure from:
/// `read_5p` is the walk offset unchanged and `read_3p` is `None`.
#[inline]
fn as_sequenced(seq_len: i64, mirror: bool, walk_off: i64, hard_5p: i64, hard_3p: i64) -> (i64, Option<i64>) {
    if seq_len <= 0 {
        return (walk_off, None);
    }
    let five_in_seq = if mirror { seq_len - 1 - walk_off } else { walk_off };
    let three_in_seq = seq_len - 1 - five_in_seq;
    (five_in_seq + hard_5p, Some(three_in_seq + hard_3p))
}

/// `Column::qual` carries -1 for "no read base here". Null says the same thing
/// and, unlike -1, drops out of a `qual >= 30` filter on its own.
#[inline]
fn append_qual(b: &mut UInt8Builder, c: &Column) {
    if c.qual >= 0 {
        b.append_value(c.qual as u8);
    } else {
        b.append_null();
    }
}

impl ArrowBuilder for HitBuilder {
    fn finish(&mut self) -> Vec<ArrayRef> {
        HitBuilder::finish(self)
    }
    fn len(&self) -> usize {
        HitBuilder::len(self)
    }
    fn schema(&self) -> SchemaRef {
        self.arrow_schema()
    }
    fn column_encodings(&self) -> Vec<(String, ColumnEncoding)> {
        HitBuilder::column_encodings(self)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::Column as AlnColumn;
    use crate::column_ring::ColumnRing;
    use crate::query::QuerySet;
    use crate::seq::Seq;
    use arrow::array::{
        Array, Int32Array, Int64Array, ListArray, StringArray, UInt32Array, UInt64Array,
        UInt8Array,
    };

    fn queries(toml: &str) -> QuerySet {
        let specs = crate::test_support::queries(toml);
        QuerySet::compile(&specs).unwrap()
    }

    fn col(read: Seq, refr: Seq, pos: i64, off: i64, qual: i32) -> AlnColumn {
        AlnColumn { read, refr, refr_pos: pos, read_off: off, qual }
    }

    /// A ring holding `cols`, stored in order, as the scanner would leave it
    /// having just read the last one.
    fn ring_of(span: usize, cols: &[AlnColumn]) -> ColumnRing {
        let mut r = ColumnRing::new(span);
        for c in cols {
            r.store(*c);
        }
        r
    }

    fn field_names(b: &HitBuilder) -> Vec<String> {
        b.arrow_fields().iter().map(|f| f.name().clone()).collect()
    }

    /// The schema `finish` fills must match the one `arrow_schema` declares, or
    /// `RecordBatch::try_new` fails at the first flush rather than at build
    /// time. Nothing but a test compares them.
    #[test]
    fn the_declared_schema_matches_the_arrays_produced() {
        for flat in [true, false] {
            let mut b = HitBuilder::new(RecordField::CORE, 4, flat, 0);
            let arrays = b.finish();
            let fields = b.arrow_fields();
            assert_eq!(
                arrays.len(),
                fields.len(),
                "flat={flat}: {} arrays for {} fields",
                arrays.len(),
                fields.len()
            );
            for (a, f) in arrays.iter().zip(fields.iter()) {
                assert_eq!(
                    a.data_type(),
                    f.data_type(),
                    "flat={flat}: column '{}' has the wrong type",
                    f.name()
                );
            }
        }
    }

    /// Anchor-only queries get scalar base columns; anything capturing more
    /// gets lists, because a shard has one schema.
    #[test]
    fn the_capture_shape_decides_the_columns() {
        let flat = field_names(&HitBuilder::new(RecordField::CORE, 4, true, 0));
        assert!(flat.contains(&"read_base".to_string()), "{flat:?}");
        assert!(!flat.iter().any(|n| n.starts_with("capture_")), "{flat:?}");

        let lists = field_names(&HitBuilder::new(RecordField::CORE, 4, false, 0));
        assert!(lists.contains(&"capture_col".to_string()), "{lists:?}");
        assert!(!lists.contains(&"read_base".to_string()), "{lists:?}");
    }

    /// The row reports the anchor column, reached through the back-index the
    /// query resolved at compile time -- not the last column stored.
    #[test]
    fn a_hit_row_reports_the_anchor_not_the_last_column() {
        // Anchor on column 0 of a two-column span: the C, not the G.
        let qs = queries("[query.m]\nread = \"Y~\"\nrefr = \"CG\"\n");
        let q = &qs.queries[0];
        assert_eq!(q.anchor, 0);

        let mut b = HitBuilder::new(RecordField::CORE, 4, qs.flat_captures, 0);
        let ring = ring_of(
            q.span,
            &[
                col(Seq::C, Seq::C, 100, 10, 30),
                col(Seq::G, Seq::G, 101, 11, 31),
            ],
        );
        b.push_hit(q, &ring);
        let arrays = b.finish();
        let names = field_names(&HitBuilder::new(RecordField::CORE, 4, qs.flat_captures, 0));
        let at = |n: &str| arrays[names.iter().position(|x| x == n).unwrap()].clone();

        let pos = at("refr_pos");
        let pos = pos.as_any().downcast_ref::<Int64Array>().unwrap();
        assert_eq!(pos.value(0), 100, "the anchor's coordinate, not the last column's");

        let read = at("read_base");
        let read = read.as_any().downcast_ref::<StringArray>().unwrap();
        assert_eq!(read.value(0), "C");
    }

    /// Captures are written in ascending column order with their own
    /// coordinates, since across an indel they do not track the anchor.
    #[test]
    fn list_captures_carry_their_own_coordinates() {
        let qs = queries("[query.m]\nread = \"Y~\"\nrefr = \"CG\"\nmark = \"+^\"\n");
        assert!(!qs.flat_captures, "two captures means the list form");
        let q = &qs.queries[0];
        assert_eq!(&*q.captures, &[0, 1]);

        let mut b = HitBuilder::new(RecordField::CORE, 4, false, 0);
        let ring = ring_of(
            q.span,
            &[
                col(Seq::C, Seq::C, 100, 10, 30),
                // A deletion: reference base, no read base, so no quality.
                col(Seq::GAP, Seq::G, 101, 11, -1),
            ],
        );
        b.push_hit(q, &ring);
        let arrays = b.finish();
        let names = field_names(&HitBuilder::new(RecordField::CORE, 4, false, 0));
        let at = |n: &str| arrays[names.iter().position(|x| x == n).unwrap()].clone();

        let cols = at("capture_col");
        let cols = cols.as_any().downcast_ref::<ListArray>().unwrap();
        let inner = cols.value(0);
        let inner = inner.as_any().downcast_ref::<Int32Array>().unwrap();
        assert_eq!(inner.values(), &[0, 1]);

        // A column with no read base has a null quality, not a sentinel: -1
        // would survive a `qual >= 30` filter.
        let qual = at("capture_qual");
        let qual = qual.as_any().downcast_ref::<ListArray>().unwrap();
        let inner = qual.value(0);
        let inner = inner.as_any().downcast_ref::<UInt8Array>().unwrap();
        assert!(!inner.is_null(0));
        assert!(inner.is_null(1), "the deletion's quality must be null");
    }

    /// The two offsets are measured from opposite ends of the same read, so
    /// they sum to `read_len - 1` -- which is the whole of what `read_read_3p`
    /// means, and the only thing a consumer can rely on.
    #[test]
    fn the_two_offsets_are_measured_from_opposite_ends() {
        let qs = queries("[query.m]\nread = \"Y~\"\nrefr = \"CG\"\n");
        let q = &qs.queries[0];
        let mut b = HitBuilder::new(RecordField::CORE, 4, qs.flat_captures, 0);
        b.read_len = 100;

        // Anchor at offset 10 of a 100-base read: 89 bases follow it.
        let ring = ring_of(
            q.span,
            &[
                col(Seq::C, Seq::C, 100, 10, 30),
                col(Seq::G, Seq::G, 101, 11, 31),
            ],
        );
        b.push_hit(q, &ring);
        let arrays = b.finish();
        let names = field_names(&HitBuilder::new(RecordField::CORE, 4, qs.flat_captures, 0));
        let at = |n: &str| arrays[names.iter().position(|x| x == n).unwrap()].clone();

        let five = at("read_5p");
        let five = five.as_any().downcast_ref::<Int64Array>().unwrap();
        let three = at("read_3p");
        let three = three.as_any().downcast_ref::<Int64Array>().unwrap();
        assert_eq!(five.value(0), 10);
        assert_eq!(three.value(0), 89);
        assert_eq!(five.value(0) + three.value(0), 99, "the two must mirror");
    }

    /// Offsets are as sequenced: read 1 keeps its walk offset, read 2's is
    /// mirrored, and a record with no SEQ has no 3' end to measure from.
    #[test]
    fn offsets_are_measured_from_the_ends_as_sequenced() {
        assert_eq!(as_sequenced(0, false, 4, 0, 0), (4, None));
        assert_eq!(as_sequenced(100, false, 0, 0, 0), (0, Some(99)), "read 1: walk start is its 5' end");
        assert_eq!(as_sequenced(100, false, 99, 0, 0), (99, Some(0)));
        assert_eq!(as_sequenced(100, true, 0, 0, 0), (99, Some(0)), "read 2: walk start is its 3' end");
        assert_eq!(as_sequenced(100, true, 99, 0, 0), (0, Some(99)));
        // Flank columns lie outside the read at either end, and the mirror
        // carries that through rather than clamping it away.
        assert_eq!(as_sequenced(100, false, -2, 0, 0), (-2, Some(101)));
        assert_eq!(as_sequenced(100, true, -2, 0, 0), (101, Some(-2)));
        // Hard clips were sequenced: 60 bases clipped off the 5' end and 5 off
        // the 3' end put SEQ's first base 60 from the read's start.
        assert_eq!(as_sequenced(20, false, 0, 60, 5), (60, Some(24)));
        assert_eq!(as_sequenced(20, true, 0, 60, 5), (79, Some(5)));
    }

    /// Captures null the 3' offset exactly where they null the 5' one: a
    /// deletion has no read base, so neither offset means anything there.
    #[test]
    fn a_capture_without_a_read_base_has_neither_offset() {
        let qs = queries("[query.m]\nread = \"Y~\"\nrefr = \"CG\"\nmark = \"+^\"\n");
        let q = &qs.queries[0];
        let mut b = HitBuilder::new(RecordField::CORE, 4, false, 0);
        b.read_len = 50;
        let ring = ring_of(
            q.span,
            &[
                col(Seq::C, Seq::C, 100, 10, 30),
                col(Seq::GAP, Seq::G, 101, 11, -1),
            ],
        );
        b.push_hit(q, &ring);
        let arrays = b.finish();
        let names = field_names(&HitBuilder::new(RecordField::CORE, 4, false, 0));
        let at = |n: &str| arrays[names.iter().position(|x| x == n).unwrap()].clone();

        for c in ["capture_read_5p", "capture_read_3p"] {
            let a = at(c);
            let a = a.as_any().downcast_ref::<ListArray>().unwrap();
            let inner = a.value(0);
            let inner = inner.as_any().downcast_ref::<Int64Array>().unwrap();
            assert!(!inner.is_null(0), "{c}: the matched base has an offset");
            assert!(inner.is_null(1), "{c}: the deletion must be null");
        }

        let a = at("capture_read_3p");
        let a = a.as_any().downcast_ref::<ListArray>().unwrap();
        let inner = a.value(0);
        let inner = inner.as_any().downcast_ref::<Int64Array>().unwrap();
        assert_eq!(inner.value(0), 39, "50 - 1 - 10");
    }

    /// A pad is not a base of the read, so a hit anchored on one has no offset,
    /// as a scalar or as a capture, whatever the walk carried for it.
    #[test]
    fn a_hit_anchored_on_a_pad_has_no_offset() {
        for flat in [true, false] {
            let src = if flat { "[query.p]\nread = \"_\"\nrefr = \"~\"\n" } else { "[query.p]\nread = \"_\"\nrefr = \"~\"\n[query.w]\nread = \"_~\"\nrefr = \"~~\"\nmark = \"+^\"\n" };
            let qs = queries(src);
            let q = &qs.queries[0];
            let mut b = HitBuilder::new(RecordField::CORE, 4, flat, 0);
            b.read_len = 50;
            b.hard_5p = 3;
            let ring = ring_of(q.span, &[col(Seq::PAD, Seq::G, 99, crate::column::PAD_READ_OFF, -1)]);
            b.push_hit(q, &ring);
            let names = field_names(&HitBuilder::new(RecordField::CORE, 4, flat, 0));
            let arrays = b.finish();
            let at = |n: &str| arrays[names.iter().position(|x| x == n).unwrap()].clone();
            for c in ["read_5p", "read_3p"] {
                assert!(at(c).is_null(0), "flat {flat}: {c} of a pad");
            }
            assert!(!at("refr_pos").is_null(0), "flat {flat}: a pad keeps its reference position");
            if !flat {
                let a = at("capture_read_5p");
                let a = a.as_any().downcast_ref::<ListArray>().unwrap().value(0);
                assert!(a.is_null(0), "the pad capture has no offset either");
            }
        }
    }

    /// A clip is a sequenced base, so it keeps its offset, as the anchor and as
    /// a capture, though it has no quality. Capture offsets follow the symbol,
    /// not the quality: an aligned base in a record with QUAL `*` has an offset.
    #[test]
    fn a_clip_and_a_base_without_quality_keep_their_offsets() {
        let qs = queries("[query.m]\nread = \":~\"\nrefr = \"~~\"\nmark = \"+^\"\n");
        let q = &qs.queries[0];
        let mut b = HitBuilder::new(RecordField::CORE, 4, false, 0);
        b.read_len = 50;
        let ring = ring_of(q.span, &[col(Seq::CLIP, Seq::G, 99, 1, -1), col(Seq::C, Seq::C, 100, 2, -1)]);
        b.push_hit(q, &ring);
        let names = field_names(&HitBuilder::new(RecordField::CORE, 4, false, 0));
        let arrays = b.finish();
        let at = |n: &str| arrays[names.iter().position(|x| x == n).unwrap()].clone();
        let off = at("read_5p");
        assert_eq!(off.as_any().downcast_ref::<Int64Array>().unwrap().value(0), 1, "the clip anchor");
        assert!(at("qual").is_null(0));
        let list = at("capture_read_5p");
        let list = list.as_any().downcast_ref::<ListArray>().unwrap().value(0);
        let list = list.as_any().downcast_ref::<Int64Array>().unwrap();
        assert_eq!((list.value(0), list.value(1)), (1, 2), "neither capture is null");
    }

    /// Every column the builder declares has a storage policy, and every policy
    /// names a column that exists. The writer errors on the second of those, so
    /// this is what keeps that error out of a real run.
    #[test]
    fn every_column_has_a_workable_encoding() {
        for flat in [true, false] {
            let b = HitBuilder::new(RecordField::ALL, 4, flat, 0);
            let fields = b.arrow_fields();
            let encs = b.column_encodings();
            assert_eq!(encs.len(), fields.len(), "flat={flat}: a policy per column");
            for ((name, enc), f) in encs.iter().zip(fields.iter()) {
                assert_eq!(name, f.name(), "flat={flat}: policies are in schema order");
                // A list column's policy is about the item inside it.
                let dt = match f.data_type() {
                    DataType::List(item) => item.data_type().clone(),
                    dt => dt.clone(),
                };
                assert!(
                    enc.supports(&dt),
                    "flat={flat}: `{name}` is {dt:?}, which cannot be stored as {}",
                    enc.name()
                );
            }
        }
    }

    /// A hitless row is null in every query column but still carries the record
    /// id, which is what makes the file an inventory rather than a hit list.
    #[test]
    fn a_hitless_row_is_null_where_a_query_would_be() {
        let mut b = HitBuilder::new(RecordField::CORE, 4, true, 0);
        b.rid = 12;
        b.push_hitless();
        assert_eq!(b.len(), 1);

        let arrays = b.finish();
        let names = field_names(&HitBuilder::new(RecordField::CORE, 4, true, 0));
        let at = |n: &str| arrays[names.iter().position(|x| x == n).unwrap()].clone();

        let rid = at("record_id");
        let rid = rid.as_any().downcast_ref::<UInt64Array>().unwrap();
        assert_eq!(rid.value(0), 12);
        for c in ["name", "read_base"] {
            assert!(at(c).is_null(0), "{c} should be null on a hitless row");
        }
    }

    /// `finish` resets the row count, so a builder is reusable across batches.
    #[test]
    fn finishing_empties_the_builder() {
        let mut b = HitBuilder::new(RecordField::CORE, 4, true, 0);
        b.push_hitless();
        assert_eq!(b.len(), 1);
        let _ = b.finish();
        assert_eq!(b.len(), 0);
    }

    /// The shard id is per row, not per file, so a consumer reading several
    /// shards off one stream can still tell which worker a row came from.
    #[test]
    fn every_row_carries_the_shard_id() {
        let qs = queries("[query.m]\nread = \"Y~\"\nrefr = \"CG\"\n");
        let q = &qs.queries[0];
        let mut b = HitBuilder::new(RecordField::CORE, 4, qs.flat_captures, 7);
        assert!(field_names(&b).contains(&"shard".to_string()));

        let ring = ring_of(
            q.span,
            &[
                col(Seq::C, Seq::C, 100, 10, 30),
                col(Seq::G, Seq::G, 101, 11, 31),
            ],
        );
        b.push_hit(q, &ring);
        b.push_hit(q, &ring);
        let arrays = b.finish();

        let shard = arrays[0]
            .as_any()
            .downcast_ref::<UInt32Array>()
            .expect("shard is the first column and is UInt32");
        assert_eq!(shard.len(), 2);
        assert_eq!(shard.value(0), 7);
        assert_eq!(shard.value(1), 7);
    }
}