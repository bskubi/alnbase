//! Extract the record columns once for each record, and append them to each row.
//!
//! A [`Column`] owns its field, its scratch cell and its Arrow builder
//! together, and the set of them is a `Vec`. The alternative -- parallel
//! structures with one `Option` per possible column, kept in lockstep by hand
//! -- fails by producing a silently all-null column rather than an error, so
//! the three parts are deliberately one thing.
//!
//! # Why extraction and appending are still separate
//!
//! The hits file is denormalized: one record produces one row per pattern that
//! fired on it, and every one of those rows repeats the record's columns.
//! Formatting the CIGAR or the SEQ once per row would redo that work a dozen
//! times per read. So [`RecordColumns::extract`] runs once per record and fills
//! each column's cell, and [`RecordColumns::append_row`] runs once per row and
//! copies cells into builders.
//!
//! # Why aux tags are read in one pass
//!
//! `Record::aux` is `bam_aux_get`, which walks the aux block from its first
//! tag every time it is called. A selection naming several tags therefore
//! walked the block once per tag, and on a CRAM of bisulfite reads -- a dozen
//! tags or more, the interesting ones often last -- that walk was the largest
//! single cost on the reader thread, which extracts the partition key for
//! every record before anything can run in parallel. So the aux columns are
//! filled together, by one walk that stops as soon as every wanted tag has
//! been seen. [`AuxWalk`] is that walk; it skips the tags nobody asked for by
//! their size alone, without parsing or validating them.

use std::fmt::Write as _;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, BooleanBuilder, Float64Builder, Int32Builder, Int64Builder, StringBuilder,
    UInt8Builder, UInt16Builder,
};
use arrow::datatypes::{DataType, Field};
use rust_htslib::bam::record::{Cigar, CigarStringView};
use rust_htslib::bam::{HeaderView, Record};

use crate::aux::{int_value, le, AuxWalk, RawAux};
use crate::encoding::ColumnEncoding;
use crate::partition::Fnv;
use crate::record_field::{AuxKind, RecordField};

/// One extracted value, or its absence.
///
/// `Str` carries no payload: the text lives in the column's own `String`, which
/// keeps its allocation across records.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Cell {
    Null,
    Bool(bool),
    U8(u8),
    U16(u16),
    I32(i32),
    I64(i64),
    F64(f64),
    Str,
}

impl Cell {
    /// A distinct byte per variant, mixed into the digest ahead of the value.
    ///
    /// Without it an absent aux tag and an empty one would hash alike, and so
    /// would a false and a zero — the values are distinguishable and their
    /// partitions should be too.
    fn tag(self) -> u8 {
        match self {
            Cell::Null => 0,
            Cell::Bool(_) => 1,
            Cell::U8(_) => 2,
            Cell::U16(_) => 3,
            Cell::I32(_) => 4,
            Cell::I64(_) => 5,
            Cell::F64(_) => 6,
            Cell::Str => 7,
        }
    }

    fn as_bool(self) -> Option<bool> {
        match self {
            Cell::Bool(v) => Some(v),
            _ => None,
        }
    }
    fn as_u8(self) -> Option<u8> {
        match self {
            Cell::U8(v) => Some(v),
            _ => None,
        }
    }
    fn as_u16(self) -> Option<u16> {
        match self {
            Cell::U16(v) => Some(v),
            _ => None,
        }
    }
    fn as_i32(self) -> Option<i32> {
        match self {
            Cell::I32(v) => Some(v),
            _ => None,
        }
    }
    fn as_i64(self) -> Option<i64> {
        match self {
            Cell::I64(v) => Some(v),
            _ => None,
        }
    }
    fn as_f64(self) -> Option<f64> {
        match self {
            Cell::F64(v) => Some(v),
            _ => None,
        }
    }
}

/// Storage, keyed by Arrow type rather than by field — a dozen fields share
/// each arm.
enum ColBuilder {
    Bool(BooleanBuilder),
    U8(UInt8Builder),
    U16(UInt16Builder),
    I32(Int32Builder),
    I64(Int64Builder),
    F64(Float64Builder),
    Str(StringBuilder),
}

impl ColBuilder {
    fn new(dt: &DataType, cap: usize) -> Self {
        match dt {
            DataType::Boolean => ColBuilder::Bool(BooleanBuilder::with_capacity(cap)),
            DataType::UInt8 => ColBuilder::U8(UInt8Builder::with_capacity(cap)),
            DataType::UInt16 => ColBuilder::U16(UInt16Builder::with_capacity(cap)),
            DataType::Int32 => ColBuilder::I32(Int32Builder::with_capacity(cap)),
            DataType::Int64 => ColBuilder::I64(Int64Builder::with_capacity(cap)),
            DataType::Float64 => ColBuilder::F64(Float64Builder::with_capacity(cap)),
            DataType::Utf8 => ColBuilder::Str(StringBuilder::with_capacity(cap, cap * 16)),
            other => unreachable!("no builder for {other:?}; RecordField::data_type is the only source of these"),
        }
    }

    fn finish(&mut self) -> ArrayRef {
        match self {
            ColBuilder::Bool(b) => Arc::new(b.finish()),
            ColBuilder::U8(b) => Arc::new(b.finish()),
            ColBuilder::U16(b) => Arc::new(b.finish()),
            ColBuilder::I32(b) => Arc::new(b.finish()),
            ColBuilder::I64(b) => Arc::new(b.finish()),
            ColBuilder::F64(b) => Arc::new(b.finish()),
            ColBuilder::Str(b) => Arc::new(b.finish()),
        }
    }
}

/// What the whole record needs computed once, shared by every column.
struct Ctx {
    cigar: Option<CigarStringView>,
    /// What the run's strand rule made of this record. Passed in rather than
    /// derived here so that the strand columns and the coordinates the walk
    /// emitted are the same decision written down twice, not two decisions
    /// that happen to agree today.
    /// `None` when the run's rule declined to call this record, which makes
    /// every strand column null: the input does not say, and a column that
    /// said anyway would read like a measurement.
    strand: Option<crate::strand_rule::StrandCall>,
}

impl Ctx {
    /// Only called by columns whose `needs_cigar` is true, which is what put
    /// the `Some` here in the first place.
    fn cigar(&self) -> &CigarStringView {
        self.cigar
            .as_ref()
            .expect("a cigar column was selected but no cigar was parsed")
    }
}

struct Column {
    field: RecordField,
    name: String,
    dtype: DataType,
    cell: Cell,
    /// Scratch for Utf8 columns; keeps its allocation across records.
    text: String,
    builder: ColBuilder,
}

impl Column {
    fn new(field: RecordField, cap: usize) -> Self {
        let dtype = field.data_type();
        Self {
            field,
            name: field.name(),
            cell: Cell::Null,
            text: String::new(),
            builder: ColBuilder::new(&dtype, cap),
            dtype,
        }
    }

    #[inline]
    fn set_text(&mut self, bytes: &[u8]) {
        self.text.clear();
        match std::str::from_utf8(bytes) {
            Ok(s) => self.text.push_str(s),
            Err(_) => self.text.push_str(&String::from_utf8_lossy(bytes)),
        }
        self.cell = Cell::Str;
    }

    /// Reference name for a tid, or null when the tid is unset.
    #[inline]
    fn set_ref_name(&mut self, header: &HeaderView, tid: i32) {
        if tid < 0 {
            self.cell = Cell::Null;
        } else {
            // Borrows the header, not self, so no copy is needed.
            self.set_text(header.tid2name(tid as u32));
        }
    }

    fn extract(&mut self, r: &Record, h: &HeaderView, ctx: &Ctx) {
        use RecordField as F;
        self.cell = match self.field {
            F::Qname => {
                self.set_text(r.qname());
                return;
            }
            F::Flags => Cell::U16(r.flags()),
            F::Tid => opt_tid(r.tid()),
            F::RefName => {
                self.set_ref_name(h, r.tid());
                return;
            }
            F::Pos => Cell::I64(r.pos()),
            F::EndPos => Cell::I64(ctx.cigar().end_pos()),
            F::UnclippedStart => Cell::I64(r.pos() - clip_len(ctx.cigar(), Side::Lead)),
            F::UnclippedEnd => Cell::I64(ctx.cigar().end_pos() + clip_len(ctx.cigar(), Side::Trail)),
            F::MapQ => Cell::U8(r.mapq()),
            F::Cigar => {
                self.text.clear();
                let _ = write!(self.text, "{}", ctx.cigar());
                Cell::Str
            }
            F::ReadLen => Cell::I64(r.seq_len() as i64),
            F::HardClip5p => Cell::I64(crate::alignment::hard_clips_as_sequenced(r).0),
            F::HardClip3p => Cell::I64(crate::alignment::hard_clips_as_sequenced(r).1),
            F::SoftClip5p => Cell::I64(crate::alignment::soft_clips_as_sequenced(r).0),
            F::SoftClip3p => Cell::I64(crate::alignment::soft_clips_as_sequenced(r).1),

            F::IsPaired => Cell::Bool(r.is_paired()),
            F::IsProperPair => Cell::Bool(r.is_proper_pair()),
            F::IsUnmapped => Cell::Bool(r.is_unmapped()),
            F::IsMateUnmapped => Cell::Bool(r.is_mate_unmapped()),
            F::IsReverse => Cell::Bool(r.is_reverse()),
            F::IsMateReverse => Cell::Bool(r.is_mate_reverse()),
            F::IsFirstInTemplate => Cell::Bool(r.is_first_in_template()),
            F::IsLastInTemplate => Cell::Bool(r.is_last_in_template()),
            F::IsSecondary => Cell::Bool(r.is_secondary()),
            F::IsSupplementary => Cell::Bool(r.is_supplementary()),
            F::IsDuplicate => Cell::Bool(r.is_duplicate()),
            F::IsQcFail => Cell::Bool(r.is_quality_check_failed()),
            // Not flag bits — the three of them are the strand call, written
            // out. `walk_alignment` was handed the same call, so these columns
            // and the emitted coordinates cannot disagree.
            F::RefrStrand => match ctx.strand {
                Some(s) => {
                    self.text.clear();
                    self.text.push_str(s.conversion.name());
                    Cell::Str
                }
                None => Cell::Null,
            },

            // Null when the rule named only the conversion strand: the input
            // does not distinguish OT from CTOT, and picking one would be a
            // guess that reads like a measurement.
            F::ConvStrand => match ctx.strand.and_then(|s| s.origin) {
                Some(origin) => {
                    self.text.clear();
                    self.text.push_str(origin.name());
                    Cell::Str
                }
                None => Cell::Null,
            },

            F::ReadReverse => match ctx.strand {
                Some(s) => Cell::Bool(s.sequenced.is_reverse()),
                None => Cell::Null,
            },

            F::InsertSize => Cell::I64(r.insert_size()),
            F::MateTid => opt_tid(r.mtid()),
            F::MateRefName => {
                self.set_ref_name(h, r.mtid());
                return;
            }
            F::MatePos => {
                if r.mtid() < 0 {
                    Cell::Null
                } else {
                    Cell::I64(r.mpos())
                }
            }

            F::SeqAscii => {
                // Indexing decodes one 4-bit nibble to its ASCII base, so this
                // reads straight out of the record with no intermediate Vec.
                self.text.clear();
                let s = r.seq();
                self.text.reserve(s.len());
                for i in 0..s.len() {
                    self.text.push(s[i] as char);
                }
                Cell::Str
            }
            F::QualPhred => {
                // htslib hands back raw phred; SAM writes phred+33. A record
                // with '*' for QUAL has 0xff in every position: no quality.
                self.text.clear();
                let q = r.qual();
                if q.first() == Some(&0xff) {
                    Cell::Null
                } else {
                    self.text.reserve(q.len());
                    for &p in q {
                        self.text.push(p.saturating_add(33) as char);
                    }
                    Cell::Str
                }
            }

            // Filled afterwards, by `RecordColumns::extract_aux`, in one pass
            // over the aux block shared by every aux column. Null is what a
            // column keeps when its tag is not found.
            F::Aux(..) => Cell::Null,
        };
    }

    /// Fill this column from its tag's raw value, as found by [`AuxWalk`].
    ///
    /// A type that does not suit the column's kind means null, as does a
    /// string that is not UTF-8 -- the same answers `Record::aux` gives, which
    /// `tests::raw_aux_reads_agree_with_rust_htslib` holds this to.
    /// Returns whether the value was dropped because its text is not UTF-8,
    /// which is a property of the data rather than of the selection and so is
    /// worth counting; every other mismatch is just a null.
    fn set_aux(&mut self, kind: AuxKind, field: &RawAux<'_>) -> bool {
        let v = field.value;
        let mut bad_text = false;
        self.cell = match (kind, field.ty) {
            (AuxKind::Int, t) => match int_value(t, v) {
                Some(x) => Cell::I64(x),
                None => Cell::Null,
            },
            (AuxKind::Float, b'f') => Cell::F64(f32::from_le_bytes(le(v)) as f64),
            (AuxKind::Float, b'd') => Cell::F64(f64::from_le_bytes(le(v))),
            // An integer tag read as a float is a widening, not a mismatch.
            (AuxKind::Float, t) => match int_value(t, v) {
                Some(x) => Cell::F64(x as f64),
                None => Cell::Null,
            },
            // `H` is hex digits, which htslib hands back as a string too.
            (AuxKind::Str, b'Z' | b'H') => match std::str::from_utf8(v) {
                Ok(s) => {
                    self.text.clear();
                    self.text.push_str(s);
                    Cell::Str
                }
                Err(_) => {
                    bad_text = true;
                    Cell::Null
                }
            },
            (AuxKind::Str, b'A') => {
                self.text.clear();
                self.text.push(v[0] as char);
                Cell::Str
            }
            // `B` arrays, comma-joined. `value` is the subtype byte, the
            // count, then the elements; `AuxWalk` has already checked that
            // the elements fit.
            (AuxKind::Str, b'B') => {
                let data = &v[5..];
                match v[0] {
                    b'c' => self.join(data.iter().map(|&b| b as i8)),
                    b'C' => self.join(data.iter()),
                    b's' => self.join(data.chunks_exact(2).map(|c| i16::from_le_bytes(le(c)))),
                    b'S' => self.join(data.chunks_exact(2).map(|c| u16::from_le_bytes(le(c)))),
                    b'i' => self.join(data.chunks_exact(4).map(|c| i32::from_le_bytes(le(c)))),
                    b'I' => self.join(data.chunks_exact(4).map(|c| u32::from_le_bytes(le(c)))),
                    b'f' => self.join(data.chunks_exact(4).map(|c| f32::from_le_bytes(le(c)))),
                    _ => Cell::Null,
                }
            }
            _ => Cell::Null,
        };
        bad_text
    }

    fn join<T: std::fmt::Display>(&mut self, it: impl Iterator<Item = T>) -> Cell {
        self.text.clear();
        for (i, v) in it.enumerate() {
            if i > 0 {
                self.text.push(',');
            }
            let _ = write!(self.text, "{v}");
        }
        Cell::Str
    }

    #[inline]
    fn append_row(&mut self) {
        // Disjoint fields: the builder borrow and the cell/text reads do not
        // overlap, so this needs no dance to satisfy borrowck.
        match &mut self.builder {
            ColBuilder::Bool(b) => b.append_option(self.cell.as_bool()),
            ColBuilder::U8(b) => b.append_option(self.cell.as_u8()),
            ColBuilder::U16(b) => b.append_option(self.cell.as_u16()),
            ColBuilder::I32(b) => b.append_option(self.cell.as_i32()),
            ColBuilder::I64(b) => b.append_option(self.cell.as_i64()),
            ColBuilder::F64(b) => b.append_option(self.cell.as_f64()),
            ColBuilder::Str(b) => match self.cell {
                Cell::Str => b.append_value(&self.text),
                _ => b.append_null(),
            },
        }
    }

    fn arrow_field(&self) -> Field {
        // Every record column is nullable; see RecordField::data_type.
        Field::new(&self.name, self.dtype.clone(), true)
    }

    /// Fold the extracted value into `h`.
    ///
    /// Text is length-prefixed, so the boundary between two fields of a
    /// composite key is itself hashed: `("ab", "c")` and `("a", "bc")` are the
    /// same bytes in the same order and must not collide.
    #[inline]
    fn digest(&self, h: &mut Fnv) {
        h.byte(self.cell.tag());
        match self.cell {
            Cell::Null => {}
            Cell::Bool(v) => h.byte(v as u8),
            Cell::U8(v) => h.u64(v as u64),
            Cell::U16(v) => h.u64(v as u64),
            Cell::I32(v) => h.u64(v as i64 as u64),
            Cell::I64(v) => h.u64(v as u64),
            // By bits, so the hash is exact rather than rounded. NaN keys are
            // not a concern: no aux tag this reads produces one.
            Cell::F64(v) => h.u64(v.to_bits()),
            Cell::Str => {
                h.u64(self.text.len() as u64);
                h.bytes(self.text.as_bytes());
            }
        }
    }
}

/// The columns reading one aux tag.
///
/// Usually one column, but kept a list so that nothing here depends on
/// [`RecordField::resolve`] having dropped same-tag duplicates.
struct AuxWant {
    tag: [u8; 2],
    /// Indices into `RecordColumns::cols`, with the kind each reads as.
    cols: Vec<(usize, AuxKind)>,
    /// Whether this record's walk has already filled these columns. The
    /// first occurrence of a tag wins, as it does for `bam_aux_get`.
    seen: bool,
}

#[inline]
fn opt_tid(tid: i32) -> Cell {
    if tid < 0 {
        Cell::Null
    } else {
        Cell::I32(tid)
    }
}

enum Side {
    Lead,
    Trail,
}

/// Soft + hard clip at one end of the alignment, in reference order.
fn clip_len(cigar: &CigarStringView, side: Side) -> i64 {
    let mut n = 0i64;
    let mut seen_aligned = false;
    for op in cigar.iter() {
        match op {
            Cigar::SoftClip(l) | Cigar::HardClip(l) => {
                let at_lead = !seen_aligned;
                let want_lead = matches!(side, Side::Lead);
                if at_lead == want_lead {
                    n += *l as i64;
                }
            }
            _ => seen_aligned = true,
        }
    }
    n
}

// ------------------------------------------------------------------ public --

pub struct RecordColumns {
    cols: Vec<Column>,
    needs_cigar: bool,
    /// One entry per distinct aux tag among `cols`, in column order. Empty
    /// when no aux column is selected, which skips the walk entirely.
    aux: Vec<AuxWant>,
    rows: usize,
    /// What the last extraction found wrong with the record's aux block.
    trouble: AuxTrouble,
}

/// What was wrong with one record's aux block, if anything.
///
/// Both cases leave a column null for a reason the selection cannot explain,
/// so a run either stops on them or counts them; see
/// [`crate::batch::OnBadData`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AuxTrouble {
    /// The walk stopped at a field it could not size, so every tag after it
    /// read as absent.
    pub damaged: bool,
    /// A `Z` or `H` value was not valid UTF-8 and became null rather than
    /// being lossily converted into a different, plausible value.
    pub non_utf8: bool,
}

impl AuxTrouble {
    pub fn any(self) -> bool {
        self.damaged || self.non_utf8
    }

    /// What to tell the user, as a sentence fragment.
    pub fn describe(self) -> &'static str {
        match (self.damaged, self.non_utf8) {
            (true, true) => "a damaged aux field, and text that is not UTF-8",
            (true, false) => "a damaged aux field, so tags after it read as absent",
            _ => "text in an aux tag that is not UTF-8",
        }
    }
}

impl RecordColumns {
    /// `fields` is already resolved: groups expanded, duplicates dropped,
    /// column order fixed. See [`RecordField::resolve`].
    pub fn new(fields: &[RecordField], cap: usize) -> Self {
        let mut aux: Vec<AuxWant> = Vec::new();
        for (i, f) in fields.iter().enumerate() {
            if let RecordField::Aux(tag, kind) = *f {
                match aux.iter_mut().find(|w| w.tag == tag) {
                    Some(w) => w.cols.push((i, kind)),
                    None => aux.push(AuxWant { tag, cols: vec![(i, kind)], seen: false }),
                }
            }
        }
        Self {
            needs_cigar: fields.iter().any(|f| f.needs_cigar()),
            cols: fields.iter().map(|f| Column::new(*f, cap)).collect(),
            aux,
            rows: 0,
            trouble: AuxTrouble::default(),
        }
    }

    /// Read every column off `record`. Call once per record.
    ///
    /// Infallible by design: a missing aux tag, an unplaced mate and invalid
    /// UTF-8 in a name are all nulls, not errors.
    pub fn extract(
        &mut self,
        record: &Record,
        header: &HeaderView,
        strand: Option<crate::strand_rule::StrandCall>,
    ) {
        // `Record::cigar()` allocates, so it is built once here rather than
        // once per column that wants it — and not at all when none does.
        let ctx = Ctx {
            cigar: self.needs_cigar.then(|| record.cigar()),
            strand,
        };
        for c in &mut self.cols {
            c.extract(record, header, &ctx);
        }
        self.trouble = AuxTrouble::default();
        if !self.aux.is_empty() {
            self.extract_aux(record);
        }
    }

    /// Fill every aux column from one walk over the record's aux block.
    ///
    /// `Column::extract` has already set each of them to null, so a tag the
    /// walk never reaches stays null. The walk stops once every wanted tag has
    /// been seen, which for a key of tags near the front of the block is a
    /// few fields in.
    fn extract_aux(&mut self, record: &Record) {
        let Self { cols, aux, trouble, .. } = self;
        for w in aux.iter_mut() {
            w.seen = false;
        }
        let mut unseen = aux.len();
        let mut walk = AuxWalk::new(record);
        for field in walk.by_ref() {
            let Some(w) = aux.iter_mut().find(|w| w.tag == field.tag) else {
                continue;
            };
            if w.seen {
                continue;
            }
            w.seen = true;
            for &(i, kind) in &w.cols {
                trouble.non_utf8 |= cols[i].set_aux(kind, &field);
            }
            unseen -= 1;
            if unseen == 0 {
                // Every wanted tag is filled, but the rest of the block still
                // decides whether it is damaged -- and on a damaged block a
                // tag that was *not* found may have been unreachable rather
                // than absent, which is the thing worth reporting.
                for _ in walk.by_ref() {}
                break;
            }
        }
        trouble.damaged = walk.stopped_early();
    }

    /// What the last [`extract`](Self::extract) found wrong with the record's
    /// aux block.
    pub fn aux_trouble(&self) -> AuxTrouble {
        self.trouble
    }

    /// Hash of the current extraction, over every column in selection order.
    ///
    /// Call after [`extract`](Self::extract) and never mind the builders: this
    /// is how a partition key is computed, and computing it here rather than
    /// from a second extractor is what makes `--partition-by ref_name` group
    /// by exactly what the `ref_name` column would say.
    ///
    /// Order matters, so `--partition-by qname,mapq` and `--partition-by
    /// mapq,qname` are different keys over the same values. They partition the
    /// same records the same way regardless, since both are deterministic.
    pub fn digest(&self) -> u64 {
        let mut h = Fnv::new();
        for c in &self.cols {
            c.digest(&mut h);
        }
        h.finish()
    }

    /// Append one row from the current extraction. Call once per emitted row.
    pub fn append_row(&mut self) {
        for c in &mut self.cols {
            c.append_row();
        }
        self.rows += 1;
    }

    /// Drain every column, in schema order.
    pub fn finish(&mut self) -> Vec<ArrayRef> {
        self.rows = 0;
        self.cols.iter_mut().map(|c| c.finish_col()).collect()
    }



    /// Number of columns, i.e. arrays `finish` will return.
    pub fn width(&self) -> usize {
        self.cols.len()
    }

    /// Schema order. Same order as [`finish`](Self::finish), by construction —
    /// both walk `self.cols`.
    pub fn arrow_fields(&self) -> Vec<Field> {
        self.cols.iter().map(|c| c.arrow_field()).collect()
    }

    /// Each column's storage policy, paired with the name it goes into the
    /// schema under.
    ///
    /// By name rather than by index because the writer merges these with the
    /// hit columns' policies and with the user's `--column-encoding`
    /// overrides, and a name is the only thing all three share. Aux tags name
    /// themselves, so `-F ML:B` is overridden as `ML=plain`.
    pub fn column_encodings(&self) -> Vec<(String, ColumnEncoding)> {
        self.cols.iter().map(|c| (c.name.clone(), c.field.encoding())).collect()
    }
}

impl Column {
    fn finish_col(&mut self) -> ArrayRef {
        self.builder.finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn schema_matches_the_selection() {
        let fields = RecordField::CORE;
        let cols = RecordColumns::new(fields, 8);
        let arrow = cols.arrow_fields();
        assert_eq!(arrow.len(), fields.len());
        for (a, f) in arrow.iter().zip(fields) {
            assert_eq!(a.name(), &f.name());
            assert_eq!(a.data_type(), &f.data_type());
            assert!(a.is_nullable());
        }
    }

    #[test]
    fn every_field_builds_a_column() {
        // Catches a DataType with no ColBuilder arm.
        let mut all = RecordField::ALL.to_vec();
        all.push(RecordField::Aux(*b"NM", AuxKind::Int));
        all.push(RecordField::Aux(*b"AS", AuxKind::Float));
        all.push(RecordField::Aux(*b"RG", AuxKind::Str));
        let cols = RecordColumns::new(&all, 4);
        assert_eq!(cols.width(), all.len());
    }

    #[test]
    fn finish_returns_one_array_per_field() {
        let mut cols = RecordColumns::new(RecordField::ALL, 4);
        let arrays = cols.finish();
        assert_eq!(arrays.len(), cols.arrow_fields().len());
    }

    // ---- extraction ------------------------------------------------------

    use crate::test_support::{built, header, LAST_IN_TEMPLATE, REVERSE};
    use arrow::array::{
        Array, BooleanArray, Int32Array, Int64Array, StringArray, UInt8Array,
    };
    use rust_htslib::bam::record::{Aux, Cigar};

    const CONTIG: &str = "chr1";

    /// The call the run makes for a record, which is what the strand columns
    /// are written from. Directional is the only rule these tests need: the
    /// point here is that the columns say what the call said, not how a
    /// particular rule reaches it.
    fn call(rec: &Record) -> crate::strand_rule::StrandCall {
        crate::tags::Library::Directional.call(rec)
    }

    /// Extract one record and read every column back as text, so a test can
    /// assert on values without knowing each column's Arrow type.
    fn extract_row(fields: &[RecordField], rec: &Record) -> Vec<Option<String>> {
        extract_row_as(fields, rec, call(rec))
    }

    /// The same, for a run whose rule reached a call the FLAG does not imply.
    fn extract_row_as(
        fields: &[RecordField],
        rec: &Record,
        strand: crate::strand_rule::StrandCall,
    ) -> Vec<Option<String>> {
        let view = header(CONTIG, 1000);
        let mut cols = RecordColumns::new(fields, 4);
        cols.extract(rec, &view, Some(strand));
        cols.append_row();
        let arrays = cols.finish();
        arrays
            .iter()
            .map(|a| {
                if a.is_null(0) {
                    return None;
                }
                let d = |s: String| Some(s);
                if let Some(x) = a.as_any().downcast_ref::<StringArray>() {
                    d(x.value(0).to_string())
                } else if let Some(x) = a.as_any().downcast_ref::<Int64Array>() {
                    d(x.value(0).to_string())
                } else if let Some(x) = a.as_any().downcast_ref::<Int32Array>() {
                    d(x.value(0).to_string())
                } else if let Some(x) = a.as_any().downcast_ref::<UInt8Array>() {
                    d(x.value(0).to_string())
                } else if let Some(x) = a.as_any().downcast_ref::<BooleanArray>() {
                    d(x.value(0).to_string())
                } else {
                    Some("<unhandled type>".to_string())
                }
            })
            .collect()
    }

    fn one(field: RecordField, rec: &Record) -> Option<String> {
        extract_row(&[field], rec).remove(0)
    }

    /// Hard clips are attributed to the read's sequenced ends: a reverse
    /// record's 5' end is its CIGAR's last operation. A '*' QUAL is null.
    #[test]
    fn hard_clips_follow_the_sequenced_ends_and_missing_quality_is_null() {
        let cig = [Cigar::HardClip(3), Cigar::Match(4), Cigar::HardClip(2)];
        let fwd = built(b"r", b"ACGT", &cig, 0, 0);
        let rev = built(b"r", b"ACGT", &cig, 0, crate::test_support::REVERSE);
        assert_eq!(one(RecordField::HardClip5p, &fwd).as_deref(), Some("3"));
        assert_eq!(one(RecordField::HardClip3p, &fwd).as_deref(), Some("2"));
        assert_eq!(one(RecordField::HardClip5p, &rev).as_deref(), Some("2"));
        assert_eq!(one(RecordField::HardClip3p, &rev).as_deref(), Some("3"));

        // Soft clips inside hard clips, attributed the same way.
        let cig = [Cigar::HardClip(3), Cigar::SoftClip(1), Cigar::Match(2), Cigar::SoftClip(4)];
        let fwd = built(b"r", b"ACGTACG", &cig, 0, 0);
        let rev = built(b"r", b"ACGTACG", &cig, 0, crate::test_support::REVERSE);
        assert_eq!(one(RecordField::SoftClip5p, &fwd).as_deref(), Some("1"));
        assert_eq!(one(RecordField::SoftClip3p, &fwd).as_deref(), Some("4"));
        assert_eq!(one(RecordField::SoftClip5p, &rev).as_deref(), Some("4"));
        assert_eq!(one(RecordField::SoftClip3p, &rev).as_deref(), Some("1"));
        let none = built(b"r", b"ACGT", &[Cigar::HardClip(2), Cigar::Match(4)], 0, 0);
        assert_eq!(one(RecordField::SoftClip5p, &none).as_deref(), Some("0"));
        assert_eq!(one(RecordField::SoftClip3p, &none).as_deref(), Some("0"));

        let mut no_qual = built(b"r", b"ACGT", &[Cigar::Match(4)], 0, 0);
        no_qual.set(b"r", Some(&rust_htslib::bam::record::CigarString(vec![Cigar::Match(4)])), b"ACGT", &[0xff; 4]);
        assert_eq!(one(RecordField::QualPhred, &no_qual), None);
    }

    #[test]
    fn the_basic_columns_read_off_the_record() {
        let cig = [Cigar::SoftClip(2), Cigar::Match(6), Cigar::SoftClip(2)];
        let rec = built(b"read42", b"ACGTACGTAC", &cig, 100, 0);

        assert_eq!(one(RecordField::Qname, &rec).as_deref(), Some("read42"));
        assert_eq!(one(RecordField::RefName, &rec).as_deref(), Some(CONTIG));
        assert_eq!(one(RecordField::Pos, &rec).as_deref(), Some("100"));
        assert_eq!(one(RecordField::EndPos, &rec).as_deref(), Some("106"));
        assert_eq!(one(RecordField::MapQ, &rec).as_deref(), Some("60"));
        assert_eq!(one(RecordField::Cigar, &rec).as_deref(), Some("2S6M2S"));
        assert_eq!(one(RecordField::ReadLen, &rec).as_deref(), Some("10"));
        assert_eq!(one(RecordField::SeqAscii, &rec).as_deref(), Some("ACGTACGTAC"));
    }

    /// Unclipped coordinates extend the alignment back over the clips, which is
    /// what `clip_len` is for.
    #[test]
    fn unclipped_coordinates_span_the_clips() {
        let cig = [Cigar::HardClip(3), Cigar::SoftClip(2), Cigar::Match(6), Cigar::SoftClip(4)];
        let rec = built(b"r", b"ACGTACGTACGT", &cig, 100, 0);
        assert_eq!(one(RecordField::Pos, &rec).as_deref(), Some("100"));
        assert_eq!(
            one(RecordField::UnclippedStart, &rec).as_deref(),
            Some("95"),
            "3H + 2S before the first aligned base"
        );
        assert_eq!(one(RecordField::EndPos, &rec).as_deref(), Some("106"));
        assert_eq!(one(RecordField::UnclippedEnd, &rec).as_deref(), Some("110"));
    }

    /// `clip_len` splits on whether an aligned op has been seen, so clips at
    /// each end are attributed to that end and nothing else counts.
    #[test]
    fn clip_len_attributes_each_clip_to_its_own_end() {
        let rec = built(
            b"r",
            b"ACGTACGTACGT",
            &[Cigar::HardClip(3), Cigar::SoftClip(2), Cigar::Match(6), Cigar::SoftClip(4)],
            0,
            0,
        );
        let c = rec.cigar();
        assert_eq!(clip_len(&c, Side::Lead), 5, "hard and soft both count");
        assert_eq!(clip_len(&c, Side::Trail), 4);

        // No clips at all.
        let plain = built(b"r", b"ACGTAC", &[Cigar::Match(6)], 0, 0);
        assert_eq!(clip_len(&plain.cigar(), Side::Lead), 0);
        assert_eq!(clip_len(&plain.cigar(), Side::Trail), 0);
    }

    /// `strand` is the fragment's strand, not FLAG 0x10, and uses the same
    /// derivation as the walk so the column and the coordinates agree.
    #[test]
    fn strand_is_the_fragment_strand_not_the_flag_bit() {
        let m = [Cigar::Match(4)];
        let plain = built(b"r", b"ACGT", &m, 0, 0);
        assert_eq!(one(RecordField::RefrStrand, &plain).as_deref(), Some("+"));
        assert_eq!(one(RecordField::IsReverse, &plain).as_deref(), Some("false"));

        // R1 reverse: the fragment reads bottom strand.
        let rev = built(b"r", b"ACGT", &m, 0, REVERSE);
        assert_eq!(one(RecordField::RefrStrand, &rev).as_deref(), Some("-"));
        assert_eq!(one(RecordField::IsReverse, &rev).as_deref(), Some("true"));

        // R2 reverse: the same fragment strand as an unreversed R1, which is
        // exactly where `strand` and `is_reverse` part company.
        let r2rev = built(b"r", b"ACGT", &m, 0, REVERSE | LAST_IN_TEMPLATE);
        assert_eq!(one(RecordField::RefrStrand, &r2rev).as_deref(), Some("+"));
        assert_eq!(one(RecordField::IsReverse, &r2rev).as_deref(), Some("true"));
    }

    /// `conv_strand` names all four strands, and splits each `strand` value by
    /// mate: OT and CTOT are `+`, OB and CTOB are `-`.
    #[test]
    fn conv_strand_is_the_directional_conversion_strand() {
        let m = [Cigar::Match(4)];
        for (flags, conv, strand) in [
            (0, "OT", "+"),
            (REVERSE, "OB", "-"),
            (REVERSE | LAST_IN_TEMPLATE, "CTOT", "+"),
            (LAST_IN_TEMPLATE, "CTOB", "-"),
        ] {
            let rec = built(b"r", b"ACGT", &m, 0, flags);
            assert_eq!(one(RecordField::ConvStrand, &rec).as_deref(), Some(conv), "flags {flags:#x}");
            assert_eq!(one(RecordField::RefrStrand, &rec).as_deref(), Some(strand), "flags {flags:#x}");
        }
    }

    /// `read_reverse` is the direction the read was sequenced in, which for a
    /// directional library is FLAG 0x10 -- and the call is what says so, not
    /// the bit.
    #[test]
    fn read_reverse_follows_the_call() {
        let m = [Cigar::Match(4)];
        for (flags, reverse) in [
            (0, "false"),
            (REVERSE, "true"),
            (REVERSE | LAST_IN_TEMPLATE, "true"),
            (LAST_IN_TEMPLATE, "false"),
        ] {
            let rec = built(b"r", b"ACGT", &m, 0, flags);
            assert_eq!(
                one(RecordField::ReadReverse, &rec).as_deref(),
                Some(reverse),
                "flags {flags:#x}"
            );
        }
    }

    /// The three strand columns say what the rule concluded, even when the
    /// FLAG would have said something else.
    ///
    /// This is the shape a two-way aligner produces: the reads are the two
    /// orientations of the top strand's converted copy, so the conversion
    /// strand is `+` for both and there is no OT/CTOT to report. A run against
    /// such input used to get `-` here from a flag bit that means something
    /// else, and an invented `OB`.
    #[test]
    fn the_strand_columns_are_the_call_not_the_flag() {
        use crate::strand_rule::{ConvStrand, SeqDir, StrandCall};
        let rec = built(b"r", b"ACGT", &[Cigar::Match(4)], 0, REVERSE);
        let two_way =
            StrandCall { conversion: ConvStrand::Plus, sequenced: SeqDir::Reverse, origin: None };
        let fields = [RecordField::RefrStrand, RecordField::ConvStrand, RecordField::ReadReverse];
        let row = extract_row_as(&fields, &rec, two_way);
        assert_eq!(row[0].as_deref(), Some("+"), "0x10 is set but the rule says plus");
        assert_eq!(row[1], None, "two strands named, so there is no strand of origin to give");
        assert_eq!(row[2].as_deref(), Some("true"));
    }

    /// An unmapped mate has no reference, and a null says so rather than -1.
    #[test]
    fn absent_mate_coordinates_are_null() {
        let rec = built(b"r", b"ACGT", &[Cigar::Match(4)], 0, 0);
        assert_eq!(one(RecordField::MateTid, &rec), None, "mtid is -1");
        assert_eq!(one(RecordField::MateRefName, &rec), None);
    }

    /// A tag that is absent, and one whose type does not match, both read as
    /// null -- an error would kill the worker on the first read without it.
    #[test]
    fn aux_tags_read_when_present_and_null_when_not() {
        let mut rec = built(b"r", b"ACGT", &[Cigar::Match(4)], 0, 0);
        rec.push_aux(b"NM", Aux::I32(4)).unwrap();
        rec.push_aux(b"RG", Aux::String("grp1")).unwrap();

        assert_eq!(
            one(RecordField::Aux(*b"NM", AuxKind::Int), &rec).as_deref(),
            Some("4")
        );
        assert_eq!(
            one(RecordField::Aux(*b"RG", AuxKind::Str), &rec).as_deref(),
            Some("grp1")
        );
        // Never pushed.
        assert_eq!(one(RecordField::Aux(*b"ZZ", AuxKind::Int), &rec), None);
        // Present, but not a string.
        assert_eq!(one(RecordField::Aux(*b"NM", AuxKind::Str), &rec), None);
    }

    /// The reason extraction and appending are separate: one extract feeds many
    /// rows, and every row must carry the same values.
    #[test]
    fn one_extraction_feeds_many_identical_rows() {
        let view = header(CONTIG, 1000);
        let fields = [RecordField::Qname, RecordField::Pos, RecordField::Cigar];
        let mut cols = RecordColumns::new(&fields, 8);

        let rec = built(b"readA", b"ACGTAC", &[Cigar::Match(6)], 7, 0);
        cols.extract(&rec, &view, Some(call(&rec)));
        for _ in 0..3 {
            cols.append_row();
        }
        assert_eq!(cols.arrow_fields().len(), cols.width(), "one field per column");

        // A second record appends one more row with its own values.
        let rec2 = built(b"readB", b"ACGTAC", &[Cigar::Match(6)], 11, 0);
        cols.extract(&rec2, &view, Some(call(&rec2)));
        cols.append_row();

        let arrays = cols.finish();
        let qname = arrays[0].as_any().downcast_ref::<StringArray>().unwrap();
        let pos = arrays[1].as_any().downcast_ref::<Int64Array>().unwrap();
        assert_eq!(
            (0..4).map(|i| qname.value(i)).collect::<Vec<_>>(),
            vec!["readA", "readA", "readA", "readB"]
        );
        assert_eq!((0..4).map(|i| pos.value(i)).collect::<Vec<_>>(), vec![7, 7, 7, 11]);
        assert_eq!(cols.finish().iter().map(|a| a.len()).max(), Some(0), "finish empties the builders");
    }

    // ---- aux: the single walk against rust-htslib's own reading ----------

    use rust_htslib::bam::record::AuxArray;
    use rust_htslib::htslib;

    impl Column {
    /// The per-tag `Record::aux` reading this module used before `AuxWalk`,
    /// kept as the reference the walk is checked against.
    ///
    /// Absent, and type-mismatched, both mean null.
    ///
    /// `Record::aux` returns `Err` when the tag is not on the record,
    /// which is the common case — propagating that would kill the worker on the
    /// first read without it.
    fn extract_aux_via_htslib(&mut self, r: &Record, tag: &[u8; 2], kind: AuxKind) {
        let Ok(v) = r.aux(tag) else {
            self.cell = Cell::Null;
            return;
        };
        self.cell = match (kind, v) {
            (AuxKind::Int, Aux::I8(x)) => Cell::I64(x as i64),
            (AuxKind::Int, Aux::U8(x)) => Cell::I64(x as i64),
            (AuxKind::Int, Aux::I16(x)) => Cell::I64(x as i64),
            (AuxKind::Int, Aux::U16(x)) => Cell::I64(x as i64),
            (AuxKind::Int, Aux::I32(x)) => Cell::I64(x as i64),
            (AuxKind::Int, Aux::U32(x)) => Cell::I64(x as i64),
            (AuxKind::Float, Aux::Float(x)) => Cell::F64(x as f64),
            (AuxKind::Float, Aux::Double(x)) => Cell::F64(x),
            // An integer tag read as a float is a widening, not a mismatch.
            (AuxKind::Float, Aux::I8(x)) => Cell::F64(x as f64),
            (AuxKind::Float, Aux::U8(x)) => Cell::F64(x as f64),
            (AuxKind::Float, Aux::I16(x)) => Cell::F64(x as f64),
            (AuxKind::Float, Aux::U16(x)) => Cell::F64(x as f64),
            (AuxKind::Float, Aux::I32(x)) => Cell::F64(x as f64),
            (AuxKind::Float, Aux::U32(x)) => Cell::F64(x as f64),
            (AuxKind::Str, Aux::String(s)) => {
                self.text.clear();
                self.text.push_str(s);
                Cell::Str
            }
            (AuxKind::Str, Aux::Char(c)) => {
                self.text.clear();
                self.text.push(c as char);
                Cell::Str
            }
            (AuxKind::Str, Aux::HexByteArray(s)) => {
                self.text.clear();
                self.text.push_str(s);
                Cell::Str
            }
            // `B` arrays, comma-joined.
            (AuxKind::Str, Aux::ArrayI8(a)) => self.join(a.iter()),
            (AuxKind::Str, Aux::ArrayU8(a)) => self.join(a.iter()),
            (AuxKind::Str, Aux::ArrayI16(a)) => self.join(a.iter()),
            (AuxKind::Str, Aux::ArrayU16(a)) => self.join(a.iter()),
            (AuxKind::Str, Aux::ArrayI32(a)) => self.join(a.iter()),
            (AuxKind::Str, Aux::ArrayU32(a)) => self.join(a.iter()),
            (AuxKind::Str, Aux::ArrayFloat(a)) => self.join(a.iter()),
            _ => Cell::Null,
        };
    }

    }

    /// Every aux column of `fields`, as `(cell, text)`: once through the walk
    /// and once through `Record::aux`, tag by tag.
    fn both_ways(fields: &[RecordField], rec: &Record) -> (Vec<(Cell, String)>, Vec<(Cell, String)>) {
        let view = header(CONTIG, 1000);
        let mut cols = RecordColumns::new(fields, 0);
        cols.extract(rec, &view, Some(call(rec)));
        let walked = cols
            .cols
            .iter()
            .map(|c| (c.cell, if c.cell == Cell::Str { c.text.clone() } else { String::new() }))
            .collect();
        let direct = fields
            .iter()
            .map(|f| {
                let mut c = Column::new(*f, 0);
                if let RecordField::Aux(tag, kind) = *f {
                    c.extract_aux_via_htslib(rec, &tag, kind);
                } else {
                    c.extract(rec, &view, &Ctx { cigar: Some(rec.cigar()), strand: Some(call(rec)) });
                }
                (c.cell, if c.cell == Cell::Str { c.text.clone() } else { String::new() })
            })
            .collect();
        (walked, direct)
    }

    fn every_kind(tags: &[&[u8; 2]]) -> Vec<RecordField> {
        let mut out = Vec::new();
        for t in tags {
            for k in [AuxKind::Int, AuxKind::Float, AuxKind::Str] {
                out.push(RecordField::Aux(**t, k));
            }
        }
        out
    }

    /// Append a raw aux field, bypassing rust-htslib's checks -- the only way
    /// to build a duplicate tag or a damaged value.
    fn raw_append(rec: &mut Record, tag: &[u8; 2], ty: u8, bytes: &[u8]) {
        let r = unsafe {
            htslib::bam_aux_append(
                rec.inner_mut() as *mut htslib::bam1_t,
                tag.as_ptr() as *const std::os::raw::c_char,
                ty as std::os::raw::c_char,
                bytes.len() as i32,
                bytes.as_ptr() as *const u8,
            )
        };
        assert_eq!(r, 0, "bam_aux_append failed");
    }

    /// Every SAM aux type, read as every kind, gives exactly what
    /// `Record::aux` gives -- values, widenings, joins and nulls alike.
    #[test]
    fn raw_aux_reads_agree_with_rust_htslib() {
        let mut rec = built(b"r", b"ACGTACGT", &[Cigar::Match(8)], 0, 0);
        let i8s: Vec<i8> = vec![-3, 0, 7];
        let u8s: Vec<u8> = vec![0, 200];
        let i16s: Vec<i16> = vec![-300, 12];
        let u16s: Vec<u16> = vec![65000];
        let i32s: Vec<i32> = vec![-70000, 1, 2];
        let u32s: Vec<u32> = vec![4_000_000_000];
        let f32s: Vec<f32> = vec![1.5, -0.25, 3.0];
        let empty: Vec<u8> = vec![];
        rec.push_aux(b"XA", Aux::Char(b'q')).unwrap();
        rec.push_aux(b"Xc", Aux::I8(-5)).unwrap();
        rec.push_aux(b"XC", Aux::U8(250)).unwrap();
        rec.push_aux(b"Xs", Aux::I16(-1234)).unwrap();
        rec.push_aux(b"XS", Aux::U16(54321)).unwrap();
        rec.push_aux(b"Xi", Aux::I32(-99999)).unwrap();
        rec.push_aux(b"XI", Aux::U32(3_000_000_000)).unwrap();
        rec.push_aux(b"Xf", Aux::Float(2.5)).unwrap();
        rec.push_aux(b"Xd", Aux::Double(-1e10)).unwrap();
        rec.push_aux(b"XZ", Aux::String("C.zZ.hH")).unwrap();
        rec.push_aux(b"Xe", Aux::String("")).unwrap();
        rec.push_aux(b"XH", Aux::HexByteArray("1AE3")).unwrap();
        rec.push_aux(b"B1", Aux::ArrayI8(AuxArray::from(&i8s))).unwrap();
        rec.push_aux(b"B2", Aux::ArrayU8(AuxArray::from(&u8s))).unwrap();
        rec.push_aux(b"B3", Aux::ArrayI16(AuxArray::from(&i16s))).unwrap();
        rec.push_aux(b"B4", Aux::ArrayU16(AuxArray::from(&u16s))).unwrap();
        rec.push_aux(b"B5", Aux::ArrayI32(AuxArray::from(&i32s))).unwrap();
        rec.push_aux(b"B6", Aux::ArrayU32(AuxArray::from(&u32s))).unwrap();
        rec.push_aux(b"B7", Aux::ArrayFloat(AuxArray::from(&f32s))).unwrap();
        rec.push_aux(b"B8", Aux::ArrayU8(AuxArray::from(&empty))).unwrap();

        let tags: Vec<&[u8; 2]> = vec![
            b"XA", b"Xc", b"XC", b"Xs", b"XS", b"Xi", b"XI", b"Xf", b"Xd", b"XZ", b"Xe", b"XH",
            b"B1", b"B2", b"B3", b"B4", b"B5", b"B6", b"B7", b"B8", b"ZZ",
        ];
        // All at once, so the walk serves many columns from one pass...
        let fields = every_kind(&tags);
        let (walked, direct) = both_ways(&fields, &rec);
        for ((f, w), d) in fields.iter().zip(&walked).zip(&direct) {
            assert_eq!(w, d, "{} disagrees", f.name());
        }
        // ...and one at a time, so a lone column is not relying on company.
        for f in &fields {
            let (w, d) = both_ways(std::slice::from_ref(f), &rec);
            assert_eq!(w, d, "{} disagrees on its own", f.name());
        }
        // A check that the comparison is not vacuous, counted by hand: six
        // integer types read as Int; those six plus `f` and `d` as Float; and
        // `A`, both `Z`s, `H` and all eight arrays as Str. The absent tag and
        // every other pairing is null.
        let filled = walked.iter().filter(|(c, _)| *c != Cell::Null).count();
        assert_eq!(filled, 6 + 8 + 12, "{walked:?}");
    }

    /// `bam_aux_get` returns a tag's first occurrence, and so does the walk.
    #[test]
    fn a_repeated_tag_reads_its_first_occurrence() {
        let mut rec = built(b"r", b"ACGT", &[Cigar::Match(4)], 0, 0);
        raw_append(&mut rec, b"RG", b'Z', b"first\0");
        raw_append(&mut rec, b"NM", b'i', &7i32.to_le_bytes());
        raw_append(&mut rec, b"RG", b'Z', b"second\0");
        let fields = [
            RecordField::Aux(*b"RG", AuxKind::Str),
            RecordField::Aux(*b"NM", AuxKind::Int),
        ];
        let (walked, direct) = both_ways(&fields, &rec);
        assert_eq!(walked, direct);
        assert_eq!(walked[0], (Cell::Str, "first".to_string()));
    }

    /// Text that is not UTF-8 is null, as `Record::aux` makes it, rather than
    /// lossily converted.
    #[test]
    fn non_utf8_aux_text_is_null() {
        let mut rec = built(b"r", b"ACGT", &[Cigar::Match(4)], 0, 0);
        raw_append(&mut rec, b"XZ", b'Z', &[b'a', 0xff, b'b', 0]);
        raw_append(&mut rec, b"RG", b'Z', b"ok\0");
        let fields = every_kind(&[b"XZ", b"RG"]);
        let (walked, direct) = both_ways(&fields, &rec);
        assert_eq!(walked, direct);
        assert_eq!(walked[2].0, Cell::Null, "XZ as text");
        assert_eq!(walked[5], (Cell::Str, "ok".to_string()), "and the walk carried on past it");
    }

    /// A field the walk cannot size ends it: tags before are read, tags after
    /// are null -- which is also where `bam_aux_get` gives up.
    #[test]
    fn a_damaged_field_ends_the_walk_where_htslib_ends_it() {
        let mut rec = built(b"r", b"ACGT", &[Cigar::Match(4)], 0, 0);
        raw_append(&mut rec, b"NM", b'i', &3i32.to_le_bytes());
        raw_append(&mut rec, b"XX", b'?', &[1, 2, 3]);
        raw_append(&mut rec, b"RG", b'Z', b"late\0");
        let fields = [
            RecordField::Aux(*b"NM", AuxKind::Int),
            RecordField::Aux(*b"XX", AuxKind::Str),
            RecordField::Aux(*b"RG", AuxKind::Str),
        ];
        let (walked, direct) = both_ways(&fields, &rec);
        assert_eq!(walked, direct);
        assert_eq!(walked[0].0, Cell::I64(3));
        assert_eq!(walked[2].0, Cell::Null);
    }

    /// The walk reuses its per-record state, so a tag seen on one record must
    /// not leave its column filled on the next record, which lacks it.
    #[test]
    fn aux_state_does_not_leak_between_records() {
        let view = header(CONTIG, 1000);
        let fields = [
            RecordField::Aux(*b"RG", AuxKind::Str),
            RecordField::Aux(*b"NM", AuxKind::Int),
        ];
        let mut cols = RecordColumns::new(&fields, 4);

        let mut with = built(b"a", b"ACGT", &[Cigar::Match(4)], 0, 0);
        with.push_aux(b"RG", Aux::String("g1")).unwrap();
        with.push_aux(b"NM", Aux::I32(2)).unwrap();
        cols.extract(&with, &view, Some(call(&with)));
        cols.append_row();

        let mut without = built(b"b", b"ACGT", &[Cigar::Match(4)], 0, 0);
        without.push_aux(b"NM", Aux::I32(5)).unwrap();
        cols.extract(&without, &view, Some(call(&without)));
        cols.append_row();

        let arrays = cols.finish();
        let rg = arrays[0].as_any().downcast_ref::<StringArray>().unwrap();
        let nm = arrays[1].as_any().downcast_ref::<Int64Array>().unwrap();
        assert_eq!(rg.value(0), "g1");
        assert!(rg.is_null(1), "RG from the first record leaked into the second");
        assert_eq!((nm.value(0), nm.value(1)), (2, 5));
    }
}