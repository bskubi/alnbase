//! Define each record column in one place.
//!
//! A column's name, its Arrow type, its aliases, which groups it belongs to,
//! its help text and how it is read off a
//! [`Record`](rust_htslib::bam::Record) all live here. Nothing else in the
//! crate keeps a list of fields — `cli_query.rs` parses through
//! [`RecordField::parse`], and `record_columns.rs` builds through
//! [`RecordField::data_type`] and `Column::extract`.
//!
//! # Adding a field
//!
//! 1. Add the variant.
//! 2. Add a row to [`RecordField::ALL`] and an arm to `name`/`data_type`.
//! 3. Add an arm to `Column::extract` in `record_columns.rs`.
//! 4. Add an arm to [`RecordField::sam_fields`], naming what the column reads.
//!
//! Steps 2 to 4 are exhaustive matches, so forgetting any of them is a compile
//! error rather than a silently missing column. Step 4 matters more than it
//! looks: on a CRAM input, a part of the record no selected column names is
//! not decoded at all, so an arm that under-reports reads back as a column of
//! placeholder values rather than an error. [`RecordField::encoding`] has a
//! catch-all arm and so does not force the question, which is deliberate: a new
//! column with no opinion about its storage should get the writer's default
//! rather than a guess it was made to write down.

use std::collections::HashMap;

use arrow::datatypes::DataType;

use crate::encoding::ColumnEncoding;

/// How an aux tag's value is read. Matches SAM's own type codes: the letter is
/// what the user writes after the colon, and several letters map to one kind
/// because SAM distinguishes integer widths where Arrow does not.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum AuxKind {
    /// `i c C s S I` — any width, widened to Int64.
    Int,
    /// `f d` — widened to Float64.
    Float,
    /// `Z A H B` — text; `B` arrays are comma-joined.
    Str,
}

impl AuxKind {
    pub fn data_type(self) -> DataType {
        match self {
            AuxKind::Int => DataType::Int64,
            AuxKind::Float => DataType::Float64,
            AuxKind::Str => DataType::Utf8,
        }
    }

    /// The canonical code, for round-tripping a spec back to text.
    pub fn code(self) -> char {
        match self {
            AuxKind::Int => 'i',
            AuxKind::Float => 'f',
            AuxKind::Str => 'Z',
        }
    }

    fn from_code(c: &str) -> Option<Self> {
        match c {
            "i" | "c" | "C" | "s" | "S" | "I" => Some(AuxKind::Int),
            "f" | "d" => Some(AuxKind::Float),
            "Z" | "A" | "H" | "B" => Some(AuxKind::Str),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RecordField {
    // ---- identity and placement
    Qname,
    Flags,
    Tid,
    RefName,
    Pos,
    EndPos,
    UnclippedStart,
    UnclippedEnd,
    MapQ,
    Cigar,
    ReadLen,
    HardClip5p,
    HardClip3p,
    SoftClip5p,
    SoftClip3p,

    // ---- flag bits, one column each
    IsPaired,
    IsProperPair,
    IsUnmapped,
    IsMateUnmapped,
    IsReverse,
    IsMateReverse,
    IsFirstInTemplate,
    IsLastInTemplate,
    IsSecondary,
    IsSupplementary,
    IsDuplicate,
    IsQcFail,

    // ---- strand
    /// Reference strand the pattern matched: `+` when the walk ran forward
    /// along the reference, `-` when it ran backward with both sides reverse
    /// complemented. Emitted as text, not a boolean, because `true` does not
    /// say true of what.
    ///
    /// The conversion strand the run's strand rule named, which is the same
    /// value `walk_alignment` keys on, so it agrees with the emitted
    /// `refr_pos` and `read_off` by construction.
    ///
    /// This is a property of the *fragment*, not the alignment: both mates of a
    /// pair carry the same value, and it differs from `is_reverse` for exactly
    /// the mate that is reverse complemented relative to its fragment. Do not
    /// read it as FLAG 0x10 — that is `is_reverse`, which is also available.
    /// For bisulfite data this is MethylDackel's OT/OB, and it is the strand
    /// the cytosine sits on.
    RefrStrand,
    /// Bisulfite conversion strand: `OT`, `OB`, `CTOT` or `CTOB`.
    ///
    /// The strand of origin the run's strand rule named, so it agrees with the
    /// tagged BAM's XR/XG-style tags and with the walk: `OB` and `CTOB` are
    /// exactly the records whose `strand` is `-`. Splits what `strand` merges
    /// -- read 1 and read 2 of the same fragment -- which per-strand files and
    /// M-bias need.
    ///
    /// Null when the rule names only the conversion strand. Some aligners
    /// record two strands rather than four, and an invented `OT`/`CTOT` split
    /// would be a guess presented as a reading.
    ConvStrand,
    /// Whether the read as sequenced runs backwards along the reference.
    ///
    /// Not FLAG 0x10: `is_reverse` is the bit, and three of the surveyed
    /// aligners write something else in it -- BSBolt puts the converted
    /// reference strand there, so its reverse-flagged records are not the
    /// reverse-sequenced ones. This column is what the run's strand rule
    /// concluded, and it is the direction `off_5p` and `off_3p` count from.
    ReadReverse,

    // ---- mate
    InsertSize,
    MateTid,
    MateRefName,
    MatePos,

    // ---- sequence
    SeqAscii,
    QualPhred,

    /// A two-character SAM tag read as `kind`.
    Aux([u8; 2], AuxKind),
}

impl RecordField {
    /// Every named column, in a sensible default order.
    ///
    /// Aux tags are absent because they are unbounded — they must be named
    /// explicitly. This is also the column order the writer emits, so the
    /// order here is user-visible.
    pub const ALL: &'static [RecordField] = &[
        RecordField::Qname,
        RecordField::Flags,
        RecordField::Tid,
        RecordField::RefName,
        RecordField::Pos,
        RecordField::EndPos,
        RecordField::UnclippedStart,
        RecordField::UnclippedEnd,
        RecordField::MapQ,
        RecordField::Cigar,
        RecordField::ReadLen,
        RecordField::HardClip5p,
        RecordField::HardClip3p,
        RecordField::SoftClip5p,
        RecordField::SoftClip3p,
        RecordField::IsPaired,
        RecordField::IsProperPair,
        RecordField::IsUnmapped,
        RecordField::IsMateUnmapped,
        RecordField::IsReverse,
        RecordField::IsMateReverse,
        RecordField::IsFirstInTemplate,
        RecordField::IsLastInTemplate,
        RecordField::IsSecondary,
        RecordField::IsSupplementary,
        RecordField::IsDuplicate,
        RecordField::IsQcFail,
        RecordField::RefrStrand,
        RecordField::ConvStrand,
        RecordField::ReadReverse,
        RecordField::InsertSize,
        RecordField::MateTid,
        RecordField::MateRefName,
        RecordField::MatePos,
        RecordField::SeqAscii,
        RecordField::QualPhred,
    ];

    /// The default extra columns when `--field` is not given.
    ///
    /// Deliberately small. Every row already carries the hit columns —
    /// `shard`, `record_id`, `name`, `off_5p`, `off_3p`, `refr_pos`, `qual`,
    /// and the captured bases (`read_base`/`refr_base`, or the `capture_*`
    /// lists) — so the default only adds what those cannot express:
    ///
    /// - `ref_name`, because `refr_pos` alone does not say which contig.
    /// - `strand`, because a methylation call is meaningless without knowing
    ///   which reference strand it was read from. This is the same derivation
    ///   `walk_alignment` keys on, so it agrees with the emitted coordinates by
    ///   construction. `is_reverse` is the raw flag bit and a different
    ///   question; it is one `-F` away.
    ///
    /// Everything else is a filter someone can apply upstream — `mapq` and the
    /// QC flags are decisions better made once, on the BAM, than carried on
    /// every hit row.
    pub const CORE: &'static [RecordField] = &[
        RecordField::RefName,
        RecordField::RefrStrand,
    ];

    /// The boolean columns.
    ///
    /// Flag bits, and `read_reverse`, which reads like one but is the strand
    /// rule's conclusion rather than a bit off the record.
    pub const BOOLS: &'static [RecordField] = &[
        RecordField::IsPaired,
        RecordField::IsProperPair,
        RecordField::IsUnmapped,
        RecordField::IsMateUnmapped,
        RecordField::IsReverse,
        RecordField::IsMateReverse,
        RecordField::IsFirstInTemplate,
        RecordField::IsLastInTemplate,
        RecordField::IsSecondary,
        RecordField::IsSupplementary,
        RecordField::IsDuplicate,
        RecordField::IsQcFail,
        RecordField::ReadReverse,
    ];

    /// Everything about the mate.
    ///
    /// Deliberately not in `ALL` order: a selection is emitted in the order it
    /// was given, so this puts the mate's placement first and its flag bits
    /// after, which is how someone reading the columns wants them. `bools` and
    /// `coords` happen to match `ALL`; this one does not, and that is a choice
    /// rather than an oversight.
    pub const MATE: &'static [RecordField] = &[
        RecordField::InsertSize,
        RecordField::MateTid,
        RecordField::MateRefName,
        RecordField::MatePos,
        RecordField::IsMateUnmapped,
        RecordField::IsMateReverse,
    ];

    /// Reference placement, including the clip-aware bounds.
    pub const COORDS: &'static [RecordField] = &[
        RecordField::Tid,
        RecordField::RefName,
        RecordField::Pos,
        RecordField::EndPos,
        RecordField::UnclippedStart,
        RecordField::UnclippedEnd,
    ];

    /// Canonical column name. Aux tags name themselves.
    pub fn name(self) -> String {
        match self.static_name() {
            Some(n) => n.to_string(),
            None => match self {
                RecordField::Aux(tag, _) => String::from_utf8_lossy(&tag).into_owned(),
                _ => unreachable!("every non-aux field has a static name"),
            },
        }
    }

    /// The name as a `&'static str`, or `None` for aux tags, whose names are
    /// only known at runtime.
    pub fn static_name(self) -> Option<&'static str> {
        use RecordField as F;
        Some(match self {
            F::Qname => "qname",
            F::Flags => "flags",
            F::Tid => "tid",
            F::RefName => "ref_name",
            F::Pos => "pos",
            F::EndPos => "end_pos",
            F::UnclippedStart => "unclipped_start",
            F::UnclippedEnd => "unclipped_end",
            F::MapQ => "mapq",
            F::Cigar => "cigar",
            F::ReadLen => "read_len",
            F::HardClip5p => "hard_clip_5p",
            F::HardClip3p => "hard_clip_3p",
            F::SoftClip5p => "soft_clip_5p",
            F::SoftClip3p => "soft_clip_3p",
            F::IsPaired => "is_paired",
            F::IsProperPair => "is_proper_pair",
            F::IsUnmapped => "is_unmapped",
            F::IsMateUnmapped => "is_mate_unmapped",
            F::IsReverse => "is_reverse",
            F::IsMateReverse => "is_mate_reverse",
            F::IsFirstInTemplate => "is_first_in_template",
            F::IsLastInTemplate => "is_last_in_template",
            F::IsSecondary => "is_secondary",
            F::IsSupplementary => "is_supplementary",
            F::IsDuplicate => "is_duplicate",
            F::IsQcFail => "is_qc_fail",
            F::RefrStrand => "strand",
            F::ConvStrand => "conv_strand",
            F::ReadReverse => "read_reverse",
            F::InsertSize => "insert_size",
            F::MateTid => "mate_tid",
            F::MateRefName => "mate_ref_name",
            F::MatePos => "mate_pos",
            F::SeqAscii => "seq_ascii",
            F::QualPhred => "qual_phred",
            F::Aux(..) => return None,
        })
    }

    /// The Arrow type of the column.
    ///
    /// Every record column is nullable. Several genuinely can be absent (an
    /// unmapped mate has no `mate_ref_name`, `qual_phred` is empty when QUAL is
    /// `*`, an aux tag may not be on the record), and declaring the rest
    /// non-nullable would only buy a panic at `RecordBatch::try_new` the first
    /// time a file violated the assumption.
    pub fn data_type(self) -> DataType {
        use RecordField as F;
        match self {
            F::Qname
            | F::RefName
            | F::Cigar
            | F::MateRefName
            | F::SeqAscii
            | F::QualPhred
            | F::RefrStrand
            | F::ConvStrand => DataType::Utf8,
            F::Flags => DataType::UInt16,
            F::MapQ => DataType::UInt8,
            F::Tid | F::MateTid => DataType::Int32,
            F::Pos
            | F::EndPos
            | F::UnclippedStart
            | F::UnclippedEnd
            | F::ReadLen
            | F::HardClip5p
            | F::HardClip3p
            | F::SoftClip5p
            | F::SoftClip3p
            | F::InsertSize
            | F::MatePos => DataType::Int64,
            F::IsPaired
            | F::IsProperPair
            | F::IsUnmapped
            | F::IsMateUnmapped
            | F::IsReverse
            | F::IsMateReverse
            | F::IsFirstInTemplate
            | F::IsLastInTemplate
            | F::IsSecondary
            | F::IsSupplementary
            | F::IsDuplicate
            | F::IsQcFail
            | F::ReadReverse => DataType::Boolean,
            F::Aux(_, kind) => kind.data_type(),
        }
    }

    /// How the column is stored in parquet.
    ///
    /// The hits file is denormalized, so a record that fires several queries
    /// writes several contiguous rows carrying identical record columns. Every
    /// choice below is about that shape plus the column's own cardinality:
    ///
    /// - `qname` is long, repeats in runs, and consecutive names off one
    ///   instrument share everything but their coordinates, so it is the one
    ///   column where prefix-delta beats both alternatives. A dictionary over
    ///   it would be tens of thousands of long strings and would overflow its
    ///   page part way through a row group; plain would spell each name out
    ///   once per row that mentions it.
    /// - `seq_ascii` and `qual_phred` are a distinct long string per record
    ///   with no structure across records. Nothing to build a dictionary from,
    ///   nothing to delta against: hand them to the codec.
    /// - Coordinates move in small steps, and in a coordinate-sorted BAM they
    ///   move monotonically, so delta stores them in a few bits each. It costs
    ///   nothing when the BAM is unsorted -- delta of an arbitrary jump is the
    ///   same width as the value.
    /// - Everything else -- names, flags, strand, CIGAR, mapq -- has few
    ///   distinct values relative to the rows, which is what a dictionary is.
    ///
    /// Aux tags are guesses about a tag this code has never seen, so they go by
    /// kind: text and integer tags are usually low-cardinality labels and codes
    /// (`RG`, `NM`, `AS`), floats usually are not. The guess is wrong for the
    /// per-base array tags -- `MM`, `ML`, and anything else that is really a
    /// vector serialised into a string, which is distinct per record and long.
    /// `--column-encoding ML=plain` is the escape hatch, and the help text for
    /// that flag names those tags.
    pub fn encoding(self) -> ColumnEncoding {
        use ColumnEncoding as E;
        use RecordField as F;
        match self {
            F::Qname => E::DeltaText,
            F::SeqAscii | F::QualPhred => E::Plain,
            F::Pos
            | F::EndPos
            | F::UnclippedStart
            | F::UnclippedEnd
            | F::MatePos
            | F::InsertSize => E::Delta,
            F::Aux(_, AuxKind::Float) => E::Plain,
            _ => E::Dictionary,
        }
    }

    /// Whether reading this column needs the record's CIGAR parsed.
    ///
    /// `Record::cigar()` allocates, so `RecordColumns` builds it once per
    /// record and only when some selected column asks for it.
    pub fn needs_cigar(self) -> bool {
        matches!(
            self,
            RecordField::Cigar
                | RecordField::EndPos
                | RecordField::UnclippedStart
                | RecordField::UnclippedEnd
        )
    }

    /// The parts of the record this column reads, as htslib `SAM_*` bits.
    ///
    /// A CRAM reader is told the union of these over every column it serves,
    /// and skips decoding the data series nothing asked for -- read names and
    /// mate fields, most often, which a methylation scan rarely wants. A BAM
    /// reader ignores the setting, since BAM decodes whole records regardless.
    ///
    /// Exhaustive with no catch-all, so a new field has to answer. When in
    /// doubt an arm should over-report: naming a field that is not needed
    /// costs decode time, while omitting one that is silently yields a
    /// placeholder value in place of the real one.
    ///
    /// The mate-dependent flag bits (`0x8`, `0x20`) ask for the mate fields
    /// too. CRAM stores them either explicitly or by reference to the mate
    /// record in the same slice, and which data series that needs depends on
    /// how the file was written, so a flag column asks for all of it.
    pub fn sam_fields(self) -> u32 {
        use rust_htslib::htslib as h;
        use RecordField as F;
        let mate = h::sam_fields_SAM_RNEXT | h::sam_fields_SAM_PNEXT | h::sam_fields_SAM_TLEN;
        match self {
            F::Qname => h::sam_fields_SAM_QNAME,
            F::Flags | F::IsMateUnmapped | F::IsMateReverse => h::sam_fields_SAM_FLAG | mate,
            F::IsPaired
            | F::IsProperPair
            | F::IsUnmapped
            | F::IsReverse
            | F::IsFirstInTemplate
            | F::IsLastInTemplate
            | F::IsSecondary
            | F::IsSupplementary
            | F::IsDuplicate
            | F::IsQcFail
            | F::RefrStrand
            | F::ConvStrand
            | F::ReadReverse => h::sam_fields_SAM_FLAG,
            F::Tid | F::RefName => h::sam_fields_SAM_RNAME,
            F::Pos => h::sam_fields_SAM_POS,
            F::EndPos | F::UnclippedStart | F::UnclippedEnd => {
                h::sam_fields_SAM_POS | h::sam_fields_SAM_CIGAR
            }
            F::MapQ => h::sam_fields_SAM_MAPQ,
            F::Cigar => h::sam_fields_SAM_CIGAR,
            F::ReadLen | F::SeqAscii => h::sam_fields_SAM_SEQ,
            F::HardClip5p | F::HardClip3p | F::SoftClip5p | F::SoftClip3p => {
                h::sam_fields_SAM_FLAG | h::sam_fields_SAM_CIGAR
            }
            F::QualPhred => h::sam_fields_SAM_QUAL,
            F::InsertSize => h::sam_fields_SAM_TLEN,
            F::MateTid | F::MateRefName => h::sam_fields_SAM_RNEXT,
            F::MatePos => h::sam_fields_SAM_RNEXT | h::sam_fields_SAM_PNEXT,
            F::Aux(tag, _) if &tag == b"RG" => h::sam_fields_SAM_RGAUX,
            F::Aux(..) => h::sam_fields_SAM_AUX,
        }
    }

    /// Whether this column reads a tag CRAM may leave out of the file and
    /// regenerate from the reference on decode.
    ///
    /// `MD` and `NM` are the two. Regenerating them means comparing every
    /// aligned base against the reference as the record is decoded, so a run
    /// that reads neither turns that off -- see
    /// [`crate::bam_io::limit_cram_decoding`].
    pub fn needs_md_nm(self) -> bool {
        matches!(self, RecordField::Aux(tag, _) if &tag == b"MD" || &tag == b"NM")
    }

    /// One-line description, for `--help`.
    pub fn help(self) -> &'static str {
        use RecordField as F;
        match self {
            F::Qname => "read name",
            F::Flags => "the FLAG bitfield as a number",
            F::Tid => "reference index; null when unplaced",
            F::RefName => "reference name; null when unplaced",
            F::Pos => "0-based leftmost aligned position",
            F::EndPos => "0-based exclusive end of the alignment",
            F::UnclippedStart => "pos minus any leading soft/hard clip",
            F::UnclippedEnd => "end_pos plus any trailing soft/hard clip",
            F::MapQ => "mapping quality",
            F::Cigar => "CIGAR string",
            F::ReadLen => "length of SEQ (excludes hard-clipped bases)",
            F::HardClip5p => "hard-clipped bases at the read's sequenced 5' end (off_5p already includes them)",
            F::HardClip3p => "hard-clipped bases at the read's sequenced 3' end (off_3p already includes them)",
            F::SoftClip5p => "soft-clipped bases at the read's sequenced 5' end (off_5p already includes them)",
            F::SoftClip3p => "soft-clipped bases at the read's sequenced 3' end (off_3p already includes them)",
            F::RefrStrand => "reference strand the pattern matched, '+' or '-'",
            F::ConvStrand => "strand of origin: OT, OB, CTOT or CTOB; null when the rule names only two strands",
            F::ReadReverse => "true when the read was sequenced backwards along the reference (not FLAG 0x10)",
            F::InsertSize => "TLEN",
            F::MateTid => "mate reference index; null when the mate is unplaced",
            F::MateRefName => "mate reference name; null when the mate is unplaced",
            F::MatePos => "mate position; null when the mate is unplaced",
            F::SeqAscii => "SEQ as ASCII bases",
            F::QualPhred => "QUAL as phred+33 text; null when QUAL is '*'",
            F::Aux(..) => "aux tag",
            _ => "flag bit",
        }
    }
}

// ------------------------------------------------------------------ groups --

/// Named groups, built from the tables above so they cannot drift from them.
pub const GROUPS: &[(&str, &[RecordField])] = &[
    ("core", RecordField::CORE),
    ("all", RecordField::ALL),
    ("bools", RecordField::BOOLS),
    ("mate", RecordField::MATE),
    ("coords", RecordField::COORDS),
];

/// Accepted alternative spellings, mostly the SAM column names.
const ALIASES: &[(&str, RecordField)] = &[
    ("rname", RecordField::RefName),
    ("refname", RecordField::RefName),
    ("rnext", RecordField::MateRefName),
    ("mrnm", RecordField::MateRefName),
    ("pnext", RecordField::MatePos),
    ("mpos", RecordField::MatePos),
    ("tlen", RecordField::InsertSize),
    ("isize", RecordField::InsertSize),
    ("seq", RecordField::SeqAscii),
    ("qual", RecordField::QualPhred),
    ("mtid", RecordField::MateTid),
    ("is_qcfail", RecordField::IsQcFail),
    // Kept so scripts written against the old boolean column still select
    // something, though the values are now '+'/'-' rather than true/false.
    ("is_ot", RecordField::RefrStrand),
    ("refr_strand", RecordField::RefrStrand),
    ("is_read1", RecordField::IsFirstInTemplate),
    ("is_read2", RecordField::IsLastInTemplate),
    ("qname", RecordField::Qname),
];

// ------------------------------------------------------------------- parse --

/// One `--field` token: either a single column or a named group.
#[derive(Clone, Debug)]
pub enum FieldSpec {
    One(RecordField),
    Group(&'static [RecordField]),
}

impl FieldSpec {
    pub fn expand(&self) -> &[RecordField] {
        match self {
            FieldSpec::One(f) => std::slice::from_ref(f),
            FieldSpec::Group(g) => g,
        }
    }
}

impl RecordField {
    /// Parse one `--field` token.
    ///
    /// Accepts a canonical name, an alias, a group name, or `TAG:type` using
    /// SAM's own type codes. Case- and hyphen-insensitive for names; aux tags
    /// and their type codes are case-sensitive, because SAM's codes are.
    pub fn parse(s: &str) -> Result<FieldSpec, String> {
        let t = s.trim();
        if t.is_empty() {
            return Err("empty field name".to_string());
        }

        // `TAG:type`, using SAM's own type codes.
        if let Some((tag, ty)) = t.split_once(':') {
            let b = tag.as_bytes();
            if b.len() != 2 || !b.iter().all(|c| c.is_ascii_alphanumeric()) {
                return Err(format!(
                    "`{tag}` is not a valid aux tag: expected two alphanumeric characters, e.g. NM:i"
                ));
            }
            let kind = AuxKind::from_code(ty).ok_or_else(|| {
                format!(
                    "unknown aux type `{ty}` in `{t}`: expected one of \
                     i/c/C/s/S/I (integer), f/d (float), Z/A/H/B (text)"
                )
            })?;
            return Ok(FieldSpec::One(RecordField::Aux([b[0], b[1]], kind)));
        }

        let norm = t.to_ascii_lowercase().replace('-', "_");

        if let Some((_, g)) = GROUPS.iter().find(|(n, _)| *n == norm) {
            return Ok(FieldSpec::Group(g));
        }
        if let Some(f) = RecordField::ALL
            .iter()
            .copied()
            .find(|f| f.static_name() == Some(norm.as_str()))
        {
            return Ok(FieldSpec::One(f));
        }
        if let Some((_, f)) = ALIASES.iter().find(|(n, _)| *n == norm) {
            return Ok(FieldSpec::One(*f));
        }

        Err(format!(
            "unknown field `{t}`\n  named columns: {}\n  groups: {}\n  aux tags: TAG:type, e.g. NM:i, AS:f, RG:Z",
            RecordField::ALL
                .iter()
                .map(|f| f.name())
                .collect::<Vec<_>>()
                .join(", "),
            GROUPS.iter().map(|(n, _)| *n).collect::<Vec<_>>().join(", "),
        ))
    }

    /// Resolve a parsed selection: groups expanded, duplicates dropped, order
    /// as given.
    ///
    /// `field` is the base set and *replaces* the default; when it is empty,
    /// [`CORE`](Self::CORE) stands. `add` is appended to whichever of those
    /// applies, so `--add-field` extends a selection instead of replacing it,
    /// and extending the default costs no restating of it.
    ///
    /// Duplicates are dropped by *column name*, not by variant, so `NM:i` and
    /// `NM:Z` collide — two columns cannot share a name in one schema, and the
    /// first spelling wins. That also makes adding a column already in the base
    /// set a no-op rather than an error, which is what you want when `-F` is
    /// scripted against a `-f` chosen elsewhere.
    pub fn resolve(field: &[FieldSpec], add: &[FieldSpec]) -> Vec<RecordField> {
        let mut out: Vec<RecordField> = Vec::new();
        let mut seen: HashMap<String, RecordField> = HashMap::new();

        let mut push = |f: RecordField, out: &mut Vec<RecordField>| {
            if seen.insert(f.name(), f).is_none() {
                out.push(f);
            }
        };

        if field.is_empty() {
            for f in RecordField::CORE {
                push(*f, &mut out);
            }
        } else {
            for spec in field {
                for f in spec.expand() {
                    push(*f, &mut out);
                }
            }
        }
        for spec in add {
            for f in spec.expand() {
                push(*f, &mut out);
            }
        }
        out
    }

    /// `--help` text, generated from the tables so the list cannot go stale.
    pub fn long_help() -> String {
        let mut s = String::from(
            "Columns to write alongside each hit, comma-separated and/or repeated.\n\n\
             Order is preserved and duplicates are dropped, so this also sets column order.\n\
             Defaults to `core` if not given.\n\n\
             `--field` REPLACES the default; `--add-field`/`-F` ADDS to it, so\n\
             `-F NM:i` is `core` plus that tag and needs no restating of `core`. The two\n\
             combine: `-f bools -F NM:i` is the flag columns plus the tag.\n\n\
             These are in addition to the hit columns, which are on every row whatever\n\
             you pass here: shard, record_id, name, off_5p, off_3p,\n\
             refr_pos, qual, and the matched bases. off_5p counts from the 5' end of\n\
             the read as sequenced and off_3p from the 3' end, so the last base\n\
             sequenced is off_3p 0.\n\nGroups:\n",
        );
        for (name, members) in GROUPS {
            let list = if *name == "all" {
                "every named column below".to_string()
            } else {
                members.iter().map(|f| f.name()).collect::<Vec<_>>().join(",")
            };
            s.push_str(&format!("  {name:<8}{list}\n"));
        }
        s.push_str("\nNamed columns:\n");
        for f in RecordField::ALL {
            s.push_str(&format!("  {:<22}{}\n", f.name(), f.help()));
        }
        s.push_str(
            "\nAux tags are given as TAG:type using SAM's type codes:\n  \
             i c C s S I  integer   (e.g. NM:i)   -> Int64\n  \
             f d          float     (e.g. AS:f)   -> Float64\n  \
             Z A H B      text      (e.g. RG:Z)   -> Utf8, B arrays comma-joined\n\n\
             Absent aux tags are null, as are type-mismatched ones.\n",
        );
        s
    }
}

// ------------------------------------------------------------------- tests --

#[cfg(test)]
mod tests {
    use super::*;

    fn specs(args: &[&str]) -> Vec<FieldSpec> {
        args.iter().map(|a| RecordField::parse(a).unwrap()).collect()
    }

    /// A base selection, no additions.
    fn fields(args: &[&str]) -> Vec<RecordField> {
        RecordField::resolve(&specs(args), &[])
    }

    /// A base selection plus `--add-field` entries.
    fn fields_plus(base: &[&str], add: &[&str]) -> Vec<RecordField> {
        RecordField::resolve(&specs(base), &specs(add))
    }

    #[test]
    fn every_name_round_trips() {
        for f in RecordField::ALL {
            assert_eq!(fields(&[&f.name()]), vec![*f], "{}", f.name());
        }
    }

    #[test]
    fn the_default_is_small_and_real() {
        assert_eq!(fields(&[]), RecordField::CORE.to_vec());
        for f in RecordField::CORE {
            assert!(RecordField::ALL.contains(f), "{f:?} is not in ALL");
        }
    }

    /// The record columns and the hit columns share one schema, so a record
    /// column may not take a name the hit side already uses.
    ///
    /// Asked of the builder rather than of a list written down here. A hand
    /// list is a second statement of the hit schema, and this one had already
    /// drifted from it -- it still named `pattern` and `qual_min`, columns the
    /// file has not had for some time, and would have said nothing about a new
    /// one. Both capture shapes, because they declare different columns.
    #[test]
    fn no_record_column_collides_with_a_hit_column() {
        use crate::hits::HitBuilder;

        for flat in [true, false] {
            let hit: Vec<String> = HitBuilder::new(&[], 1, flat, 0)
                .arrow_fields()
                .iter()
                .map(|f| f.name().clone())
                .collect();
            for f in RecordField::ALL {
                assert!(
                    !hit.contains(&f.name()),
                    "record column `{}` collides with a hit column (flat={flat})",
                    f.name()
                );
            }
        }
    }

    #[test]
    fn every_field_has_a_name_and_a_type() {
        for f in RecordField::ALL {
            assert!(f.static_name().is_some(), "{f:?} has no static name");
            let _ = f.data_type();
        }
    }

    #[test]
    fn strand_is_text_not_a_flag_bit() {
        // It reads like a flag column but is not one: BOOLS must not carry it,
        // and it must not be typed Boolean.
        assert!(!RecordField::BOOLS.contains(&RecordField::RefrStrand));
        assert_eq!(RecordField::RefrStrand.data_type(), DataType::Utf8);
        assert_eq!(RecordField::RefrStrand.name(), "strand");
        // The old spelling still selects it.
        assert_eq!(fields(&["is_ot"]), vec![RecordField::RefrStrand]);
        // And is distinct from the raw flag bit, which is still available.
        assert_ne!(fields(&["strand"]), fields(&["is_reverse"]));
    }

    #[test]
    fn bools_is_exactly_the_boolean_columns() {
        // `BOOLS` is written out separately from `ALL`, so nothing but this
        // holds them in step. The direction that matters is the one a reviewer
        // would not notice: a Boolean column added to `ALL` and forgotten here
        // makes `-F bools` quietly incomplete, with no error anywhere.
        let from_all: Vec<RecordField> = RecordField::ALL
            .iter()
            .copied()
            .filter(|f| f.data_type() == DataType::Boolean)
            .collect();
        assert_eq!(
            RecordField::BOOLS.to_vec(),
            from_all,
            "BOOLS and the Boolean columns of ALL have diverged"
        );
    }

    /// Every mate-related column belongs to the `mate` group.
    ///
    /// Unlike `bools` there is no type to derive this from, so the membership
    /// rule is the name: a column about the mate says so.
    #[test]
    fn mate_holds_every_mate_column() {
        let by_name: Vec<RecordField> = RecordField::ALL
            .iter()
            .copied()
            .filter(|f| {
                let n = f.name();
                n.starts_with("mate_") || n == "insert_size" || n == "is_mate_unmapped"
                    || n == "is_mate_reverse"
            })
            .collect();
        for f in by_name {
            assert!(
                RecordField::MATE.contains(&f),
                "{:?} looks like a mate column but the `mate` group omits it",
                f
            );
        }
    }

    #[test]
    fn names_are_unique() {
        let mut seen = std::collections::HashSet::new();
        for f in RecordField::ALL {
            assert!(seen.insert(f.name()), "duplicate column name {}", f.name());
        }
    }

    #[test]
    fn groups_contain_only_known_fields() {
        for (name, members) in GROUPS {
            for m in *members {
                assert!(
                    RecordField::ALL.contains(m),
                    "group `{name}` has {m:?}, which is not in ALL"
                );
            }
        }
    }

    #[test]
    fn aux_type_codes_route_to_the_right_kind() {
        assert_eq!(fields(&["NM:i"]), vec![RecordField::Aux(*b"NM", AuxKind::Int)]);
        assert_eq!(fields(&["AS:f"]), vec![RecordField::Aux(*b"AS", AuxKind::Float)]);
        assert_eq!(fields(&["RG:Z"]), vec![RecordField::Aux(*b"RG", AuxKind::Str)]);
        assert_eq!(fields(&["NM:C"]), vec![RecordField::Aux(*b"NM", AuxKind::Int)]);
    }

    #[test]
    fn aliases_resolve() {
        assert_eq!(fields(&["rname"]), vec![RecordField::RefName]);
        assert_eq!(fields(&["tlen"]), vec![RecordField::InsertSize]);
        assert_eq!(fields(&["SEQ"]), vec![RecordField::SeqAscii]);
        assert_eq!(fields(&["is-qcfail"]), vec![RecordField::IsQcFail]);
    }

    #[test]
    fn order_is_preserved_and_dupes_dropped() {
        assert_eq!(
            fields(&["mapq", "cigar", "mapq", "rname", "ref_name"]),
            vec![RecordField::MapQ, RecordField::Cigar, RecordField::RefName]
        );
    }

    #[test]
    fn group_then_member_does_not_reorder() {
        // `mapq` is not in CORE, so it lands after it rather than moving.
        let mut want = RecordField::CORE.to_vec();
        want.push(RecordField::MapQ);
        assert_eq!(fields(&["core", "mapq"]), want);
        // And a member already in the group does not move to the front.
        assert_eq!(fields(&["core", "ref_name"]), RecordField::CORE.to_vec());
    }

    #[test]
    fn flags_is_the_column_not_the_group() {
        assert_eq!(fields(&["flags"]), vec![RecordField::Flags]);
        assert_eq!(fields(&["bools"]), RecordField::BOOLS.to_vec());
    }

    #[test]
    fn same_tag_two_types_does_not_duplicate_a_column() {
        assert_eq!(fields(&["NM:i", "NM:Z"]), vec![RecordField::Aux(*b"NM", AuxKind::Int)]);
    }

    #[test]
    fn add_extends_the_default_without_restating_it() {
        let got = fields_plus(&[], &["NM:i", "is_reverse"]);
        let mut want = RecordField::CORE.to_vec();
        want.push(RecordField::Aux(*b"NM", AuxKind::Int));
        want.push(RecordField::IsReverse);
        assert_eq!(got, want);
    }

    #[test]
    fn add_extends_an_explicit_selection_and_appends_after_it() {
        assert_eq!(
            fields_plus(&["mapq"], &["cigar"]),
            vec![RecordField::MapQ, RecordField::Cigar]
        );
    }

    #[test]
    fn adding_something_already_present_is_a_noop() {
        assert_eq!(fields_plus(&["mapq", "cigar"], &["mapq"]), fields(&["mapq", "cigar"]));
        // Including when it is present by way of the default.
        assert_eq!(fields_plus(&[], &["ref_name"]), fields(&[]));
    }

    #[test]
    fn bad_input_is_rejected() {
        assert!(RecordField::parse("nope").is_err());
        assert!(RecordField::parse("NMX:i").is_err());
        assert!(RecordField::parse("NM:q").is_err());
        assert!(RecordField::parse("").is_err());
    }
}