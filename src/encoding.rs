//! Say how one output column is stored. The statement is made for each column,
//! and not for the file as a whole.
//!
//! A parquet writer's default is to dictionary-encode everything and fall back
//! to plain once the dictionary outgrows its page. That default is right for
//! most of a hits file and wrong for the columns that dominate its size, and
//! the two failure modes point in opposite directions:
//!
//! - A dictionary built over a column that is nearly all distinct values --
//!   `seq_ascii`, `qual_phred`, an `ML:B` array -- costs a hash per row, grows
//!   to the 1 MB page limit, and is then abandoned mid-chunk. The work is
//!   wasted, and the file pays for the bytes twice. It pays once in the dead
//!   dictionary, and once in the plain values that follow it.
//! - Plain-encoding a column that repeats -- `ref_name`, `strand`, `expr` --
//!   hands the codec megabytes of identical bytes to rediscover, per page.
//!
//! Which side a column falls on is a property of the column, known where the
//! column is defined, so that is where it is stated: [`RecordField::encoding`]
//! for record columns and `HitBuilder::column_encodings` for the hit columns.
//! This module is only the vocabulary and the one place that translates it into
//! writer properties.
//!
//! The table below measures what this is worth. Against the same run with every
//! column forced to `dict`, it wrote 30.7 MB instead of 49.7 MB under LZ4, and
//! 19.2 MB instead of 30.1 MB under zstd. The wall time was the same either way.
//! The dictionaries it removes were the expensive ones to build.
//!
//! # The shape of a hits file
//!
//! Rows are denormalized. A record that fires n queries writes n rows, those
//! rows are contiguous, and they carry the same record columns. A record column
//! is therefore not one value per row. It is a run, and a run is what the
//! dictionary and delta encodings are for.
//!
//! `qname` is the clearest case. Its values are long, they repeat in runs, and
//! consecutive names off one instrument share a long prefix. [`DeltaText`]
//! encodes all three of those properties. A dictionary of 30 000 distinct names
//! encodes none of them.
//!
//! [`RecordField::encoding`]: crate::record_field::RecordField::encoding
//! [`DeltaText`]: ColumnEncoding::DeltaText

use arrow::datatypes::DataType;
use parquet::basic::Encoding;
use parquet::file::properties::WriterPropertiesBuilder;
use parquet::schema::types::ColumnPath;

/// The storage policy for one column.
///
/// Every variant is legal for every column *as far as the file format goes*,
/// except the two delta encodings, which parquet defines only over particular
/// physical types. [`ColumnEncoding::supports`] is the check, and the writer
/// runs it against the declared schema before opening the file, so a bad
/// `--column-encoding` fails at the command line rather than at the first
/// flush.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum ColumnEncoding {
    /// Few distinct values, or many repeats. Dictionary + RLE, which collapses
    /// a run to a run length whatever the values are.
    ///
    /// The default because it is the parquet writer's default: a column with no
    /// opinion gets what it would have got anyway.
    #[default]
    Dictionary,
    /// Mostly distinct, and nothing structural to exploit. No dictionary, no
    /// per-value bookkeeping -- the block codec does what it can and nothing
    /// pays for a hash table that will be thrown away.
    Plain,
    /// Integers that move in small steps: coordinates, offsets, monotonic ids.
    /// Delta-packed, so a column of 8-byte positions costs a few bits each.
    ///
    /// Integer columns only.
    Delta,
    /// Text whose neighbours share a prefix, repeat outright, or both. Each
    /// value is stored as "how much of the previous one to keep" plus the
    /// remainder, so an exact repeat costs a length and no bytes.
    ///
    /// Byte-array columns only.
    DeltaText,
}

impl ColumnEncoding {
    /// Whether parquet defines this encoding over `dt`'s physical type.
    ///
    /// The delta encodings are the only ones that can be wrong: parquet defines
    /// `DELTA_BINARY_PACKED` over INT32/INT64 and `DELTA_BYTE_ARRAY` over
    /// BYTE_ARRAY, and a writer handed either over the wrong type fails when it
    /// first tries to flush a page -- long after the run started.
    pub fn supports(self, dt: &DataType) -> bool {
        use DataType as D;
        match self {
            ColumnEncoding::Dictionary | ColumnEncoding::Plain => true,
            ColumnEncoding::Delta => matches!(
                dt,
                D::Int8 | D::Int16 | D::Int32 | D::Int64
                    | D::UInt8 | D::UInt16 | D::UInt32 | D::UInt64
            ),
            ColumnEncoding::DeltaText => {
                matches!(dt, D::Utf8 | D::LargeUtf8 | D::Binary | D::LargeBinary)
            }
        }
    }

    /// What to accept on the command line, and what to print when refusing one.
    pub const NAMES: &'static [(&'static str, ColumnEncoding)] = &[
        ("dict", ColumnEncoding::Dictionary),
        ("plain", ColumnEncoding::Plain),
        ("delta", ColumnEncoding::Delta),
        ("delta-text", ColumnEncoding::DeltaText),
    ];

    pub fn parse(s: &str) -> Result<Self, String> {
        let norm = s.trim().to_ascii_lowercase().replace('_', "-");
        Self::NAMES
            .iter()
            .find(|(n, _)| *n == norm)
            .map(|(_, e)| *e)
            .ok_or_else(|| {
                format!(
                    "unknown encoding `{s}`: expected one of {}",
                    Self::NAMES.iter().map(|(n, _)| *n).collect::<Vec<_>>().join(", ")
                )
            })
    }

    pub fn name(self) -> &'static str {
        Self::NAMES
            .iter()
            .find(|(_, e)| *e == self)
            .map(|(n, _)| *n)
            .expect("every variant is in NAMES")
    }

    /// Set this policy on one column path.
    ///
    /// Both halves are set together on purpose. A dictionary takes precedence
    /// over the fallback encoding, so `set_column_encoding` alone is silently
    /// inert on a column that still has its dictionary -- the encoding named
    /// here would only appear after the dictionary overflowed, which is the
    /// case we are trying to avoid.
    pub fn apply(
        self,
        props: WriterPropertiesBuilder,
        path: ColumnPath,
    ) -> WriterPropertiesBuilder {
        let encoding = match self {
            ColumnEncoding::Dictionary => {
                return props.set_column_dictionary_enabled(path, true)
            }
            ColumnEncoding::Plain => Encoding::PLAIN,
            ColumnEncoding::Delta => Encoding::DELTA_BINARY_PACKED,
            ColumnEncoding::DeltaText => Encoding::DELTA_BYTE_ARRAY,
        };
        props
            .set_column_dictionary_enabled(path.clone(), false)
            .set_column_encoding(path, encoding)
    }
}

/// A `COLUMN=ENCODING` override, as `--column-encoding` takes them.
///
/// Kept as a name rather than a resolved column because the two sides are
/// resolved at different times: the flag is parsed before the field selection
/// is known, and the name may be a hit column, a record column or an aux tag.
/// The writer matches it against the schema it is about to declare, which is
/// the only point at which all three are one list.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EncodingOverride {
    pub column: String,
    pub encoding: ColumnEncoding,
}

impl EncodingOverride {
    /// Parse one `COLUMN=ENCODING` token.
    pub fn parse(s: &str) -> Result<Self, String> {
        let (col, enc) = s.split_once('=').ok_or_else(|| {
            format!("`{s}` is not a COLUMN=ENCODING pair, e.g. qname=plain or ML=plain")
        })?;
        let column = col.trim().to_string();
        if column.is_empty() {
            return Err(format!("`{s}` names no column"));
        }
        Ok(Self { column, encoding: ColumnEncoding::parse(enc)? })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The delta encodings are the ones a user can get wrong, and the check has
    /// to agree with what parquet actually defines or the error surfaces at the
    /// first flush instead of at the command line.
    #[test]
    fn only_the_delta_encodings_restrict_the_type() {
        assert!(ColumnEncoding::Delta.supports(&DataType::Int64));
        assert!(ColumnEncoding::Delta.supports(&DataType::UInt64));
        assert!(!ColumnEncoding::Delta.supports(&DataType::Utf8));
        assert!(!ColumnEncoding::Delta.supports(&DataType::Boolean));

        assert!(ColumnEncoding::DeltaText.supports(&DataType::Utf8));
        assert!(!ColumnEncoding::DeltaText.supports(&DataType::Int64));

        for dt in [DataType::Utf8, DataType::Int64, DataType::Boolean, DataType::Float64] {
            assert!(ColumnEncoding::Dictionary.supports(&dt), "{dt:?}");
            assert!(ColumnEncoding::Plain.supports(&dt), "{dt:?}");
        }
    }

    /// Every variant round-trips through the name the CLI accepts, so the help
    /// text and the parser cannot list different sets.
    #[test]
    fn every_encoding_names_itself() {
        for (n, e) in ColumnEncoding::NAMES {
            assert_eq!(ColumnEncoding::parse(n).unwrap(), *e);
            assert_eq!(e.name(), *n);
        }
        assert_eq!(ColumnEncoding::parse("DELTA_TEXT").unwrap(), ColumnEncoding::DeltaText);
        assert!(ColumnEncoding::parse("rle").is_err());
    }

    #[test]
    fn an_override_needs_both_halves() {
        let o = EncodingOverride::parse("qname=plain").unwrap();
        assert_eq!(o.column, "qname");
        assert_eq!(o.encoding, ColumnEncoding::Plain);
        // Aux tags are named by the tag, and the tag is case-sensitive.
        assert_eq!(EncodingOverride::parse("ML=plain").unwrap().column, "ML");
        assert!(EncodingOverride::parse("qname").is_err());
        assert!(EncodingOverride::parse("=plain").is_err());
        assert!(EncodingOverride::parse("qname=zstd").is_err());
    }
}
