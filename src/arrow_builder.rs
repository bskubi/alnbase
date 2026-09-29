//! Build one Arrow column. This is the interface that
//! [`crate::hit_writer::HitWriter`] drives.

use arrow::array::ArrayRef;
use arrow::datatypes::SchemaRef;

use crate::encoding::ColumnEncoding;

/// A set of Arrow column builders that can be drained into one record batch.
///
/// Implementors must keep three things in step, and nothing here can check
/// them for you:
///
/// - `finish` returns one array per field in `schema`, in the same order;
/// - all of those arrays have the same length;
/// - `finish` resets the row count, so `len` reads zero straight afterwards.
///
/// A field-count mismatch trips a `debug_assert` in the writer, but a swap of
/// two same-typed columns will not, so round-trip a record in a test. See
/// `hits::tests::the_declared_schema_matches_the_arrays_produced` for the
/// first two, and `hit_writer::tests::the_file_schema_is_the_builders`
/// for what reaches the file.
pub trait ArrowBuilder {
    fn finish(&mut self) -> Vec<ArrayRef>;

    /// Rows appended since the last `finish`.
    fn len(&self) -> usize;

    /// Constant for the lifetime of the builder.
    fn schema(&self) -> SchemaRef;

    /// How each column would like to be stored, by schema name.
    ///
    /// Advisory and partial: a name the schema does not have is an error the
    /// writer reports, but a column left out of the list gets the
    /// writer's default. Parquet-only — an IPC stream has no such concept, and
    /// the writer ignores this entirely when it is writing one.
    ///
    /// It belongs on the builder for the same reason `schema` does. The builder
    /// is the one place that knows what the columns are; a table of per-column
    /// policies kept anywhere else would be a second list to hold in step, and
    /// nothing would notice when it drifted.
    fn column_encodings(&self) -> Vec<(String, ColumnEncoding)> {
        Vec::new()
    }

}