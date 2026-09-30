//! Take hit rows and write them out, through column builders that this sink
//! owns.
//!
//! Two output formats, one type. Parquet is the default: seekable, indexed,
//! and what a later query engine wants. Arrow IPC is for the streaming case --
//! a consumer reading batches as they are produced rather than a file read
//! back afterwards. The distinction that forces two writers rather than one is
//! narrow: `parquet::ArrowWriter` needs `Write + Seek` because it rewrites a
//! footer, while an IPC stream is append-only and needs only `Write`. That is
//! also why IPC can be pointed at a pipe and parquet cannot.
//!
//! An IPC stream ends with an explicit end-of-stream marker, so a truncated
//! stream is detectable as such by the reader rather than looking like a short
//! one. `close` is what writes it, for the same reason it writes the parquet
//! footer.
//!
//! `HitWriter<B>` is generic over the builder, and it does not hold a
//! `Box<dyn ArrowBuilder>`. This gives two things:
//!
//! - The schema is taken from the builder in `new` and can never disagree with
//!   the arrays `finish` hands back. There is no second place to state it.
//! - `len()` is called once per record on the hot path. Through a generic it
//!   monomorphises to a field read the optimiser can inline; through `dyn` it
//!   is a vtable call.
//!
//! The bound lives on the `impl` block, not on the struct. That is the usual
//! Rust convention: putting `<B: ArrowBuilder>` in the type definition forces
//! every other mention of the type to repeat it. Every field and every `where`
//! clause would carry it. That gains nothing, because the methods that need the
//! bound state it themselves.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use anyhow::{anyhow, Context, Result};
use arrow::datatypes::{DataType, Schema, SchemaRef};
use arrow::ipc::writer::StreamWriter;
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::metadata::KeyValue;
use parquet::file::properties::{EnabledStatistics, WriterProperties, WriterPropertiesBuilder};
use parquet::schema::types::ColumnPath;

use crate::arrow_builder::ArrowBuilder;
use crate::encoding::{ColumnEncoding, EncodingOverride};

/// Which of the two formats a scan writes.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum OutputFormat {
    #[default]
    Parquet,
    /// Arrow IPC stream. Append-only, so it can be written to a pipe, and each
    /// batch is readable as soon as it is written.
    Ipc,
}

/// The parquet-only half of the output settings.
///
/// Grouped because they are exactly the settings an IPC stream has no answer
/// for, and because passing four more positional arguments to
/// [`HitWriter::new`] -- two of them `Option`s and two of them enums -- is how
/// a caller ends up swapping a pair of them and finding out at run time.
///
/// [`Default`] is the fast profile, which is also the default: LZ4 and
/// chunk-level statistics. Writing a hits file is part of a pipeline, not the
/// end of one, and the pipeline is usually waiting on it.
#[derive(Clone, Debug)]
pub struct ParquetOpts {
    pub compression: Compression,
    /// Page-level statistics cost a comparison per value per column and buy
    /// page-index pruning; chunk-level keeps row-group pruning, which is what
    /// a row group of 50 000 rows is for in the first place.
    pub statistics: EnabledStatistics,
    pub max_row_group_rows: Option<usize>,
    /// `--column-encoding` overrides, applied after the builder's own table so
    /// the user's answer wins. Names are matched against the declared schema.
    pub encodings: Vec<EncodingOverride>,
    /// Key-value metadata every file carries: the run manifest (see
    /// [`crate::manifest`]). Written to the parquet footer, and for an IPC
    /// stream, which has no footer, into the schema's metadata. Applies to both
    /// formats despite the struct's name.
    pub metadata: Vec<(String, String)>,
}

impl Default for ParquetOpts {
    fn default() -> Self {
        Self {
            compression: Compression::LZ4_RAW,
            statistics: EnabledStatistics::Chunk,
            max_row_group_rows: None,
            encodings: Vec::new(),
            metadata: Vec::new(),
        }
    }
}

/// The parquet column path for one arrow field, and the type stored at it.
///
/// A list column has no leaf of its own: arrow-rs writes `List(item)` as a
/// group holding a repeated group `list` holding the item, so the path to the
/// values a policy is about is `name.list.item` and the type is the item's, not
/// the list's. Deriving both from the schema rather than spelling them out is
/// what keeps a policy on `capture_read_5p` pointed at the column it names --
/// a wrong path is not an error, it is a setting that silently does nothing.
fn leaf(schema: &Schema, name: &str) -> Result<(ColumnPath, DataType)> {
    let field = schema
        .field_with_name(name)
        .map_err(|_| anyhow!(
            "no column named `{name}` in this run's output; the columns are: {}",
            schema.fields().iter().map(|f| f.name().as_str()).collect::<Vec<_>>().join(", ")
        ))?;
    Ok(match field.data_type() {
        DataType::List(item) => (
            ColumnPath::new(vec![name.to_string(), "list".to_string(), item.name().clone()]),
            item.data_type().clone(),
        ),
        dt => (ColumnPath::new(vec![name.to_string()]), dt.clone()),
    })
}

/// Set one column's policy, refusing an encoding parquet does not define over
/// that column's type.
///
/// Checked here, before the file is created, because the alternative is a
/// writer that accepts the setting and fails on the first flush -- minutes into
/// a scan, with the reason a page-level error a long way from the flag that
/// caused it.
fn apply_encoding(
    props: WriterPropertiesBuilder,
    schema: &Schema,
    name: &str,
    enc: ColumnEncoding,
) -> Result<WriterPropertiesBuilder> {
    let (path, dt) = leaf(schema, name)?;
    if !enc.supports(&dt) {
        return Err(anyhow!(
            "column `{name}` is {dt:?}, which cannot be stored as `{}`; \
             `delta` needs an integer column and `delta-text` a text one",
            enc.name()
        ));
    }
    Ok(enc.apply(props, path))
}

enum Sink {
    Parquet(Box<ArrowWriter<BufWriter<File>>>),
    Ipc(Box<StreamWriter<BufWriter<File>>>),
}

pub struct HitWriter<B> {
    builder: B,
    schema: SchemaRef,
    batch_rows: usize,
    writer: Sink,
    /// Rows handed to the file so far.
    rows: u64,
}

impl<B: ArrowBuilder> HitWriter<B> {
    /// Take ownership of `builder` and open `path` with its schema.
    ///
    /// Everything in `opts` is parquet-only and ignored for IPC, which has
    /// none of those concepts.
    ///
    /// The per-column policies are resolved before the file exists, so a
    /// mistyped `--column-encoding` or one naming a column this run does not
    /// write fails here rather than part way through the scan.
    pub fn new(
        path: &Path,
        builder: B,
        batch_rows: usize,
        format: OutputFormat,
        opts: &ParquetOpts,
    ) -> Result<Self> {
        let mut schema = builder.schema();
        let writer = match format {
            OutputFormat::Parquet => {
                // Plain footer keys, which any parquet reader sees (DuckDB's
                // `parquet_kv_metadata`); `ArrowWriter` puts the schema's own
                // metadata only inside the serialized `ARROW:schema`. The
                // schema's keys first, then the run's, once each.
                let mut kv: Vec<KeyValue> = Vec::new();
                let mut schema_keys: Vec<(&String, &String)> = schema.metadata().iter().collect();
                schema_keys.sort();
                for (k, v) in schema_keys.into_iter().chain(opts.metadata.iter().map(|(k, v)| (k, v))) {
                    if !kv.iter().any(|e| &e.key == k) {
                        kv.push(KeyValue::new(k.clone(), v.clone()));
                    }
                }
                let mut props = WriterProperties::builder()
                    .set_compression(opts.compression)
                    .set_statistics_enabled(opts.statistics)
                    .set_max_row_group_row_count(opts.max_row_group_rows)
                    .set_key_value_metadata((!kv.is_empty()).then_some(kv));
                // The builder's own table first, the user's overrides second:
                // later settings win per column path, so this is what makes
                // `--column-encoding` an override rather than a suggestion.
                for (name, enc) in builder.column_encodings() {
                    props = apply_encoding(props, &schema, &name, enc)?;
                }
                for o in &opts.encodings {
                    props = apply_encoding(props, &schema, &o.column, o.encoding)
                        .with_context(|| format!("--column-encoding {}={}", o.column, o.encoding.name()))?;
                }
                let file = File::create(path)
                    .with_context(|| format!("creating {}", path.display()))?;
                let buf_writer = BufWriter::with_capacity(4 * 1024 * 1024, file);
                Sink::Parquet(Box::new(ArrowWriter::try_new(
                    buf_writer,
                    schema.clone(),
                    Some(props.build()),
                )?))
            }
            OutputFormat::Ipc => {
                if !opts.metadata.is_empty() {
                    let mut meta = schema.metadata().clone();
                    meta.extend(opts.metadata.iter().cloned());
                    schema = std::sync::Arc::new(schema.as_ref().clone().with_metadata(meta));
                }
                let file = File::create(path)
                    .with_context(|| format!("creating {}", path.display()))?;
                let buf_writer = BufWriter::with_capacity(4 * 1024 * 1024, file);
                Sink::Ipc(Box::new(
                    StreamWriter::try_new(buf_writer, &schema)
                        .with_context(|| format!("starting an IPC stream at {}", path.display()))?,
                ))
            }
        };
        Ok(Self { builder, schema, batch_rows, writer, rows: 0 })
    }

    /// The builder, to append a row through.
    ///
    /// The returned borrow can be held across a whole record's worth of
    /// pushes; it ends at its last use, which is what lets `maybe_flush` be
    /// called straight afterwards.
    #[inline]
    pub fn builder_mut(&mut self) -> &mut B {
        &mut self.builder
    }




    #[inline]
    #[cfg_attr(not(test), allow(dead_code))] // read by the batching tests, which pin the flush boundary
    pub fn buffered(&self) -> usize {
        self.builder.len()
    }

    pub fn maybe_flush(&mut self) -> Result<()> {
        if self.builder.len() >= self.batch_rows {
            self.flush()
        } else {
            Ok(())
        }
    }

    /// Drain the builder into one record batch.
    ///
    /// A no-op on an empty builder: `RecordBatch::try_new` would accept the
    /// zero-row batch, but writing it costs a page per column for nothing.
    pub fn flush(&mut self) -> Result<()> {
        if self.builder.len() == 0 {
            return Ok(());
        }
        // Disjoint fields, so borrowck is happy with the mutable builder borrow
        // alongside the shared schema one.
        let cols = self.builder.finish();
        debug_assert_eq!(
            cols.len(),
            self.schema.fields().len(),
            "builder produced {} arrays for a {}-column schema",
            cols.len(),
            self.schema.fields().len()
        );
        let batch = RecordBatch::try_new(self.schema.clone(), cols)?;
        self.rows += batch.num_rows() as u64;
        match &mut self.writer {
            Sink::Parquet(w) => w.write(&batch)?,
            Sink::Ipc(w) => w.write(&batch)?,
        }
        Ok(())
    }

    /// Flush the tail and write the footer, returning the rows written.
    ///
    /// Consuming `self` makes it a compile error to keep using the writer, and
    /// skipping this call leaves a parquet file with no footer — unreadable,
    /// not merely short — or an IPC stream with no end-of-stream marker, which
    /// a reader reports as truncation. It cannot be done in `Drop`, because
    /// `Drop::drop` has nowhere to report an I/O failure.
    pub fn close(mut self) -> Result<u64> {
        self.flush()?;
        match self.writer {
            Sink::Parquet(w) => {
                w.close()?;
            }
            Sink::Ipc(mut w) => {
                w.finish()?;
                // StreamWriter::finish writes the marker but leaves the inner
                // BufWriter unflushed; dropping it would swallow an I/O error.
                w.into_inner()?;
            }
        }
        Ok(self.rows)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::hits::HitBuilder;
    use crate::record_field::RecordField;
    use crate::test_support::{built, header, temp_path};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use parquet::basic::Encoding;
    use rust_htslib::bam::record::Cigar;

    fn writer(tag: &str, batch_rows: usize, row_group: Option<usize>) -> HitWriter<HitBuilder> {
        let b = HitBuilder::new(RecordField::CORE, batch_rows, true, 0);
        let opts = ParquetOpts { max_row_group_rows: row_group, ..Default::default() };
        HitWriter::new(
            &temp_path(&format!("pw_{tag}"), "parquet"),
            b,
            batch_rows,
            OutputFormat::Parquet,
            &opts,
        )
        .unwrap()
    }

    /// Rows written to the file, read back.
    fn rows(tag: &str) -> usize {
        let file = std::fs::File::open(temp_path(&format!("pw_{tag}"), "parquet")).unwrap();
        let reader = ParquetRecordBatchReaderBuilder::try_new(file).unwrap().build().unwrap();
        reader.map(|b| b.unwrap().num_rows()).sum()
    }

    /// `maybe_flush` drains only once the batch is full, so a record's rows
    /// never straddle two batches.
    #[test]
    fn flushing_waits_for_a_full_batch() {
        let mut w = writer("partial", 4, None);
        for _ in 0..3 {
            w.builder_mut().push_hitless();
        }
        assert_eq!(w.buffered(), 3);
        w.maybe_flush().unwrap();
        assert_eq!(w.buffered(), 3, "three rows is not a full batch of four");

        w.builder_mut().push_hitless();
        w.maybe_flush().unwrap();
        assert_eq!(w.buffered(), 0, "a full batch drains");
    }

    /// The tail is written by `close`, so rows buffered when the scan ends are
    /// not lost.
    #[test]
    fn close_writes_the_partial_tail() {
        let mut w = writer("tail", 100, None);
        for _ in 0..7 {
            w.builder_mut().push_hitless();
        }
        assert_eq!(w.buffered(), 7, "well short of the batch size");
        w.close().unwrap();
        assert_eq!(rows("tail"), 7, "the tail must survive");
    }

    /// Flushing an empty builder writes nothing rather than a zero-row batch,
    /// which would cost a page per column for no rows.
    #[test]
    fn an_empty_flush_writes_nothing() {
        let mut w = writer("empty", 4, None);
        w.flush().unwrap();
        w.flush().unwrap();
        w.close().unwrap();
        assert_eq!(rows("empty"), 0);
        // And the file is still a readable parquet file, not a truncated one.
        let file = std::fs::File::open(temp_path("pw_empty", "parquet")).unwrap();
        assert!(ParquetRecordBatchReaderBuilder::try_new(file).is_ok());
    }

    /// Every row pushed reaches the file across several flushes.
    #[test]
    fn all_rows_survive_repeated_flushes() {
        let mut w = writer("many", 4, None);
        for _ in 0..10 {
            w.builder_mut().push_hitless();
            w.maybe_flush().unwrap();
        }
        w.close().unwrap();
        assert_eq!(rows("many"), 10);
    }

    /// The schema on the file comes from the builder, so the two cannot
    /// disagree -- there is no second place to state it.
    #[test]
    fn the_file_schema_is_the_builders() {
        let mut w = writer("schema", 4, None);
        let declared: Vec<String> =
            w.builder_mut().schema().fields().iter().map(|f| f.name().clone()).collect();
        w.builder_mut().push_hitless();
        w.close().unwrap();

        let file = std::fs::File::open(temp_path("pw_schema", "parquet")).unwrap();
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
        let on_disk: Vec<String> =
            builder.schema().fields().iter().map(|f| f.name().clone()).collect();
        assert_eq!(declared, on_disk);
    }

    /// The encodings the first row group actually used, by dotted column path.
    ///
    /// Read back off the file rather than off the properties, because the two
    /// can disagree in exactly the way that matters: a policy set on a path the
    /// writer does not have is not an error, it is a setting that does nothing.
    fn encodings_on_disk(path: &Path) -> Vec<(String, Vec<Encoding>)> {
        let file = std::fs::File::open(path).unwrap();
        let md = ParquetRecordBatchReaderBuilder::try_new(file).unwrap().metadata().clone();
        md.row_group(0)
            .columns()
            .iter()
            .map(|c| (c.column_path().string(), c.encodings().collect::<Vec<_>>()))
            .collect()
    }

    /// Write one real record's worth of rows and hand back what the file says.
    fn written_with(tag: &str, flat: bool, opts: &ParquetOpts) -> Vec<(String, Vec<Encoding>)> {
        let path = temp_path(&format!("pw_{tag}"), "parquet");
        let head = header("chr1", 1000);
        let rec = built(b"A00123:45:HXYZ:1:1101:1000:2000", b"ACGTACGTAC", &[Cigar::Match(10)], 100, 0);
        let b = HitBuilder::new(RecordField::ALL, 8, flat, 0);
        let mut w =
            HitWriter::new(&path, b, 8, OutputFormat::Parquet, opts).unwrap();
        for i in 0..4 {
            w.builder_mut().begin_record(i, &rec, &head, Some(crate::tags::Library::Directional.call(&rec)));
            w.builder_mut().push_hitless();
        }
        w.close().unwrap();
        encodings_on_disk(&path)
    }

    fn used(cols: &[(String, Vec<Encoding>)], path: &str) -> Vec<Encoding> {
        cols.iter()
            .find(|(p, _)| p == path)
            .unwrap_or_else(|| panic!("no column at `{path}`; the file has {:?}",
                cols.iter().map(|(p, _)| p.as_str()).collect::<Vec<_>>()))
            .1
            .clone()
    }

    /// Each column's declared policy is the encoding the file actually uses.
    ///
    /// The point of reading it back is the negative case: parquet takes a
    /// setting on any column path you like and silently ignores one that names
    /// nothing, so "the properties were set" proves nothing about the file.
    #[test]
    fn the_declared_policies_are_what_the_file_uses() {
        let cols = written_with("enc_flat", true, &ParquetOpts::default());

        // Prefix-delta for read names, which repeat in runs and share prefixes.
        assert!(used(&cols, "qname").contains(&Encoding::DELTA_BYTE_ARRAY));
        assert!(!used(&cols, "qname").contains(&Encoding::RLE_DICTIONARY));
        // Delta for the monotonic id and the coordinates.
        for c in ["record_id", "read_5p", "read_3p", "refr_pos", "pos"] {
            assert!(
                used(&cols, c).contains(&Encoding::DELTA_BINARY_PACKED),
                "{c} should be delta-packed, got {:?}",
                used(&cols, c)
            );
        }
        // Dictionary for the low-cardinality ones.
        for c in ["ref_name", "strand", "name"] {
            assert!(
                used(&cols, c).contains(&Encoding::RLE_DICTIONARY),
                "{c} should be dictionary-encoded, got {:?}",
                used(&cols, c)
            );
        }
        // Plain for the per-record strings, which have nothing to share.
        assert!(!used(&cols, "seq_ascii").contains(&Encoding::RLE_DICTIONARY));
    }

    /// A list column's policy is about the leaf inside the list, and the path to
    /// it is `name.list.item` — spelled by the writer, from the schema. Getting
    /// it wrong would leave the capture columns on the writer's defaults with
    /// nothing to say so.
    #[test]
    fn a_list_columns_policy_reaches_the_item_inside_it() {
        let cols = written_with("enc_list", false, &ParquetOpts::default());
        assert!(
            used(&cols, "capture_read_5p.list.item").contains(&Encoding::DELTA_BINARY_PACKED),
            "the leaf inside the list is what carries the policy"
        );
        assert!(used(&cols, "capture_read.list.item").contains(&Encoding::RLE_DICTIONARY));
    }

    /// `--column-encoding` is applied after the built-in table, so it wins.
    /// This is the `ML:B`/`MM:Z` case: a text tag that is really a per-record
    /// array, which the by-kind default would dictionary-encode.
    #[test]
    fn an_override_beats_the_builtin_policy() {
        let opts = ParquetOpts {
            encodings: vec![EncodingOverride {
                column: "qname".to_string(),
                encoding: ColumnEncoding::Plain,
            }],
            ..Default::default()
        };
        let cols = written_with("enc_over", true, &opts);
        let qname = used(&cols, "qname");
        assert!(qname.contains(&Encoding::PLAIN));
        assert!(!qname.contains(&Encoding::DELTA_BYTE_ARRAY), "the override wins");
    }

    /// Both ways an override can be wrong are refused before the file is
    /// created, rather than at the first flush minutes into a scan.
    #[test]
    fn a_bad_override_is_refused_at_open() {
        let bad = |column: &str, encoding| ParquetOpts {
            encodings: vec![EncodingOverride { column: column.to_string(), encoding }],
            ..Default::default()
        };
        let open = |tag: &str, opts: &ParquetOpts| {
            HitWriter::new(
                &temp_path(&format!("pw_{tag}"), "parquet"),
                HitBuilder::new(RecordField::CORE, 4, true, 0),
                4,
                OutputFormat::Parquet,
                opts,
            )
            .map(|_| ())
        };

        // `{:#}` because the flag it came from is the outer context and the
        // reason is the cause underneath it; both belong in the message.
        let e = open("enc_typo", &bad("qnaem", ColumnEncoding::Plain)).unwrap_err();
        let e = format!("{e:#}");
        assert!(e.contains("no column named `qnaem`"), "{e}");
        assert!(e.contains("--column-encoding qnaem=plain"), "{e}");

        // Delta is defined over integers; `strand` is text.
        let e = open("enc_type", &bad("strand", ColumnEncoding::Delta)).unwrap_err();
        let e = format!("{e:#}");
        assert!(e.contains("strand"), "{e}");
        assert!(e.contains("--column-encoding strand=delta"), "{e}");
    }

    /// The codec is per file and reaches every column chunk. `Default` is the
    /// fast profile, so this is also what a plain parquet run writes.
    #[test]
    fn the_default_profile_writes_lz4() {
        let path = temp_path("pw_codec", "parquet");
        let mut w = HitWriter::new(
            &path,
            HitBuilder::new(RecordField::CORE, 4, true, 0),
            4,
            OutputFormat::Parquet,
            &ParquetOpts::default(),
        )
        .unwrap();
        w.builder_mut().push_hitless();
        w.close().unwrap();

        let file = std::fs::File::open(&path).unwrap();
        let md = ParquetRecordBatchReaderBuilder::try_new(file).unwrap().metadata().clone();
        for c in md.row_group(0).columns() {
            assert_eq!(c.compression(), Compression::LZ4_RAW, "{}", c.column_path().string());
        }
    }

    /// The run's metadata, and the schema's own `format_version` and
    /// `coordinate_base`, are plain parquet footer keys that any reader sees, and
    /// in an IPC stream they are schema metadata. `close` reports the rows.
    #[test]
    fn metadata_reaches_the_footer_and_the_stream_schema() {
        let opts = ParquetOpts {
            metadata: vec![("alnbase_manifest".to_string(), "{\"run_id\":\"x\"}".to_string())],
            ..Default::default()
        };
        let path = temp_path("pw_meta", "parquet");
        let mut w =
            HitWriter::new(&path, HitBuilder::new(RecordField::CORE, 4, true, 0), 4, OutputFormat::Parquet, &opts)
                .unwrap();
        for _ in 0..5 {
            w.builder_mut().push_hitless();
        }
        assert_eq!(w.close().unwrap(), 5);
        let file = std::fs::File::open(&path).unwrap();
        let md = ParquetRecordBatchReaderBuilder::try_new(file).unwrap().metadata().clone();
        let kv: std::collections::HashMap<String, Option<String>> = md
            .file_metadata()
            .key_value_metadata()
            .unwrap()
            .iter()
            .map(|e| (e.key.clone(), e.value.clone()))
            .collect();
        let keys = md.file_metadata().key_value_metadata().unwrap().len();
        assert_eq!(keys, kv.len(), "no key is written twice");
        assert_eq!(kv["format_version"].as_deref(), Some(crate::hits::FORMAT_VERSION));
        assert_eq!(kv["coordinate_base"].as_deref(), Some("0"));
        assert_eq!(kv["alnbase_manifest"].as_deref(), Some("{\"run_id\":\"x\"}"));

        let path = temp_path("pw_meta", "arrows");
        let w = HitWriter::new(&path, HitBuilder::new(RecordField::CORE, 4, true, 0), 4, OutputFormat::Ipc, &opts)
            .unwrap();
        assert_eq!(w.close().unwrap(), 0);
        let r = arrow::ipc::reader::StreamReader::try_new(std::fs::File::open(&path).unwrap(), None).unwrap();
        let meta = r.schema().metadata().clone();
        assert_eq!(meta["coordinate_base"], "0");
        assert_eq!(meta["alnbase_manifest"], "{\"run_id\":\"x\"}");
    }

    /// An IPC stream carries the same schema and rows as the parquet form, and
    /// `close` is what terminates it: a reader that reaches the end without the
    /// marker reports truncation rather than a short file.
    #[test]
    fn the_ipc_stream_round_trips() {
        use arrow::ipc::reader::StreamReader;

        let path = temp_path("pw_ipc", "arrows");
        let mut w = HitWriter::new(
            &path,
            HitBuilder::new(RecordField::CORE, 2, true, 3),
            2,
            OutputFormat::Ipc,
            &ParquetOpts::default(),
        )
        .unwrap();
        let want = w.builder_mut().schema();
        for _ in 0..5 {
            w.builder_mut().push_hitless();
            w.maybe_flush().unwrap();
        }
        w.close().unwrap();

        let file = std::fs::File::open(&path).unwrap();
        let reader = StreamReader::try_new(file, None).unwrap();
        assert_eq!(reader.schema(), want, "the stream states the builder's schema");
        let n: usize = reader.map(|b| b.unwrap().num_rows()).sum();
        assert_eq!(n, 5);
    }
}