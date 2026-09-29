//! Parse the arguments of the `query` subcommand. This module does nothing else.
//!
//! The overlap tool has its own file. The two files share only the top-level
//! `Cli` enum, which lives in `main.rs`. They are independent tools that ship in
//! one binary.
//!
//! This file defines no vocabulary of alnbase's own. The record columns come
//! from [`crate::record_field`], which also builds them. The walk policy enums
//! come from [`crate::alignment`], which also acts on them. The query syntax
//! comes from [`crate::query_toml`], [`crate::lower`] and [`crate::dsl`]. Each
//! of these types has its parser at the place where the type is defined. The
//! values that the CLI accepts and the values that the code handles are
//! therefore always the same set.
//!
//! The output options are the exception, for three different reasons.
//! [`CompressionArg`] and [`StatisticsArg`] name values that belong to the
//! `parquet` crate. You cannot derive clap's `ValueEnum` for a type from another
//! crate, so this file writes the accepted names once and maps them across.
//! [`OutputFormatArg`] copies an alnbase type, [`hit_writer::OutputFormat`].
//! That copy is a choice and not a necessity: it keeps clap out of the writer.
//! [`OutputProfile`] copies nothing at all. `fast` and `small` are only a
//! convenience for the person who types the command.
//!
//! [`ParquetArgs::parquet_opts`] resolves all three of the parquet options in
//! one place, and not next to each enum. The profile and the single flags must
//! be resolved against each other. The profile selects a pair of codec and
//! statistics, and each single flag that the user gives replaces its half of
//! that pair. If you split the mapping up, this rule has no single place to
//! live.
//!
//! This file also keeps the text of each query file, next to its parsed contents
//! ([`QuerySource`]). A run copies its query files into the header of the BAM
//! that it writes, so that a tagged BAM carries the definitions of its own tags.
//! See [`crate::header_query`].
//!
//! # Why queries are resolved after parsing and not in a `value_parser`
//!
//! Queries come only from `--query-file` options, in TOML. A run checks the
//! files of the run against each other, for aliases, query names and tags. A
//! clap `value_parser` sees one value at a time and cannot do this check.
//! [`QueryArgs::resolve`] runs once, after clap collects everything.

use std::path::PathBuf;

use clap::{Args, ValueEnum};
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::EnabledStatistics;

use crate::alignment::Insertions;
use crate::dsl::{self, Aliases, QuerySpec};
use crate::encoding::EncodingOverride;
use crate::hit_writer::{self, ParquetOpts};
use crate::partition;
use crate::record_field::{FieldSpec, RecordField};
use crate::scanner::WalkConfig;
use crate::header_query::QuerySource;
use crate::query_toml;
use crate::tags::TagConfig;
use crate::strand_rule::StrandRule;

/// Which format the query subcommand writes its hits in.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum, Default)]
pub enum OutputFormatArg {
    /// Seekable and indexed, for reading back later. The default.
    #[default]
    Parquet,
    /// Arrow IPC stream. Append-only, so it can be written to a FIFO and read
    /// batch by batch as the scan proceeds, without the hits landing on disk.
    Ipc,
}

impl From<OutputFormatArg> for hit_writer::OutputFormat {
    fn from(v: OutputFormatArg) -> Self {
        match v {
            OutputFormatArg::Parquet => hit_writer::OutputFormat::Parquet,
            OutputFormatArg::Ipc => hit_writer::OutputFormat::Ipc,
        }
    }
}

/// How hard to work at making the output small.
///
/// There are two settings, and not a number to tune. The choice is nearly
/// always which of the two resources is scarce: processor time, or disk space.
/// The granular flags below override either setting.
///
/// On 400k 100bp reads with every column selected, `fast` wrote 30.7 MB and
/// `small` 19.2 MB, `small` taking about 15% longer end to end -- of which the
/// write is only a part, since both runs walk the same alignments. The gap
/// widens with thread count, where each worker writes its own shard.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum, Default)]
pub enum OutputProfile {
    /// LZ4 and chunk-level statistics. The default: a hits file is usually an
    /// intermediate, and the next step is waiting on it.
    #[default]
    Fast,
    /// Zstd at `--compression-level` and page-level statistics. For a file that
    /// will be kept, shipped, or read many more times than it is written.
    Small,
}

/// Which block codec compresses the pages.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum CompressionArg {
    /// No compression. This is usually the *slowest* option, and not the
    /// fastest. The hits file is roughly 17x larger, and writing those extra
    /// bytes costs more than LZ4 costs to avoid them. Use it only when the
    /// output goes to a pipe or to a tmpfs, where no bytes reach a disk.
    None,
    /// Very fast, decent ratio. What `--output-profile fast` uses.
    Lz4,
    /// Fast, and larger than LZ4 on this data rather than smaller.
    Snappy,
    /// Slowest of the four and much the smallest. Honours
    /// `--compression-level`.
    Zstd,
}

/// How much statistics the writer computes.
///
/// Measurably the least important of these knobs -- page-level statistics cost
/// a few percent of a run, not tens. It is here because it is free to offer and
/// because a reader that never prunes has no use for them at all.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum StatisticsArg {
    /// None at all. Nothing can be pruned without reading it.
    None,
    /// Min/max per column per row group, which is what a reader needs to skip a
    /// row group. What `--output-profile fast` uses.
    Chunk,
    /// Per page as well, for finer pruning through the page index, at a
    /// comparison per value per column while writing.
    Page,
}

/// Appended to every query parse failure.
///
/// `--help` for `--partition-by`. Built from the field tables, so the list of
/// what it accepts cannot drift from what `--field` accepts.
fn partition_help() -> String {
    format!(
        "Record columns whose combined value decides which file a record's rows go to.\n\n\
         Comma-separated and/or repeated: record column names, groups and aux tags.\n\
         Defaults to `qname` if not given.\n\n\
         Within a run, records sharing a key always land in the same file, whatever\n\
         --threads and --shards-per-worker are set to. The file is chosen by\n\n    \
         slot   = hash(key) % (threads * shards-per-worker)\n    \
         worker = slot / shards-per-worker\n    \
         shard  = slot % shards-per-worker\n\n\
         and written beside OUT as {{stem}}_{{worker}}_{{shard}}.{{ext}}: calls.parquet\n\
         becomes calls_0_0.parquet, calls_0_1.parquet, ... The hash is the same on\n\
         every run and machine, but a key is only guaranteed the same file in another\n\
         run when --threads and --shards-per-worker are both unchanged.\n\n\
         The default exists for a reason: routing by qname is what puts every alignment\n\
         of a fragment in one file, which is what lets mate-overlap resolution see both\n\
         sides. Any other key gives that up unless it is constant across a fragment --\n\
         a read group is, a reference name is not -- so a run that resolves mates per\n\
         file should keep the default. A run whose output is aggregated per site does\n\
         not care, and gets a dataset already partitioned by something it can use.\n\n\
         Examples:\n  \
         --partition-by ref_name          one contig per file\n  \
         --partition-by RG:Z              one read group per file\n  \
         --partition-by qname --shards-per-worker 8\n                                   \
         fragments kept together, 8 files per worker\n\n\
         Groups: {}\n",
        crate::record_field::GROUPS.iter().map(|(n, _)| *n).collect::<Vec<_>>().join(", "),
    )
}

#[derive(Args, Debug)]
pub struct QueryArgs {
    /// TOML file of queries and tags. Repeatable.
    ///
    /// The format and a worked example live with the parser, in
    /// [`crate::query_toml`], so `--help` cannot drift from what it accepts.
    #[arg(
        long = "query-file",
        value_name = "FILE",
        help = "TOML file of queries and tags; see --help for the format",
        long_help = query_toml::help()
    )]
    pub query_file: Vec<PathBuf>,

    /// Render each query column-aligned and exit without scanning.
    ///
    /// Shows the expression as a tree, the query in words, every operand
    /// stacked column by column, and the normal form the predicate compiled to.
    /// Worth running on any query before committing a full scan to it.
    #[arg(long)]
    pub explain: bool,

    /// Trace the queries over a synthetic record and exit.
    ///
    /// Takes an aligned READ@REFR column pair, e.g. --trace 'TTTCGTTT@TTTCGTTT'.
    /// Shows how far each pattern matched at every column, which query fired
    /// where, and the reasoning behind each result. Needs no BAM.
    #[arg(long, value_name = "READ@REFR", conflicts_with = "trace_records")]
    pub trace: Option<String>,

    /// Trace the queries over the first N records of the BAM and exit.
    ///
    /// The same report as --trace, on real alignment columns. Start with 1:
    /// the output is a screenful per record.
    #[arg(long, value_name = "N")]
    pub trace_records: Option<usize>,

    /// Show the per-column progress grid in a trace.
    ///
    /// This is the part that says where a match broke, so it is on unless a
    /// record is too wide to read.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub trace_grid: bool,

    /// Print every code a pattern can contain, and exit.
    ///
    /// Includes which characters are still free to use as an alias, which is
    /// otherwise only discoverable by guessing one and reading the error.
    #[arg(long)]
    pub list_codes: bool,

    /// What to do with inserted bases.
    ///
    /// `skip` (the default) emits no column for them, which keeps the two
    /// reference bases flanking an insertion adjacent so a pattern can span it.
    ///
    /// `emit` pairs them against a gap on the reference side. Needed to write
    /// patterns *about* insertions, at the cost of breaking patterns that span
    /// one.
    #[arg(long, value_enum, value_name = "POLICY", default_value = "skip")]
    pub insertions: Insertions,

    /// Flank columns to show past each end of the aligned read.
    ///
    /// Each carries the real reference base (or a pad past the contig end), so a
    /// pattern can see reference context beyond the read. On the read side, the
    /// columns nearest a soft-clipped end are clips (`:`), one per clipped base,
    /// with that base's offsets; the rest are pads (`_`), which have no read
    /// offset: off_5p and off_3p are null on a hit anchored on one.
    ///
    /// Defaults to k, the span of the widest query in the run: far enough past a
    /// read end for every query, the widest included, to fit a whole window
    /// into the pads.
    ///
    /// A query fires on every window of these columns it matches, including
    /// windows made only of pads. With the default, adding a wider query adds
    /// pads, so a query whose read row accepts pads (`~`, `_`) can gain hits.
    /// Set this explicitly to fix the pads whatever else is queried, or exclude
    /// all-pad windows in the query with `and not` over a pattern of pads.
    #[arg(long, value_name = "N")]
    pub end_context: Option<usize>,

    /// Reference bases to show at each end of an intron.
    ///
    /// A `CIGAR N` is not emitted base by base -- that would be tens of
    /// thousands of columns the read never covered. It emits this many
    /// reference bases at each end against `,` on the read side, and one `,@,`
    /// marker for everything between them.
    ///
    /// Defaults to k, the span of the widest query, so that a junction edge
    /// behaves like the end of a read: a pattern reaches as far into an intron
    /// as it reaches past the end of the read, and a CG whose C the read covers
    /// still matches when the G is the intron's first base.
    ///
    /// Set it explicitly when matching a splice motif. The default moves with
    /// the widest query in the run, which slides the intron's edges away from
    /// the marker; a fixed value pins the layout, so `--splice-context 2` makes
    /// `~~,@GT,` a donor site whatever else is being queried. As with
    /// `--end-context`, a query that can match windows made only of intron
    /// columns can gain hits when a wider query raises the default.
    #[arg(long, value_name = "N")]
    pub splice_context: Option<usize>,

    /// Stop when more than this fraction of read bases disagree with the reference.
    ///
    /// Counts observed columns where the read and the reference each show one
    /// unambiguous base, not counting a reference C read as T (conversion). A
    /// correct alignment disagrees at a few percent at most; a wrong reference, a
    /// BAM aligned to another assembly, or shifted coordinates disagree at about
    /// 75%. The run is checked once it has compared 100,000 bases (all threads
    /// together) and again when it ends. The run summary
    /// reports the rate. 1 turns the check off.
    #[arg(long, value_name = "FRACTION", default_value_t = 0.25, value_parser = fraction)]
    pub max_discordance: f64,

    /// Refuse records on a contig whose @SQ line in the BAM header has no M5.
    ///
    /// The reference is always checked against the header by contig name and
    /// length, and by MD5 wherever the header gives an M5. A masked analysis set
    /// or a patched assembly has the same names and lengths, so without an M5
    /// nothing tells it apart. With this flag, a record on a contig with no M5
    /// stops the run (whatever --permissive says); contigs no record is on are
    /// not checked. `samtools dict` prints M5 values for a FASTA.
    #[arg(long)]
    pub require_m5: bool,

    /// Tagging threads.
    #[arg(short = '@', long, default_value_t = 1)]
    pub threads: usize,

    /// BGZF decompression threads for reading the BAM.
    ///
    /// Separate from `--threads` because they are not the same resource. One
    /// reader feeds every worker, so decompression is serial with respect to
    /// the scan and becomes the ceiling long before the workers are busy:
    /// a single BGZF thread tops out around 150-250 MB/s regardless of how
    /// many workers are waiting on the channel. These threads also compete
    /// with the workers for cores, so more is not always better -- raise it
    /// until reads/s stops improving.
    ///
    /// 0 means min(4, --threads).
    #[arg(long, default_value_t = 0)]
    pub reader_threads: usize,

    /// BGZF compression threads for writing the output BAM.
    ///
    /// The mirror of `--reader-threads`: one writer takes every tagged record
    /// in input order, so compressing the output is serial with respect to the
    /// taggers unless it has threads of its own.
    ///
    /// 0 means min(4, --threads).
    #[arg(long, default_value_t = 0)]
    pub writer_threads: usize,

    /// Input BAM, or '-' for standard input.
    ///
    /// A single forward pass, so a pipe is as good as a file:
    ///
    ///     samtools collate -O in.bam | alnbase overlap - - | alnbase query --query-file q.toml - ref.aref out.bam
    ///
    /// Not required with --explain, --trace or --list-codes, which never open one.
    pub bam: Option<PathBuf>,
    pub refr: Option<PathBuf>,
    /// Output BAM, or '-' for standard output. Records come out in the order
    /// they went in.
    #[arg(value_name = "OUT")]
    pub out: Option<PathBuf>,

    /// Write parquet hit rows instead of a tagged BAM.
    ///
    /// OUT is then a path prefix: rows are sharded by `--partition-by` into
    /// `{stem}_{worker}_{shard}.{ext}` beside it, so `calls.parquet` becomes
    /// `calls_0_0.parquet`, ... The options under "Parquet output" apply. The
    /// query files' tags are not used.
    #[arg(long = "parquet")]
    pub parquet_output: bool,

    #[command(flatten, next_help_heading = "Parquet output (with --parquet)")]
    pub parquet: ParquetArgs,

    /// Warn and count instead of stopping on data the run cannot handle.
    ///
    /// Without it, a record on a contig the reference does not have stops the
    /// run. With it, such records are skipped and counted. Nothing about what
    /// is written changes either way.
    #[arg(long)]
    pub permissive: bool,

    /// Replace a tag an input record already has, instead of stopping.
    ///
    /// Off by default because a record that already carries, say, XM was
    /// probably tagged by another tool, and quietly replacing its calls with
    /// this run's is the kind of mistake that surfaces much later.
    #[arg(long)]
    pub overwrite_tags: bool,
}

/// The parquet output's options: shared by `query --parquet` and `extract`.
#[derive(Args, Debug)]
pub struct ParquetArgs {
    /// Record columns to emit. Names, groups and aux tags all come from
    /// `RecordField`; see `--help` for the generated list.
    #[arg(
        short = 'f',
        long = "field",
        value_name = "FIELD",
        value_parser = RecordField::parse,
        value_delimiter = ',',
        long_help = RecordField::long_help(),
        help = "Record columns to emit, replacing the default; see --help [default: core]"
    )]
    pub field: Vec<FieldSpec>,

    /// Columns to add to the default set rather than replace it.
    #[arg(
        short = 'F',
        long = "add-field",
        value_name = "FIELD",
        value_parser = RecordField::parse,
        value_delimiter = ',',
        help = "Record columns to add to the default set; see --help"
    )]
    pub add_field: Vec<FieldSpec>,

    /// Records that match no query still get a row, with the query columns
    /// null. Pass this to write only rows that matched.
    #[arg(long)]
    pub only_hits: bool,

    /// Output format for the hits.
    #[arg(long, value_enum, default_value_t = OutputFormatArg::Parquet)]
    pub output_format: OutputFormatArg,

    /// Rows held in memory before they are written, in total across every
    /// output file.
    ///
    /// Each output file buffers its rows until it holds its share of this
    /// number, so this is divided evenly over the files the run writes
    /// (`--threads` x `--shards-per-worker`). Buffered rows are most of what a
    /// run holds in memory, so this sizes memory directly, whatever the file
    /// count.
    #[arg(long, default_value_t = 500_000, value_parser = positive)]
    pub batch_rows: usize,

    /// Number of rows per row group in the output parquet
    #[arg(long, default_value_t = 50_000, value_parser = positive)]
    pub row_group_rows: usize,

    /// How hard to work at making the output small.
    ///
    /// `fast` (the default) writes LZ4 with chunk-level statistics; `small`
    /// writes zstd with page-level statistics. Either can be overridden piece
    /// by piece with --compression and --statistics.
    #[arg(long, value_enum, value_name = "PROFILE", default_value_t = OutputProfile::Fast)]
    pub output_profile: OutputProfile,

    /// Block codec for the parquet pages. Overrides --output-profile.
    #[arg(long, value_enum, value_name = "CODEC")]
    pub compression: Option<CompressionArg>,

    /// Compression level, for zstd only.
    ///
    /// Ignored unless the codec in force is zstd, whether that came from
    /// --compression or from --output-profile small.
    #[arg(long, default_value_t = 3)]
    pub compression_level: i32,

    /// How much statistics to compute. Overrides --output-profile.
    #[arg(long, value_enum, value_name = "LEVEL")]
    pub statistics: Option<StatisticsArg>,

    /// Override how one column is stored, as COLUMN=ENCODING. Repeatable.
    ///
    /// Encodings are `dict`, `plain`, `delta` (integer columns) and
    /// `delta-text` (text columns). Every column already has a default chosen
    /// for what it holds -- dictionaries for the low-cardinality ones, delta
    /// for coordinates and offsets, prefix-delta for `qname`, plain for `seq`
    /// and `qual` -- so this is for the cases that table cannot know about.
    ///
    /// Aux tags are named by the tag and default by type: text and integer tags
    /// are dictionary-encoded, float tags plain. That is right for labels and
    /// codes like RG:Z and NM:i and wrong for the per-base array tags, which
    /// are a long distinct string per record: prefer
    /// `--column-encoding ML=plain --column-encoding MM=plain` when selecting
    /// those.
    #[arg(
        long,
        value_name = "COLUMN=ENCODING",
        value_parser = EncodingOverride::parse,
        action = clap::ArgAction::Append,
        value_delimiter = ','
    )]
    pub column_encoding: Vec<EncodingOverride>,

    /// Record attributes whose combined value decides which file a record's
    /// rows are written to.
    ///
    /// The same vocabulary as `--field`, read by the same code, so a key names
    /// exactly what the column of that name would say.
    #[arg(
        long = "partition-by",
        value_name = "FIELD",
        value_parser = RecordField::parse,
        value_delimiter = ',',
        long_help = partition_help(),
        help = "Record columns that pick the output file; see --help [default: qname]"
    )]
    pub partition_by: Vec<FieldSpec>,

    /// Output files per worker.
    ///
    /// Separate from `--threads` because they answer different questions:
    /// threads are how fast the scan runs, files are how the output is
    /// partitioned for whatever reads it next. A run writes
    /// `--threads * --shards-per-worker` files in total, named
    /// `{stem}_{worker}_{shard}.{ext}` after OUT: `calls.parquet` becomes
    /// `calls_0_0.parquet`, ...
    ///
    /// Every file holds its own row builders and write buffer, so memory
    /// scales with the file count: budget roughly `--batch-rows` rows per
    /// file, and lower `--batch-rows` when asking for many files.
    #[arg(
        long,
        value_name = "N",
        default_value_t = 1,
        help = "Output files per worker; the total is --threads times this"
    )]
    pub shards_per_worker: usize,
}

impl ParquetArgs {
    /// The resolved partition key, defaulting to the qname.
    ///
    /// The only place [`partition::DEFAULT_KEY`] is applied: everything below
    /// this takes the key as given and refuses an empty one, so there is one
    /// answer to what a run with no `--partition-by` partitions by.
    pub fn partition_fields(&self) -> Vec<RecordField> {
        if self.partition_by.is_empty() {
            partition::DEFAULT_KEY.to_vec()
        } else {
            RecordField::resolve(&self.partition_by, &[])
        }
    }

    /// The resolved column set. Delegates: the replace-vs-add rule, the group
    /// expansion and the dedupe all belong with the fields, not with the parser.
    pub fn record_fields(&self) -> Vec<RecordField> {
        RecordField::resolve(&self.field, &self.add_field)
    }

    /// Each file's share of `--batch-rows`, rounded up so the shares cover the
    /// total and never reach zero.
    ///
    /// `--batch-rows` is a total across every output file, as the memory knobs
    /// are throughout, so that raising the shard count does not quietly raise
    /// the run's memory.
    pub fn batch_rows_per_file(&self, n_files: usize) -> usize {
        self.batch_rows.div_ceil(n_files.max(1)).max(1)
    }

    /// The parquet settings, profile resolved and flags applied over it.
    ///
    /// The only place the profile is turned into settings, so `fast` means one
    /// thing. Per-column encodings are passed through unresolved: they are
    /// matched against the schema by the writer, which is the first point at
    /// which the hit columns, the record columns and the aux tags are one list.
    pub fn parquet_opts(&self) -> Result<ParquetOpts, String> {
        let zstd = |lvl: i32| -> Result<Compression, String> {
            ZstdLevel::try_new(lvl)
                .map(Compression::ZSTD)
                .map_err(|_| format!("--compression-level {lvl} is not a zstd level; expected 1 to 22"))
        };
        let (profile_codec, profile_stats) = match self.output_profile {
            OutputProfile::Fast => (Compression::LZ4_RAW, EnabledStatistics::Chunk),
            OutputProfile::Small => (zstd(self.compression_level)?, EnabledStatistics::Page),
        };
        let compression = match self.compression {
            None => profile_codec,
            Some(CompressionArg::None) => Compression::UNCOMPRESSED,
            Some(CompressionArg::Lz4) => Compression::LZ4_RAW,
            Some(CompressionArg::Snappy) => Compression::SNAPPY,
            Some(CompressionArg::Zstd) => zstd(self.compression_level)?,
        };
        let statistics = match self.statistics {
            None => profile_stats,
            Some(StatisticsArg::None) => EnabledStatistics::None,
            Some(StatisticsArg::Chunk) => EnabledStatistics::Chunk,
            Some(StatisticsArg::Page) => EnabledStatistics::Page,
        };
        Ok(ParquetOpts {
            compression,
            statistics,
            max_row_group_rows: Some(self.row_group_rows),
            encodings: self.column_encoding.clone(),
            metadata: Vec::new(),
        })
    }
}

/// What a list of query files declares, merged.
pub struct QueryFiles {
    pub aliases: Aliases,
    pub queries: Vec<QuerySpec>,
    pub tags: TagConfig,
    /// The one `[strand.*]` rule the run's files declare, named after the file
    /// it came from. `None` when none of them declares one.
    pub strand: Option<crate::strand_rule::StrandRule>,
    /// Every file as read, in order.
    pub sources: Vec<QuerySource>,
}

/// Read and parse query files, in order, as TOML whatever their extension.
/// Each file's aliases merge into one table, so a letter meaning two things
/// across files is an error. Parse warnings are printed here.
pub fn read_query_files(paths: &[PathBuf]) -> Result<QueryFiles, String> {
    let mut out = QueryFiles {
        aliases: Aliases::new(),
        queries: Vec::new(),
        tags: TagConfig::default(),
        strand: None,
        sources: Vec::new(),
    };
    for path in paths {
        let src = std::fs::read_to_string(path)
            .map_err(|e| format!("reading {}: {e}", path.display()))?;
        let file = query_toml::parse_file(&src).map_err(|e| {
            let is_toml = path.extension().is_some_and(|x| x.eq_ignore_ascii_case("toml"));
            // A file in the old line syntax fails on its first line; say why.
            let hint = if is_toml { "" } else { " (query files are TOML; see --help on --query-file)" };
            format!("{}: {e}{hint}", path.display())
        })?;
        for w in &file.warnings {
            eprintln!("warning: {}: {w}", path.display());
        }
        out.aliases.merge(&file.aliases).map_err(|e| format!("{}: {e}", path.display()))?;
        out.tags.merge(file.tags).map_err(|e| format!("{}: {e}", path.display()))?;
        if let Some(rule) = file.strand {
            take_strand(&mut out.strand, rule, &path.display().to_string())?;
        }
        out.queries.extend(file.queries);
        out.sources.push(QuerySource { path: path.display().to_string(), text: src });
    }
    Ok(out)
}

/// The command-line options, among those `cmd` defines, that `m` says were
/// given explicitly rather than left at their defaults -- by long name.
pub fn given_options(m: &clap::ArgMatches, cmd: clap::Command) -> Vec<String> {
    cmd.get_arguments()
        .filter(|a| m.value_source(a.get_id().as_str()) == Some(clap::parser::ValueSource::CommandLine))
        .map(|a| a.get_long().map(|l| format!("--{l}")).unwrap_or_else(|| a.get_id().to_string()))
        .collect()
}

/// Options that do not apply to the output `query` was asked for: parquet
/// options without `--parquet`, BAM options with it. Refused rather than
/// ignored, since an option that silently does nothing reads as one that
/// worked.
pub fn misplaced_query_options(m: &clap::ArgMatches, parquet: bool) -> Vec<String> {
    use clap::Args as _;
    if parquet {
        let bam_only = clap::Command::new("bam")
            .arg(clap::Arg::new("overwrite_tags").long("overwrite-tags"))
            .arg(clap::Arg::new("writer_threads").long("writer-threads"));
        given_options(m, bam_only)
    } else {
        given_options(m, ParquetArgs::augment_args(clap::Command::new("parquet")))
    }
}

/// Queries and the alias table they were parsed against.
#[derive(Debug, Clone)]
pub struct Resolved {
    pub aliases: Aliases,
    pub queries: Vec<QuerySpec>,
    /// Tags declared by TOML query files, merged across files. Validated here
    /// so a bad declaration fails before any BAM is opened.
    pub tags: TagConfig,
    /// How this run recovers each record's strand, from the `[strand.*]`
    /// tables of one of its query files. `None` only for `--explain` and
    /// `--list-codes`, which read no records and so call no strands;
    /// [`resolve`](QueryArgs::resolve) refuses a run that would.
    pub strand: Option<crate::strand_rule::StrandRule>,
    /// Every query file as read, in command-line order, for the output header.
    pub sources: Vec<QuerySource>,
}

/// What a run with no declared strand rule is told.
///
/// It names a file rather than describing one, because the fix is to pass a
/// file: `query/strand/` ships one per surveyed aligner, and the right one
/// is a property of the BAM, which alnbase will not guess at.
pub const NO_STRAND_RULE: &str = "no strand rule: one --query-file must declare \
     [strand.original] and [strand.aligned], which say how this aligner records the strand a \
     read came from. alnbase ships one file per surveyed aligner in query/strand/ -- pass the \
     one that matches the BAM, for example --query-file query/strand/directional.toml";

/// The strand rule declared by query files that have already been read --
/// notably the ones stored in a tagged BAM's header, which is where `extract`
/// gets its definition when the command line gives none.
///
/// A re-extraction has to orient offsets the way the tagging run did, and the
/// tagging run's rule travelled with its query files, so this reads it back
/// out rather than making the user remember to pass it again.
pub fn strand_of_sources(sources: &[QuerySource]) -> Result<Option<StrandRule>, String> {
    let mut out: Option<StrandRule> = None;
    for s in sources {
        let file = query_toml::parse_file(&s.text).map_err(|e| format!("{}: {e}", s.path))?;
        if let Some(rule) = file.strand {
            take_strand(&mut out, rule, &s.path)?;
        }
    }
    Ok(out)
}

/// One run, one rule. Two of them could only disagree about the same record,
/// and there is no principled way to pick a winner -- so this is not a merge,
/// it is a check. The rule learns its file's name here because only the loader
/// knows it, and a record no rule covers must say which rule failed to cover
/// it.
fn take_strand(out: &mut Option<StrandRule>, rule: StrandRule, name: &str) -> Result<(), String> {
    if let Some(first) = out {
        return Err(format!(
            "{name}: declares a strand rule, but {} already did; a run has exactly one",
            first.source()
        ));
    }
    *out = Some(rule.named(name));
    Ok(())
}

impl QueryArgs {
    /// The walk policy. The pad and intron-context counts stay `None` unless
    /// given, because their default is the widest query's span, and only the
    /// query set knows that.
    pub fn walk_config(&self) -> WalkConfig {
        WalkConfig {
            insertions: self.insertions,
            end_context: self.end_context,
            splice_context: self.splice_context,
            max_discordance: self.max_discordance,
        }
    }

    /// True when no BAM will be opened. `--trace-records` reads one, so it is
    /// not a dry run even though it never writes output.
    pub fn is_dry_run(&self) -> bool {
        self.explain || self.list_codes || self.trace.is_some()
    }

    /// Read every query file and check the run's queries against one another.
    pub fn resolve(&self) -> Result<Resolved, String> {
        let QueryFiles { aliases, queries, tags, strand, sources } =
            read_query_files(&self.query_file)?;

        // --list-codes is a reference lookup, not a query operation.
        if queries.is_empty() && tags.is_empty() && !self.list_codes {
            return Err("no queries given; pass --query-file".to_string());
        }

        // A single file cannot contain a name conflict -- the parser rejects
        // duplicates as it reads -- but a run may combine several files, and
        // nothing has compared them until now.
        //
        // Duplicate *patterns* cost nothing. The automaton keeps one copy of
        // a pattern, so two queries that use one pattern share a single set of
        // bits. It is a duplicate *name* that is the problem here.
        dsl::validate_set(&queries)?;

        // The strand is a declaration like any other, so a run that reads
        // records has to have been given one. There is no default, because a
        // default would be a rule compiled into the program that quietly
        // disagrees with the aligner that wrote the BAM. The dry runs call no
        // strands, so they are exempt; `--trace-records` is not one of them.
        if strand.is_none() && !self.is_dry_run() {
            return Err(NO_STRAND_RULE.to_string());
        }

        Ok(Resolved { aliases, queries, tags, strand, sources })
    }

    /// The three positional paths, required unless this is a dry run.
    pub fn paths(&self) -> Result<(&PathBuf, &PathBuf, &PathBuf), String> {
        match (&self.bam, &self.refr, &self.out) {
            (Some(b), Some(r), Some(h)) => Ok((b, r, h)),
            _ => Err("BAM, reference index and output path are all required unless \
                      --explain, --trace or --list-codes is given"
                .to_string()),
        }
    }
}

/// A fraction between 0 and 1 inclusive.
fn fraction(s: &str) -> Result<f64, String> {
    match s.parse::<f64>() {
        Ok(f) if (0.0..=1.0).contains(&f) => Ok(f),
        Ok(_) => Err("must be between 0 and 1".into()),
        Err(e) => Err(e.to_string()),
    }
}

/// A count that must be at least 1.
fn positive(s: &str) -> Result<usize, String> {
    match s.parse::<usize>() {
        Ok(0) => Err("must be at least 1".into()),
        Ok(n) => Ok(n),
        Err(e) => Err(e.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::{CommandFactory, Parser};
    use crate::encoding::ColumnEncoding;

    /// A standalone parser, so these tests do not need the top-level `Cli`.
    #[derive(Parser, Debug)]
    struct Harness {
        #[command(flatten)]
        args: QueryArgs,
    }

    /// A temporary TOML query file holding `text`, named for `tag`.
    ///
    /// Tests run in parallel, and `fs::write` truncates before it writes, so a
    /// file two tests rewrite can be read empty in between. Each file is
    /// therefore written to a scratch name and renamed into place, which
    /// readers see as one step; the process id keeps two `cargo test` runs
    /// from sharing names at all.
    fn query_file(tag: &str, text: &str) -> String {
        let dir = std::env::temp_dir();
        let pid = std::process::id();
        let path = dir.join(format!("alnbase_cli_{tag}_{pid}.toml"));
        let scratch = dir.join(format!("alnbase_cli_{tag}_{pid}_{:?}.part", std::thread::current().id()));
        std::fs::write(&scratch, text).unwrap();
        std::fs::rename(&scratch, &path).unwrap();
        path.to_str().unwrap().to_string()
    }

    /// The query file every parsed command line starts with: one query.
    /// Written once per test process and only read after that.
    fn default_file() -> String {
        static PATH: std::sync::OnceLock<String> = std::sync::OnceLock::new();
        PATH.get_or_init(|| query_file("default", "[query.cg]\n[pat.p]\nread = \"C~\"\nrefr = \"CG\"\n"))
            .clone()
    }

    fn parse(extra: &[&str]) -> QueryArgs {
        let default = default_file();
        let mut argv = vec!["alnbase", "--query-file", default.as_str()];
        argv.extend_from_slice(extra);
        argv.extend(["a.bam", "r.mm", "out.bam"]);
        Harness::try_parse_from(argv).expect("parse").args
    }

    /// The shipped directional rule's path, for tests whose run has to resolve:
    /// `resolve` refuses a run that reads records without a strand rule, and
    /// most of these tests are about something else entirely.
    fn strand_arg() -> String {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("query/strand/directional.toml")
            .display()
            .to_string()
    }

    /// A `strand` tag resolves under a rule that reads the original strand from
    /// a tag and the aligned direction from the FLAG, since the pair names the
    /// strand of origin.
    #[test]
    fn a_strand_tag_resolves_under_any_rule() {
        let tag = query_file(
            "xrtag",
            "[query.cg]\n[pat.p]\nread = \"C~\"\nrefr = \"CG\"\n\n\
             [tag.XR.strand]\nCT = [\"OT\", \"CTOB\"]\nGA = [\"CTOT\", \"OB\"]\n",
        );
        let yd = query_file(
            "ydrule",
            "[strand.original]\nforward = 'YD == \"f\"'\nreverse = 'YD == \"r\"'\n\n\
             [strand.aligned]\nforward = \"not is_reverse\"\nreverse = \"is_reverse\"\n",
        );
        for rule in [yd, strand_arg()] {
            let resolved = Harness::try_parse_from([
                "alnbase", "--query-file", &tag, "--query-file", &rule, "a.bam", "r.mm", "out.bam",
            ])
            .expect("parse")
            .args
            .resolve();
            assert!(resolved.is_ok(), "{rule}: {:?}", resolved.err());
        }
    }

    /// The parquet options, parsed on their own.
    #[derive(Parser, Debug)]
    struct ParquetHarness {
        #[command(flatten)]
        args: ParquetArgs,
    }

    fn parquet(extra: &[&str]) -> ParquetArgs {
        let mut argv = vec!["alnbase"];
        argv.extend_from_slice(extra);
        ParquetHarness::try_parse_from(argv).expect("parse").args
    }

    #[test]
    fn args_are_well_formed() {
        Harness::command().debug_assert();
    }

    #[test]
    fn queries_come_from_toml_files() {
        let more = query_file(
            "more",
            "[query.mCpG]\nmark = \"...+........\"\nwhere = \"wide\"\n[pat.wide]\nread = \"~~~Y~~~~~~~~\"\nrefr = \"~~~CG~~~~~~~\"\n",
        );
        let strand = strand_arg();
        let r = parse(&["--query-file", &more, "--query-file", &strand]).resolve().unwrap();
        let names: Vec<&str> = r.queries.iter().map(|q| q.name.as_str()).collect();
        assert_eq!(names, ["cg", "mCpG"]);
        assert_eq!(r.queries[1].anchor, 3);
        assert_eq!(r.sources.len(), 3, "every file is kept for the output header");
    }

    #[test]
    fn there_is_no_command_line_query_input() {
        for flags in [&["-q", "C~@CG"][..], &["--query", "C~@CG"], &["--alias", "j={C.}"]] {
            let mut argv = vec!["alnbase"];
            argv.extend_from_slice(flags);
            assert!(Harness::try_parse_from(argv).is_err(), "{flags:?} should be rejected");
        }
        let e = Harness::try_parse_from(["alnbase", "--explain"]).unwrap().args.resolve().unwrap_err();
        assert!(e.contains("pass --query-file"), "{e}");
    }

    /// Every query file is read as TOML. One in the old line syntax fails on
    /// its first line, with a pointer to the format.
    #[test]
    fn a_line_syntax_file_is_refused_with_a_hint() {
        let path = std::env::temp_dir().join(format!("alnbase_cli_line_syntax_{}.txt", std::process::id()));
        std::fs::write(&path, "pat p\n  read C~\n  refr CG\nquery cg\n").unwrap();
        let e = parse(&["--query-file", path.to_str().unwrap()]).resolve().unwrap_err();
        assert!(e.contains("query files are TOML"), "{e}");
    }

    #[test]
    fn duplicate_names_across_files_are_rejected() {
        let again = query_file("again", "[query.cg]\n[pat.q]\nread = \"T~\"\nrefr = \"CG\"\n");
        let e = parse(&["--query-file", &again]).resolve().unwrap_err();
        assert!(e.contains("named 'cg'"), "{e}");
    }

    #[test]
    fn the_query_file_help_shows_a_valid_example() {
        let h = query_toml::help();
        assert!(h.contains("[tag.XX.bases]") && h.contains(query_toml::EXAMPLE), "{h}");
    }

    #[test]
    fn dry_runs_do_not_need_paths() {
        let default = default_file();
        let a = Harness::try_parse_from(["alnbase", "--query-file", default.as_str(), "--explain"])
            .expect("parse")
            .args;
        assert!(a.is_dry_run());
        assert!(a.paths().is_err());
        assert!(a.resolve().is_ok());
    }

    #[test]
    fn a_dash_names_standard_input() {
        let default = default_file();
        let a = Harness::try_parse_from(["alnbase", "--query-file", default.as_str(), "-", "r.mm", "out.bam"])
            .expect("parse")
            .args;
        assert_eq!(a.bam.as_ref().unwrap().as_os_str(), "-");
        assert!(a.paths().is_ok());
    }

    #[test]
    fn no_fields_means_the_default() {
        assert_eq!(parquet(&[]).record_fields(), RecordField::resolve(&[], &[]));
    }

    #[test]
    fn fields_reach_the_resolver_in_order() {
        let args = parquet(&["-f", "mapq,cigar", "-f", "NM:i"]);
        let names: Vec<String> = args.record_fields().iter().map(|f| f.name()).collect();
        assert_eq!(names, ["mapq", "cigar", "NM"]);
    }

    #[test]
    fn add_field_extends_rather_than_replaces() {
        let base = parquet(&[]).record_fields();
        let got = parquet(&["-F", "NM:i"]).record_fields();
        assert_eq!(got.len(), base.len() + 1);
        assert!(got.starts_with(&base), "additions must come after the default");
        assert_eq!(got.last().unwrap().name(), "NM");
    }

    #[test]
    fn field_and_add_field_combine() {
        let got = parquet(&["-f", "mapq", "-F", "cigar"]).record_fields();
        assert_eq!(got, vec![RecordField::MapQ, RecordField::Cigar]);
    }

    #[test]
    fn the_partition_key_defaults_to_the_qname() {
        // The default that keeps a fragment's alignments in one file.
        assert_eq!(parquet(&[]).partition_fields(), partition::DEFAULT_KEY.to_vec());
        assert_eq!(parquet(&[]).shards_per_worker, 1);
    }

    #[test]
    fn a_partition_key_resolves_like_a_field_list() {
        // Same vocabulary as --field: names, aliases, aux tags, order kept,
        // duplicates dropped.
        let args = parquet(&["--partition-by", "rname,RG:Z", "--partition-by", "ref_name"]);
        let names: Vec<String> =
            args.partition_fields().iter().map(|f| f.name()).collect();
        assert_eq!(names, ["ref_name", "RG"]);
    }

    /// The key replaces the default rather than extending it, unlike -F. A key
    /// is a grouping, not a column set, so silently adding the qname to it
    /// would put every record in its own group.
    #[test]
    fn a_partition_key_replaces_the_default() {
        let fields = parquet(&["--partition-by", "ref_name"]).partition_fields();
        assert_eq!(fields, vec![RecordField::RefName]);
        assert!(!fields.contains(&RecordField::Qname));
    }

    #[test]
    fn the_files_per_worker_is_settable() {
        assert_eq!(parquet(&["--shards-per-worker", "8"]).shards_per_worker, 8);
        // --batch-rows is a total, shared out over the files.
        let p = parquet(&["--batch-rows", "1000"]);
        assert_eq!(p.batch_rows_per_file(1), 1000);
        assert_eq!(p.batch_rows_per_file(8), 125);
        assert_eq!(p.batch_rows_per_file(3), 334, "rounded up so the shares cover the total");
        assert_eq!(parquet(&["--batch-rows", "2"]).batch_rows_per_file(8), 1, "never zero");
    }

    /// The key is independent of the emitted columns: partitioning by
    /// something does not put it in the file, and neither implies the other.
    #[test]
    fn the_partition_key_and_the_column_selection_are_separate() {
        let args = parquet(&["--partition-by", "mapq", "-f", "cigar"]);
        assert_eq!(args.partition_fields(), vec![RecordField::MapQ]);
        assert_eq!(args.record_fields(), vec![RecordField::Cigar]);
    }

    /// Parquet options are on the query command line, but only apply with
    /// --parquet; BAM options only without it.
    #[test]
    fn options_for_the_other_output_are_refused() {
        use clap::CommandFactory;
        let check = |extra: &[&str]| {
            let default = default_file();
            let mut argv = vec!["alnbase", "--query-file", default.as_str()];
            argv.extend_from_slice(extra);
            argv.extend(["a.bam", "r.mm", "out"]);
            let m = Harness::command().try_get_matches_from(argv).expect("parse");
            let parquet = m.get_flag("parquet_output");
            misplaced_query_options(&m, parquet)
        };
        assert_eq!(check(&["-F", "NM:i", "--only-hits"]), ["--add-field", "--only-hits"]);
        assert_eq!(check(&["--shards-per-worker", "4"]), ["--shards-per-worker"]);
        assert!(check(&["--parquet", "-F", "NM:i", "--shards-per-worker", "4"]).is_empty());
        assert_eq!(check(&["--parquet", "--overwrite-tags"]), ["--overwrite-tags"]);
        assert!(check(&["--overwrite-tags", "--writer-threads", "2"]).is_empty());
    }

    #[test]
    fn the_bam_output_flags() {
        let a = parse(&[]);
        assert!(!a.overwrite_tags);
        assert!(parse(&["--overwrite-tags"]).overwrite_tags);
        assert_eq!(a.writer_threads, 0);
        assert_eq!(parse(&["--writer-threads", "6"]).writer_threads, 6);
        assert_eq!(a.out.as_deref(), Some(std::path::Path::new("out.bam")));
    }

    /// Fast is the default, and fast means LZ4 with chunk statistics. The
    /// profile is the only place that mapping is written down, so this is the
    /// test that keeps the help text honest.
    #[test]
    fn the_output_defaults_to_the_fast_profile() {
        let o = parquet(&[]).parquet_opts().unwrap();
        assert_eq!(o.compression, Compression::LZ4_RAW);
        assert_eq!(o.statistics, EnabledStatistics::Chunk);
        assert_eq!(o.max_row_group_rows, Some(50_000));
        assert!(o.encodings.is_empty());
    }

    #[test]
    fn the_small_profile_is_zstd_at_the_given_level() {
        let o = parquet(&["--output-profile", "small", "--compression-level", "9"])
            .parquet_opts()
            .unwrap();
        assert_eq!(o.compression, Compression::ZSTD(ZstdLevel::try_new(9).unwrap()));
        assert_eq!(o.statistics, EnabledStatistics::Page);
    }

    /// The granular flags are overrides, so either half of a profile can be
    /// replaced without restating the other.
    #[test]
    fn a_granular_flag_beats_the_profile() {
        let o = parquet(&["--output-profile", "small", "--compression", "none"])
            .parquet_opts()
            .unwrap();
        assert_eq!(o.compression, Compression::UNCOMPRESSED, "the flag wins");
        assert_eq!(o.statistics, EnabledStatistics::Page, "the profile still sets the rest");

        let o = parquet(&["--statistics", "none"]).parquet_opts().unwrap();
        assert_eq!(o.compression, Compression::LZ4_RAW);
        assert_eq!(o.statistics, EnabledStatistics::None);
    }

    /// The level is a zstd setting, so it is only validated when zstd is what
    /// is in force -- a level left over in a script does not fail a fast run.
    #[test]
    fn the_compression_level_only_matters_to_zstd() {
        assert!(parquet(&["--compression-level", "99"]).parquet_opts().is_ok());
        let e = parquet(&["--output-profile", "small", "--compression-level", "99"])
            .parquet_opts()
            .unwrap_err();
        assert!(e.contains("1 to 22"), "{e}");
    }

    /// Overrides are passed through by name and resolved against the schema by
    /// the writer, which is the first place all three kinds of column are one
    /// list. Repeated and comma-separated both work, as everywhere else here.
    #[test]
    fn column_encodings_are_collected_by_name() {
        let o = parquet(&["--column-encoding", "ML=plain,MM=plain", "--column-encoding", "qname=dict"])
            .parquet_opts()
            .unwrap();
        assert_eq!(
            o.encodings,
            vec![
                EncodingOverride { column: "ML".into(), encoding: ColumnEncoding::Plain },
                EncodingOverride { column: "MM".into(), encoding: ColumnEncoding::Plain },
                EncodingOverride { column: "qname".into(), encoding: ColumnEncoding::Dictionary },
            ]
        );
    }

    #[test]
    fn walk_policy_defaults_to_methyldackel() {
        assert_eq!(parse(&[]).walk_config(), WalkConfig::default());
    }

    /// The intron window is derived unless asked for, and an explicit value is
    /// what pins the column layout a motif pattern is written against.
    #[test]
    fn the_splice_context_is_derived_unless_given() {
        assert_eq!(parse(&[]).walk_config().splice_context, None);
        let args = parse(&["--splice-context", "2"]);
        assert_eq!(args.walk_config().splice_context, Some(2));
        // Derived, it tracks the widest query; fixed, it does not.
        assert_eq!(WalkConfig::default().to_opts(5).splice_context, 5);
        assert_eq!(args.walk_config().to_opts(5).splice_context, 2);
    }

    #[test]
    fn walk_policy_is_settable() {
        let args = parse(&["--insertions", "emit", "--max-discordance", "0.5"]);
        assert_eq!(
            args.walk_config(),
            WalkConfig { insertions: Insertions::Emit, end_context: None, splice_context: None, max_discordance: 0.5 }
        );
        assert_eq!(parse(&[]).walk_config().max_discordance, 0.25, "on by default");
    }

    /// Tags from several files end up in one config -- once each -- and a
    /// TOML error names the file, then the line.
    #[test]
    fn tags_from_several_files_merge() {
        let queries = query_file(
            "tags_queries",
            "[pat.p]\nread = \"C~\"\nrefr = \"CG\"\n[query.cgx]\n[tag.XM.bases]\nfill = \".\"\nZ = \"cgx\"\n",
        );
        let strands = query_file("tags_strands", "[tag.XG.strand]\nCT = [\"OT\", \"CTOT\"]\nGA = [\"OB\", \"CTOB\"]\n");
        let rule = strand_arg();
        let r = parse(&["--query-file", &queries, "--query-file", &strands, "--query-file", &rule])
            .resolve()
            .unwrap();
        let tags: Vec<String> = r.tags.tags.iter().map(|t| t.name.to_string()).collect();
        assert_eq!(tags, ["XM", "XG"]);

        let e = parse(&["--query-file", &strands, "--query-file", &strands]).resolve().unwrap_err();
        assert!(e.contains("tag XG is declared in more than one query file"), "{e}");

        let broken = query_file("tags_broken", "[tag.XG.strand]\nCT = [\"OT\"]\nGA = OB\n");
        let e = parse(&["--query-file", &broken]).resolve().unwrap_err();
        assert!(e.starts_with(&format!("{broken}: line 3:")), "{e}");
    }

    /// A strand rule is the one declaration a run may hold only one of, so it
    /// is checked rather than merged, and it learns the name of the file it
    /// came from on the way through.
    #[test]
    fn one_strand_rule_comes_through_and_a_second_is_refused() {
        let rule = "[strand.original]\nforward = 'XG == \"CT\"'\nreverse = 'XG == \"GA\"'\n\
                    [strand.aligned]\nforward = 'not is_reverse'\nreverse = 'is_reverse'\n";
        let a = query_file("strand_rule_a", rule);
        let b = query_file("strand_rule_b", rule);

        let r = parse(&["--query-file", &a]).resolve().unwrap();
        assert_eq!(r.strand.as_ref().map(|s| s.source()), Some(a.as_str()));

        let e = parse(&["--query-file", &a, "--query-file", &b]).resolve().unwrap_err();
        assert!(e.starts_with(&format!("{b}: declares a strand rule")), "{e}");
        assert!(e.contains(&a), "the error names the file that got there first: {e}");
    }

    /// `extract` with no --query-file reads its definition out of the BAM's
    /// header, and the rule has to come back with it: the tags it is reading
    /// were oriented by that rule, so a re-extraction under another one would
    /// mirror offsets the tagging run did not.
    #[test]
    fn a_rule_is_recovered_from_stored_query_files() {
        let text = "[strand.original]\nforward = 'true'\n\
                    [strand.aligned]\nforward = 'not is_reverse'\nreverse = 'is_reverse'\n";
        let plain = QuerySource { path: "q.toml".into(), text: "[query.c]\nread = \"C\"\nrefr = \"C\"\n".into() };
        let ruled = QuerySource { path: "strand.toml".into(), text: text.into() };

        assert!(strand_of_sources(&[plain.clone()]).unwrap().is_none());
        let rule = strand_of_sources(&[plain.clone(), ruled.clone()]).unwrap().unwrap();
        assert_eq!(rule.source(), "strand.toml", "the rule names the stored file it came from");

        let e = strand_of_sources(&[ruled.clone(), ruled]).unwrap_err();
        assert!(e.contains("a run has exactly one"), "{e}");
    }

    /// A run that reads records and was given no rule is refused, and the
    /// refusal names a file to pass: there is no default strand rule, because
    /// a default is a rule compiled into the program that can disagree with
    /// the aligner that wrote the BAM.
    #[test]
    fn a_run_with_no_strand_rule_is_refused() {
        let e = parse(&[]).resolve().unwrap_err();
        assert!(e.contains("query/strand/directional.toml"), "{e}");
    }

    /// The dry runs that read no records are exempt, since a rule they would
    /// never consult is not something to demand. `--trace-records` reads a
    /// BAM, so it is not one of them.
    #[test]
    fn explain_and_list_codes_need_no_strand_rule() {
        assert!(parse(&["--explain"]).resolve().unwrap().strand.is_none());
        assert!(parse(&["--list-codes"]).resolve().unwrap().strand.is_none());
        assert!(parse(&["--trace", "C@C"]).resolve().unwrap().strand.is_none());
    }
}
