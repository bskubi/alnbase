//! Parse the arguments of `extract` and `dump-query`. Both commands read a
//! tagged BAM instead of writing one.

use std::path::PathBuf;

use clap::Args;

use crate::cli_query::ParquetArgs;

/// `alnbase extract`: parquet hit rows from a tag.
#[derive(Args, Debug)]
pub struct ExtractArgs {
    /// Tagged BAM, or '-' for standard input.
    pub bam: PathBuf,

    /// Output prefix: `out.parquet` writes `out_{worker}_{shard}.parquet`.
    pub out: PathBuf,

    /// The tag to extract, e.g. XM.
    ///
    /// Needed only when more than one bases tag is defined. Taken from the
    /// most recent alnbase run that stored a definition of it, unless
    /// --query-file is given.
    #[arg(long, value_name = "TAG")]
    pub tag: Option<String>,

    /// Read the tag with these query files instead of the ones stored in the
    /// BAM's header. Repeatable.
    ///
    /// Also how a tag written by another tool is read: a file whose bases tag
    /// maps that tool's characters to queries is its definition.
    #[arg(long = "query-file", value_name = "FILE")]
    pub query_file: Vec<PathBuf>,

    /// Warn and count instead of stopping on data the run cannot handle;
    /// see `alnbase query --help`.
    #[arg(long)]
    pub permissive: bool,

    /// Worker threads.
    #[arg(short = '@', long, default_value_t = 1)]
    pub threads: usize,

    /// BGZF decompression threads for reading the BAM. 0 means min(4, --threads).
    #[arg(long, default_value_t = 0)]
    pub reader_threads: usize,

    #[command(flatten, next_help_heading = "Parquet output")]
    pub parquet: ParquetArgs,
}

/// `alnbase dump-query`: print a query file stored in a BAM's header.
#[derive(Args, Debug)]
pub struct DumpQueryArgs {
    /// BAM written by `alnbase query`, or '-' for standard input.
    pub bam: PathBuf,

    /// The run to dump, by its @PG ID. Defaults to the most recent.
    #[arg(long, value_name = "ID")]
    pub pg: Option<String>,

    /// Which of the run's files to dump, counting from 0. Needed only when the
    /// run stored more than one.
    #[arg(long, value_name = "N")]
    pub file: Option<usize>,
}
