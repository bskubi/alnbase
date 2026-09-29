//! Query aligned read and reference columns, and report every position where a
//! query matches.
//!
//! This page is for a person who reads or changes the source. It explains the
//! model that the whole program is built on. It then follows one record through
//! the program, and says which module owns which decision. The manual for users
//! is separate. Read `docs/alnbase.md` for the tour, and `docs/reference/` for
//! the query language, the walk, the outputs and `overlap`.
//!
//! Rebuild these pages with `cargo docs --open`. `alnbase` is a binary, and
//! almost every item in it is private, so a plain `cargo doc` renders almost
//! nothing. `cargo docs` is the alias in `.cargo/config.toml` that supplies
//! `--document-private-items`.
//!
//! # What the program is for
//!
//! A methylation caller has its rules built in. Bismark decides what a CG is,
//! and it decides what to do when a deletion falls inside the context. You
//! cannot see these decisions, and you cannot correct one that is wrong.
//! alnbase puts the decisions into a file that you write:
//!
//! ```toml
//! [query.CG]
//! read = "C~"
//! refr = "CG"
//! ```
//!
//! This query finds a protected CG cytosine on the top strand. The read must
//! have a `C` where the reference has the `C` of a `CG`. The second column has
//! no constraint on the read side. alnbase scans a BAM for that pattern and
//! reports each place where it matches. It writes the report as tags on the
//! reads, or as one parquet row per hit.
//!
//! No part of the program is specific to methylation. A query is a pattern over
//! aligned columns. The user gives the pattern its meaning.
//!
//! A query *fires* at a position when it matches there. The source, the help
//! text and the outputs all use the word in that sense.
//!
//! # Three ideas the code assumes
//!
//! Almost every module assumes these three ideas, and the code is difficult to
//! follow without them.
//!
//! **1. A base is a set of possible values, not a letter.** A [`seq::Seq`] is
//! one byte. It holds four nucleotide bits and four flag bits.
//!
//! | bits | meaning |
//! |---|---|
//! | `A` `C` `G` `T` | the four nucleotides, one bit each |
//! | `GAP` | the other side has a base here and this side does not: a deletion on the read side, or an insertion on the reference side |
//! | `PAD` | this side has ended: the position is past the end of the read, or past the end of the contig |
//! | `CLIP` | a read base that the aligner soft-clipped. The sequencer read it, but the aligner did not align it |
//! | `SKIP` | a `CIGAR N`: reference that the alignment spans but the read does not. For RNA this is an intron |
//!
//! Each flag sits on the side that is missing something, so the side is part of
//! the meaning of the flag. The four are easy to confuse, and [`seq`] describes
//! all of them in full.
//!
//! A base is a set, so an ambiguity code is several bits set at the same time.
//! `N` is `A|C|G|T`. A pattern column matches when the observed set is a subset
//! of the set of the pattern. This also lets [`query`] resolve `{...}` groups
//! and relational codes into lookup tables when it compiles the query, so they
//! add no cost to the scan.
//!
//! **2. The column is the unit of work.** A [`column::Column`] is one aligned
//! position. It holds the read base and the reference base at that position. It
//! also holds what an output needs in order to report the position and its
//! quality: the reference coordinate, the offset into the `SEQ` of the read, and
//! the base quality. The program moves these columns past a set of patterns.
//!
//! The walk makes the columns from the CIGAR. An indel is therefore an explicit
//! column and not a hidden adjustment. The walk also continues a short distance
//! past the alignment. It adds a flank of reference at each end of the read, and
//! reference context at each edge of an intron. A pattern whose reference half
//! reaches past the part that the read covers can therefore still match. A
//! cytosine at the last base of a read keeps its `CG` context.
//!
//! **3. A query is a boolean predicate, not a pattern.** A query names one or
//! more patterns and combines them with `and`, `or` and `not`. All of the
//! operands span the same window, so "column 3" means the same column in each
//! operand. The anchor and the captures therefore belong to the query, and not
//! to one pattern. The matching engine reports which patterns fired. A separate
//! boolean layer then decides which queries fired. See [`dsl`] for the language
//! and [`predicate`] for the boolean layer.
//!
//! # How one record moves through the program
//!
//! Each step below is a module.
//!
//! 1. **Read the record.** [`bam_io`] opens the BAM, and a pipe also works.
//!    [`batch`] moves the records between threads in batches. A handoff of one
//!    record costs more than the work on that record.
//! 2. **Find the strand that the record came from.** [`strand_rule`] applies a
//!    rule that the user declares over flags and tags. Aligners do not agree
//!    about how to express the original strand, so this rule is a query file and
//!    not a built-in reading of the SAM flag. If the rule cannot classify a
//!    record, alnbase skips the record and counts it. It never guesses.
//! 3. **Walk the alignment.** [`alignment`] steps through the CIGAR against the
//!    memory-mapped reference. `alnbase index` builds that reference ([`aref`]),
//!    and [`contig_map`] matches it to the contigs of the BAM. The walk turns
//!    each CIGAR step into [`column::Column`] values. It also adds the flank
//!    columns and the intron-edge columns that are described above. There is
//!    usually one column per aligned position, but not always. An inserted base
//!    is its own column, and the interior of a long intron becomes one marker
//!    instead of one column per skipped base.
//! 4. **Match all of the patterns at the same time.** [`byg_search`] is a
//!    bit-parallel shift-and automaton, from Baeza-Yates and Gonnet's "A New
//!    Approach to Text Searching" (CACM 1992). alnbase extends it to hold all of
//!    the patterns in one bitmap, with one bit per *trace*. A trace is a prefix
//!    of a pattern that still matches the columns that came before. Advancing
//!    one column costs the same few bitwise operations for any number of
//!    patterns, so the total length of the patterns sets the cost and their
//!    number does not. [`query`] therefore keeps one copy of a pattern that two
//!    queries share, and the second query adds no cost.
//! 5. **Keep the earlier columns.** A pattern reports on its *last* column, but
//!    the coordinates that it needs belong to earlier columns. [`column_ring`]
//!    therefore keeps the last `max_span` columns, so that the anchor and the
//!    captures can point back into them.
//! 6. **Decide which queries fired.** [`predicate`] evaluates the boolean
//!    expression of each query over the bitmap of pattern hits from the
//!    automaton.
//! 7. **Write the answer.** [`bam_out`] writes tags onto the record. Under
//!    `--parquet`, [`hits`] builds one row for each hit instead. That row holds
//!    the columns of the query next to the columns of the record itself. The
//!    columns of the record come from [`record_columns`]. [`hit_writer`] and
//!    [`arrow_builder`] then write the row into a parquet file.
//!
//! [`scanner`] is the driver. It does steps 3 to 7 for one record. It also
//! fixes the order of the ring and the automaton. [`column_ring`] must store
//! the column before [`byg_search`] advances, because the back-index of the ring
//! counts from the column that the ring stored last. [`scanner`] states what the
//! wrong order costs.
//!
//! # The two scan pipelines
//!
//! Both pipelines feed [`scanner`]. They differ in what the output must
//! guarantee.
//!
//! | pipeline | used by | shape | why |
//! |---|---|---|---|
//! | [`ordered`] | `query` when it writes a BAM | one reader, a pool of taggers, one writer | a tagged BAM must come out in input order. The pipeline puts finished batches back in sequence before it writes them |
//! | [`parallel`] | `query --parquet`, `extract` | sharded: each worker owns its own output files | parquet rows need no order. The pipeline therefore has no single writer to wait on. [`partition`] selects the file for each record |
//!
//! # Where each subcommand starts
//!
//! | subcommand | entry point | what it does |
//! |---|---|---|
//! | `index` | [`aref::write_from_fasta`] | build the `.aref` reference store from a FASTA |
//! | `info` | [`run_info`] | report what an index holds, and print a region to compare against `samtools faidx` |
//! | `query` | [`run_query`] | the scan: tags on reads, or parquet hit rows. See [`cli_query`] for the arguments |
//! | `extract` | [`run_extract`] | read hits again out of a `bases` tag that an earlier `query` wrote, through [`tag_extract`] |
//! | `dump-query` | [`run_dump_query`] | print the query files that are stored in the header of a BAM ([`header_query`]) |
//! | `overlap` | [`overlap_apply::run`] | resolve bases that a fragment observes twice. See [`cli_overlap`] for the arguments |
//!
//! `overlap` is a separate tool that shares this binary. It runs upstream of any
//! caller, including a caller that is not alnbase.
//!
//! # Map of the modules
//!
//! **The query language**, in the order in which a file passes through it:
//! [`query_toml`] reads the TOML. [`lower`] resolves aliases, parses patterns,
//! and attaches anchors and captures. [`dsl`] is the language itself: columns,
//! patterns and boolean expressions. [`query`] compiles the language into an
//! automaton plan and one predicate per query. [`byg_search`] and [`predicate`]
//! then run at scan time. [`tags`] declares which BAM tags a run writes.
//!
//! **The walk:** [`seq`], [`mod@column`], [`column_ring`], [`alignment`],
//! [`aref`], [`contig_map`].
//!
//! **Output:** [`bam_out`] writes tags. [`hits`], [`record_columns`],
//! [`record_field`], [`hit_writer`] and [`arrow_builder`] write parquet rows.
//! [`encoding`] says how each output column is stored.
//! [`manifest`] and [`header_query`] write the provenance for both kinds of
//! output. A file from an earlier run therefore still records what its tags
//! mean.
//!
//! **Mate overlap:** [`overlap`] decides which bases of a fragment were observed
//! twice. It reads no records and writes none. It works on the position along
//! the molecule, and not on reference intervals, so a chimeric fragment also
//! works. [`overlap_apply`] is the part that touches real records. The two
//! modules are separate, because the reasoning is then easy to test on its own.
//!
//! **Explaining a run to its user:** [`explain`] shows what a query says, column
//! by column. [`prose`] describes a query in English. [`trace`] runs queries
//! against a single record and shows the work. [`codes`] prints the code table.
//!
//! **Plumbing:** [`bam_io`], [`batch`], [`aux`], [`ordered`], [`parallel`],
//! [`partition`], [`scanner`]. Two more modules are `#[cfg(test)]` and do not
//! appear in these pages at all. `golden` pins the output of the whole pipeline
//! byte for byte, and `test_support` holds the shared record and genome
//! fixtures.
//!
//! # Suggested reading order
//!
//! Read [`seq`] and [`mod@column`] first. They are small. The rest of the
//! program is difficult to read without them. Then read [`dsl`] to see what a
//! query is, [`byg_search`] to see how matching works, and [`scanner`] to see
//! the two joined together. [`alignment`] holds the most detail, and it is the
//! module where correctness is most fragile. Read it after you know what a
//! `Column` is for. [`overlap`] is self-contained, and you can read it at any
//! point.
//!
//! # Conventions to know before you change anything
//!
//! - **Coordinates.** Reference positions are 0-based inside the program.
//!   alnbase reports them in the form that each output format requires. Read
//!   offsets (`off_5p`, `off_3p`) are defined against the original read as the
//!   sequencer read it, and include hard-clipped bases. They are not defined
//!   against the BAM record. A column carries the plain `SEQ` offset, which is a
//!   different number, and [`hits`] adds the hard clip back on when it builds a
//!   row. See `docs/design/off-by-one-safeguards.md`.
//! - **Filtering belongs elsewhere.** alnbase has no filter for MAPQ, base
//!   quality or flags. `samtools view` does that upstream, and the parquet
//!   output carries the columns that you filter on downstream. alnbase does not
//!   add a second implementation of a filter that another program already has.
//! - **Knowledge ships as query files, not as code.** The rules for each aligner
//!   and each protocol live in `query/`, which users can read and fork.
//!   Support for a new aligner therefore never needs a recompile.

mod alignment;
mod aux;
mod bam_out;
mod batch;
mod arrow_builder;
mod bam_io;
mod byg_search;
mod cli_extract;
mod cli_query;
mod cli_overlap;
mod codes;
mod column;
mod column_ring;
mod contig_map;
mod dsl;
mod encoding;
mod explain;
mod golden;
mod header_query;
mod hits;
mod parallel;
mod partition;
mod ordered;
mod hit_writer;
mod predicate;
mod prose;
mod query;
mod query_toml;
mod record_columns;
mod record_field;
mod aref;
mod overlap;
mod overlap_apply;
mod scanner;
#[cfg(test)]
mod test_support;
mod trace;
mod seq;
mod tag_extract;
mod strand_rule;
mod tags;
mod lower;
mod manifest;

use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{anyhow, Result};
use clap::Parser;
use mimalloc::MiMalloc;

use cli_extract::{DumpQueryArgs, ExtractArgs};
use cli_query::QueryArgs;
use cli_overlap::OverlapArgs;
use ordered::OrderedConfig;
use query::QuerySet;
use aref::{write_from_fasta, Aref};
use bam_out::{BamOptions, BamStream};

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

/// Query aligned read and reference columns.
///
/// You write a pattern over aligned columns, and alnbase records where that
/// pattern fires. It writes the record as tags on the reads, or as a table with
/// one row per hit.
///
/// The subcommands are independent tools that share one binary, and each
/// subcommand owns its own argument file. `overlap` in particular is made to run
/// upstream of any caller, including a caller that is not alnbase.
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(clap::Subcommand, Debug)]
pub enum Commands {
    /// Make reference genome index
    Index(IndexArgs),

    /// Report what is in a reference index
    Info(InfoArgs),

    /// Tag a BAM from boolean queries over read/reference column patterns
    Query(QueryArgs),

    /// Write parquet hit rows from a tag in a BAM
    Extract(ExtractArgs),

    /// Print a query file stored in a BAM's header by `alnbase query`
    DumpQuery(DumpQueryArgs),

    /// Resolve duplicated observations within a sequenced fragment
    Overlap(OverlapArgs),
}

#[derive(clap::Args, Debug)]
pub struct IndexArgs {
    pub input_fasta: PathBuf,
    pub output_index: PathBuf,
}

#[derive(clap::Args, Debug)]
pub struct InfoArgs {
    /// Reference index written by `alnbase index`.
    pub index: PathBuf,

    /// List every contig, rather than only the first and last few.
    #[arg(long)]
    pub all: bool,

    /// Print the sequence over a region, as alnbase reads it, for example
    /// `chr19:1000-1100`.
    ///
    /// The region is 1-based and inclusive, as it is in `samtools faidx`. You
    /// can therefore compare the two outputs directly, which is the purpose of
    /// this option. If you build an index from one FASTA and give a downstream
    /// tool another FASTA, each tool makes calls that agree with themselves and
    /// disagree with each other. This is how you see that.
    ///
    /// The output is uppercase. The index stores which base is at a position,
    /// and not whether a masker made it lowercase. Compare the two sides without
    /// case by piping both through `tr a-z A-Z`.
    #[arg(long, value_name = "CHR:FROM-TO")]
    pub seq: Option<String>,
}

/// Print `region` of `aref` as FASTA. It uses the coordinates and the line
/// width of `samtools faidx`, so that a `diff` of the two outputs tells you
/// whether the two references agree.
fn print_region(aref: &Aref, region: &str) -> Result<()> {
    let (name, span) = region
        .split_once(':')
        .ok_or_else(|| anyhow!("expected CHR:FROM-TO, got {region}"))?;
    let (from, to) = span
        .split_once('-')
        .ok_or_else(|| anyhow!("expected CHR:FROM-TO, got {region}"))?;
    let from: i64 = from.replace(',', "").parse().map_err(|_| anyhow!("{region}: FROM is not a number"))?;
    let to: i64 = to.replace(',', "").parse().map_err(|_| anyhow!("{region}: TO is not a number"))?;
    if from < 1 || to < from {
        return Err(anyhow!("{region}: FROM must be at least 1 and TO at least FROM"));
    }
    let tid = aref
        .tid(name)
        .ok_or_else(|| anyhow!("the index has no contig {name}"))?;
    let len = aref.len_of(tid).unwrap_or(0) as i64;
    if from > len {
        return Err(anyhow!("{name} is {len} bases, so {from} is past its end"));
    }
    // faidx is 1-based inclusive; the store is 0-based half-open.
    let seq = aref
        .slice(tid, from - 1, to.min(len))
        .ok_or_else(|| anyhow!("could not read {region}"))?;
    let text: String = seq.iter().map(|s| s.to_ascii() as char).collect();
    println!(">{name}:{from}-{}", to.min(len));
    for line in text.as_bytes().chunks(60) {
        println!("{}", String::from_utf8_lossy(line));
    }
    Ok(())
}

/// `alnbase info`: report what an index holds, without opening a BAM.
///
/// The main reason for this subcommand is the question "is this the reference
/// that the BAM was aligned against?". You can answer it here, instead of
/// running a scan and reading the error that it gives you.
fn run_info(args: InfoArgs) -> Result<()> {
    let aref = Aref::open(&args.index)?;
    if let Some(region) = &args.seq {
        return print_region(&aref, region);
    }
    let size = std::fs::metadata(&args.index).map(|m| m.len()).unwrap_or(0);
    print!("{}", info_text(&aref, &args.index.display().to_string(), size, args.all));
    Ok(())
}

/// Build the report as text, so that a test can hold the whole report and not
/// the pieces that it is made of.
fn info_text(aref: &Aref, path: &str, size: u64, all: bool) -> String {
    use std::fmt::Write as _;

    let n = aref.n_contigs();
    let total: usize = (0..n).filter_map(|i| aref.len_of(i)).sum();
    let mut out = String::new();
    let _ = writeln!(out, "{path}");
    let _ = writeln!(out, "  format version  {}", aref.version());
    let _ = writeln!(out, "  contigs         {n}");
    let _ = writeln!(out, "  total bases     {total}");
    let _ = writeln!(out, "  file size       {size} bytes");
    if n == 0 {
        return out;
    }

    let _ = writeln!(out, "\ncontigs");
    // A whole-genome index has thousands of contigs, nearly all of them
    // unplaced scaffolds; the first and last few answer "which assembly is
    // this" without burying it.
    const SHOWN: usize = 5;
    let width = (0..n).filter_map(|i| aref.name(i)).map(str::len).max().unwrap_or(0);
    let show = |out: &mut String, i: usize| {
        if let (Some(name), Some(len), Some(md5)) = (aref.name(i), aref.len_of(i), aref.md5_hex(i)) {
            let _ = writeln!(out, "  {name:<width$}  {len:>12}  {md5}");
        }
    };
    if all || n <= SHOWN * 2 + 1 {
        for i in 0..n {
            show(&mut out, i);
        }
    } else {
        for i in 0..SHOWN {
            show(&mut out, i);
        }
        let _ = writeln!(out, "  ... {} more, --all to list them", n - SHOWN * 2);
        for i in n - SHOWN..n {
            show(&mut out, i);
        }
    }
    out
}

#[cfg(test)]
mod info_tests {
    use super::*;

    /// The region printer agrees with `samtools faidx`. It is 1-based and
    /// inclusive, and it clips at the end of the contig. With any other
    /// behaviour, a comparison of two references would compare two different
    /// windows and would not say so.
    #[test]
    fn info_seq_matches_faidx_coordinates() {
        let dir = std::env::temp_dir();
        let fa = dir.join(format!("alnbase_seq_{}.fa", std::process::id()));
        let idx = dir.join(format!("alnbase_seq_{}.aref", std::process::id()));
        let seq = "ACGTTGCAAGGCTTACGATC";
        std::fs::write(&fa, format!(">c1\n{seq}\n")).unwrap();
        write_from_fasta(&fa, &idx).unwrap();
        let aref = Aref::open(&idx).unwrap();

        let at = |from: i64, to: i64| {
            let tid = aref.tid("c1").unwrap();
            let s = aref.slice(tid, from - 1, to.min(seq.len() as i64)).unwrap();
            s.iter().map(|b| b.to_ascii() as char).collect::<String>()
        };
        assert_eq!(at(1, 5), &seq[0..5]);
        assert_eq!(at(11, 20), &seq[10..20]);
        assert_eq!(at(20, 20), &seq[19..20], "a single base");
        assert_eq!(at(18, 40), &seq[17..], "clipped at the contig end");
        assert!(print_region(&aref, "c1:0-5").is_err(), "0 is not a faidx coordinate");
        assert!(print_region(&aref, "cX:1-5").is_err());
        assert!(print_region(&aref, "c1:5").is_err());
        assert!(print_region(&aref, "c1:9-3").is_err(), "TO before FROM");
    }

    #[test]
    fn info_reports_the_index() {
        let dir = std::env::temp_dir();
        let fa = dir.join(format!("alnbase_info_{}.fa", std::process::id()));
        let idx = dir.join(format!("alnbase_info_{}.aref", std::process::id()));
        std::fs::write(&fa, ">chr1\nACGTACGTAC\n>chr2\nACGT\n").unwrap();
        write_from_fasta(&fa, &idx).unwrap();
        let aref = Aref::open(&idx).unwrap();

        assert_eq!(
            info_text(&aref, "ref.aref", 123, false),
            "ref.aref\n  \
             format version  3\n  \
             contigs         2\n  \
             total bases     14\n  \
             file size       123 bytes\n\ncontigs\n  \
             chr1            10  45aff2fecf7615d56bc0567dffab9fa8\n  \
             chr2             4  f1f8f4bf413b16ad135722aa4591043e\n"
        );
        // The digests are the SAM @SQ M5 values; `samtools dict` on the same
        // FASTA writes exactly these.
    }

    /// Most of a whole-genome index is unplaced scaffolds, so the command
    /// shortens the list unless you ask for the full list.
    #[test]
    fn a_long_contig_list_is_abridged() {
        let dir = std::env::temp_dir();
        let fa = dir.join(format!("alnbase_info_many_{}.fa", std::process::id()));
        let idx = dir.join(format!("alnbase_info_many_{}.aref", std::process::id()));
        let body: String = (0..14).map(|i| format!(">c{i}\nACGT\n")).collect();
        std::fs::write(&fa, body).unwrap();
        write_from_fasta(&fa, &idx).unwrap();
        let aref = Aref::open(&idx).unwrap();

        let short = info_text(&aref, "x", 0, false);
        assert!(short.contains("... 4 more, --all to list them"), "{short}");
        assert!(short.contains("  c0 ") && short.contains("  c13 "), "both ends are shown");
        assert!(!short.contains("  c7 "), "the middle is not");

        let full = info_text(&aref, "x", 0, true);
        assert!(full.contains("  c7 ") && !full.contains("more, --all"), "{full}");
    }
}

fn run_index(args: IndexArgs) -> Result<()> {
    write_from_fasta(args.input_fasta, args.output_index)?;
    Ok(())
}

/// Scan a BAM for boolean queries over read and reference columns. Write the
/// BAM back out, in order, with the tags that the query files declare.
///
/// There is one code path for any number of threads. With `--threads 1` the
/// path is a reader that feeds one tagger and a writer. That costs a few channel
/// sends per 512 records, and it keeps the code that runs in production the same
/// as the code that runs in a test.
fn run_query(args: QueryArgs) -> Result<()> {
    // Queries come only from --query-file TOML files. See `QueryArgs::resolve`
    // for why this happens after clap has parsed everything.
    let mut resolved = args.resolve().map_err(|e| anyhow!(e))?;
    // Every path below this line takes the run's one strand rule; `resolve`
    // has already refused a run that reads records without one.
    let strand = resolved.strand.take().map(Arc::new);

    // A reference lookup: prints whatever aliases this invocation declared
    // alongside the fixed table, so it answers "what is still free" too.
    if args.list_codes {
        print!("{}", codes::table(&resolved.aliases));
        return Ok(());
    }


    let queries = QuerySet::compile(&resolved.queries).map_err(|e| anyhow!("{e}"))?;

    // Columns that parsed but contribute nothing. Not fatal, but almost always
    // a typo, and silence here means a query quietly does less than it says.
    for w in &queries.warnings {
        eprintln!("{w}");
    }

    // `--explain` renders every query column-aligned and stops. The DSL is dense
    // enough that seeing the operands stacked is the difference between trusting
    // a query and hoping about it — in the methyl-Hi-C case the whole point is
    // that the exclusion's chimeric CG lands on the same columns as the
    // cytosine, which no amount of squinting at the one-line form reveals.
    if args.explain {
        for (spec, q) in resolved.queries.iter().zip(queries.queries.iter()) {
            let (form, groups) = q.predicate.shape();
            print!(
                "{}",
                explain::explain_with(spec, Some((form, groups)), &resolved.aliases)
            );
            println!();
        }
        // The set, after the detail: one line each, which is where an odd
        // query stands out against the ordinary ones.
        if queries.queries.len() > 1 {
            println!("compiled shapes");
            for line in queries.shapes() {
                println!("  {line}");
            }
        }
        return Ok(());
    }

    // A synthetic trace needs no alignment, so it runs before any path is
    // required. Both trace paths hand `trace::trace` the same `&[Column]`, so
    // the two cannot disagree about what a trace means.
    if let Some(pair) = &args.trace {
        let cols = trace::columns_from_pair(pair).map_err(|e| anyhow!(e))?;
        print!(
            "{}",
            trace::trace(&queries, &resolved.queries, &cols, "synthetic", args.trace_grid)
        );
        return Ok(());
    }

    // Past every dry run, so the rule `resolve` insisted on is here.
    let strand = strand.expect("resolve() requires a strand rule for a run that reads records");

    let (bam, refr_path, out_path) = args.paths().map_err(|e| anyhow!(e))?;

    // Tracing real records reuses the scanner's own walk, so the columns it
    // reports are exactly the columns the scan would see -- flank and
    // insertion policy included. A trace that built its own columns would be
    // able to disagree with the run it is supposed to explain. Checked before
    // `--parquet`, like every other dry run, so it never starts a full scan.
    if let Some(n) = args.trace_records {
        let refr = Arc::new(Aref::open(refr_path)?);
        return trace_bam(bam, refr, &queries, &resolved.queries, &args, n, &strand);
    }

    if args.parquet_output {
        return run_query_parquet(&args, queries, &resolved.sources, bam, refr_path, out_path, strand);
    }

    // Resolved once, here, before any file is opened. The BAM output is one
    // file, in input order; `--parquet` took the sharded path above.
    let out = bam_out::output_path(out_path)?;
    let output = BamStream::new(
        &resolved.tags,
        &queries,
        BamOptions {
            walk: args.walk_config(),
            strand,
            overwrite_tags: args.overwrite_tags,
            bad_data: batch::OnBadData::from_flag(args.permissive),
            command_line: std::env::args().collect::<Vec<_>>().join(" "),
            query_sources: resolved.sources.clone(),
        },
    )?;
    for w in output.warnings() {
        eprintln!("{w}");
    }

    // Shared, immutable, and read from every worker.
    let queries = Arc::new(queries);
    let refr = Arc::new(Aref::open(refr_path)?);

    let n_workers = args.threads.max(1);
    let auto = |given: usize| if given == 0 { n_workers.min(4) } else { given };
    let cfg = OrderedConfig {
        n_workers,
        reader_threads: auto(args.reader_threads),
        writer_threads: auto(args.writer_threads),
        out,
        bad_data: batch::OnBadData::from_flag(args.permissive),
        require_m5: args.require_m5,
    };

    let stats = ordered::scan_ordered(bam, refr, queries, &cfg, &output)?;
    report(&stats, "tagged");

    Ok(())
}

/// Trace the first `n` mapped records of a BAM.
fn trace_bam(
    bam: &std::path::Path,
    refr: Arc<Aref>,
    queries: &QuerySet,
    specs: &[dsl::QuerySpec],
    args: &QueryArgs,
    n: usize,
    strand: &Arc<strand_rule::StrandRule>,
) -> Result<()> {
    use rust_htslib::bam::{Read as _, Reader};

    let mut reader = Reader::from_path(bam)?;
    let header = reader.header().clone();
    let mut contigs = contig_map::ContigMap::build(&header, &refr)?;
    if args.require_m5 {
        contigs = contigs.require_m5(&header);
    }
    if let Some(w) = contigs.warning() {
        eprintln!("{w}");
    }
    let walk = args.walk_config().to_opts(queries.max_span);

    let mut read_seq = Vec::with_capacity(300);
    let mut seen = 0usize;
    // Checked as a scan checks it, so a trace of a BAM against the wrong
    // reference says so rather than explaining calls that mean nothing.
    let mut concordance = scanner::Concordance::default();

    for rec in reader.records() {
        let rec = rec?;
        match scanner::walkability(&rec, &contigs) {
            scanner::Walkability::Walkable => {}
            scanner::Walkability::MissingM5 => {
                let name = String::from_utf8_lossy(header.tid2name(rec.tid() as u32));
                return Err(batch::missing_m5_error(&name, &String::from_utf8_lossy(rec.qname())));
            }
            _ => continue,
        }
        let mut cols = Vec::new();
        alignment::walk_alignment(
            &rec,
            &refr,
            &contigs,
            &mut read_seq,
            walk,
            match strand.call(&rec)? {
                Some(call) => call,
                // The trace shows the walk, and there is no walk without a
                // call. Named rather than skipped silently: a trace is asked
                // for to find out what happens to a particular read.
                None => {
                    eprintln!(
                        "{}: this run's strand rule declines this record, so there is no walk \
                         to trace",
                        String::from_utf8_lossy(rec.qname())
                    );
                    continue;
                }
            },
            &mut |c| cols.push(c),
        )?;

        let name = String::from_utf8_lossy(rec.qname()).into_owned();
        let contig = if rec.tid() >= 0 {
            String::from_utf8_lossy(header.tid2name(rec.tid() as u32)).into_owned()
        } else {
            "*".to_string()
        };
        let label = format!("{name}  {contig}:{} (1-based POS)", rec.pos() + 1);
        print!("{}", trace::trace(queries, specs, &cols, &label, args.trace_grid));
        let mut this = scanner::Concordance::default();
        for c in &cols {
            this.observe(c);
        }
        println!(
            "reference concordance: {} of {} compared read bases differ ({:.1}%; a reference C read as T is not counted)",
            this.discordant,
            this.compared,
            this.percent()
        );
        println!();
        concordance.add(this);

        seen += 1;
        if seen >= n {
            break;
        }
    }
    if seen == 0 {
        eprintln!("no mapped records to trace");
    }
    let max = args.max_discordance;
    concordance.check(max, scanner::Concordance::FINAL_MIN_COMPARED, "in the traced records")?;
    if concordance.exceeds(max, 1) {
        eprintln!(
            "warning: {:.1}% of the {} read bases traced disagree with the reference, above \
             --max-discordance {max}; too few to judge (a scan judges from {} bases), but a wrong \
             reference, a BAM aligned to another assembly, or shifted coordinates look like this.",
            concordance.percent(),
            concordance.compared,
            scanner::Concordance::FINAL_MIN_COMPARED
        );
    }
    Ok(())
}

fn main() -> Result<()> {
    use clap::{CommandFactory, FromArgMatches};
    // Parsed in two steps so `query` can see which options were given, not
    // just their values; see `cli_query::misplaced_query_options`.
    let matches = Cli::command().get_matches();
    let cli = Cli::from_arg_matches(&matches).map_err(|e| e.exit()).unwrap();

    match cli.command {
        Commands::Index(args) => run_index(args)?,
        Commands::Info(args) => run_info(args)?,
        Commands::Query(args) => {
            let m = matches.subcommand_matches("query").expect("the query subcommand");
            let wrong = cli_query::misplaced_query_options(m, args.parquet_output);
            if !wrong.is_empty() {
                let why = if args.parquet_output { "not used with --parquet" } else { "only used with --parquet" };
                return Err(anyhow!("{why}: {}", wrong.join(", ")));
            }
            run_query(args)?
        }
        Commands::Extract(args) => run_extract(args)?,
        Commands::DumpQuery(args) => run_dump_query(args)?,
        Commands::Overlap(args) => overlap_apply::run(args)?,
    }
    Ok(())
}

/// `alnbase query --parquet`: write hit rows, in shards, directly from the
/// walk.
fn run_query_parquet(
    args: &QueryArgs,
    queries: QuerySet,
    sources: &[header_query::QuerySource],
    bam: &std::path::Path,
    refr_path: &std::path::Path,
    out_path: &std::path::Path,
    strand: Arc<strand_rule::StrandRule>,
) -> Result<()> {
    let p = &args.parquet;
    let mut cfg = shard_config(p, args.threads, args.reader_threads, out_path, args.permissive, strand.clone());
    cfg.require_m5 = args.require_m5;
    let refr = Arc::new(Aref::open(refr_path)?);
    // One reader for the manifest's contig list and the scan, so standard
    // input works.
    let reader = bam_io::open_reader(bam, cfg.reader_threads)?;

    let walk = args.walk_config();
    let opts = walk.to_opts(queries.max_span);
    let manifest = manifest::Manifest {
        alnbase_version: env!("CARGO_PKG_VERSION"),
        format_version: hits::FORMAT_VERSION,
        coordinate_base: manifest::COORDINATE_BASE,
        run_id: manifest::run_id(),
        command: "query",
        command_line: std::env::args().collect(),
        input: bam.display().to_string(),
        reference: Some(refr_path.display().to_string()),
        strand: strand.source().to_string(),
        output: manifest_output(p, &cfg, queries.flat_captures)?,
        query_files: sources
            .iter()
            .map(|s| manifest::QueryFile { path: s.path.clone(), text: s.text.clone() })
            .collect(),
        extract: None,
        walk: Some(manifest::Walk {
            insertions: match opts.insertions {
                alignment::Insertions::Skip => "skip",
                alignment::Insertions::Emit => "emit",
            },
            end_context: opts.flank,
            end_context_given: walk.end_context.is_some(),
            splice_context: opts.splice_context,
            splice_context_given: walk.splice_context.is_some(),
            widest_span: queries.max_span,
            max_discordance: walk.max_discordance,
            require_m5: args.require_m5,
        }),
        contigs: {
            use rust_htslib::bam::Read;
            manifest::contigs(reader.header(), Some(&refr))
        },
    };

    let mut parquet = p.parquet_opts().map_err(|e| anyhow!(e))?;
    parquet.metadata = manifest.footer();
    let out = scanner::OutputConfig {
        bad_data: batch::OnBadData::from_flag(args.permissive),
        fields: p.record_fields(),
        batch_rows: p.batch_rows_per_file(args.threads.max(1) * p.shards_per_worker),
        parquet,
        keep_hitless: !p.only_hits,
        format: p.output_format.into(),
        walk,
        strand,
    };
    manifest.remove_stale()?;
    let queries = Arc::new(queries);
    let stats = parallel::scan_sharded_reader(reader, Some(refr), queries, &cfg, &scanner::ParquetOutput::new(out))?;
    report(&stats, "scanned");
    finish_manifest(&manifest, &stats)
}

/// The output half of the manifest of a hit table.
fn manifest_output(
    p: &cli_query::ParquetArgs,
    cfg: &parallel::ShardConfig,
    flat_captures: bool,
) -> Result<manifest::Output> {
    let router = cfg.router()?;
    let files = (0..router.slots())
        .map(|slot| {
            let (w, s) = router.split(slot);
            cfg.out_path(w, s).file_name().unwrap_or_default().to_string_lossy().into_owned()
        })
        .collect();
    Ok(manifest::Output {
        format: match p.output_format {
            cli_query::OutputFormatArg::Parquet => "parquet",
            cli_query::OutputFormatArg::Ipc => "ipc",
        },
        prefix: cfg.out_prefix.display().to_string(),
        n_workers: router.n_workers(),
        shards_per_worker: router.shards_per_worker(),
        partition_by: cfg.partition_by.iter().map(|f| f.name()).collect(),
        fields: p.record_fields().iter().map(|f| f.name()).collect(),
        only_hits: p.only_hits,
        capture_layout: if flat_captures { "flat" } else { "lists" },
        files,
    })
}

/// Write the completeness marker for a run that succeeded. Also report the
/// files that look like files of this run but are not.
fn finish_manifest(manifest: &manifest::Manifest, stats: &batch::Stats) -> Result<()> {
    let path = manifest.write_complete(stats)?;
    let foreign = manifest.foreign_files();
    if !foreign.is_empty() {
        let shown: Vec<&str> = foreign.iter().take(5).map(String::as_str).collect();
        eprintln!(
            "warning: {} file(s) beside the output follow its naming but are not part of this run \
             (left by an earlier run that wrote more files?): {}{}. {} lists this run's files.",
            foreign.len(),
            shown.join(", "),
            if foreign.len() > shown.len() { ", ..." } else { "" },
            path.display()
        );
    }
    Ok(())
}

/// How many decompression threads to use. This is `--reader-threads` when you
/// give it. If you do not, it is one thread per worker, up to four threads. The
/// scaling of htslib itself becomes flat at about four.
fn reader_threads(threads: usize, reader_threads: usize) -> usize {
    if reader_threads == 0 { threads.max(1).min(4) } else { reader_threads }
}

/// Work out the sharding of a parquet output from its options.
fn shard_config(
    p: &cli_query::ParquetArgs,
    threads: usize,
    reader_threads: usize,
    out: &std::path::Path,
    permissive: bool,
    strand: Arc<strand_rule::StrandRule>,
) -> parallel::ShardConfig {
    let n_workers = threads.max(1);
    // Said once, here, because it is a property of the arguments rather than of
    // anything the scan discovers: a key that can split a fragment across files
    // is legitimate, and silently wrong only for a run that then resolves mate
    // overlaps file by file.
    let partition_by = p.partition_fields();
    if let Some(w) = partition::fragment_warning(&partition_by) {
        eprintln!("{w}");
    }
    parallel::ShardConfig {
        n_workers,
        shards_per_worker: p.shards_per_worker,
        partition_by,
        reader_threads: self::reader_threads(threads, reader_threads),
        out_prefix: out.to_path_buf(),
        bad_data: batch::OnBadData::from_flag(permissive),
        require_m5: false,
        strand,
    }
}

fn report(stats: &batch::Stats, verb: &str) {
    eprintln!(
        "read {} records, {verb} {}, skipped {}{}",
        stats.records_read,
        stats.records_scanned,
        stats.records_skipped,
        match (stats.records_skipped > 0, stats.skipped_written, verb) {
            (true, true, "tagged") => " (written untagged)",
            (true, true, _) => " (given a null row)",
            _ => "",
        }
    );
    if stats.unplaced_hits > 0 {
        eprintln!(
            "  {} hits were anchored on a pad or a soft-clipped base, so had no aligned base to tag",
            stats.unplaced_hits
        );
    }
    if stats.non_utf8_aux > 0 {
        eprintln!(
            "  {} records had text in an aux tag that is not UTF-8, read as null",
            stats.non_utf8_aux
        );
    }
    if stats.damaged_aux > 0 {
        eprintln!(
            "  {} records had a damaged aux block; tags after the damage could not be read",
            stats.damaged_aux
        );
    }
    if stats.unknown_strand > 0 {
        // Said separately from the line above because these records are not
        // written the way the others are: a tagged BAM keeps them, untagged,
        // but a hit table gives them no row at all rather than a null one --
        // a null row would say the record was read and matched nothing, which
        // is a different claim from "this run cannot say which strand it is".
        eprintln!(
            "  of those, {} were declined by this run's strand rule: no strand, so no walk{}. \
             Extend its [strand.*] tables to cover them.",
            stats.unknown_strand,
            if verb == "tagged" { ", and they were written untagged" } else { " and no row" }
        );
    }
    if stats.records_off_reference > 0 {
        eprintln!("  of those, {} were on contigs the reference does not have", stats.records_off_reference);
    }
    if stats.concordance.compared > 0 {
        let c = stats.concordance;
        eprintln!(
            "  {} of {} compared read bases differ from the reference ({:.2}%; a reference C read as T is not counted)",
            c.discordant,
            c.compared,
            100.0 * c.discordant as f64 / c.compared as f64
        );
    }
    if stats.records_without_seq > 0 {
        eprintln!("  of those, {} were mapped but had no SEQ ('*')", stats.records_without_seq);
    }
}

fn header_text(bam: &std::path::Path) -> Result<String> {
    let reader = bam_io::open_reader(bam, 1)?;
    use rust_htslib::bam::Read;
    Ok(String::from_utf8_lossy(reader.header().as_bytes()).into_owned())
}

/// `alnbase extract`.
fn run_extract(args: ExtractArgs) -> Result<()> {
    let want = args.tag.as_deref().map(tags::TagName::parse).transpose().map_err(|e| anyhow!(e))?;
    let mut given_sources = Vec::new();
    let mut rule = None;
    let given = if args.query_file.is_empty() {
        None
    } else {
        let files = cli_query::read_query_files(&args.query_file).map_err(|e| anyhow!(e))?;
        given_sources = files.sources;
        rule = files.strand;
        Some((files.queries, files.tags))
    };
    // One reader for the header and the records, so standard input works: the
    // header is read, the definition chosen from it, and the same reader
    // scans on.
    let p = &args.parquet;
    let reader = bam_io::open_reader(&args.bam, reader_threads(args.threads, args.reader_threads))?;
    let (header, contigs) = {
        use rust_htslib::bam::Read;
        (
            String::from_utf8_lossy(reader.header().as_bytes()).into_owned(),
            manifest::contigs(reader.header(), None),
        )
    };
    let choice = tag_extract::choose(&header, given, want).map_err(|e| anyhow!(e))?;
    eprintln!("extracting tag {} with {}", choice.tag, choice.source);

    let queries = QuerySet::compile(&choice.queries).map_err(|e| anyhow!("{e}"))?;
    let sources = if given_sources.is_empty() { &choice.stored_files } else { &given_sources };

    // With no --query-file, the definition came out of the BAM's header, and so
    // does the rule that oriented the offsets when the tags were written: a
    // re-extraction that mirrored them differently would silently disagree with
    // the run it is reading.
    let strand = match rule {
        Some(r) => Arc::new(r),
        None => Arc::new(
            cli_query::strand_of_sources(sources)
                .map_err(|e| anyhow!(e))?
                .ok_or_else(|| anyhow!(cli_query::NO_STRAND_RULE))?,
        ),
    };
    let cfg = shard_config(p, args.threads, args.reader_threads, &args.out, args.permissive, strand.clone());
    let manifest = manifest::Manifest {
        alnbase_version: env!("CARGO_PKG_VERSION"),
        format_version: hits::FORMAT_VERSION,
        coordinate_base: manifest::COORDINATE_BASE,
        run_id: manifest::run_id(),
        command: "extract",
        command_line: std::env::args().collect(),
        input: args.bam.display().to_string(),
        reference: None,
        strand: strand.source().to_string(),
        // `extract` always writes flat captures; see `TagExtract`.
        output: manifest_output(p, &cfg, true)?,
        query_files: sources
            .iter()
            .map(|s| manifest::QueryFile { path: s.path.clone(), text: s.text.clone() })
            .collect(),
        extract: Some(manifest::Extract { tag: choice.tag.to_string(), definition: choice.source.clone() }),
        walk: None,
        contigs,
    };
    let mut parquet = p.parquet_opts().map_err(|e| anyhow!(e))?;
    parquet.metadata = manifest.footer();
    let out = scanner::OutputConfig {
        bad_data: batch::OnBadData::from_flag(args.permissive),
        fields: p.record_fields(),
        batch_rows: p.batch_rows_per_file(args.threads.max(1) * p.shards_per_worker),
        parquet,
        keep_hitless: !p.only_hits,
        format: p.output_format.into(),
        walk: scanner::WalkConfig::default(),
        strand: strand.clone(),
    };
    let extract = tag_extract::TagExtract::new(&choice, &queries, out, strand)?;
    manifest.remove_stale()?;
    let stats = parallel::scan_sharded_reader(reader, None, Arc::new(queries), &cfg, &extract)?;
    report(&stats, "extracted from");
    finish_manifest(&manifest, &stats)
}

/// `alnbase dump-query`.
fn run_dump_query(args: DumpQueryArgs) -> Result<()> {
    let text = header_query::dump(&header_text(&args.bam)?, args.pg.as_deref(), args.file)
        .map_err(|e| anyhow!(e))?;
    print!("{text}");
    Ok(())
}