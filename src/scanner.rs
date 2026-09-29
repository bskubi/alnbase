//! Drive the whole scan for one record. The scanner extracts the record columns,
//! walks the alignment, and writes the rows.
//!
//! # Two outputs
//!
//! [`Matcher`] holds the walk and the automaton. Both outputs of `query` use it:
//! the tagged BAM ([`crate::bam_out`]) and, with `--parquet`, the rows that
//! [`Scanner`] and [`ParquetOutput`] write here.
//!
//! The row files are denormalized. There is one [`Scanner`] per file, and each
//! worker writes `--shards-per-worker` files. [`crate::parallel`] routes the
//! records to them. Every row holds the columns of the record and the query that
//! fired. You do not need a join to read a hit.
//!
//! A record that fires no query still gets a row. The query columns of that row
//! are null. A record that the walk does not touch also gets such a row. An
//! unmapped read is an example: it goes to [`Scanner::pass_through`] and not to
//! the walk. The file therefore lists every record that the scan saw, and not
//! only the records that matched. You can count a denominator from it.
//! `--only-hits` clears `keep_hitless` and drops both kinds of row.
//!
//! # Order of work in the walk
//!
//! The walk does three things for each column, in this order. First it stores
//! the column in the ring. Second it advances the automaton. Third it evaluates
//! the queries.
//!
//! The first two steps do not depend on each other, but both must run before the
//! evaluation. The store is the more important of the two. A query that ends on
//! this column reads its anchor back out of the ring, and the back-index counts
//! from the column that the walk stored last. If you evaluate before you store,
//! every coordinate in the output moves by one. The walk keeps this order in one
//! place, so that a caller cannot get it wrong.
//!
//! A query cannot fire until the ring holds its full span. The test
//! `q.span <= ring.stored()` makes sure of this. Shorter queries continue to
//! fire at the same time. The first columns of a record are therefore not dead.
//! They are only too few for the widest query in the run.
//!
//! # Why the `any()` gate is safe
//!
//! `QuerySet::compile` rejects each query whose predicate is true against an
//! all-zero hit bitmap. No query can therefore fire at a column that has no
//! pattern match. The `any()` test costs little, and it can therefore skip the
//! full set of predicates. The same rule makes a bare `not X` illegal.
//!
//! The code borrows `QuerySet` and does not own it, because `State` and `Hits`
//! borrow from its plan and because the workers share one automaton. It shares
//! `Aref` through an `Arc`, so that it maps the reference once for any number of
//! threads.

use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use anyhow::{bail, Result};
use rust_htslib::bam::{HeaderView, Record};

use crate::alignment::{walk_alignment, Insertions, WalkOpts, WALK_SAM_FIELDS};
use crate::column::Column;
use crate::column_ring::ColumnRing;
use crate::contig_map::ContigMap;
use crate::hits::HitBuilder;
use crate::hit_writer::{HitWriter, OutputFormat, ParquetOpts};
use crate::parallel::{ShardSink, SinkCounts, SinkFactory};
use crate::predicate::Predicate;
use crate::query::{Query, QuerySet};
use crate::record_field::RecordField;
use crate::aref::Aref;
use crate::seq::Seq;
use crate::byg_search::{Hits, State};
use crate::strand_rule::{StrandCall, StrandRule};

/// Everything about the output file that is not its path.
///
/// Every worker gets a clone of it, so two shards cannot have different
/// schemas or different flush policies.
#[derive(Clone, Debug)]
pub struct OutputConfig {
    /// What to do about a record whose aux block cannot be read in full.
    pub bad_data: crate::batch::OnBadData,
    /// Already resolved: groups expanded, duplicates dropped, order fixed.
    pub fields: Vec<RecordField>,
    pub batch_rows: usize,
    /// Compression, statistics and per-column encodings. The writer ignores all
    /// of these when `format` is IPC, because IPC has none of them.
    pub parquet: ParquetOpts,
    /// Write a row with null query columns for a record that matches nothing.
    pub keep_hitless: bool,
    /// Parquet is the default. Use IPC when you stream the hits to a consumer
    /// instead of writing them for a later pass.
    pub format: OutputFormat,
    pub walk: WalkConfig,
    /// The rule that calls the strand of each record. Every part of the run
    /// that needs a strand takes it from here. The walk direction, the written
    /// offsets and the strand columns therefore cannot disagree.
    pub strand: Arc<StrandRule>,
}

/// The number of rows that the builders of a file hold when the code opens it.
///
/// This value is not `batch_rows`. If the code reserved a full batch at the
/// start, it would reserve every column of every file at once. At the default
/// of 500,000 rows that is tens of megabytes per file, times
/// `--threads * --shards-per-worker` files. The kernel must then fault in and
/// zero all of those buffers before the first row arrives. Profiling measured
/// this at 7-8% of the CPU time of a scan.
///
/// The large reservation also gains nothing after the first batch. An Arrow
/// builder gives its buffers to the array when it does `finish`, so every later
/// batch grows from empty in any case. Growth by doubling from this small size
/// costs a few reallocations per batch, and it touches only the memory that the
/// rows really use. For a file that receives few rows, that is much less
/// memory.
const INITIAL_ROWS: usize = 4096;

/// How often the read disagrees with the reference, over the bases that a scan
/// compared.
///
/// This is a check for the wrong reference. It catches a BAM that was aligned
/// to another assembly, and it catches shifted coordinates. Both of these
/// produce output that looks like data.
///
/// Only an observed column counts, and only when both sides are one
/// unambiguous base. A reference C that the read shows as T is not a
/// disagreement. It is the conversion that the assays alnbase was built for.
/// The walk runs along the conversion strand, so this is the only form that
/// conversion takes on any strand. On a correct alignment the remaining
/// disagreements are sequencing errors and variants, which are usually a few
/// percent at most. A wrong reference or a shifted coordinate gives about 75%.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Concordance {
    pub compared: u64,
    pub discordant: u64,
}

impl Concordance {
    /// How many bases the whole run must compare before `--max-discordance`
    /// can stop it early. A few unusual reads at the start therefore cannot
    /// stop a run.
    pub const MIN_COMPARED: u64 = 100_000;
    /// How many bases a finished run must have compared before the check
    /// judges it at all.
    pub const FINAL_MIN_COMPARED: u64 = 1_000;

    #[inline]
    pub fn observe(&mut self, c: &Column) {
        let single = |s: Seq| matches!(s, Seq::A | Seq::C | Seq::G | Seq::T);
        if c.is_observed() && single(c.read) && single(c.refr) {
            self.compared += 1;
            if c.read != c.refr && !(c.refr == Seq::C && c.read == Seq::T) {
                self.discordant += 1;
            }
        }
    }

    /// Stop with an error when the rate is above `max` and the run compared at
    /// least `min` bases.
    pub fn check(&self, max: f64, min: u64, when: &str) -> Result<()> {
        if self.exceeds(max, min) {
            bail!(
                "{:.1}% of the {} read bases compared {when} disagree with the reference (a \
                 reference C read as T is not counted), above --max-discordance {max}. A wrong \
                 reference, a BAM aligned to a different assembly, or shifted coordinates cause \
                 this; pass --max-discordance 1 to turn the check off.",
                self.percent(),
                self.compared,
            );
        }
        Ok(())
    }

    /// The discordant bases as a percentage of the compared bases. This is 0
    /// when the run compared no bases.
    pub fn percent(&self) -> f64 {
        if self.compared == 0 { 0.0 } else { 100.0 * self.discordant as f64 / self.compared as f64 }
    }

    pub fn add(&mut self, other: Concordance) {
        self.compared += other.compared;
        self.discordant += other.discordant;
    }

    /// Whether the rate is above `max`, which is a fraction, and the run
    /// compared at least `min` bases.
    pub fn exceeds(&self, max: f64, min: u64) -> bool {
        max < 1.0 && self.compared >= min && self.discordant as f64 > max * self.compared as f64
    }
}

/// The concordance of the run. Every worker shares one of these, so the check
/// judges the whole run and not the share of one worker. A check per worker
/// would let a wrong reference pass whenever the input is split into small
/// enough pieces that no single worker reaches the minimum.
#[derive(Debug)]
pub struct SharedConcordance {
    compared: AtomicU64,
    discordant: AtomicU64,
    max: f64,
}

impl SharedConcordance {
    /// `max` is the `--max-discordance` fraction. A value of 1 turns the check
    /// off.
    pub fn new(max: f64) -> Arc<Self> {
        Arc::new(Self { compared: AtomicU64::new(0), discordant: AtomicU64::new(0), max })
    }

    pub fn total(&self) -> Concordance {
        Concordance { compared: self.compared.load(Ordering::Relaxed), discordant: self.discordant.load(Ordering::Relaxed) }
    }

    /// Add the counts of one record. Stop the run when the total of the run is
    /// above the threshold and the run has compared at least
    /// [`Concordance::MIN_COMPARED`] bases.
    fn add(&self, record: Concordance) -> Result<()> {
        self.compared.fetch_add(record.compared, Ordering::Relaxed);
        self.discordant.fetch_add(record.discordant, Ordering::Relaxed);
        self.check(Concordance::MIN_COMPARED, "so far")
    }

    /// The check at the end of the run, over every base that the run compared.
    pub fn check_final(&self) -> Result<()> {
        self.check(Concordance::FINAL_MIN_COMPARED, "in total")
    }

    fn check(&self, min: u64, when: &str) -> Result<()> {
        self.total().check(self.max, min, when)
    }
}

/// The walk and the automaton for one record. This is the part of a scan that
/// decides which queries fire, and where. It does not depend on what the run
/// writes about them.
///
/// Both outputs drive a record through this code, so a hit means the same thing
/// in a parquet row and in a BAM tag.
pub struct Matcher<'a> {
    queries: &'a QuerySet,
    ring: ColumnRing,
    refr: Arc<Aref>,
    contigs: Arc<ContigMap>,
    state: State<'a>,
    hits: Hits,
    read_seq: Vec<Seq>,
    walk_opts: WalkOpts,
    strand: Arc<StrandRule>,
    concordance: Arc<SharedConcordance>,
}

impl<'a> Matcher<'a> {
    pub fn new(
        queries: &'a QuerySet,
        refr: Arc<Aref>,
        contigs: Arc<ContigMap>,
        walk: WalkConfig,
        strand: Arc<StrandRule>,
        concordance: Arc<SharedConcordance>,
    ) -> Self {
        // One derivation. A query of span k reaches k-1 columns past the
        // alignment, so the walk emits at least that many flanking columns at
        // each end (k by default; see `WalkConfig::to_opts`) — without them a
        // cytosine at the last aligned base has no second column and cannot
        // match — and the ring must hold k columns for the anchor lookup. Both
        // come from `max_span`.
        let span = queries.max_span;
        Self {
            state: queries.plan.new_state(),
            hits: queries.plan.new_hits(),
            ring: ColumnRing::new(span),
            queries,
            refr,
            contigs,
            read_seq: Vec::with_capacity(300),
            walk_opts: walk.to_opts(span),
            strand,
            concordance,
        }
    }

    /// How this run calls the strand of `record`. The method is here so that
    /// the walk and the consumer of its columns ask one thing, one time, for
    /// each record.
    ///
    /// It calls [`crate::strand_rule::StrandRule::call`]. `Ok(None)` means that
    /// the rule of the user does not classify this record. The caller then
    /// skips the record and counts the skip. It does not guess a default, so
    /// that the count can tell the user that their rule is incomplete.
    pub fn strand(&self, record: &Record) -> Result<Option<StrandCall>> {
        self.strand.call(record)
    }

    /// Walk `record` along `strand`. Call `on_hit(index, query, ring)` for
    /// every query that fires, in column order. The anchor of the query is
    /// `ring.get(q.anchor_back)`. Return whether any query fired.
    ///
    /// The caller passes `strand` in, and this method does not call for it,
    /// because the caller needs the same value for the row that it writes. See
    /// [`Matcher::strand`].
    pub fn run<F>(&mut self, record: &Record, strand: StrandCall, on_hit: &mut F) -> Result<bool>
    where
        F: FnMut(usize, &Query, &ColumnRing),
    {
        debug_assert_eq!(
            walkability(record, &self.contigs),
            Walkability::Walkable,
            "unwalkable records must be filtered before they reach the scanner"
        );

        // Destructured so the closure below borrows individual fields rather
        // than all of `self`.
        let Self { queries, ring, refr, contigs, state, hits, read_seq, walk_opts, concordance, strand: _ } = self;
        let mut this_record = Concordance::default();

        ring.clear();
        state.reset();

        let mut fired = false;
        walk_alignment(record, refr, contigs, read_seq, *walk_opts, strand, &mut |c: Column| {
            // The ring must hold this column before a query ending on it can
            // read its own anchor back out.
            ring.store(c);
            this_record.observe(&c);
            state.step(QuerySet::compat_idx(c.read, c.refr), &mut *hits);

            // No pattern matched, so no query can be true — see the module
            // comment. One word scan instead of N predicate evaluations, at the
            // overwhelming majority of columns.
            if !hits.inner.any() {
                return;
            }

            let words = hits.inner.words();
            for (i, q) in queries.queries.iter().enumerate() {
                // A query cannot fire before its own span has been stored;
                // shorter queries can still fire in the meantime. It fires on
                // any window of the columns the walk emitted, pads and intron
                // context included; see `WalkConfig` for how many there are.
                if q.span <= ring.stored() && q.predicate.eval(words) {
                    fired = true;
                    on_hit(i, q, ring);
                }
            }
        })?;
        concordance.add(this_record)?;
        Ok(fired)
    }

}

pub struct Scanner<'a> {
    matcher: Matcher<'a>,
    writer: HitWriter<HitBuilder>,
    keep_hitless: bool,
    bad_data: crate::batch::OnBadData,
    damaged_aux: u64,
    non_utf8_aux: u64,
    /// The records that the strand rule of the run did not call. The scan
    /// skips these records.
    unknown_strand: u64,

    n_scanned: u64,
}

impl<'a> Scanner<'a> {
    pub fn open(
        queries: &'a QuerySet,
        refr: Arc<Aref>,
        contigs: Arc<ContigMap>,
        path: &Path,
        cfg: &OutputConfig,
        shard: u32,
        concordance: Arc<SharedConcordance>,
    ) -> Result<Self> {
        // `flat_captures` is a property of the whole query set, so every
        // worker derives the same schema from the same place.
        let cap = cfg.batch_rows.min(INITIAL_ROWS);
        let sink = HitBuilder::new(&cfg.fields, cap, queries.flat_captures, shard);
        let writer = HitWriter::new(path, sink, cfg.batch_rows, cfg.format, &cfg.parquet)?;
        Ok(Self::new(queries, refr, contigs, writer, cfg, concordance))
    }

    pub fn new(
        queries: &'a QuerySet,
        refr: Arc<Aref>,
        contigs: Arc<ContigMap>,
        writer: HitWriter<HitBuilder>,
        cfg: &OutputConfig,
        concordance: Arc<SharedConcordance>,
    ) -> Self {
        Self {
            matcher: Matcher::new(queries, refr, contigs, cfg.walk, cfg.strand.clone(), concordance),
            writer,
            keep_hitless: cfg.keep_hitless,
            bad_data: cfg.bad_data,
            damaged_aux: 0,
            non_utf8_aux: 0,
            unknown_strand: 0,
            n_scanned: 0,
        }
    }


    pub fn scan_with_id(&mut self, rid: u64, record: &Record, header: &HeaderView) -> Result<()> {
        let Self {
            matcher, writer, keep_hitless, bad_data, damaged_aux, non_utf8_aux, unknown_strand,
            n_scanned,
        } = self;

        // Borrows the builder out of the writer for the length of the record;
        // the borrow ends at its last use, which is what lets
        // `writer.maybe_flush()` run straight afterwards.
        let Some(strand) = matcher.strand(record)? else {
            // The run's rule declines this record -- its `unknown` key claims
            // it. Skipping is the only honest answer: the walk has a direction
            // only because the call gives it one, and a row written on a
            // guessed strand reads exactly like a measurement. The count is
            // what turns that into something a user can act on, by extending
            // the rule until it covers their input.
            *unknown_strand += 1;
            return Ok(());
        };
        let sink = writer.builder_mut();
        sink.begin_record(rid, record, header, Some(strand));
        let trouble = sink.aux_trouble();
        if trouble.any() {
            if *bad_data == crate::batch::OnBadData::Stop {
                anyhow::bail!(
                    "record {}: {}. Pass --permissive to read what can be read and count \
                     these records instead.",
                    String::from_utf8_lossy(record.qname()),
                    trouble.describe()
                );
            }
            *damaged_aux += u64::from(trouble.damaged);
            *non_utf8_aux += u64::from(trouble.non_utf8);
        }

        let fired = matcher.run(record, strand, &mut |_, q, ring| sink.push_hit(q, ring))?;

        if !fired && *keep_hitless {
            sink.push_hitless();
        }

        // Flushing is deliberately out here rather than inside the emit loop, so
        // a record's rows never straddle two batches — which is also what makes
        // row groups line up on record boundaries.
        writer.maybe_flush()?;
        *n_scanned += 1;
        Ok(())
    }


    /// Handle a record that the scan cannot walk. Write one row whose values
    /// are all null, unless the run wants only the rows that matched.
    pub fn pass_through(
        &mut self,
        rid: u64,
        record: &Record,
        header: &HeaderView,
    ) -> Result<()> {
        if !self.keep_hitless {
            return Ok(());
        }
        // A record that is not walked is not held to the rule: an unmapped
        // read carries none of the evidence a rule reads, and failing the run
        // over records it was never going to scan would make every shipped
        // rule unusable on a real BAM. Its strand columns are null instead.
        let strand = self.matcher.strand(record).unwrap_or(None);
        let sink = self.writer.builder_mut();
        sink.begin_record(rid, record, header, strand);
        sink.push_hitless();
        self.writer.maybe_flush()?;
        Ok(())
    }

    #[cfg(test)]
    pub fn finish(self) -> Result<u64> {
        Ok(self.close()?.0)
    }

    /// Close the file. Return the number of records scanned and the number of
    /// rows written.
    fn close(self) -> Result<(u64, u64)> {
        let n = self.n_scanned;
        let rows = self.writer.close()?;
        Ok((n, rows))
    }
}

/// The parquet output of `query --parquet`, as a [`SinkFactory`]. It makes one
/// [`Scanner`] per file.
pub struct ParquetOutput {
    cfg: OutputConfig,
    concordance: Arc<SharedConcordance>,
}

impl ParquetOutput {
    pub fn new(cfg: OutputConfig) -> Self {
        let concordance = SharedConcordance::new(cfg.walk.max_discordance);
        Self { cfg, concordance }
    }
}

impl SinkFactory for ParquetOutput {
    type Sink<'q> = Scanner<'q>;

    /// An unmapped record still gets a row, and the query columns of that row
    /// are null. The output therefore accounts for every record in the input.
    fn keeps_unscanned(&self) -> bool {
        true
    }

    /// The fields that the walk needs, and the columns that this output
    /// writes. See [`crate::bam_io::limit_cram_decoding`].
    fn cram_decoding(&self) -> Option<(u32, bool)> {
        let fields = &self.cfg.fields;
        let required = fields.iter().fold(WALK_SAM_FIELDS, |m, f| m | f.sam_fields());
        Some((required, fields.iter().any(|f| f.needs_md_nm())))
    }

    fn open<'q>(
        &self,
        queries: &'q QuerySet,
        refr: Option<Arc<Aref>>,
        contigs: Option<Arc<ContigMap>>,
        _header: &HeaderView,
        path: &Path,
        slot: u32,
    ) -> Result<Scanner<'q>> {
        let (Some(refr), Some(contigs)) = (refr, contigs) else {
            bail!("the parquet query output walks alignments, so it needs a reference");
        };
        Scanner::open(queries, refr, contigs, path, &self.cfg, slot, Arc::clone(&self.concordance))
    }

    fn concordance(&self) -> Option<&SharedConcordance> {
        Some(&self.concordance)
    }
}

impl ShardSink for Scanner<'_> {
    fn scan(&mut self, rid: u64, record: &mut Record, header: &HeaderView) -> Result<()> {
        self.scan_with_id(rid, record, header)
    }

    fn pass(&mut self, rid: u64, record: &mut Record, header: &HeaderView) -> Result<()> {
        self.pass_through(rid, record, header)
    }

    fn finish(self) -> Result<SinkCounts> {
        let counts = SinkCounts {
            damaged_aux: self.damaged_aux,
            non_utf8_aux: self.non_utf8_aux,
            unknown_strand: self.unknown_strand,
            unplaced_hits: 0,
            scanned: 0,
            rows: 0,
        };
        let (scanned, rows) = self.close()?;
        Ok(SinkCounts { scanned, rows, ..counts })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WalkConfig {
    pub insertions: Insertions,
    /// How many pad columns the walk shows past each end of the aligned read.
    /// `None` derives the number, as `k`, where `k` is the span of the widest
    /// query in the run. The widest query can then reach past the end of a read,
    /// and a whole window also fits inside the pads.
    ///
    /// A query fires on any window of the columns that the walk emits, and this
    /// includes a window that holds only pads. With the derived default, a wider
    /// query that you add to a run therefore adds pads, and a query that can
    /// match an all-pad window can gain hits. A read row of `~` or `_` is such a
    /// query. To keep the columns of a run fixed for any set of queries, set
    /// this value yourself. To exclude the all-pad windows instead, write that
    /// in the query, with `and not` over an all-pad pattern.
    pub end_context: Option<usize>,
    /// Stop the run when more than this fraction of the compared read bases
    /// disagree with the reference. A value of 1 turns the check off. See
    /// [`Concordance`].
    pub max_discordance: f64,
    /// How many reference bases the walk shows at each end of an intron.
    /// `None` derives the number.
    ///
    /// The derived number follows `flank`, so the edge of a junction behaves
    /// like the end of a read. A pattern then reaches as far into an intron as
    /// it reaches past the end of a read. This is the correct default for a
    /// pattern that wants reference context. It is the wrong default for a motif
    /// pattern, because the layout moves. The distance from the marker to the
    /// first base of the intron is `context - p` for a pattern at offset `p`,
    /// and the *widest* query in the run sets `context`. An unrelated wider
    /// query that you add to the run therefore slides the intron edges away from
    /// the marker.
    ///
    /// Set the value yourself and the layout stops moving, which is what a motif
    /// pattern needs. With `--splice-context 2`, `~~,@GT,` is a donor site in
    /// every run.
    pub splice_context: Option<usize>,
}

impl Default for WalkConfig {
    fn default() -> Self {
        Self { insertions: Insertions::Skip, end_context: None, splice_context: None, max_discordance: 0.25 }
    }
}

impl WalkConfig {
    pub fn to_opts(self, max_span: usize) -> WalkOpts {
        // One derivation for both boundaries: the widest query's span, past the
        // end of the read *and* into an intron. Deriving the two from the same
        // number is what makes an intron edge behave like a read end. A span-k
        // pattern needs only k-1 columns to reach past its last read column;
        // the k-th means every query, the widest included, also meets windows
        // made only of pads, so no query is the one exception to that rule.
        let reach = max_span;
        WalkOpts {
            flank: self.end_context.unwrap_or(reach),
            splice_context: self.splice_context.unwrap_or(reach),
            insertions: self.insertions,
        }
    }
}
/// Why the scan can or cannot walk a record against the reference.
///
/// The reader decides this one time. The tagged BAM, `query --parquet` and
/// `--trace-records` therefore skip and count exactly the same records.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Walkability {
    Walkable,
    /// FLAG 0x4. There is no alignment to walk.
    Unmapped,
    /// The record is aligned to a contig that the reference index does not
    /// have.
    OffReference,
    /// The record is mapped, but its SEQ is `*`. There are no read bases to
    /// compare.
    NoSeq,
    /// The record is on a contig whose `@SQ` line has no M5, and the run uses
    /// `--require-m5`.
    MissingM5,
}

pub fn walkability(record: &Record, contigs: &ContigMap) -> Walkability {
    if record.is_unmapped() {
        Walkability::Unmapped
    } else if contigs.refr_tid(record.tid()).is_none() {
        Walkability::OffReference
    } else if contigs.lacks_required_m5(record.tid()) {
        Walkability::MissingM5
    } else if record.seq_len() == 0 {
        Walkability::NoSeq
    } else {
        Walkability::Walkable
    }
}

#[cfg(test)]
mod tests {
    use crate::tags::Library;
    use super::*;
    use crate::test_support::{built, contig_map, header, shared_reference, temp_path, LAST_IN_TEMPLATE, REVERSE};
    use arrow::array::{Array, Int64Array, StringArray, UInt64Array};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use rust_htslib::bam::record::Cigar;

    const CONTIG: &str = "chr1";
    //                      0123456789...
    const GENOME: &str = "TTTCGTTTCGTTTT";

    fn queries(toml: &str) -> QuerySet {
        QuerySet::compile(&crate::test_support::queries(toml)).unwrap_or_else(|e| panic!("{e}"))
    }

    fn config(keep_hitless: bool) -> OutputConfig {
        OutputConfig {
            bad_data: crate::batch::OnBadData::Stop,
            fields: RecordField::CORE.to_vec(),
            batch_rows: 16,
            parquet: ParquetOpts::default(),
            keep_hitless,
            format: OutputFormat::Parquet,
            walk: WalkConfig::default(),
            strand: crate::strand_rule::directional(),
        }
    }

    /// Scan `records`, then read the parquet file back as columns.
    ///
    /// It is important that this helper reads the file and does not inspect the
    /// builder. Only a read of the file can catch a schema that disagrees with
    /// the arrays, and only a read exercises the flush path and the footer path
    /// that `finish` exists for.
    fn scan_to_rows(
        tag: &str,
        qs: &QuerySet,
        cfg: &OutputConfig,
        records: &[Record],
    ) -> Vec<(String, i64, String, Option<String>)> {
        let refr = shared_reference(&format!("scan_{tag}"), CONTIG, GENOME);
        let head = header(CONTIG, GENOME.len());
        let path = temp_path(&format!("scan_{tag}"), "parquet");
        {
            let map = contig_map(&refr, CONTIG, GENOME.len());
            let concordance = SharedConcordance::new(cfg.walk.max_discordance);
            let mut sc = Scanner::open(qs, refr, map, &path, cfg, 0, concordance).unwrap();
            for (i, r) in records.iter().enumerate() {
                sc.scan_with_id(i as u64, r, &head).unwrap();
            }
            assert_eq!(sc.finish().unwrap(), records.len() as u64);
        }

        let file = std::fs::File::open(&path).unwrap();
        let reader = ParquetRecordBatchReaderBuilder::try_new(file).unwrap().build().unwrap();
        let mut out = Vec::new();
        for batch in reader {
            let batch = batch.unwrap();
            let name = |c: &str| batch.column_by_name(c).unwrap().clone();
            let names = name("name");
            let names = names.as_any().downcast_ref::<StringArray>().unwrap();
            let pos = name("refr_pos");
            let pos = pos.as_any().downcast_ref::<Int64Array>().unwrap();
            let rid = name("record_id");
            let rid = rid.as_any().downcast_ref::<UInt64Array>().unwrap();
            let rb = batch.column_by_name("read_base").map(|c| {
                c.as_any().downcast_ref::<StringArray>().unwrap().clone()
            });
            for i in 0..batch.num_rows() {
                out.push((
                    if names.is_null(i) { String::new() } else { names.value(i).to_string() },
                    if pos.is_null(i) { -1 } else { pos.value(i) },
                    rid.value(i).to_string(),
                    rb.as_ref().and_then(|a| {
                        if a.is_null(i) { None } else { Some(a.value(i).to_string()) }
                    }),
                ));
            }
        }
        out
    }

    /// A CG in the reference that is methylated in the read. The scan reports
    /// it at the reference coordinate of the anchor.
    #[test]
    fn a_hit_lands_on_the_anchor_coordinate() {
        // GENOME has CG at 3-4 and at 8-9. The read matches the reference.
        let qs = queries("[query.mCpG]\nread = \"Y~\"\nrefr = \"CG\"\n");
        let rec = built(b"r1", b"TTTCGTTTCGTTTT", &[Cigar::Match(14)], 0, 0);
        let rows = scan_to_rows("anchor", &qs, &config(false), &[rec]);

        assert_eq!(rows.len(), 2, "two CpGs in the reference: {rows:?}");
        assert_eq!(rows[0].0, "mCpG");
        // The anchor is column 0 of the pattern, so it reports the C.
        assert_eq!(rows[0].1, 3);
        assert_eq!(rows[1].1, 8);
        assert_eq!(rows[0].3.as_deref(), Some("C"), "read base at the anchor");
    }

    /// Only an observed, unambiguous pair of bases counts. A reference C that
    /// the read shows as T is conversion and not disagreement. The check also
    /// waits until the run has compared enough bases.
    #[test]
    fn concordance_counts_what_it_should() {
        let col = |read, refr| Column { read, refr, refr_pos: 0, read_off: 0, qual: 30 };
        let mut c = Concordance::default();
        c.observe(&col(Seq::A, Seq::A));
        c.observe(&col(Seq::T, Seq::C)); // conversion
        c.observe(&col(Seq::G, Seq::A)); // disagreement
        c.observe(&col(Seq::N, Seq::A)); // ambiguous: not compared
        c.observe(&col(Seq::GAP, Seq::A)); // deletion: not compared
        c.observe(&col(Seq::PAD, Seq::A)); // flank: not observed
        assert_eq!(c, Concordance { compared: 3, discordant: 1 });
        assert!(!c.exceeds(0.25, Concordance::MIN_COMPARED), "too few bases to judge");
        let many = Concordance { compared: Concordance::MIN_COMPARED, discordant: Concordance::MIN_COMPARED / 2 };
        assert!(many.exceeds(0.25, Concordance::MIN_COMPARED));
        assert!(!many.exceeds(1.0, Concordance::MIN_COMPARED), "1 turns the check off");

        // The shared total is what is judged, however the run is split across
        // workers: two workers of half the minimum each still stop the run.
        let shared = SharedConcordance::new(0.25);
        let half = Concordance { compared: Concordance::MIN_COMPARED / 2, discordant: Concordance::MIN_COMPARED / 4 };
        assert!(shared.add(half).is_ok(), "one worker's share alone is below the minimum");
        assert!(shared.add(half).is_err(), "the run's total is not");
        // A small run is still judged when it ends.
        let small = SharedConcordance::new(0.25);
        small.add(Concordance { compared: 2_000, discordant: 1_400 }).unwrap();
        assert!(small.check_final().is_err());
    }

    /// The widest query sizes the pads unless you give a pad count, and a query
    /// also fires on an all-pad window. With the default, a wider query can
    /// therefore add hits to a narrower one. With a fixed pad count it cannot.
    /// A query that must not match an all-pad window says so with `and not`.
    #[test]
    fn pads_follow_the_widest_query_unless_fixed() {
        // Read covers 4..6 (GTT). GENOME's CG at 3-4 straddles the read's first
        // base; the CG at 8-9 starts two columns past the read's end.
        let rec = || built(b"r1", b"GTT", &[Cigar::Match(3)], 4, 0);
        let narrow = "[query.cg]\nread = \"~~\"\nrefr = \"CG\"\n";
        let wide = "[query.wide]\nread = \"~~~~\"\nrefr = \"AAAA\"\n";
        let cg = |tag: &str, toml: &str, cfg: &OutputConfig| {
            let rows = scan_to_rows(tag, &queries(toml), cfg, &[rec()]);
            rows.iter().filter(|r| r.0 == "cg").map(|r| r.1).collect::<Vec<_>>()
        };
        let both = format!("{narrow}{wide}");

        // Derived: two pads per end alone (7-8), four beside the span-4 query (7-10).
        assert_eq!(cg("pads_alone", narrow, &config(false)), vec![3]);
        assert_eq!(cg("pads_beside", &both, &config(false)), vec![3, 8], "the all-pad window at 8-9 fires");

        // Fixed: the same columns, so the same hits, whatever else is queried.
        let mut fixed = config(false);
        fixed.walk.end_context = Some(3);
        assert_eq!(cg("pads_fixed_alone", narrow, &fixed), cg("pads_fixed_beside", &both, &fixed));

        // Excluding all-pad windows in the query itself.
        let no_pads = "[pat.cg_ref]\nread = \"~~\"\nrefr = \"CG\"\n\
                       [pat.all_pad]\nread = \"__\"\nrefr = \"~~\"\n\
                       [query.cg]\nwhere = \"cg_ref and not all_pad\"\n";
        assert_eq!(cg("pads_not", &format!("{no_pads}{wide}"), &config(false)), vec![3]);
    }

    /// With a fixed pad count and a fixed intron-context count, the queries are
    /// independent of each other. This test uses random alignments, which cover
    /// every flag case and the CIGAR operations M, I, D, N, S and H, and it uses
    /// both insertion policies. For each query, the hits of that query alone
    /// equal its hits in a run that also holds every other query and one much
    /// wider query. The test compares hits by record, reference position and
    /// walk offset.
    #[test]
    fn with_fixed_context_every_query_fires_identically_alone_and_in_company() {
        const PANEL: &[(&str, &str, &str)] = &[
            ("ref_cg", "~~", "CG"),
            ("meth_cg", "C~", "CG"),
            ("pad_before", "_N", "~N"),
            ("pad_after", "N_", "N~"),
            ("ref_a_off_read", "_", "A"),
            ("junction_edge", "~,", "~~"),
            ("intron_context", ",", "C"),
            ("marker", "~,~", "~,~"),
            ("deletion", ".", "T"),
            ("inserted", "A", "."),
            ("mismatch", "/", "N"),
            ("off_contig", "~", "_"),
        ];
        const WIDE: &str = "[query.wide]\nread = \"~~~~~~~~~\"\nrefr = \"~~~~~~~~~\"\n";
        let toml = |(name, read, refr): &(&str, &str, &str)| {
            format!("[query.{name}]\nread = \"{read}\"\nrefr = \"{refr}\"\n")
        };

        // A small deterministic generator: the test must be reproducible.
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut rand = move |n: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % n
        };
        // Mostly ACGT, with the occasional N.
        let genome: String =
            (0..60).map(|_| if rand(20) == 0 { 'N' } else { b"ACGT"[rand(4) as usize] as char }).collect();
        let refr = shared_reference("independence", CONTIG, &genome);
        let map = contig_map(&refr, CONTIG, genome.len());

        let mut records = Vec::new();
        for i in 0..300 {
            let mut cigar = Vec::new();
            let (mut qlen, mut rlen) = (0u32, 0i64);
            if rand(4) == 0 { cigar.push(Cigar::HardClip(1 + rand(3) as u32)); }
            if rand(3) == 0 { let n = 1 + rand(3) as u32; cigar.push(Cigar::SoftClip(n)); qlen += n; }
            for k in 0..(1 + rand(3)) {
                if k > 0 {
                    let long = rand(3) == 0;
                    let n = 1 + rand(if long { 12 } else { 3 }) as u32;
                    match rand(3) {
                        0 => { cigar.push(Cigar::Ins(n)); qlen += n; }
                        1 => { cigar.push(Cigar::Del(n)); rlen += n as i64; }
                        _ => { cigar.push(Cigar::RefSkip(n)); rlen += n as i64; }
                    }
                }
                let n = 1 + rand(6) as u32;
                cigar.push(Cigar::Match(n));
                qlen += n;
                rlen += n as i64;
            }
            if rand(3) == 0 { let n = 1 + rand(3) as u32; cigar.push(Cigar::SoftClip(n)); qlen += n; }
            if rand(4) == 0 { cigar.push(Cigar::HardClip(1 + rand(3) as u32)); }
            if rlen >= genome.len() as i64 { continue; }
            let pos = rand((genome.len() as i64 - rlen + 1) as u64) as i64;
            let seq: Vec<u8> = (0..qlen).map(|_| b"ACGT"[rand(4) as usize]).collect();
            let flags = [0u16, REVERSE, LAST_IN_TEMPLATE | 1, LAST_IN_TEMPLATE | REVERSE | 1][i % 4];
            records.push(built(format!("r{i}").as_bytes(), &seq, &cigar, pos, flags));
        }

        let hits = |qs: &QuerySet, walk: WalkConfig, only: &str| {
            let mut m = Matcher::new(qs, Arc::clone(&refr), Arc::clone(&map), walk, crate::strand_rule::directional(), SharedConcordance::new(1.0));
            let mut out = Vec::new();
            for (r, rec) in records.iter().enumerate() {
                m.run(rec, Library::Directional.call(rec), &mut |_, q: &Query, ring: &ColumnRing| {
                    if q.name == only {
                        let a = ring.get(q.anchor_back);
                        out.push((r, a.refr_pos, a.read_off));
                    }
                })
                .unwrap();
            }
            out
        };

        let company = queries(&(PANEL.iter().map(toml).collect::<String>() + WIDE));
        for insertions in [Insertions::Skip, Insertions::Emit] {
            // Pads and intron context fixed, as a run that wants its queries
            // independent of each other sets them.
            let walk = WalkConfig { insertions, end_context: Some(8), splice_context: Some(3), ..WalkConfig::default() };
            let mut fired = Vec::new();
            for entry in PANEL {
                let alone = hits(&queries(&toml(entry)), walk, entry.0);
                let together = hits(&company, walk, entry.0);
                assert_eq!(alone, together, "{} ({insertions:?}) depends on the other queries", entry.0);
                if !alone.is_empty() {
                    fired.push(entry.0);
                }
            }
            // Guard against a vacuous pass. Nothing here runs off a contig end,
            // and skipped insertions produce no columns.
            let mut expected = vec![
                "ref_cg", "meth_cg", "pad_before", "pad_after", "ref_a_off_read",
                "junction_edge", "intron_context", "marker", "deletion", "mismatch",
            ];
            if insertions == Insertions::Emit {
                expected.push("inserted");
            }
            let missing: Vec<_> = expected.iter().filter(|n| !fired.contains(n)).collect();
            assert!(missing.is_empty(), "never fired ({insertions:?}): {missing:?}");
        }
    }

    /// An unmethylated cytosine reads as T. It still matches `Y`, and the scan
    /// reports it as T. This is what makes the read base the methylation call.
    #[test]
    fn the_read_base_distinguishes_converted_from_not() {
        let qs = queries("[query.mCpG]\nread = \"Y~\"\nrefr = \"CG\"\n");
        // Position 3 converted to T, position 8 left as C.
        let rec = built(b"r1", b"TTTTGTTTCGTTTT", &[Cigar::Match(14)], 0, 0);
        let rows = scan_to_rows("convert", &qs, &config(false), &[rec]);
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].3.as_deref(), Some("T"), "converted");
        assert_eq!(rows[1].3.as_deref(), Some("C"), "not converted");
    }

    /// `~` accepts a junction, so a pattern that wants reference context
    /// reaches into an intron exactly as it reaches past the end of a read.
    ///
    /// The CG of the genome at 3-4 crosses the boundary. The read covers the C,
    /// and the G is the first base of the intron. That is a real CG in the
    /// genome, and the evidence of the read about the C is real, so the query
    /// fires. This is the same rule that makes the flank exist. Write the read
    /// side as `{N._}` instead and the query does not fire. That spelling is
    /// what `~` meant before a junction became its own value.
    #[test]
    fn reference_context_reaches_into_an_intron() {
        // 4M 6N 4M over TTTCGTTTCGTTTT: the read covers 0-3 and 10-13.
        let cig = [Cigar::Match(4), Cigar::RefSkip(6), Cigar::Match(4)];
        let rec = || built(b"r1", b"TTTCGTTT", &cig, 0, 0);

        let loose = queries("[query.loose]\nread = \"~~\"\nrefr = \"CG\"\n");
        let rows = scan_to_rows("splice_loose", &loose, &config(false), &[rec()]);
        assert_eq!(rows.len(), 2, "the boundary CG and the one inside the context: {rows:?}");
        assert_eq!(rows[0].1, 3, "anchored on the observed C");
        assert_eq!(rows[0].3.as_deref(), Some("C"), "the read base is real evidence");
        // The derived context is the span, 2 here, so the intron's last two
        // bases (8-9) are shown and hold a CG of their own: a window made only
        // of intron context, which fires like an all-pad window does.
        assert_eq!((rows[1].1, rows[1].3.as_deref()), (8, Some(",")));

        let strict = queries("[alias]\ng = \"{N._}\"\n[query.strict]\nread = \"gg\"\nrefr = \"CG\"\n");
        let rows = scan_to_rows("splice_strict", &strict, &config(false), &[rec()]);
        assert!(rows.is_empty(), "excluding junctions excludes the boundary: {rows:?}");
    }

    /// A pattern can match the elision marker, and only a pattern that names
    /// the marker can match it.
    #[test]
    fn a_query_can_ask_for_a_junction() {
        let cig = [Cigar::Match(4), Cigar::RefSkip(6), Cigar::Match(4)];
        let qs = queries("[query.j]\nread = \"~,~\"\nrefr = \"~,~\"\n");
        // A motif pins the context. The derived one (3, the span) would cover
        // this 6-base intron completely, leaving no marker to match.
        let mut cfg = config(false);
        cfg.walk.splice_context = Some(2);

        let rows = scan_to_rows(
            "splice_named",
            &qs,
            &cfg,
            &[built(b"r1", b"TTTCGTTT", &cig, 0, 0)],
        );
        assert_eq!(rows.len(), 1, "one junction, one row: {rows:?}");

        // The same query on an unspliced read finds nothing: `,` is specific,
        // not a wildcard that happens to match gaps.
        let rows = scan_to_rows(
            "splice_none",
            &qs,
            &config(false),
            &[built(b"r1", b"TTTCGTTTCGTT", &[Cigar::Match(12)], 0, 0)],
        );
        assert!(rows.is_empty(), "no junction, no row: {rows:?}");
    }

    /// A deletion and a junction are different values, so a deletion pattern
    /// does not fire on an intron. Earlier versions gave the two the same
    /// column, which made every spliced read look like one very large deletion.
    #[test]
    fn a_deletion_pattern_does_not_match_an_intron() {
        let del = [Cigar::Match(4), Cigar::Del(6), Cigar::Match(4)];
        let skip = [Cigar::Match(4), Cigar::RefSkip(6), Cigar::Match(4)];
        let qs = queries("[query.d]\nread = \"~.~\"\nrefr = \"~C~\"\n");

        let rows = scan_to_rows(
            "gap_del",
            &qs,
            &config(false),
            &[built(b"r1", b"TTTTTTTT", &del, 0, 0)],
        );
        assert!(!rows.is_empty(), "a deletion still matches a deletion pattern");

        let rows = scan_to_rows(
            "gap_skip",
            &qs,
            &config(false),
            &[built(b"r1", b"TTTTTTTT", &skip, 0, 0)],
        );
        assert!(rows.is_empty(), "an intron is not a deletion: {rows:?}");
    }

    /// A read base outside the set of the pattern stops the query from firing
    /// at all.
    #[test]
    fn a_mismatching_read_base_does_not_fire() {
        let qs = queries("[query.mCpG]\nread = \"Y~\"\nrefr = \"CG\"\n");
        // Both cytosines read A, which is not in Y.
        let rec = built(b"r1", b"TTTAGTTTAGTTTT", &[Cigar::Match(14)], 0, 0);
        let rows = scan_to_rows("mismatch", &qs, &config(false), &[rec]);
        assert!(rows.is_empty(), "{rows:?}");
    }

    /// With `keep_hitless`, a record that fires nothing still appears in the
    /// file. The file therefore lists everything that the scan saw.
    #[test]
    fn hitless_records_are_kept_or_dropped_on_request() {
        let qs = queries("[query.mCpG]\nread = \"Y~\"\nrefr = \"CG\"\n");
        let miss = || built(b"r1", b"TTTAGTTTAGTTTT", &[Cigar::Match(14)], 0, 0);

        let kept = scan_to_rows("keep", &qs, &config(true), &[miss()]);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].0, "", "the query columns are null");
        assert_eq!(kept[0].1, -1, "so is the coordinate");

        let dropped = scan_to_rows("drop", &qs, &config(false), &[miss()]);
        assert!(dropped.is_empty());
    }

    /// The record ids come from the caller, and the sink writes them through
    /// unchanged. Overlap resolution breaks ties on them across the shards.
    #[test]
    fn record_ids_are_the_callers() {
        let qs = queries("[query.mCpG]\nread = \"Y~\"\nrefr = \"CG\"\n");
        // 'A' at position 8 keeps the second CpG from matching, so each record
        // contributes exactly one row and the ids line up one to one.
        let a = built(b"r1", b"TTTCGTTTATTTTT", &[Cigar::Match(14)], 0, 0);
        let b = built(b"r2", b"TTTCGTTTATTTTT", &[Cigar::Match(14)], 0, 0);
        let refr = shared_reference("scan_ids", CONTIG, GENOME);
        let head = header(CONTIG, GENOME.len());
        let path = temp_path("scan_ids", "parquet");
        {
            let cfg = config(false);
            let map = contig_map(&refr, CONTIG, GENOME.len());
            let mut sc = Scanner::open(&qs, refr, map, &path, &cfg, 0, SharedConcordance::new(1.0)).unwrap();
            sc.scan_with_id(41, &a, &head).unwrap();
            sc.scan_with_id(42, &b, &head).unwrap();
            sc.finish().unwrap();
        }
        let file = std::fs::File::open(&path).unwrap();
        let reader = ParquetRecordBatchReaderBuilder::try_new(file).unwrap().build().unwrap();
        let mut ids = Vec::new();
        for batch in reader {
            let batch = batch.unwrap();
            let c = batch.column_by_name("record_id").unwrap().clone();
            let c = c.as_any().downcast_ref::<UInt64Array>().unwrap();
            for i in 0..batch.num_rows() {
                ids.push(c.value(i));
            }
        }
        assert_eq!(ids, vec![41, 42]);
    }

    /// On the bottom strand the scan reports the coordinate of the cytosine as
    /// it sits on the reference. This is the position that MethylDackel calls on
    /// that strand.
    #[test]
    fn a_bottom_strand_read_reports_reference_coordinates() {
        let qs = queries("[query.mCpG]\nread = \"Y~\"\nrefr = \"CG\"\n");
        // Same bases; the strand is carried by the flag alone.
        let rec = built(b"r1", b"TTTCGTTTCGTTTT", &[Cigar::Match(14)], 0, REVERSE);
        let rows = scan_to_rows("ob", &qs, &config(false), &[rec]);
        // On the bottom strand the walk sees the reverse complement, so the CpGs
        // it finds are the ones whose C sits at the G of a forward CpG: 4 and 9.
        let mut got: Vec<i64> = rows.iter().map(|r| r.1).collect();
        got.sort_unstable();
        assert_eq!(got, vec![4, 9], "rows: {rows:?}");
    }

    /// Unless you fix it, the flank follows the set of queries. It is the span
    /// of the widest query, which is one more than a query of span `k` needs in
    /// order to reach past the alignment. The widest query therefore also meets
    /// windows that hold only pads.
    #[test]
    fn the_walk_flank_follows_the_longest_query() {
        assert_eq!(WalkConfig::default().to_opts(1).flank, 1);
        assert_eq!(WalkConfig::default().to_opts(2).flank, 2);
        assert_eq!(WalkConfig::default().to_opts(12).flank, 12);
        assert_eq!(WalkConfig::default().to_opts(0).flank, 0);
        let fixed = WalkConfig { end_context: Some(3), splice_context: Some(1), ..WalkConfig::default() };
        assert_eq!((fixed.to_opts(12).flank, fixed.to_opts(12).splice_context), (3, 1));
    }

    /// A CG whose G is the first base past the alignment still fires, because
    /// the flank supplied that column.
    #[test]
    fn a_cytosine_at_the_last_aligned_base_still_matches() {
        let qs = queries("[query.mCpG]\nread = \"Y~\"\nrefr = \"CG\"\n");
        // Align only through position 3, the C of the CpG at 3-4. The G at 4 has
        // to come from the flank.
        let rec = built(b"r1", b"TTTC", &[Cigar::Match(4)], 0, 0);
        let rows = scan_to_rows("flank", &qs, &config(false), &[rec]);
        assert_eq!(rows.len(), 1, "{rows:?}");
        assert_eq!(rows[0].1, 3);
    }
}