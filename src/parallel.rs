//! Run the scan in parallel, and give each worker its own output files.
//!
//! This is the pipeline for the parquet outputs. `query --parquet` gets its rows
//! from walking alignments ([`crate::scanner::ParquetOutput`]), and `extract`
//! gets its rows from a tag ([`crate::tag_extract::TagExtract`]). The tagged BAM
//! keeps input order instead, and it uses [`crate::ordered`].
//!
//! A [`SinkFactory`] owns what a worker does with a record and which file the
//! worker writes. The reader, the routing and the workers in this module are
//! therefore the same for every kind of output.
//!
//! One reader thread pulls the records from the BAM. It marks the records that
//! the scan cannot walk, which are the unmapped records and the records on a
//! contig that the reference does not have. It drops those records unless the
//! output [keeps them](SinkFactory::keeps_unscanned). It gives each record that
//! it keeps a `record_id` that increases over the whole run. It then routes the
//! record to an output file. To select the file it hashes a **partition key**,
//! which is a combination of record attributes that the user names with
//! `--partition-by`. The default key is the qname. The hash selects one slot out
//! of `n_workers * shards_per_worker`. The quotient names the worker that owns
//! the slot, and the remainder names which file of that worker it is. See
//! [`crate::partition`] for the arithmetic.
//!
//! Each worker owns one sink per file that it writes. For `query --parquet` that
//! sink is a [`Scanner`](crate::scanner::Scanner). Each file is named
//! `{stem}_{worker}_{shard}.{ext}` and sits beside the output path, so
//! `calls.parquet` becomes `calls_0_0.parquet` and so on. See
//! [`ShardConfig::out_path`]. Each file holds denormalized rows: the columns of
//! the record and the query that fired, side by side.
//!
//! # Why the reader routes by a record key and not by chunk
//!
//! Under the default key, all of the alignments of one fragment land in the same
//! file. Mate-overlap resolution needs this, because it must see both sides. A
//! hash of the qname guarantees it without any coordination between the workers.
//! The alignments of one fragment can be anywhere in the genome, as they often
//! are in Hi-C data, and they still meet in one file.
//!
//! If you name a different key, you lose that guarantee. In exchange you group
//! the output by something that a consumer cares about, and `calls_*.parquet`
//! then reads as a partitioned dataset: one read group per file, or one contig
//! per file, with no shuffle afterwards. Whether the exchange is safe depends on
//! what you run next, so alnbase gives a warning at the call site and does not
//! make a rule here.
//!
//! # Why a worker writes more than one file
//!
//! The number of files and the number of threads are different questions. The
//! files are how the output is partitioned for whatever reads it. The threads
//! are how fast the scan runs. `--shards-per-worker` separates the two. A scan
//! with 4 threads can still write 64 files, and more threads do not change how
//! many pieces the output arrives in. They change only which thread writes each
//! piece.
//!
//! The cost is memory. Every file is a separate sink with its own Arrow builders
//! and its own write buffer, and each sink fills to `batch_rows` rows before it
//! drains. The number of rows in memory therefore grows with the number of
//! files, and not with the number of threads. Use a smaller `--batch-rows` when
//! you write many files.
//!
//! # Two invariants that the downstream resolver needs
//!
//! **The reader assigns `record_id` before it routes the record.** The value
//! therefore increases in BAM order across all of the files. The overlap policy
//! breaks ties by `record_id`, so a counter per worker would quietly change
//! which alignment wins.
//!
//! **Each file stays sorted by position.** A subsequence of a coordinate-sorted
//! stream is still coordinate-sorted. Every `calls_{w}_{s}.parquet` therefore
//! meets the ordering requirement of the resolver on its own.
//!
//! # How to aggregate afterwards
//!
//! Under the default key, overlap resolution is correct within each file,
//! because a fragment never spans two files. Per-site *counts* are not correct
//! within each file. Several fragments can cover one position, and the reads of
//! those fragments spread across the files. You therefore have two options:
//!
//! - resolve each file, then add up the per-site counts, or
//! - give all of the files to one resolver run, which merges them by position.
//!
//! Both options are correct. The first is simple, and the second avoids a second
//! pass. In both cases the files share one schema, because the schema is a
//! function of the resolved field list alone. `calls_*.parquet` therefore reads
//! as a single dataset.

use std::path::{Path, PathBuf};
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::sync::Arc;
use std::thread;

use anyhow::{Context, Result};
use rust_htslib::bam::{HeaderView, Read, Record};

use crate::bam_io::limit_cram_decoding;
use crate::batch::{count_walkability, thread_start_error, OnBadData, Recs, Recycler, Slot, Stats, BATCH};
use crate::contig_map::ContigMap;
use crate::partition::{PartitionKey, Router};
use crate::query::QuerySet;
use crate::scanner::{walkability, Walkability};
use crate::record_field::RecordField;
use crate::aref::Aref;

/// Batches in flight per worker. Bounds memory and applies backpressure when a
/// worker falls behind.
const QUEUE_DEPTH: usize = 64;

/// Makes the sinks one kind of output writes: one per file, each opened by
/// the worker thread that owns the file.
///
/// Everything about the output lives behind this trait, so the reader, the
/// router and the workers are the same code whatever is being written.
pub trait SinkFactory: Sync {
    /// A sink borrows the compiled queries, which live as long as its worker.
    type Sink<'q>: ShardSink;

    /// Whether a record the scan cannot walk -- unmapped, or on a contig the
    /// reference lacks -- still reaches a sink, to be written as it is.
    ///
    /// An output that is a copy of its input wants them, or it would silently
    /// lose reads. Both parquet outputs want them too, and write them as a row
    /// with the query columns null, so the table accounts for every record
    /// read. An output that answers `false` never sees them, and the records it
    /// does see are numbered without gaps.
    fn keeps_unscanned(&self) -> bool;

    /// What a CRAM reader must decode for this output: required `SAM_*` bits
    /// and whether MD/NM are read. `None` decodes everything, which is what an
    /// output that writes whole records needs.
    fn cram_decoding(&self) -> Option<(u32, bool)>;

    /// Open the sink for one file. `slot` is the file's global routing slot.
    /// The reference and contig map are there when the run was given a
    /// reference; an output that walks alignments needs them, one that only
    /// reads what the records already carry does not.
    fn open<'q>(
        &self,
        queries: &'q QuerySet,
        refr: Option<Arc<Aref>>,
        contigs: Option<Arc<ContigMap>>,
        header: &HeaderView,
        path: &Path,
        slot: u32,
    ) -> Result<Self::Sink<'q>>;

    /// The run's read/reference concordance, for an output that walks
    /// alignments; checked once more when the run ends. See
    /// [`crate::scanner::SharedConcordance`].
    fn concordance(&self) -> Option<&crate::scanner::SharedConcordance> {
        None
    }
}

/// One output file's writer.
pub trait ShardSink {
    /// A record the scan walks. May modify it, e.g. to add tags.
    fn scan(&mut self, rid: u64, record: &mut Record, header: &HeaderView) -> Result<()>;

    /// A record the scan cannot walk. Only called when the factory
    /// [keeps unscanned records](SinkFactory::keeps_unscanned). Takes the
    /// same arguments as [`scan`](Self::scan) because an output that writes a
    /// row per record needs them for this record too.
    fn pass(&mut self, rid: u64, record: &mut Record, header: &HeaderView) -> Result<()>;

    /// Finish the file -- flush, footer, close -- and report what it saw.
    fn finish(self) -> Result<SinkCounts>;
}

/// What one sink reports when it finishes.
#[derive(Debug, Default, Clone, Copy)]
pub struct SinkCounts {
    pub scanned: u64,
    /// Records whose aux block could not be read in full; see
    /// [`crate::record_columns::AuxTrouble`].
    pub damaged_aux: u64,
    pub non_utf8_aux: u64,
    /// Records the run's strand rule declined to call, which are skipped
    /// rather than written on a guessed strand.
    pub unknown_strand: u64,
    /// Hits a sink could not place in its output -- for a per-base tag, an
    /// anchor past either end of the read.
    pub unplaced_hits: u64,
    /// Rows written to the file.
    pub rows: u64,
}

#[derive(Clone, Debug)]
pub struct ShardConfig {
    pub n_workers: usize,
    /// Output files each worker writes, the `N` of `{stem}_{worker}_{0..N}`.
    pub shards_per_worker: usize,
    /// Record attributes whose combined value picks the file. Already
    /// resolved — see [`RecordField::resolve`] — and non-empty; the default
    /// is applied where the arguments are read, so that
    /// [`crate::partition::DEFAULT_KEY`] has exactly one caller.
    pub partition_by: Vec<RecordField>,
    /// BGZF decompression threads for the single reader. See `--reader-threads`.
    pub reader_threads: usize,
    /// `{stem}_{worker}_{shard}.{ext}` is written beside this.
    pub out_prefix: PathBuf,
    /// What to do about a record the run cannot handle; see [`OnBadData`].
    pub bad_data: OnBadData,
    /// Refuse records on contigs whose `@SQ` line has no M5 (`--require-m5`).
    pub require_m5: bool,
    /// How to call each record's strand. The reader thread needs it for a
    /// partition key naming a strand column; the workers get the same rule
    /// through [`crate::scanner::OutputConfig`].
    pub strand: std::sync::Arc<crate::strand_rule::StrandRule>,
}

impl ShardConfig {
    /// The routing arithmetic this config implies. Fails on a zero count.
    pub fn router(&self) -> Result<Router> {
        Router::new(self.n_workers, self.shards_per_worker)
    }

    /// A fresh extractor for the key. One per reader; it is stateful.
    pub fn key(&self) -> Result<PartitionKey> {
        PartitionKey::new(&self.partition_by, self.strand.clone())
    }

    pub fn out_path(&self, worker: usize, shard: usize) -> PathBuf {
        numbered(&self.out_prefix, &[worker, shard])
    }
}

/// `foo.parquet` + [3, 1] -> `foo_3_1.parquet`; `foo` + [3, 1] -> `foo_3_1`.
///
/// Underscores rather than dots, so the indices are part of the stem and the
/// only dot in the name is the one before the extension: `hits_3_1.parquet`
/// has the extension a reader expects, where `hits.3.1.parquet` presents three
/// suffixes and invites a tool to take `.1.parquet` or `.3.1.parquet` for one.
///
/// Both indices are always written, including for a one-file-per-worker run:
/// a name that changes shape with the settings is a special case every
/// downstream glob would have to know about, and `hits_*.parquet` matches
/// either way.
fn numbered(prefix: &Path, ids: &[usize]) -> PathBuf {
    let mut tail = String::new();
    for i in ids {
        tail.push('_');
        tail.push_str(&i.to_string());
    }
    match prefix.extension().and_then(|e| e.to_str()) {
        Some(ext) => {
            let stem = prefix.file_stem().unwrap_or_default().to_string_lossy().into_owned();
            prefix.with_file_name(format!("{stem}{tail}.{ext}"))
        }
        None => {
            let name = prefix.file_name().unwrap_or_default().to_string_lossy().into_owned();
            prefix.with_file_name(format!("{name}{tail}"))
        }
    }
}

/// A run of records bound for one of a worker's files.
///
/// Batching per file rather than per worker is what lets the worker hand a
/// whole batch to one scanner: the file is named once for 512 records instead
/// of being decided per record, and each file's rows arrive in runs.
struct Batch {
    /// Index within the receiving worker's scanners, not a global slot.
    shard: usize,
    recs: Recs,
}

/// Scan `bam`, writing `shards_per_worker` files per worker with the sinks
/// `output` makes.
///
/// Opens the file itself, which tests want and production does not: the real
/// callers need to read the header before deciding how to scan, so they open
/// the reader and call [`scan_sharded_reader`]. Test-only for that reason.
#[cfg(test)]
pub fn scan_sharded<F: SinkFactory>(
    bam: &Path,
    refr: Option<Arc<Aref>>,
    queries: Arc<QuerySet>,
    cfg: &ShardConfig,
    output: &F,
) -> Result<Stats> {
    // Checked before the input is opened, as well as inside, so a bad config
    // fails without touching the BAM.
    cfg.router()?;
    cfg.key()?;

    // A single forward pass, so a pipe works as well as a file. Workers get
    // their records over a channel either way.
    //
    // This one reader feeds every worker, so its decompression rate bounds the
    // whole scan. It gets its own thread count rather than `n_workers`: those
    // threads sit inside htslib, not in the scanner pool, and the two compete
    // for the same cores.
    // Cloned because `read()` needs the reader mutably.
    let reader = crate::bam_io::open_reader(bam, cfg.reader_threads.max(1))?;
    scan_sharded_reader(reader, refr, queries, cfg, output)
}

/// The sharded scan, from a reader that is already open.
///
/// This is the entry point every subcommand uses. Taking an open reader rather
/// than a path is what lets a caller read the header first and decide how to
/// scan from what it finds — which is the only option when the BAM is arriving
/// on standard input and the header cannot be read twice.
pub fn scan_sharded_reader<F: SinkFactory>(
    mut reader: rust_htslib::bam::Reader,
    refr: Option<Arc<Aref>>,
    queries: Arc<QuerySet>,
    cfg: &ShardConfig,
    output: &F,
) -> Result<Stats> {
    // Both fail on a count of zero or an empty key, before a file is created.
    let router = cfg.router()?;
    let mut key = cfg.key()?;
    let header = reader.header().clone();

    // On a CRAM, decode only what this run reads, when the output says it
    // reads less than everything. The partition key is read by the reader
    // here, whatever the output. A BAM ignores both settings.
    if let Some((required, md_nm)) = output.cram_decoding() {
        let required = cfg.partition_by.iter().fold(required, |m, f| m | f.sam_fields());
        let md_nm = md_nm || cfg.partition_by.iter().any(|f| f.needs_md_nm());
        limit_cram_decoding(&mut reader, required, md_nm)?;
    }
    let keep_unscanned = output.keeps_unscanned();

    // Resolved once, from the header that every record will carry, and shared
    // with the workers. Built here rather than inside `Scanner` because the
    // reader below needs it too: a record on a contig the reference lacks is
    // filtered out with the unmapped ones, which keeps `record_id` dense.
    // Without a reference every mapped record is walkable: there is no contig
    // list for it to be missing from.
    let contigs = match &refr {
        Some(r) => {
            let mut c = ContigMap::build(&header, r)?;
            if cfg.require_m5 {
                c = c.require_m5(&header);
            }
            let c = Arc::new(c);
            if let Some(w) = c.warning() {
                eprintln!("{w}");
            }
            Some(c)
        }
        None => None,
    };

    // One channel per worker, not per file: backpressure is about the thread
    // that would fall behind, and a worker's files are written by one thread.
    let mut senders: Vec<SyncSender<Batch>> = Vec::with_capacity(router.n_workers());
    let mut receivers: Vec<Receiver<Batch>> = Vec::with_capacity(router.n_workers());
    for _ in 0..router.n_workers() {
        let (tx, rx) = sync_channel::<Batch>(QUEUE_DEPTH);
        senders.push(tx);
        receivers.push(rx);
    }
    // One return channel shared by every worker, sized to what can be in
    // flight at once; see [`Recycler`].
    let (returns, returned) = sync_channel::<Recs>(QUEUE_DEPTH * router.n_workers());
    let mut pool = Recycler::new(returned, None);

    let stats = thread::scope(|scope| -> Result<Stats> {
        // ---- workers -------------------------------------------------------
        let mut handles = Vec::with_capacity(router.n_workers());
        for (i, rx) in receivers.into_iter().enumerate() {
            let refr = refr.clone();
            let contigs = contigs.clone();
            let queries = Arc::clone(&queries);
            let header = header.clone();
            let returns = returns.clone();
            // Resolved out here so the worker closure carries paths rather
            // than the config, and so a bad prefix fails the same way for
            // every worker.
            let files: Vec<(PathBuf, u32)> = (0..router.shards_per_worker())
                .map(|s| (cfg.out_path(i, s), router.slot_of(i, s) as u32))
                .collect();
            let worker = thread::Builder::new().name(format!("alnbase-worker-{i}"));
            handles.push(worker.spawn_scoped(scope, move || -> Result<(SinkCounts, Vec<u64>)> {
                // Every file is opened up front, so a run that writes nothing
                // to one of them still leaves a complete, readable file rather
                // than a hole in the glob.
                let mut sinks = Vec::with_capacity(files.len());
                for (path, slot) in &files {
                    sinks.push(
                        output
                            .open(
                                &queries,
                                refr.clone(),
                                contigs.clone(),
                                &header,
                                path,
                                *slot,
                            )
                            .with_context(|| format!("worker {i} opening {}", path.display()))?,
                    );
                }
                for mut batch in rx {
                    let sink = &mut sinks[batch.shard];
                    for (slot, rec) in batch.recs.iter_mut() {
                        if slot.walkable {
                            sink.scan(slot.id, rec, &header)?;
                        } else {
                            sink.pass(slot.id, rec, &header)?;
                        }
                    }
                    // Never blocks: a full channel just drops the batch.
                    let _ = returns.try_send(batch.recs);
                }
                // Flushes and closes every file; without it a file may be
                // truncated or missing its footer.
                let mut total = SinkCounts::default();
                let mut rows = Vec::with_capacity(sinks.len());
                for sink in sinks {
                    let c = sink.finish()?;
                    total.scanned += c.scanned;
                    total.damaged_aux += c.damaged_aux;
                    total.non_utf8_aux += c.non_utf8_aux;
                    total.unknown_strand += c.unknown_strand;
                    total.unplaced_hits += c.unplaced_hits;
                    rows.push(c.rows);
                }
                Ok((total, rows))
            }).map_err(|e| thread_start_error(&format!("worker {i}"), e))?);
        }

        // ---- reader --------------------------------------------------------
        // Filtering happens here so `record_id` stays dense, and assignment
        // happens here so it stays in BAM order across all files.
        //
        // Staging is per file rather than per worker, so at most
        // `slots * BATCH` records are held back from the workers.
        let mut staging: Vec<Recs> =
            (0..router.slots()).map(|_| Vec::with_capacity(BATCH)).collect();
        let mut stats = Stats { skipped_written: keep_unscanned, ..Stats::default() };
        let mut next_rid: u64 = 0;
        // Only the workers' clones should keep the return channel open.
        drop(returns);

        // `read` into a record from the pool rather than `records()`, which
        // allocates a new one per read. A skipped record is simply read over.
        let mut rec = pool.record();
        let read_result = (|| -> Result<()> {

            while let Some(read) = reader.read(&mut rec) {
                read?;
                stats.records_read += 1;
                // The reference has no sequence for an off-reference contig, so
                // there is nothing to walk against. Deciding here rather than in
                // the worker keeps `record_id` dense and the skip counted once.
                // Without a reference (`extract`) nothing is walked, so only an
                // unmapped record is set aside.
                let why = match contigs.as_ref() {
                    Some(c) => walkability(&rec, c),
                    None if rec.is_unmapped() => Walkability::Unmapped,
                    None => Walkability::Walkable,
                };
                let walkable = count_walkability(why, &rec, &header, cfg.bad_data, &mut stats)?;
                if !walkable && !keep_unscanned {
                    continue;
                }
                // Every record read gets an id, walkable or not, so the ids are
                // the input's own numbering and the output lines up with the BAM.
                let slot = Slot { id: next_rid, walkable };
                next_rid += 1;
                // Hashed before the record is moved into staging. An unscanned
                // record hashes like any other, so it lands beside the records it
                // would group with -- its mate, under the default key.
                let dest = router.slot(key.hash(&rec, &header));
                let filled = std::mem::replace(&mut rec, pool.record());
                staging[dest].push((slot, filled));

                if staging[dest].len() >= BATCH {
                    // Unbounded, so this always returns a vector.
                    let fresh = pool.vec().expect("an unbounded pool always has a vector");
                    let recs = std::mem::replace(&mut staging[dest], fresh);
                    let (w, shard) = router.split(dest);
                    // A closed channel means that worker has already failed; stop
                    // reading and let its own error be the one reported below.
                    if senders[w].send(Batch { shard, recs }).is_err() {
                        return Ok(());
                    }
                }
            }
            for (slot, recs) in staging.into_iter().enumerate() {
                if !recs.is_empty() {
                    let (w, shard) = router.split(slot);
                    let _ = senders[w].send(Batch { shard, recs });
                }
            }
            Ok(())
        })();
        drop(senders); // closes the channels, ending the worker loops

        // A worker's own error is the root cause whenever there is one: the
        // reader only notices a failed worker as a closed channel, which says
        // nothing about why. So every worker is joined before anything is
        // returned, and a worker error wins over the reader's.
        let mut worker_error = None;
        stats.file_rows = vec![0; router.slots()];
        for (i, h) in handles.into_iter().enumerate() {
            match h.join() {
                Ok(Ok((c, rows))) => {
                    for (s, n) in rows.into_iter().enumerate() {
                        stats.file_rows[router.slot_of(i, s)] = n;
                    }
                    stats.records_scanned += c.scanned;
                    stats.damaged_aux += c.damaged_aux;
                    stats.non_utf8_aux += c.non_utf8_aux;
                    stats.unknown_strand += c.unknown_strand;
                    // Skipped, like a record the walk cannot walk, so that
                    // read = scanned + skipped holds however a record was lost.
                    stats.records_skipped += c.unknown_strand;
                    stats.unplaced_hits += c.unplaced_hits;
                }
                Ok(Err(e)) => {
                    worker_error.get_or_insert(e);
                }
                Err(_) => {
                    worker_error.get_or_insert_with(|| anyhow::anyhow!("worker {i} panicked"));
                }
            }
        }
        if let Some(e) = worker_error {
            return Err(e);
        }
        read_result?;
        // The walk checks the run's total as it goes, but a run too small, or
        // split too thinly, to reach that minimum is judged here, on everything
        // it compared.
        if let Some(c) = output.concordance() {
            stats.concordance = c.total();
            c.check_final()?;
        }
        Ok(stats)
    });

    // A failed run removes every file it may have written: a partial file set,
    // some of them valid parquet with a subset of the rows, looks like a
    // result.
    if stats.is_err() {
        remove_outputs(cfg, &router);
    }
    stats
}

/// Remove every output file a run is configured to write, ignoring any that
/// were never created. Only regular files are touched.
fn remove_outputs(cfg: &ShardConfig, router: &Router) {
    for w in 0..router.n_workers() {
        for s in 0..router.shards_per_worker() {
            let path = cfg.out_path(w, s);
            if std::fs::metadata(&path).map(|m| m.is_file()).unwrap_or(false) {
                let _ = std::fs::remove_file(&path);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dsl::QuerySpec;
    use crate::hit_writer::{OutputFormat, ParquetOpts};
    use crate::partition::DEFAULT_KEY;
    use crate::scanner::{OutputConfig, ParquetOutput, WalkConfig};
    use crate::test_support::{built, shared_reference, temp_path, write_bam, REVERSE};
    use arrow::array::{Array, StringArray, UInt32Array, UInt64Array};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use rust_htslib::bam::record::Cigar;

    const CONTIG: &str = "chr1";
    const GENOME: &str = "TTTCGTTTCGTTTT";

    /// The file names are a user-visible contract: the docs tell people to
    /// glob `hits_*.parquet`.
    #[test]
    fn files_are_numbered_before_the_extension() {
        assert_eq!(
            numbered(Path::new("hits.parquet"), &[3, 0]),
            PathBuf::from("hits_3_0.parquet")
        );
        assert_eq!(numbered(Path::new("hits"), &[0, 2]), PathBuf::from("hits_0_2"));
        assert_eq!(
            numbered(Path::new("/tmp/out/hits.parquet"), &[11, 7]),
            PathBuf::from("/tmp/out/hits_11_7.parquet"),
            "the directory is preserved"
        );
    }

    /// The routing config and the parquet settings these tests scan with. The
    /// tests below exercise the sharding through the parquet output, which
    /// reads back as a table and so says exactly where every record went.
    struct TestCfg {
        shard: ShardConfig,
        out: OutputConfig,
    }

    impl std::ops::Deref for TestCfg {
        type Target = ShardConfig;
        fn deref(&self) -> &ShardConfig {
            &self.shard
        }
    }

    fn scan_sharded(
        bam: &Path,
        refr: Arc<Aref>,
        queries: Arc<QuerySet>,
        cfg: &TestCfg,
    ) -> Result<Stats> {
        super::scan_sharded(bam, Some(refr), queries, &cfg.shard, &ParquetOutput::new(cfg.out.clone()))
    }

    fn config(n_workers: usize, per_worker: usize, tag: &str) -> TestCfg {
        keyed_config(n_workers, per_worker, DEFAULT_KEY, tag)
    }

    fn keyed_config(
        n_workers: usize,
        per_worker: usize,
        partition_by: &[RecordField],
        tag: &str,
    ) -> TestCfg {
        TestCfg {
            shard: ShardConfig {
                strand: crate::strand_rule::directional(),
                n_workers,
                shards_per_worker: per_worker,
                partition_by: partition_by.to_vec(),
                reader_threads: 1,
                out_prefix: temp_path(&format!("par_{tag}"), "parquet"),
                bad_data: crate::batch::OnBadData::Stop,
                require_m5: false,
            },
            out: OutputConfig {
                bad_data: crate::batch::OnBadData::Stop,
                fields: RecordField::CORE.to_vec(),
                batch_rows: 8,
                parquet: ParquetOpts::default(),
                keep_hitless: true,
                format: OutputFormat::Parquet,
                walk: WalkConfig::default(),
                strand: crate::strand_rule::directional(),
                },
        }
    }

    fn query_set() -> Arc<QuerySet> {
        let specs: Vec<QuerySpec> = vec![crate::test_support::query("[pattern.p1]\nread = \"Y~\"\nrefr = \"CG\"\n\n[query.m]\nmark = \"+.\"\nwhere = \"p1\"\n")];
        Arc::new(QuerySet::compile(&specs).unwrap())
    }

    /// Every file the config names, whether or not anything was routed to it.
    fn files(cfg: &ShardConfig) -> Vec<PathBuf> {
        let r = cfg.router().unwrap();
        (0..r.slots())
            .map(|slot| {
                let (w, s) = r.split(slot);
                cfg.out_path(w, s)
            })
            .collect()
    }

    /// One file's rows, as `(record_id, shard, strand)`. `strand` is a CORE
    /// column and never null on a mapped record, so it reads as a plain string.
    fn rows_of(path: &Path) -> Vec<(u64, u32, String)> {
        let file = std::fs::File::open(path)
            .unwrap_or_else(|e| panic!("opening {}: {e}", path.display()));
        let reader = ParquetRecordBatchReaderBuilder::try_new(file).unwrap().build().unwrap();
        let mut out = Vec::new();
        for batch in reader {
            let batch = batch.unwrap();
            let col = |n: &str| batch.column_by_name(n).unwrap().clone();
            let ids = col("record_id");
            let ids = ids.as_any().downcast_ref::<UInt64Array>().unwrap();
            let shard = col("shard");
            let shard = shard.as_any().downcast_ref::<UInt32Array>().unwrap();
            let strand = col("strand");
            let strand = strand.as_any().downcast_ref::<StringArray>().unwrap();
            for r in 0..batch.num_rows() {
                out.push((ids.value(r), shard.value(r), strand.value(r).to_string()));
            }
        }
        out
    }

    /// Every `record_id` written across every file.
    fn ids_across_files(cfg: &ShardConfig) -> Vec<u64> {
        files(cfg).iter().flat_map(|p| rows_of(p)).map(|r| r.0).collect()
    }

    fn reads(n: u32) -> Vec<Record> {
        (0..n)
            .map(|i| {
                built(
                    format!("read{i}").as_bytes(),
                    b"TTTCGTTTCGTTTT",
                    &[Cigar::Match(14)],
                    0,
                    0,
                )
            })
            .collect()
    }

    /// `--require-m5` stops a run at the first record on a contig whose `@SQ`
    /// line has no M5 (the test header has none), even when permissive, and the
    /// failed run leaves no files.
    #[test]
    fn a_required_m5_that_is_missing_stops_the_run() {
        let bam = write_bam("par_m5", CONTIG, GENOME.len(), &reads(4));
        let refr = shared_reference("par_m5", CONTIG, GENOME);
        let mut cfg = config(2, 1, "m5");
        cfg.shard.require_m5 = true;
        cfg.shard.bad_data = crate::batch::OnBadData::WarnAndCount;
        let err = scan_sharded(&bam, Arc::clone(&refr), query_set(), &cfg).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("--require-m5") && msg.contains(CONTIG), "{msg}");
        assert!(files(&cfg).iter().all(|f| !f.exists()), "no output is left");

        cfg.shard.require_m5 = false;
        assert!(scan_sharded(&bam, refr, query_set(), &cfg).is_ok(), "names and lengths agree");
    }

    /// Ids are assigned by the reader before routing, so they are dense and in
    /// BAM order across all files taken together -- which is what lets the
    /// overlap resolver break ties on them.
    #[test]
    fn record_ids_are_dense_across_every_file() {
        let recs = reads(20);
        let bam = write_bam("par_ids", CONTIG, GENOME.len(), &recs);
        let refr = shared_reference("par_ids", CONTIG, GENOME);

        let cfg = config(3, 1, "ids");
        let stats = scan_sharded(&bam, refr, query_set(), &cfg).unwrap();
        assert_eq!(stats.records_read, 20);
        assert_eq!(stats.records_scanned, 20);
        assert_eq!(stats.records_skipped, 0);

        let mut ids = ids_across_files(&cfg);
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids, (0..20).collect::<Vec<u64>>(), "ids must be dense and complete");
    }

    /// An unmapped record is not walked, but it is still a record of the
    /// input: it gets an all-null row and its own id, so the ids are the
    /// input's own numbering and every record is accounted for.
    #[test]
    fn unmapped_records_get_a_row_and_keep_the_ids_aligned() {
        let mut recs = reads(6);
        for (i, r) in recs.iter_mut().enumerate() {
            if i % 2 == 1 {
                r.set_unmapped();
            }
        }
        let bam = write_bam("par_unmapped", CONTIG, GENOME.len(), &recs);
        let refr = shared_reference("par_unmapped", CONTIG, GENOME);

        let cfg = config(2, 1, "unmapped");
        let stats = scan_sharded(&bam, refr, query_set(), &cfg).unwrap();
        assert_eq!(stats.records_read, 6);
        assert_eq!(stats.records_skipped, 3, "three were not walked");
        assert_eq!(stats.records_scanned, 3);

        let mut ids = ids_across_files(&cfg);
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids, (0..6).collect::<Vec<u64>>(), "every record has its input ordinal");
    }

    /// A rule built from the tables a query file would declare, for the two
    /// tests below. `unknown` is the escape a rule uses to say a record is not
    /// one it can call.
    fn rule(original: &[(&str, &str)]) -> std::sync::Arc<crate::strand_rule::StrandRule> {
        let own = |t: &[(&str, &str)]| -> Vec<(String, String)> {
            t.iter().map(|(k, v)| ((*k).to_string(), (*v).to_string())).collect()
        };
        let aligned = own(&[("forward", "not is_reverse"), ("reverse", "is_reverse")]);
        let r = crate::strand_rule::StrandRule::from_tables(&own(original), &aligned).unwrap();
        std::sync::Arc::new(r.named("test.toml"))
    }

    /// The walk takes its strand from the run's rule and not from the FLAG.
    /// These records are forward and unpaired, which the directional rule
    /// calls OT and writes on the plus strand; a rule that gives them a reverse
    /// original strand has to move every row to the minus strand, or the rule
    /// is not being read.
    #[test]
    fn the_walk_takes_its_strand_from_the_rule() {
        let bam = write_bam("par_rule", CONTIG, GENOME.len(), &reads(4));
        let refr = shared_reference("par_rule", CONTIG, GENOME);

        let mut cfg = config(2, 1, "rule");
        let inverted = rule(&[("reverse", "not is_reverse"), ("forward", "is_reverse")]);
        cfg.shard.strand = inverted.clone();
        cfg.out.strand = inverted;

        let stats = scan_sharded(&bam, refr, query_set(), &cfg).unwrap();
        assert_eq!(stats.records_scanned, 4);
        let strands: Vec<String> = files(&cfg.shard).iter().flat_map(|p| rows_of(p)).map(|r| r.2).collect();
        assert!(!strands.is_empty(), "the records were written");
        assert!(strands.iter().all(|s| s == "-"), "the rule decides the strand: {strands:?}");
    }

    /// A record the rule declines -- its `unknown` key claims it -- is skipped
    /// and counted rather than walked on a guessed strand. The count is what
    /// tells a user their rule does not cover their input.
    #[test]
    fn a_record_the_rule_declines_is_skipped_and_counted() {
        let mut recs = reads(6);
        for (i, r) in recs.iter_mut().enumerate() {
            if i % 2 == 1 {
                r.set_reverse();
            }
        }
        let bam = write_bam("par_unknown", CONTIG, GENOME.len(), &recs);
        let refr = shared_reference("par_unknown", CONTIG, GENOME);

        let mut cfg = config(2, 1, "unknown");
        let declines_reverse = rule(&[("forward", "not is_reverse"), ("unknown", "is_reverse")]);
        cfg.shard.strand = declines_reverse.clone();
        cfg.out.strand = declines_reverse;

        let stats = scan_sharded(&bam, refr, query_set(), &cfg).unwrap();
        assert_eq!(stats.records_read, 6);
        assert_eq!(stats.unknown_strand, 3, "the three reverse records were declined");
        assert_eq!(stats.records_scanned, 3);
        assert_eq!(stats.records_skipped, 3, "declined records are skipped like unwalkable ones");

        let mut ids = ids_across_files(&cfg.shard);
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids, vec![0, 2, 4], "no row is written on a strand the rule would not give");
    }

    /// `--only-hits` means only rows that matched, so it drops the unmapped
    /// rows too rather than keeping a second kind of empty row.
    #[test]
    fn only_hits_drops_the_unmapped_rows() {
        let mut recs = reads(6);
        for (i, r) in recs.iter_mut().enumerate() {
            if i % 2 == 1 {
                r.set_unmapped();
            }
        }
        let bam = write_bam("par_unmapped_only", CONTIG, GENOME.len(), &recs);
        let refr = shared_reference("par_unmapped_only", CONTIG, GENOME);
        let mut cfg = config(1, 1, "unmapped_only");
        cfg.out.keep_hitless = false;
        scan_sharded(&bam, refr, query_set(), &cfg).unwrap();

        let mut ids = ids_across_files(&cfg);
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids, vec![0, 2, 4], "only the walked records that matched");
    }

    /// A worker count of zero is refused rather than silently scanning nothing,
    /// and so is a file count of zero.
    #[test]
    fn a_zero_count_is_an_error() {
        let bam = write_bam("par_zero", CONTIG, GENOME.len(), &[]);
        let run = |cfg: &TestCfg| {
            let refr = shared_reference("par_zero", CONTIG, GENOME);
            scan_sharded(&bam, refr, query_set(), cfg).unwrap_err().to_string()
        };
        assert!(run(&config(0, 1, "zero_w")).contains("at least 1"));
        assert!(run(&config(1, 0, "zero_s")).contains("at least 1"));
    }

    /// A key naming nothing would funnel every record into one file, so it is
    /// refused rather than defaulted here -- the default belongs to the CLI.
    #[test]
    fn an_empty_partition_key_is_an_error() {
        let bam = write_bam("par_nokey", CONTIG, GENOME.len(), &[]);
        let refr = shared_reference("par_nokey", CONTIG, GENOME);
        let cfg = keyed_config(2, 2, &[], "nokey");
        let err = scan_sharded(&bam, refr, query_set(), &cfg).unwrap_err();
        assert!(err.to_string().contains("at least one field"), "{err}");
    }

    /// One worker and several must scan the same records; only the file split
    /// differs. Same for one file per worker and several.
    #[test]
    fn the_file_layout_does_not_change_what_is_scanned() {
        let recs = reads(12);
        let bam = write_bam("par_count", CONTIG, GENOME.len(), &recs);

        let mut seen = Vec::new();
        let layouts = [(1usize, 1usize, "one"), (4, 1, "four"), (2, 3, "six")];
        for (workers, per_worker, tag) in layouts {
            let refr = shared_reference(&format!("par_count_{tag}"), CONTIG, GENOME);
            let cfg = config(workers, per_worker, tag);
            let stats = scan_sharded(&bam, refr, query_set(), &cfg).unwrap();
            assert_eq!(stats.records_scanned, 12);
            let mut ids = ids_across_files(&cfg);
            ids.sort_unstable();
            ids.dedup();
            seen.push(ids);
        }
        assert_eq!(seen[0], seen[1], "the worker count changed which records were scanned");
        assert_eq!(seen[1], seen[2], "the file count changed which records were scanned");
    }

    /// Every named file exists and is a readable parquet file, including one
    /// that nothing was routed to -- a hole in the glob would be a partial
    /// dataset that reads as a complete one.
    #[test]
    fn every_file_is_written_even_when_empty() {
        // Two records into six files: most of them get nothing.
        let recs = reads(2);
        let bam = write_bam("par_empty", CONTIG, GENOME.len(), &recs);
        let refr = shared_reference("par_empty", CONTIG, GENOME);
        let cfg = config(2, 3, "empty");
        scan_sharded(&bam, refr, query_set(), &cfg).unwrap();

        let paths = files(&cfg);
        assert_eq!(paths.len(), 6, "two workers x three files each");
        let total: usize = paths.iter().map(|p| rows_of(p).len()).sum();
        assert!(total > 0, "the records went somewhere");
    }

    /// The `shard` column is the slot: which of the run's files a row is in,
    /// so a consumer reading several files off one stream can demultiplex them
    /// without the filename.
    #[test]
    fn the_shard_column_is_the_files_slot() {
        let recs = reads(24);
        let bam = write_bam("par_slot", CONTIG, GENOME.len(), &recs);
        let refr = shared_reference("par_slot", CONTIG, GENOME);
        let cfg = config(3, 2, "slot");
        scan_sharded(&bam, refr, query_set(), &cfg).unwrap();

        let r = cfg.router().unwrap();
        for slot in 0..r.slots() {
            let (w, s) = r.split(slot);
            for row in rows_of(&cfg.out_path(w, s)) {
                assert_eq!(row.1 as usize, slot, "row in {w}.{s} claims shard {}", row.1);
            }
        }
    }

    /// The point of a key: every record sharing one lands in one file.
    ///
    /// Asserted as "no key value appears in two files" rather than "these two
    /// values land in files 1 and 3", which would be a test of the digest's
    /// exact byte encoding rather than of the routing. The complementary fact,
    /// that distinct keys spread rather than collapsing onto one file, belongs
    /// to the hash and is tested there.
    #[test]
    fn records_sharing_a_key_share_a_file() {
        let recs: Vec<Record> = (0..12)
            .map(|i| {
                built(
                    format!("read{i}").as_bytes(),
                    b"TTTCGTTTCGTTTT",
                    &[Cigar::Match(14)],
                    0,
                    if i % 2 == 0 { 0 } else { REVERSE },
                )
            })
            .collect();
        let bam = write_bam("par_key", CONTIG, GENOME.len(), &recs);
        let refr = shared_reference("par_key", CONTIG, GENOME);
        // Six files, but a key with two values can only ever occupy two of them.
        let cfg = keyed_config(2, 3, &[RecordField::RefrStrand], "key");
        let stats = scan_sharded(&bam, refr, query_set(), &cfg).unwrap();
        assert_eq!(stats.records_scanned, 12);

        // Where each strand's rows were found. A record emits a row per hit,
        // so the row count is a property of the query; the record ids are what
        // says every record was written.
        let mut home: std::collections::HashMap<String, PathBuf> =
            std::collections::HashMap::new();
        let mut ids = Vec::new();
        for path in files(&cfg) {
            for (rid, _, strand) in rows_of(&path) {
                ids.push(rid);
                let first = home.entry(strand.clone()).or_insert_with(|| path.clone()).clone();
                assert_eq!(
                    first,
                    path,
                    "strand {strand} is split between {} and {}",
                    first.display(),
                    path.display()
                );
            }
        }
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids, (0..12).collect::<Vec<u64>>(), "every record was written");
        assert_eq!(home.len(), 2, "both strands were seen");
    }

    /// Records are read into buffers a worker has already used, so a short
    /// record following a long one must not carry any of the long one's name,
    /// bases or tags. Enough records to recycle several batches, of varying
    /// lengths, with a tag present on only some of them; every row is then
    /// checked against the record it claims to be.
    #[test]
    fn recycled_records_carry_nothing_over() {
        use rust_htslib::bam::record::Aux;

        const N: usize = 5 * BATCH + 37;
        let bases = b"ACGT";
        let expected = |i: usize| {
            let len = 4 + (i * 7) % 11; // 4..=14, never the same twice running
            let seq: Vec<u8> = (0..len).map(|j| bases[(i + j) % 4]).collect();
            let tag = (i % 3 == 0).then(|| format!("xb{i}"));
            (format!("read{i}_{}", "n".repeat(i % 13)), seq, tag)
        };
        let recs: Vec<Record> = (0..N)
            .map(|i| {
                let (name, seq, tag) = expected(i);
                let mut r =
                    built(name.as_bytes(), &seq, &[Cigar::Match(seq.len() as u32)], 0, 0);
                r.push_aux(b"RG", Aux::String(["g0", "g1"][i % 2])).unwrap();
                if let Some(t) = &tag {
                    r.push_aux(b"XB", Aux::String(t)).unwrap();
                }
                r
            })
            .collect();
        let bam = write_bam("par_recycle", CONTIG, GENOME.len(), &recs);
        let refr = shared_reference("par_recycle", CONTIG, GENOME);

        let xb = RecordField::Aux(*b"XB", crate::record_field::AuxKind::Str);
        let rg = RecordField::Aux(*b"RG", crate::record_field::AuxKind::Str);
        let mut cfg = keyed_config(2, 2, &[rg], "recycle");
        cfg.out.fields = vec![RecordField::Qname, RecordField::SeqAscii, xb];
        cfg.out.batch_rows = 1000;
        // The reads are synthetic sequence, not drawn from GENOME.
        cfg.out.walk.max_discordance = 1.0;
        let stats = scan_sharded(&bam, refr, query_set(), &cfg).unwrap();
        assert_eq!(stats.records_scanned, N as u64);

        let mut seen = vec![false; N];
        for path in files(&cfg) {
            let file = std::fs::File::open(&path).unwrap();
            let reader =
                ParquetRecordBatchReaderBuilder::try_new(file).unwrap().build().unwrap();
            for batch in reader {
                let batch = batch.unwrap();
                let col = |n: String| batch.column_by_name(&n).unwrap().clone();
                let ids = col("record_id".into());
                let ids = ids.as_any().downcast_ref::<UInt64Array>().unwrap();
                let names = col(RecordField::Qname.name());
                let names = names.as_any().downcast_ref::<StringArray>().unwrap();
                let seqs = col(RecordField::SeqAscii.name());
                let seqs = seqs.as_any().downcast_ref::<StringArray>().unwrap();
                let tags = col(xb.name());
                let tags = tags.as_any().downcast_ref::<StringArray>().unwrap();
                for r in 0..batch.num_rows() {
                    // Ids are dense in BAM order and every record is mapped,
                    // so the id is the record's index.
                    let i = ids.value(r) as usize;
                    let (name, seq, tag) = expected(i);
                    assert_eq!(names.value(r), name, "record {i}");
                    assert_eq!(seqs.value(r).as_bytes(), &seq[..], "record {i}");
                    let got = (!tags.is_null(r)).then(|| tags.value(r).to_string());
                    assert_eq!(got, tag, "record {i}");
                    seen[i] = true;
                }
            }
        }
        assert!(seen.iter().all(|s| *s), "every record was written");
    }

    /// Every column of every file, rendered as text in file order, for
    /// comparing two runs cell by cell.
    fn rendered(cfg: &ShardConfig) -> Vec<String> {
        use arrow::util::display::{ArrayFormatter, FormatOptions};
        let mut out = Vec::new();
        for path in files(cfg) {
            let file = std::fs::File::open(&path).unwrap();
            let reader =
                ParquetRecordBatchReaderBuilder::try_new(file).unwrap().build().unwrap();
            for batch in reader {
                let batch = batch.unwrap();
                let opts = FormatOptions::default().with_null("<null>");
                for (field, col) in batch.schema().fields().iter().zip(batch.columns()) {
                    let fmt = ArrayFormatter::try_new(col.as_ref(), &opts).unwrap();
                    for r in 0..batch.num_rows() {
                        out.push(format!("{}[{r}]={}", field.name(), fmt.value(r)));
                    }
                }
            }
        }
        out
    }

    /// A CRAM decodes only the fields a run asks for, so the one thing keeping
    /// a column correct on a CRAM is `RecordField::sam_fields` naming what that
    /// column reads. Selecting one field at a time makes each field's arm the
    /// only thing standing between its column and a placeholder value, and the
    /// same records as a BAM, where nothing is skipped, are the answer.
    ///
    /// The reads are proper pairs whose mates share a slice, which is where
    /// CRAM derives flags, mate position and TLEN rather than storing them,
    /// and they carry correct `MD`/`NM` tags, which CRAM drops on write and
    /// regenerates on read only when asked.
    #[test]
    fn a_cram_scans_exactly_like_the_same_bam() {
        use crate::record_field::AuxKind;
        use rust_htslib::bam::header::HeaderRecord;
        use rust_htslib::bam::record::Aux;
        use rust_htslib::bam::{Format, Header, Writer};

        let tag = "par_cram";
        let refr = shared_reference(tag, CONTIG, GENOME);
        let fasta = crate::test_support::temp_path(tag, "fa");

        // A CRAM reader looks its reference up by MD5 along REF_PATH before
        // trying the header's UR, and REF_PATH defaults to EBI's server. An
        // empty local directory makes that lookup miss at once, so the test
        // uses the UR below and never touches the network.
        //
        // SAFETY: the variable is set before any CRAM is opened, and no other
        // test reads it or opens a CRAM, so nothing reads it concurrently.
        let no_refs = std::env::temp_dir().join("alnbase_test_no_ref_path");
        std::fs::create_dir_all(&no_refs).unwrap();
        unsafe { std::env::set_var("REF_PATH", format!("{}/%s", no_refs.display())) };

        // (name, pos, seq, flags, mate pos, tlen, NM, MD). Mate 1 matches the
        // reference; mate 2 carries one mismatch so NM and MD are not trivial.
        const P: u16 = 0x1 | 0x2;
        let mut recs = Vec::new();
        for i in 0..40usize {
            let name = format!("frag{i}");
            let fwd_first = i % 2 == 0;
            let (f1, f2) = if fwd_first {
                (P | 0x20 | 0x40, P | 0x10 | 0x80)
            } else {
                (P | 0x10 | 0x80, P | 0x20 | 0x40)
            };
            let mates: [(i64, &[u8], u16, i64, i64, i32, &str); 2] = [
                (0, b"TTTCGT", f1, 4, 10, 0, "6"),
                (4, b"GTATCG", f2, 0, -10, 1, "2T3"),
            ];
            for (pos, seq, flags, mpos, tlen, nm, md) in mates {
                let mut r = built(name.as_bytes(), seq, &[Cigar::Match(6)], pos, flags);
                r.set_mtid(0);
                r.set_mpos(mpos);
                r.set_insert_size(tlen);
                r.set_mapq(10 + i as u8);
                r.push_aux(b"NM", Aux::I32(nm)).unwrap();
                r.push_aux(b"MD", Aux::String(md)).unwrap();
                r.push_aux(b"RG", Aux::String(["g0", "g1"][i % 2])).unwrap();
                if i % 3 != 0 {
                    r.push_aux(b"XB", Aux::String(&format!("xb{i}"))).unwrap();
                }
                recs.push(r);
            }
        }

        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", CONTIG);
        sq.push_tag(b"LN", GENOME.len());
        // rust-htslib writes the header as the writer is created, before
        // `set_reference` can run, so the reference is named here instead.
        // Without it htslib embeds a consensus of the reads, which can store
        // MD/NM verbatim rather than regenerating them.
        sq.push_tag(b"UR", fasta.to_str().unwrap());
        header.push_record(&sq);
        let mut rg = HeaderRecord::new(b"RG");
        rg.push_tag(b"ID", "g0");
        header.push_record(&rg);
        let mut rg = HeaderRecord::new(b"RG");
        rg.push_tag(b"ID", "g1");
        header.push_record(&rg);

        let bam = crate::test_support::temp_path(tag, "bam");
        let cram = crate::test_support::temp_path(tag, "cram");
        {
            let mut w = Writer::from_path(&bam, &header, Format::Bam).unwrap();
            for r in &recs {
                w.write(r).unwrap();
            }
            let mut w = Writer::from_path(&cram, &header, Format::Cram).unwrap();
            w.set_reference(&fasta).unwrap();
            for r in &recs {
                w.write(r).unwrap();
            }
        }

        let mut selections: Vec<Vec<RecordField>> =
            RecordField::ALL.iter().map(|f| vec![*f]).collect();
        for (t, k) in [(b"NM", AuxKind::Int), (b"MD", AuxKind::Str), (b"RG", AuxKind::Str), (b"XB", AuxKind::Str)] {
            selections.push(vec![RecordField::Aux(*t, k)]);
        }
        selections.push(vec![RecordField::Aux(*b"XB", AuxKind::Str), RecordField::Aux(*b"RG", AuxKind::Str)]);

        let mut wrong = Vec::new();
        for (n, fields) in selections.iter().enumerate() {
            let run = |input: &Path, which: &str| {
                // A flag-only key, so the key asks for nothing the field does not.
                let mut cfg = keyed_config(1, 1, &[RecordField::RefrStrand], &format!("cram_{which}_{n}"));
                cfg.out.fields = fields.clone();
                scan_sharded(input, Arc::clone(&refr), query_set(), &cfg).unwrap();
                rendered(&cfg)
            };
            let from_bam = run(&bam, "bam");
            let from_cram = run(&cram, "cram");
            let label = fields.iter().map(|f| f.name()).collect::<Vec<_>>().join(",");
            assert!(!from_bam.is_empty(), "{label}: nothing written");
            if let Some((c, b)) = from_cram.iter().zip(&from_bam).find(|(c, b)| c != b) {
                wrong.push(format!("{label}: CRAM has {c}, BAM has {b}"));
            } else if from_cram.len() != from_bam.len() {
                wrong.push(format!("{label}: {} cells from CRAM, {} from BAM", from_cram.len(), from_bam.len()));
            }
        }
        // Every disagreeing selection at once, rather than the first.
        assert!(wrong.is_empty(), "CRAM and BAM disagree:\n{}", wrong.join("\n"));
    }
}