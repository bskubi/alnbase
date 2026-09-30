//! Run the scan with one reader, a pool of taggers, and one writer. The output
//! keeps the exact order of the input.
//!
//! ```text
//! reader ──(numbered batches)──▶ taggers ──(tagged batches)──▶ writer
//!    ▲                                                           │
//!    └───────────────────── emptied batches ─────────────────────┘
//! ```
//!
//! The reader cuts the input into consecutive batches and numbers them. Any
//! tagger takes the next batch, tags its records in place and passes it on.
//! The writer holds batches that arrive early and writes strictly in number
//! order, then hands each written batch back to the reader, whose records are
//! read into again.
//!
//! # Why nothing can block forever
//!
//! A fixed number of batch vectors exists, `threads * INFLIGHT_PER_WORKER`,
//! and the reader cannot start a batch without one. It gets a vector back only
//! when the writer has written a batch. So at most that many batches are ever
//! in flight, and every channel is sized to hold them all: no send ever waits
//! for room. The only waits are for work -- a tagger for a batch, the writer
//! for the next number, the reader for a vector -- and each is for something
//! already in the pipeline. The reader sends a full batch before waiting for
//! a vector, so the batch the writer needs next is never the one held back.
//!
//! The bound is also the memory bound: a run holds at most that many batches'
//! records, and allocates none after warming up.
//!
//! # Stopping
//!
//! A failing stage has to stop the others, and a stage blocked waiting has to
//! notice:
//!
//! - **A tagger fails** (or panics): it tells the writer, which stops. With
//!   the writer gone, the reader's wait for a vector and the other taggers'
//!   sends fail, and they stop too.
//! - **The writer fails**: its channels close, with the same effect.
//! - **The reader fails**: its channel closes, the taggers finish and stop,
//!   and the writer sees that input ended without the reader saying so.
//!
//! Every stage that stops because of another reports [`Stopped`], and the
//! error returned is the first that is not. The file written so far is then
//! removed, rather than left looking like a complete BAM of fewer records --
//! htslib writes the end-of-file marker however a file is closed.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{sync_channel, SyncSender};
use std::sync::{Arc, Mutex};
use std::thread;

use anyhow::{anyhow, Result};
use rust_htslib::bam::Read;

use crate::bam_io::open_reader;
use crate::bam_out::{BamStream, TagCounts};
use crate::batch::{count_walkability, thread_start_error, OnBadData, Recs, Recycler, Slot, Stats, BATCH};
use crate::contig_map::ContigMap;
use crate::query::QuerySet;
use crate::scanner::walkability;
use crate::aref::Aref;

/// Batches in flight per tagger. Enough slack that one slow batch does not
/// idle the others at once; small enough that memory stays a few batches per
/// thread.
const INFLIGHT_PER_WORKER: usize = 8;

#[derive(Debug, Clone)]
pub struct OrderedConfig {
    /// Tagging threads.
    pub n_workers: usize,
    /// BGZF decompression threads for the reader.
    pub reader_threads: usize,
    /// BGZF compression threads for the writer.
    pub writer_threads: usize,
    /// The output BAM, or `-` for standard output.
    pub out: PathBuf,
    /// What to do about a record the run cannot handle; see [`OnBadData`].
    pub bad_data: OnBadData,
    /// Refuse records on contigs whose `@SQ` line has no M5 (`--require-m5`).
    pub require_m5: bool,
}

/// A batch and its place in the input.
struct Numbered {
    seq: u64,
    recs: Recs,
}

/// What taggers send the writer.
enum Done {
    Batch(Numbered),
    /// A tagger failed; stop.
    Failed,
}

/// The error of a stage that stopped because another stage failed. Never the
/// error a run reports, unless nothing better is known.
#[derive(Debug)]
pub struct Stopped;

impl std::fmt::Display for Stopped {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("stopped because another stage of the scan failed")
    }
}

impl std::error::Error for Stopped {}

fn stopped() -> anyhow::Error {
    anyhow::Error::new(Stopped)
}

/// Tells the writer to stop if dropped while its thread is panicking, so a
/// tagger that panics cannot leave the writer waiting for a batch that will
/// never come.
struct PanicAlarm<'a> {
    done: &'a SyncSender<Done>,
}

impl Drop for PanicAlarm<'_> {
    fn drop(&mut self) {
        if thread::panicking() {
            // `try_send`: never block while unwinding. There is always room,
            // since the channel holds every batch plus one message per tagger.
            let _ = self.done.try_send(Done::Failed);
        }
    }
}

/// Scan `bam`, writing every record, tagged, to `cfg.out` in input order.
pub fn scan_ordered(
    bam: &Path,
    refr: Arc<Aref>,
    queries: Arc<QuerySet>,
    cfg: &OrderedConfig,
    output: &BamStream,
) -> Result<Stats> {
    let n_workers = cfg.n_workers.max(1);

    // Whole records are written, so a CRAM is decoded in full: no
    // `limit_cram_decoding` here.
    let mut reader = open_reader(bam, cfg.reader_threads.max(1))?;
    let header = reader.header().clone();
    let mut contigs = ContigMap::build(&header, &refr)?;
    if cfg.require_m5 {
        contigs = contigs.require_m5(&header);
    }
    let contigs = Arc::new(contigs);
    if let Some(w) = contigs.warning() {
        eprintln!("{w}");
    }
    // Created before any thread starts, so a bad path fails at once.
    let writer = output.writer(&cfg.out, &header, cfg.writer_threads)?;

    let limit = n_workers * INFLIGHT_PER_WORKER;
    let (work_tx, work_rx) = sync_channel::<Numbered>(limit);
    let work_rx = Arc::new(Mutex::new(work_rx));
    let (done_tx, done_rx) = sync_channel::<Done>(limit + n_workers);
    let (return_tx, return_rx) = sync_channel::<Recs>(limit);
    let mut pool = Recycler::new(return_rx, Some(limit));
    // Set by the reader once it has read to the end, so the writer can tell a
    // finished input from an abandoned one.
    let read_all = AtomicBool::new(false);

    let (read, tagged, written) = thread::scope(|scope| {
        // ---- writer --------------------------------------------------------
        let read_all = &read_all;
        let writer_thread = thread::Builder::new().name("alnbase-writer".into()).spawn_scoped(scope, move || -> Result<()> {
            let mut writer = writer;
            let mut next = 0u64;
            let mut early: BTreeMap<u64, Recs> = BTreeMap::new();
            for msg in done_rx {
                let Numbered { seq, recs } = match msg {
                    Done::Batch(b) => b,
                    Done::Failed => return Err(stopped()),
                };
                early.insert(seq, recs);
                while let Some(recs) = early.remove(&next) {
                    for (_, rec) in &recs {
                        writer.write(rec)?;
                    }
                    next += 1;
                    // Never waits: no more vectors exist than the channel
                    // holds. Fails only once the reader is done, which is fine.
                    let _ = return_tx.send(recs);
                }
            }
            // Every tagger has stopped. That is the end of the input only if
            // the reader said so and nothing is still waiting its turn.
            if !read_all.load(Ordering::SeqCst) || !early.is_empty() {
                return Err(stopped());
            }
            writer.close()
        });
        let writer_thread = match writer_thread {
            Ok(h) => h,
            // Nothing else has started, so there is nothing to stop.
            Err(e) => return (Err(thread_start_error("writer", e)), Vec::new(), Ok(Ok(()))),
        };

        // ---- taggers -------------------------------------------------------
        let mut taggers = Vec::with_capacity(n_workers);
        let mut spawn_error = None;
        for i in 0..n_workers {
            let work_rx = Arc::clone(&work_rx);
            let done = done_tx.clone();
            let queries = Arc::clone(&queries);
            let refr = Arc::clone(&refr);
            let contigs = Arc::clone(&contigs);
            let tagger = thread::Builder::new().name(format!("alnbase-tagger-{i}"));
            let spawned = tagger.spawn_scoped(scope, move || -> Result<TagCounts> {
                let _alarm = PanicAlarm { done: &done };
                let mut tagger = output.tagger(&queries, refr, contigs);
                loop {
                    let next = work_rx.lock().unwrap_or_else(|p| p.into_inner()).recv();
                    let Ok(mut batch) = next else { break };
                    for (slot, rec) in batch.recs.iter_mut() {
                        if !slot.walkable {
                            continue;
                        }
                        if let Err(e) = tagger.tag(rec) {
                            let _ = done.try_send(Done::Failed);
                            return Err(e);
                        }
                    }
                    if done.send(Done::Batch(batch)).is_err() {
                        return Err(stopped());
                    }
                }
                Ok(tagger.counts())
            });
            match spawned {
                Ok(h) => taggers.push(h),
                Err(e) => {
                    // Reading nothing ends the taggers already started, and
                    // the writer then stops without the reader's all-clear.
                    spawn_error = Some(thread_start_error(&format!("tagger {i}"), e));
                    break;
                }
            }
        }
        // Only the threads' copies may keep these channels open.
        drop(done_tx);
        drop(work_rx);

        // ---- reader --------------------------------------------------------
        let read = spawn_error.map(Err).unwrap_or_else(|| (|| -> Result<Stats> {
            let mut stats = Stats { skipped_written: true, ..Stats::default() };
            let mut seq = 0u64;
            let mut next_id = 0u64;
            let mut batch = pool.vec().ok_or_else(stopped)?;
            let mut rec = pool.record();
            while let Some(r) = reader.read(&mut rec) {
                r?;
                stats.records_read += 1;
                // The reference has no sequence for an off-reference contig,
                // so there is nothing to walk against. Such records, and
                // unmapped ones, travel untagged and keep their place.
                let why = walkability(&rec, &contigs);
                let walkable = count_walkability(why, &rec, &header, cfg.bad_data, &mut stats)?;
                let slot = Slot { id: next_id, walkable };
                next_id += 1;
                let filled = std::mem::replace(&mut rec, pool.record());
                batch.push((slot, filled));

                if batch.len() == BATCH {
                    // Send before waiting for a vector: the batch the writer
                    // needs next may be this one.
                    let full = Numbered { seq, recs: std::mem::take(&mut batch) };
                    work_tx.send(full).map_err(|_| stopped())?;
                    seq += 1;
                    batch = pool.vec().ok_or_else(stopped)?;
                }
            }
            if !batch.is_empty() {
                work_tx.send(Numbered { seq, recs: batch }).map_err(|_| stopped())?;
            }
            read_all.store(true, Ordering::SeqCst);
            Ok(stats)
        })());
        // Ends the taggers once they have drained what was sent.
        drop(work_tx);
        // Returned vectors are no longer needed; the writer's sends may fail.
        drop(pool);

        let tagged: Vec<thread::Result<Result<TagCounts>>> =
            taggers.into_iter().map(|h| h.join()).collect();
        let written = writer_thread.join();
        (read, tagged, written)
    });

    // ---- outcome -----------------------------------------------------------
    let mut errors: Vec<anyhow::Error> = Vec::new();
    let mut stats = match read {
        Ok(s) => Some(s),
        Err(e) => {
            errors.push(e);
            None
        }
    };
    let mut counts = TagCounts::default();
    for t in tagged {
        match t {
            Ok(Ok(c)) => {
                counts.scanned += c.scanned;
                counts.unplaced_hits += c.unplaced_hits;
                counts.damaged_aux += c.damaged_aux;
                counts.unknown_strand += c.unknown_strand;
            }
            Ok(Err(e)) => errors.push(e),
            Err(_) => errors.push(anyhow!("a tagging thread panicked")),
        }
    }
    match written {
        Ok(Ok(())) => {}
        Ok(Err(e)) => errors.push(e),
        Err(_) => errors.push(anyhow!("the writing thread panicked")),
    }

    // Judged on the whole run once it has ended, like the parquet outputs.
    if errors.is_empty() {
        if let Err(e) = output.concordance().check_final() {
            errors.push(e);
        }
    }
    if !errors.is_empty() {
        remove_partial(&cfg.out);
        let first = errors.iter().position(|e| !e.is::<Stopped>()).unwrap_or(0);
        return Err(errors.swap_remove(first));
    }
    let stats = stats.as_mut().expect("no errors, so the reader finished");
    stats.records_scanned = counts.scanned;
    stats.unplaced_hits = counts.unplaced_hits;
    stats.damaged_aux = counts.damaged_aux;
    // A record the rule declines is written through untagged, which is what a
    // record the walk cannot walk gets too, so it is counted the same way:
    // read = tagged + skipped still holds.
    stats.unknown_strand = counts.unknown_strand;
    stats.records_skipped += counts.unknown_strand;
    stats.concordance = output.concordance().total();
    Ok(stats.clone())
}

/// Remove an output file a failed run left behind. Only a regular file: not
/// standard output, and not a pipe or device the user named.
fn remove_partial(out: &Path) {
    if out.as_os_str() == "-" {
        return;
    }
    if std::fs::metadata(out).map(|m| m.is_file()).unwrap_or(false) {
        match std::fs::remove_file(out) {
            Ok(()) => eprintln!("removed the incomplete output {}", out.display()),
            Err(e) => eprintln!("could not remove the incomplete output {}: {e}", out.display()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam_out::BamOptions;
    use crate::query_toml;
    use crate::scanner::WalkConfig;
    use crate::test_support::{built, shared_reference, temp_path, write_bam};
    use rust_htslib::bam::record::{Aux, Cigar};
    use rust_htslib::bam::{Reader, Record};

    const CONTIG: &str = "chr1";
    const GENOME: &str = "TTTCGTTTCGTTTT";
    const UNMAPPED: u16 = 0x4;

    /// XM for CpGs.
    const TOML: &str = r#"
[query.TG]
mark = "+."
where = "tg"
[query.CG]
mark = "+."
where = "cg"
[pattern.tg]
read = "T~"
refr = "CG"
[pattern.cg]
read = "C~"
refr = "CG"
[tag.XM.bases]
fill = "."
z = "TG"
Z = "CG"
"#;

    fn stream(toml: &str) -> (Arc<QuerySet>, BamStream) {
        let file = query_toml::parse_file(toml).unwrap_or_else(|e| panic!("{e}"));
        let queries = QuerySet::compile(&file.queries).unwrap();
        let opts = BamOptions {
            walk: WalkConfig::default(),
            strand: crate::strand_rule::directional(),
            overwrite_tags: false,
            bad_data: crate::batch::OnBadData::Stop,
            command_line: "test".into(),
            query_sources: Vec::new(),
        };
        let s = BamStream::new(&file.tags, &queries, opts).unwrap();
        (Arc::new(queries), s)
    }

    fn run(tag: &str, recs: &[Record], workers: usize) -> (Result<Stats>, PathBuf) {
        let (queries, output) = stream(TOML);
        let bam = write_bam(tag, CONTIG, GENOME.len(), recs);
        let refr = shared_reference(tag, CONTIG, GENOME);
        let out = temp_path(&format!("{tag}_out"), "bam");
        let _ = std::fs::remove_file(&out);
        let cfg = OrderedConfig {
            n_workers: workers,
            reader_threads: 1,
            writer_threads: 2,
            out: out.clone(),
            bad_data: crate::batch::OnBadData::Stop,
            require_m5: false,
        };
        (scan_ordered(&bam, refr, queries, &cfg, &output), out)
    }

    fn read_back(path: &Path) -> Vec<(String, Option<String>)> {
        let mut reader = Reader::from_path(path).unwrap();
        reader
            .records()
            .map(|r| {
                let r = r.unwrap();
                let xm = match r.aux(b"XM") {
                    Ok(Aux::String(s)) => Some(s.to_string()),
                    _ => None,
                };
                (String::from_utf8_lossy(r.qname()).into_owned(), xm)
            })
            .collect()
    }

    /// Enough records for many batches across several taggers, with read
    /// lengths varying so that recycled records are reused at other sizes and
    /// unmapped records scattered through: the output is the input, in order,
    /// every mapped record tagged for its own sequence and nothing else.
    #[test]
    fn the_output_is_the_input_in_order_with_tags() {
        let n = 7 * BATCH + 123;
        let seq_for = |i: usize| -> (Vec<u8>, String) {
            // Alternate the call at the first CpG, and vary the length so the
            // second CpG's C is sometimes off the end. When the read ends on
            // that C, the walk's flank still supplies the reference G after
            // it, so the call is made.
            let first = if i % 2 == 0 { b'T' } else { b'C' };
            let len = 6 + i % 9; // 6..=14
            let mut seq: Vec<u8> = GENOME.as_bytes()[..len].to_vec();
            seq[3] = first;
            let mut xm = vec![b'.'; len];
            xm[3] = if first == b'T' { b'z' } else { b'Z' };
            if len >= 9 {
                xm[8] = b'Z';
            }
            (seq, String::from_utf8(xm).unwrap())
        };
        let recs: Vec<Record> = (0..n)
            .map(|i| {
                let (seq, _) = seq_for(i);
                let flags = if i % 11 == 5 { UNMAPPED } else { 0 };
                built(format!("r{i:05}").as_bytes(), &seq, &[Cigar::Match(seq.len() as u32)], 0, flags)
            })
            .collect();

        let (stats, out) = run("ordered_order", &recs, 3);
        let stats = stats.unwrap();
        let got = read_back(&out);
        assert_eq!(got.len(), n);
        for (i, (name, xm)) in got.iter().enumerate() {
            assert_eq!(name, &format!("r{i:05}"), "record {i} out of order");
            if i % 11 == 5 {
                assert_eq!(xm, &None, "unmapped {name} is untagged");
            } else {
                assert_eq!(xm.as_deref(), Some(seq_for(i).1.as_str()), "{name}");
            }
        }
        let unmapped = (0..n).filter(|i| i % 11 == 5).count() as u64;
        assert_eq!(stats.records_read, n as u64);
        assert_eq!(stats.records_skipped, unmapped);
        assert_eq!(stats.records_scanned, n as u64 - unmapped);
    }

    #[test]
    fn an_empty_input_writes_an_empty_bam() {
        let (stats, out) = run("ordered_empty", &[], 2);
        assert_eq!(stats.unwrap().records_read, 0);
        assert!(read_back(&out).is_empty());
    }

    /// A tagger failing far into the run -- many batches already written,
    /// more in flight on other threads -- stops everything without hanging,
    /// reports the tagger's error rather than a knock-on one, and removes the
    /// partial file.
    #[test]
    fn a_failure_mid_run_stops_every_stage_and_removes_the_output() {
        let n = 20 * BATCH;
        let bad = 13 * BATCH + 7;
        let recs: Vec<Record> = (0..n)
            .map(|i| {
                let mut r = built(format!("r{i:05}").as_bytes(), GENOME.as_bytes(), &[Cigar::Match(14)], 0, 0);
                if i == bad {
                    // Already tagged, and overwriting is off: tagging fails.
                    r.push_aux(b"XM", Aux::String("stale")).unwrap();
                }
                r
            })
            .collect();
        for workers in [1, 4] {
            let tag = format!("ordered_fail_{workers}");
            let (res, out) = run(&tag, &recs, workers);
            let err = res.expect_err("the run fails");
            assert!(!err.is::<Stopped>(), "the real error is reported: {err:#}");
            let msg = format!("{err:#}");
            assert!(msg.contains(&format!("record r{bad:05} already carries tag XM")), "{msg}");
            assert!(!out.exists(), "the partial output is removed");
        }
    }

    /// A truncated input ends the run with the reader's error, not a
    /// knock-on one, and leaves no output that looks complete.
    #[test]
    fn a_truncated_input_is_an_error_and_leaves_no_output() {
        let recs: Vec<Record> = (0..4 * BATCH)
            .map(|i| built(format!("r{i:05}").as_bytes(), GENOME.as_bytes(), &[Cigar::Match(14)], 0, 0))
            .collect();
        let (queries, output) = stream(TOML);
        let bam = write_bam("ordered_trunc", CONTIG, GENOME.len(), &recs);
        let bytes = std::fs::read(&bam).unwrap();
        let cut = temp_path("ordered_trunc_cut", "bam");
        std::fs::write(&cut, &bytes[..bytes.len() * 2 / 3]).unwrap();
        let refr = shared_reference("ordered_trunc", CONTIG, GENOME);
        let out = temp_path("ordered_trunc_out", "bam");
        let _ = std::fs::remove_file(&out);
        let cfg = OrderedConfig {
            n_workers: 2,
            reader_threads: 1,
            writer_threads: 1,
            out: out.clone(),
            bad_data: crate::batch::OnBadData::Stop,
            require_m5: false,
        };
        let err = scan_ordered(&cut, refr, queries, &cfg, &output).expect_err("truncated input");
        assert!(!err.is::<Stopped>(), "{err:#}");
        assert!(!out.exists());
    }
    /// A record on a contig the reference does not have stops the run, names
    /// the contig, and leaves no output behind.
    #[test]
    fn a_record_off_the_reference_stops_the_run() {
        let (queries, output) = stream(TOML);
        let m = [Cigar::Match(6)];
        let recs = vec![
            built(b"ok", b"TTTCGT", &m, 0, 0),
            {
                let mut r = built(b"elsewhere", b"TTTCGT", &m, 0, 0);
                r.set_tid(1);
                r
            },
        ];
        let bam = crate::test_support::write_bam_many(
            "ordered_offref",
            &[(CONTIG, GENOME.len()), ("chrUn", 100)],
            &recs,
        );
        let refr = shared_reference("ordered_offref", CONTIG, GENOME);
        let out = temp_path("ordered_offref_out", "bam");
        let _ = std::fs::remove_file(&out);
        let cfg = OrderedConfig {
            n_workers: 1,
            reader_threads: 1,
            writer_threads: 1,
            out: out.clone(),
            bad_data: crate::batch::OnBadData::Stop,
            require_m5: false,
        };
        let err = scan_ordered(&bam, Arc::clone(&refr), Arc::clone(&queries), &cfg, &output)
            .expect_err("an off-reference record stops the run");
        let msg = format!("{err:#}");
        assert!(msg.contains("record elsewhere") && msg.contains("chrUn"), "{msg}");
        assert!(msg.contains("--permissive"), "the message says how to proceed: {msg}");
        assert!(!out.exists(), "no partial output is left");

        // With the flag: written untagged, counted, and the run succeeds.
        let cfg = OrderedConfig { bad_data: crate::batch::OnBadData::WarnAndCount, ..cfg };
        let stats = scan_ordered(&bam, refr, queries, &cfg, &output).unwrap();
        assert_eq!((stats.records_read, stats.records_off_reference), (2, 1));
        let names: Vec<String> = read_back(&out).into_iter().map(|(n, _)| n).collect();
        assert_eq!(names, ["ok", "elsewhere"], "both records are still written");
    }

    /// A header listing contigs the reference lacks is fine on its own: only a
    /// record actually on one is a problem. This is the whole-genome-header,
    /// single-chromosome-reference case, which must work without the flag.
    #[test]
    fn a_header_contig_nothing_aligns_to_is_fine() {
        let (queries, output) = stream(TOML);
        let recs = vec![built(b"a", b"TTTCGT", &[Cigar::Match(6)], 0, 0)];
        let bam = crate::test_support::write_bam_many(
            "ordered_wide_header",
            &[(CONTIG, GENOME.len()), ("chr2", 1000), ("chrUn", 100)],
            &recs,
        );
        let refr = shared_reference("ordered_wide_header", CONTIG, GENOME);
        let out = temp_path("ordered_wide_header_out", "bam");
        let cfg = OrderedConfig {
            n_workers: 1,
            reader_threads: 1,
            writer_threads: 1,
            out: out.clone(),
            bad_data: crate::batch::OnBadData::Stop,
            require_m5: false,
        };
        let stats = scan_ordered(&bam, refr, queries, &cfg, &output).unwrap();
        assert_eq!((stats.records_read, stats.records_off_reference), (1, 0));
    }

}
