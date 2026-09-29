//! Move records between threads in batches. This module holds the batch itself,
//! the pool that recycles it, and the counts that a run reports.

use std::sync::mpsc::Receiver;

use rust_htslib::bam::Record;

/// What to do about a record the run cannot handle correctly.
///
/// The default is to stop. A BAM aligned against a different reference is
/// usually a mistake upstream, and so is a BAM whose aux blocks are damaged. A
/// run that quietly skips those records writes an output file that gives no
/// sign that records are missing. `--permissive` states that the caller knows
/// about the problem. It turns each such record into a counted warning.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OnBadData {
    Stop,
    WarnAndCount,
}

impl OnBadData {
    pub fn from_flag(permissive: bool) -> Self {
        if permissive { OnBadData::WarnAndCount } else { OnBadData::Stop }
    }
}

/// The error a record on a contig the reference lacks produces, unless the
/// run is permissive.
///
/// Raised at the first such record rather than from the header, because a
/// header is routinely broader than the reference -- a whole-genome `@SQ`
/// list against a single-chromosome reference is a normal way to work, and
/// rejecting it at startup would break that. A *record* there is different:
/// its calls would silently be missing from the output.
use crate::scanner::Walkability;

pub fn off_reference_error(name: &str, qname: &str) -> anyhow::Error {
    anyhow::anyhow!(
        "record {qname} is on contig {name}, which the reference does not have, so it \
         cannot be scanned. Pass --permissive to skip such records and count them \
         instead, or use the reference the BAM was aligned against."
    )
}

/// The error for a record on a contig with no `@SQ M5`, under `--require-m5`.
pub fn missing_m5_error(name: &str, qname: &str) -> anyhow::Error {
    anyhow::anyhow!(
        "record {qname} is on contig {name}, whose @SQ line has no M5 checksum, and \
         --require-m5 was given. Add M5 to the header's @SQ lines (`samtools dict REF.fa` \
         prints them for the FASTA the BAM was aligned against; `samtools reheader` \
         applies an edited header), or leave out --require-m5 to check names and lengths only."
    )
}

/// Count a record the reader has classified, and decide whether it is walked.
///
/// One place for the skip bookkeeping, shared by the tagged-BAM and hit-table
/// readers so they skip, count and refuse exactly the same records. An
/// off-reference record stops the run unless `bad_data` is permissive.
pub fn count_walkability(
    why: Walkability,
    rec: &rust_htslib::bam::Record,
    header: &rust_htslib::bam::HeaderView,
    bad_data: OnBadData,
    stats: &mut Stats,
) -> anyhow::Result<bool> {
    match why {
        Walkability::Walkable => return Ok(true),
        Walkability::OffReference if bad_data == OnBadData::Stop => {
            let name = String::from_utf8_lossy(header.tid2name(rec.tid() as u32));
            let qname = String::from_utf8_lossy(rec.qname()).into_owned();
            return Err(off_reference_error(&name, &qname));
        }
        Walkability::OffReference => stats.records_off_reference += 1,
        // Asked for explicitly, so --permissive does not relax it.
        Walkability::MissingM5 => {
            let name = String::from_utf8_lossy(header.tid2name(rec.tid() as u32));
            let qname = String::from_utf8_lossy(rec.qname()).into_owned();
            return Err(missing_m5_error(&name, &qname));
        }
        Walkability::NoSeq => stats.records_without_seq += 1,
        Walkability::Unmapped => {}
    }
    stats.records_skipped += 1;
    Ok(false)
}

/// The error for a thread the operating system refused to start.
///
/// `thread::Scope::spawn` panics in that case; callers use
/// `Builder::spawn_scoped` and this instead, so the run ends with a message
/// naming the usual cause rather than a panic.
pub fn thread_start_error(what: &str, e: std::io::Error) -> anyhow::Error {
    anyhow::anyhow!(
        "could not start the {what} thread: {e}. The system refused a new thread; the \
         per-user process/thread limit (`ulimit -u`) is the usual cause, so lower --threads \
         or raise the limit."
    )
}

/// Records are moved in batches so a channel is touched once per batch rather
/// than once per record.
pub const BATCH: usize = 512;

/// A batch's records on their way somewhere: to a worker, or back to the
/// reader to be filled again.
///
/// `id` is the record's ordinal in the input, assigned to every record read
/// whether or not the scan can walk it, so it lines up one-to-one with the
/// BAM. `walkable` is false for a record the scan cannot walk -- unmapped, or
/// on a contig the reference lacks -- which an output may still want to write.
pub type Recs = Vec<(Slot, Record)>;

/// Where a record sits in the input, and whether the scan can walk it.
#[derive(Debug, Clone, Copy)]
pub struct Slot {
    pub id: u64,
    pub walkable: bool,
}

#[derive(Debug, Default, Clone)]
pub struct Stats {
    pub records_read: u64,
    pub records_scanned: u64,
    pub records_skipped: u64,
    /// Of `records_skipped`, those on a contig the reference does not have.
    /// Broken out because it is the one skip reason that usually means
    /// something is misconfigured rather than something is unmapped.
    pub records_off_reference: u64,
    /// Of `records_skipped`, mapped records whose SEQ is `*`.
    pub records_without_seq: u64,
    /// Whether skipped records still reached the output: written unchanged in
    /// a tagged BAM, given a null row in a hit table.
    pub skipped_written: bool,
    /// Records whose aux block was damaged, so their tags were written at
    /// the front of the block rather than appended. Only counted in a
    /// permissive run; otherwise such a record stops the run.
    pub damaged_aux: u64,
    /// Records with an aux value that was not valid UTF-8, read as null
    /// rather than lossily converted into a different, plausible value.
    pub non_utf8_aux: u64,
    /// Records the run's strand rule declined to call. Skipped, and counted
    /// here so that a rule which does not cover the input says so rather than
    /// quietly shrinking the output.
    pub unknown_strand: u64,
    /// Hits the output could not place -- for a per-base tag, an anchor past
    /// either end of the read.
    pub unplaced_hits: u64,
    /// Read/reference agreement over every walked base; see
    /// [`crate::scanner::Concordance`].
    pub concordance: crate::scanner::Concordance,
    /// Rows written to each output file, by slot; empty for a tagged BAM.
    pub file_rows: Vec<u64>,
}

/// Records the workers have finished with, for the reader to read into again.
///
/// # Why records go back rather than being dropped
///
/// A fresh `Record` starts with no data buffer, so reading into one makes
/// htslib grow the buffer from nothing -- a `malloc` for the name in
/// `Record::new`, then `realloc`s as the CIGAR, sequence, qualities and aux
/// tags are copied in. Once per record, on the one thread every record passes
/// through, that was about a tenth of the whole scan's CPU, and under glibc it
/// also queued on the allocator lock the decode threads were using. A record
/// that has already held a read of about the same length needs none of it:
/// `sam_read1` overwrites it in place.
///
/// Reusing a record is safe because `sam_read1` overwrites all of it: core
/// fields, name, CIGAR, sequence, qualities and the whole aux block, so a tag
/// added last time is gone once the record is read into again. The one piece
/// of state it would not overwrite, `Record::cache_cigar`, is never called.
///
/// # Two ways to be fed
///
/// **Unbounded** (`limit` is `None`): whoever returns batches does so with
/// `try_send` and may drop them, and [`vec`](Self::vec) makes a new vector
/// whenever none is free. Returning adds no way to block.
///
/// **Bounded** (`limit` is `Some(n)`): at most `n` batch vectors ever exist.
/// Once all are out, [`vec`](Self::vec) waits for one to come back. That is
/// the backpressure of the ordered pipeline: a batch can only start when an
/// earlier one has been written, so no more than `n` batches are ever in
/// flight, and a channel sized for `n` can never be full.
pub struct Recycler {
    returned: Receiver<Recs>,
    records: Vec<Record>,
    /// Emptied batch vectors, kept for their allocation.
    vecs: Vec<Recs>,
    limit: Option<usize>,
    made: usize,
}

impl Recycler {
    pub fn new(returned: Receiver<Recs>, limit: Option<usize>) -> Self {
        Self { returned, records: Vec::new(), vecs: Vec::new(), limit, made: 0 }
    }

    /// A record to read into: a returned one when there is one, else new.
    #[inline]
    pub fn record(&mut self) -> Record {
        if self.records.is_empty() {
            // One batch at a time: enough to be going on with, and it keeps
            // the free list from holding much more than a batch.
            if let Ok(recs) = self.returned.try_recv() {
                self.take_back(recs);
            }
        }
        self.records.pop().unwrap_or_default()
    }

    /// An empty vector to stage a batch in. `None` only when bounded, every
    /// vector is out, and the side returning them has gone -- the pipeline
    /// has stopped.
    #[inline]
    pub fn vec(&mut self) -> Option<Recs> {
        if let Some(v) = self.vecs.pop() {
            return Some(v);
        }
        match self.limit {
            Some(limit) if self.made >= limit => {
                let recs = self.returned.recv().ok()?;
                self.take_back(recs);
                self.vecs.pop()
            }
            _ => {
                self.made += 1;
                Some(Vec::with_capacity(BATCH))
            }
        }
    }

    fn take_back(&mut self, mut recs: Recs) {
        self.records.extend(recs.drain(..).map(|(_, r)| r));
        self.vecs.push(recs);
    }
}
