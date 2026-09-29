//! Choose the output file that a record goes to.
//!
//! A scan given the output path `OUT` writes one file per `(worker, shard)`,
//! named `{stem}_{worker}_{shard}.{ext}` beside it: `calls.parquet` becomes
//! `calls_0_0.parquet`, `calls_0_1.parquet`, ..., and a path with no
//! extension just gains the suffix. See
//! [`ShardConfig::out_path`](crate::parallel::ShardConfig::out_path).
//!
//! Both indices come from one number: the hash of a **partition key**, a
//! combination of record-level attributes the user names with
//! `--partition-by`.
//!
//! ```text
//!     slot   = hash(key)  %  (n_workers * shards_per_worker)
//!     worker = slot       /  shards_per_worker
//!     shard  = slot       %  shards_per_worker
//! ```
//!
//! So the modulo picks a file out of the whole set and the division says which
//! worker owns it. Three consequences follow:
//!
//! - Within one run, records sharing a key always land in the same file,
//!   whatever the worker count and whatever `--shards-per-worker` is. That is
//!   the property the whole scheme exists for.
//! - Across runs, the guarantee is weaker. The hash itself is deterministic
//!   across runs and machines, but the modulo is taken over the total file
//!   count. A key therefore keeps its slot, which is the `shard` column, only
//!   while `n_workers * shards_per_worker` is unchanged. It keeps its file name
//!   only while both counts are unchanged, because one slot splits into a
//!   different `(worker, shard)` pair under a different `shards_per_worker`.
//! - With `--shards-per-worker 1` the slot *is* the worker index, so the
//!   routing is exactly what it was before per-worker sharding existed.
//!
//! # The key uses the record-field vocabulary
//!
//! A partition key names the same fields `--field` names, and the same code
//! reads them. A key of `ref_name` therefore groups by exactly what the
//! `ref_name` column would say, aux tags included:
//! [`PartitionKey::hash`] runs
//! [`RecordColumns::extract`](crate::record_columns::RecordColumns::extract)
//! and hashes the cells it filled. There is no second extractor to disagree
//! with the first, and a field added to `RecordField` can be partitioned on at
//! once.
//!
//! # What a key other than the qname costs
//!
//! Routing by qname puts every alignment of a fragment in one file, which is
//! what lets mate-overlap resolution see both sides. Any other key gives that
//! up, unless the key happens to be constant across a fragment. A read group is
//! constant; a reference name is not. [`fragment_warning`] says so at the point
//! of use. The key is still allowed, because a run that never resolves mates
//! has no reason to care.

use anyhow::{bail, Result};
use rust_htslib::bam::{HeaderView, Record};

use crate::record_columns::RecordColumns;
use crate::record_field::RecordField;

/// The key used when `--partition-by` is not given.
///
/// The qname, because it is the one key that keeps a fragment's alignments
/// together. This is the only place that default is written down.
pub const DEFAULT_KEY: &[RecordField] = &[RecordField::Qname];

// -------------------------------------------------------------------- hash --

/// FNV-1a, fed a byte at a time.
///
/// Small, dependency-free and deterministic across runs and machines, which is
/// all routing asks of a hash — it is not a digest and nothing downstream
/// stores it. The one property that matters is that equal keys hash equally,
/// so the same fragment cannot be split across files.
#[derive(Clone, Copy, Debug)]
pub struct Fnv(u64);

impl Default for Fnv {
    fn default() -> Self {
        Self::new()
    }
}

impl Fnv {
    const OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;

    #[inline]
    pub fn new() -> Self {
        Self(Self::OFFSET_BASIS)
    }

    #[inline]
    pub fn byte(&mut self, b: u8) {
        self.0 ^= b as u64;
        self.0 = self.0.wrapping_mul(Self::PRIME);
    }

    #[inline]
    pub fn bytes(&mut self, bs: &[u8]) {
        for b in bs {
            self.byte(*b);
        }
    }

    /// Little-endian, so the same value hashes the same on any target.
    #[inline]
    pub fn u64(&mut self, v: u64) {
        self.bytes(&v.to_le_bytes());
    }

    #[inline]
    pub fn finish(self) -> u64 {
        self.0
    }
}

// --------------------------------------------------------------------- key --

/// The record attributes that choose a file, and the machinery to read them.
///
/// Lives on the reader thread, which owns the one instance: extraction reuses
/// its buffers across records, so hashing a record allocates nothing.
pub struct PartitionKey {
    /// The key as named, kept for [`Debug`] and for nothing else — hashing
    /// reads the columns below.
    fields: Vec<RecordField>,
    /// Capacity zero: the Arrow builders inside are never appended to. Only
    /// `extract` is used, for its cells. Reusing `RecordColumns` here is what
    /// keeps a partition key and a column of the same name in agreement.
    cols: RecordColumns,
    /// How this run calls a record's strand, for a key naming `strand`,
    /// `conv_strand` or `read_reverse`. Held here because the key is hashed on
    /// the reader thread, before any worker has looked at the record, so it
    /// cannot borrow the call the walk will be given -- it has to make the
    /// same one, which is why both come from the same rule.
    strand: std::sync::Arc<crate::strand_rule::StrandRule>,
}

/// Written out rather than derived, and it has to stay that way.
///
/// `#[derive(Debug)]` here would require it of [`RecordColumns`], which would
/// require it of its `Column`s, which would require it of the Arrow builders
/// they own — a cascade of derives down three types to print a key that is
/// really just a list of field names. Those names are the whole of what anyone
/// wants to see, so this prints them and stops.
impl std::fmt::Debug for PartitionKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PartitionKey({})", names(&self.fields))
    }
}

impl PartitionKey {
    /// `fields` is already resolved — see [`RecordField::resolve`].
    ///
    /// An empty list is refused rather than defaulted: every record would hash
    /// alike and the whole scan would funnel into one file, which is never
    /// what someone meant. The default lives in [`DEFAULT_KEY`], applied where
    /// the arguments are read.
    pub fn new(fields: &[RecordField], strand: std::sync::Arc<crate::strand_rule::StrandRule>) -> Result<Self> {
        if fields.is_empty() {
            bail!("a partition key must name at least one field");
        }
        Ok(Self { fields: fields.to_vec(), cols: RecordColumns::new(fields, 0), strand })
    }

    /// Hash this record's key. Called once per record on the reader thread.
    #[inline]
    pub fn hash(&mut self, record: &Record, header: &HeaderView) -> u64 {
        // A rule that refuses this record, or cannot decide, gives it a null
        // strand cell and so a file of its own. Nothing is lost by that: the
        // record either fails the run when the scan reaches it, or is skipped
        // and counted there. Only the scan is in a position to say which, so
        // this does not try.
        self.cols.extract(record, header, self.strand.call(record).unwrap_or(None));
        self.cols.digest()
    }
}

/// The key as the user would write it, for messages.
pub fn names(fields: &[RecordField]) -> String {
    fields.iter().map(|f| f.name()).collect::<Vec<_>>().join(",")
}

/// What to tell the user about a key that may split a fragment, or `None` when
/// the key is the qname alone and cannot.
///
/// A warning rather than an error: splitting fragments is fine for a scan
/// whose output is aggregated per site, and it is the only way to get files
/// grouped by contig or by read group. It is silently wrong only for a run
/// that then resolves mate overlaps per file, which is why it is said out loud.
pub fn fragment_warning(fields: &[RecordField]) -> Option<String> {
    if fields == DEFAULT_KEY {
        return None;
    }
    let extra = names(
        &fields
            .iter()
            .copied()
            .filter(|f| *f != RecordField::Qname)
            .collect::<Vec<_>>(),
    );
    Some(if fields.contains(&RecordField::Qname) {
        format!(
            "warning: partitioning by {}: a fragment's alignments stay in one file only \
             if {extra} is constant across the fragment",
            names(fields)
        )
    } else {
        format!(
            "warning: partitioning by {} rather than qname: a fragment's alignments can \
             land in different files, so mate-overlap resolution cannot see both sides \
             of one fragment within a file",
            names(fields)
        )
    })
}

// ------------------------------------------------------------------ router --

/// Hash to (worker, shard), and the file count that follows from both.
///
/// Copied into the reader loop, so the arithmetic is a field read and two
/// divisions rather than a config lookup.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Router {
    n_workers: usize,
    shards_per_worker: usize,
}

impl Router {
    pub fn new(n_workers: usize, shards_per_worker: usize) -> Result<Self> {
        if n_workers == 0 {
            bail!("n_workers must be at least 1");
        }
        if shards_per_worker == 0 {
            bail!("shards per worker must be at least 1");
        }
        // `slots` multiplies these, so a pathological pair would wrap.
        if n_workers.checked_mul(shards_per_worker).is_none() {
            bail!("{n_workers} workers x {shards_per_worker} shards each overflows");
        }
        Ok(Self { n_workers, shards_per_worker })
    }

    pub fn n_workers(&self) -> usize {
        self.n_workers
    }

    pub fn shards_per_worker(&self) -> usize {
        self.shards_per_worker
    }

    /// Output files in total, across every worker.
    pub fn slots(&self) -> usize {
        self.n_workers * self.shards_per_worker
    }

    /// The file a hash lands in, as a single index over every file.
    ///
    /// Also the value written to the `shard` column, so a row says which file
    /// it came from without reference to the filename.
    #[inline]
    pub fn slot(&self, hash: u64) -> usize {
        (hash % self.slots() as u64) as usize
    }

    /// Split a slot into the worker that owns it and its index within that
    /// worker. The division is the worker; the remainder is the shard.
    #[inline]
    pub fn split(&self, slot: usize) -> (usize, usize) {
        (slot / self.shards_per_worker, slot % self.shards_per_worker)
    }

    /// The inverse of [`split`](Self::split).
    #[inline]
    pub fn slot_of(&self, worker: usize, shard: usize) -> usize {
        worker * self.shards_per_worker + shard
    }
}

// ------------------------------------------------------------------- tests --

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record_field::AuxKind;
    use crate::test_support::{built, header};
    use rust_htslib::bam::record::{Aux, Cigar};

    const CONTIG: &str = "chr1";

    fn key(fields: &[RecordField]) -> PartitionKey {
        PartitionKey::new(fields, crate::strand_rule::directional()).unwrap()
    }

    fn rec(name: &[u8]) -> Record {
        built(name, b"ACGTACGT", &[Cigar::Match(8)], 0, 0)
    }

    fn hash_of(k: &mut PartitionKey, r: &Record) -> u64 {
        k.hash(r, &header(CONTIG, 1000))
    }

    // ---- the hash --------------------------------------------------------

    #[test]
    fn the_hash_is_fnv1a() {
        // The offset basis, unchanged by an empty input.
        assert_eq!(Fnv::new().finish(), 0xcbf2_9ce4_8422_2325);

        let of = |s: &[u8]| {
            let mut h = Fnv::new();
            h.bytes(s);
            h.finish()
        };
        assert_ne!(of(b"a"), of(b"b"));
        // Order matters, so two names that are anagrams do not collide.
        assert_ne!(of(b"ab"), of(b"ba"));
    }

    // ---- the key ---------------------------------------------------------

    /// Every alignment of one fragment must reach one file, or mate-overlap
    /// resolution cannot see both sides. That rests on the hash depending on
    /// nothing but the key.
    #[test]
    fn one_qname_always_hashes_the_same() {
        let mut k = key(DEFAULT_KEY);
        let names: [&[u8]; 4] = [b"readA", b"readB", b"read_with_a_much_longer_name", b"x"];
        for name in names {
            let a = hash_of(&mut k, &rec(name));
            let b = hash_of(&mut k, &rec(name));
            assert_eq!(a, b, "hashing is not deterministic");
        }
    }

    /// The same qname on records that differ in every other way still hashes
    /// alike -- the key reads the named fields and nothing else.
    #[test]
    fn only_the_named_fields_are_read() {
        let mut k = key(DEFAULT_KEY);
        let plain = built(b"r", b"ACGTACGT", &[Cigar::Match(8)], 0, 0);
        let rev = crate::test_support::REVERSE;
        let moved = built(b"r", b"TTTTTTTT", &[Cigar::Match(8)], 500, rev);
        assert_eq!(hash_of(&mut k, &plain), hash_of(&mut k, &moved));
    }

    /// Distinct names must not all collapse onto one file, or sharding buys
    /// nothing.
    #[test]
    fn distinct_qnames_spread_across_slots() {
        let mut k = key(DEFAULT_KEY);
        let router = Router::new(4, 1).unwrap();
        let spread: std::collections::HashSet<usize> = (0..64u32)
            .map(|i| router.slot(hash_of(&mut k, &rec(format!("read{i}").as_bytes()))))
            .collect();
        assert!(spread.len() > 1, "the hash sends every qname to one file");
    }

    /// A key over several fields is more than the fields concatenated: each
    /// value is framed, so a boundary cannot move between them unnoticed.
    #[test]
    fn the_field_boundary_is_part_of_the_key() {
        let mut k = key(&[RecordField::Qname, RecordField::Aux(*b"RG", AuxKind::Str)]);
        let mut with_rg = |name: &[u8], rg: &str| {
            let mut r = rec(name);
            r.push_aux(b"RG", Aux::String(rg)).unwrap();
            hash_of(&mut k, &r)
        };
        // "ab" + "c" and "a" + "bc" are the same bytes in the same order.
        let joined = with_rg(&b"ab"[..], "c");
        let split = with_rg(&b"a"[..], "bc");
        assert_ne!(joined, split);
    }

    /// An absent value is not a value: the tag byte in the digest is what
    /// keeps a record without the aux tag from hashing like one that has it.
    #[test]
    fn an_absent_value_does_not_hash_like_a_present_one() {
        let mut k = key(&[RecordField::Aux(*b"RG", AuxKind::Str)]);
        let absent = hash_of(&mut k, &rec(b"r"));
        let present = {
            let mut r = rec(b"r");
            r.push_aux(b"RG", Aux::String("x")).unwrap();
            hash_of(&mut k, &r)
        };
        assert_ne!(absent, present);
        // And every record missing the tag still lands together, which is what
        // makes a partition by an optional tag usable at all.
        assert_eq!(absent, hash_of(&mut k, &rec(b"other")));
    }

    /// Fields with different values give different keys -- the point of
    /// partitioning by them.
    #[test]
    fn different_values_hash_differently() {
        let mut k = key(&[RecordField::Pos]);
        let a = hash_of(&mut k, &built(b"r", b"ACGT", &[Cigar::Match(4)], 10, 0));
        let b = hash_of(&mut k, &built(b"r", b"ACGT", &[Cigar::Match(4)], 11, 0));
        assert_ne!(a, b);
    }

    #[test]
    fn an_empty_key_is_refused() {
        let err = PartitionKey::new(&[], crate::strand_rule::directional()).unwrap_err();
        assert!(err.to_string().contains("at least one field"), "{err}");
    }

    /// `Debug` names the key rather than dumping the column machinery, which
    /// is also what keeps the impl hand-written -- see the note on it.
    #[test]
    fn debug_prints_the_field_names() {
        let k = key(&[RecordField::Qname, RecordField::RefName]);
        assert_eq!(format!("{k:?}"), "PartitionKey(qname,ref_name)");
    }

    // ---- the warning -----------------------------------------------------

    #[test]
    fn only_the_qname_key_passes_without_a_warning() {
        assert!(fragment_warning(DEFAULT_KEY).is_none());

        let w = fragment_warning(&[RecordField::RefName]).unwrap();
        assert!(w.contains("ref_name"), "{w}");
        assert!(w.contains("both sides"), "{w}");

        // Qname plus something else: fragments survive only if that something
        // is constant across the fragment, which is a different warning.
        let w = fragment_warning(&[RecordField::Qname, RecordField::RefName]).unwrap();
        assert!(w.contains("constant across the fragment"), "{w}");
    }

    // ---- the router ------------------------------------------------------

    #[test]
    fn a_slot_splits_into_a_worker_and_a_shard() {
        let r = Router::new(3, 4).unwrap();
        assert_eq!(r.slots(), 12);
        for slot in 0..r.slots() {
            let (w, s) = r.split(slot);
            assert!(w < 3, "worker {w} out of range");
            assert!(s < 4, "shard {s} out of range");
            assert_eq!(r.slot_of(w, s), slot, "split and slot_of disagree");
        }
    }

    /// Integer division picks the worker, the remainder picks the file within
    /// it -- so a worker's shards are contiguous slots.
    #[test]
    fn the_division_is_the_worker_and_the_remainder_the_shard() {
        let r = Router::new(2, 3).unwrap();
        assert_eq!(r.split(0), (0, 0));
        assert_eq!(r.split(2), (0, 2));
        assert_eq!(r.split(3), (1, 0));
        assert_eq!(r.split(5), (1, 2));
    }

    /// One file per worker is the old behaviour: the slot is the worker index,
    /// so routing is unchanged from before per-worker sharding existed.
    #[test]
    fn one_shard_per_worker_routes_by_worker_alone() {
        let r = Router::new(5, 1).unwrap();
        for h in [0u64, 1, 7, 12345, u64::MAX] {
            let (w, s) = r.split(r.slot(h));
            assert_eq!(w, (h % 5) as usize);
            assert_eq!(s, 0);
            assert_eq!(r.slot_of(w, s), w);
        }
    }

    #[test]
    fn every_slot_is_reachable() {
        let r = Router::new(2, 3).unwrap();
        let seen: std::collections::HashSet<usize> =
            (0..600u64).map(|h| r.slot(h)).collect();
        assert_eq!(seen.len(), r.slots(), "some file would never be written to");
    }

    #[test]
    fn a_zero_count_is_refused() {
        assert!(Router::new(0, 1).unwrap_err().to_string().contains("at least 1"));
        assert!(Router::new(1, 0).unwrap_err().to_string().contains("at least 1"));
    }
}