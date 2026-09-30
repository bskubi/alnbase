//! Turn one record's CIGAR into a stream of aligned columns.
//!
//! This module holds [`walk_alignment`]. That function reads a BAM record and
//! the reference, and calls a closure once for each [`Column`]. Everything
//! after it in the program works on those columns. Nothing after it looks at a
//! CIGAR again.
//!
//! This is the module where correctness is most fragile. Two cursors, two
//! coordinate systems and two walk directions all have to stay in step. A
//! mistake here does not crash the program. It moves a coordinate or an offset
//! by one, and the output still looks like data. Read [`crate::column`] and
//! [`crate::seq`] before you change anything here.
//!
//! # Emitted order
//!
//! The walk emits columns 5' to 3' along the *conversion strand*, and not along
//! the reference. The caller supplies that strand as a [`StrandCall`], which
//! [`crate::strand_rule`] derives from a rule that the user writes. The walk
//! never reads the SAM flag itself. An aligner can put something other than the
//! read's orientation in flag `0x10`, and BSBolt does.
//!
//! A record on the original top strand walks forward. The columns come out in
//! increasing reference order, and the bases come out as the BAM stores them.
//!
//! A record on the original bottom strand walks backward. The CIGAR runs in
//! reverse, `refr_pos` decreases across the record, and the walk reverse
//! complements both the read base and the reference base. The variable
//! `is_ob` in the code marks this case. Several small helpers exist only to
//! serve it: `orient` reverses and complements bases, `orient_gen` reverses
//! anything that travels beside the bases but does not depend on the strand,
//! and `take_q` and `take_r` move the two cursors in the right direction.
//!
//! Emitted order is the reason a pattern can be written once and still work on
//! both strands. The user writes `CG`, and a record from either strand that
//! carries a CG produces the columns `C` then `G`.
//!
//! # What each CIGAR operation produces
//!
//! | operation | columns | reference side | read side |
//! |---|---|---|---|
//! | `M` `=` `X` | one per base | the reference base | the read base, with its quality |
//! | `D` | one per deleted base | the reference base | `GAP` |
//! | `I` | one per inserted base, or none | `GAP` | the read base, with its quality |
//! | `N` | a small fixed number, see below | the reference base, or `SKIP` | `SKIP` |
//! | `S` | none | | |
//! | `H` `P` | none | | |
//!
//! A deletion and an insertion are therefore ordinary columns. They are not a
//! hidden adjustment to a coordinate. A pattern can ask about them, because
//! `GAP` is a value like any other.
//!
//! [`Insertions`] controls whether the walk emits an insertion at all. An
//! inserted base has no reference coordinate. If the walk emits it, it sits
//! between two reference bases that are still next to each other, and it breaks
//! any pattern that spans it. If the walk skips it, a pattern about an
//! insertion becomes impossible to write. Neither answer is correct for every
//! query, so the option exists. `read_off` counts the inserted bases in both
//! cases, so it stays a true offset into `SEQ`.
//!
//! # The flank
//!
//! `WalkOpts::flank` adds columns past each end of the read. These columns
//! carry a real reference base. The read side carries `CLIP` where the aligner
//! soft-clipped a base, and `PAD` past the end of the read.
//!
//! The flank is what lets a cytosine at the last base of a read keep its `CG`
//! context. Without it the `G` has no column, the pattern cannot match, and the
//! last base of every read silently reports nothing.
//!
//! [`crate::scanner::WalkConfig::to_opts`] sets the width. By default it uses
//! the span of the widest query in the run, so a narrow query can gain hits
//! when the user adds a wider one. `--end-context` sets the width directly.
//!
//! # Introns
//!
//! A `CIGAR N` run is an intron. The walk does not emit one column per skipped
//! base. It emits `WalkOpts::splice_context` reference bases at each end of the
//! intron, against `SKIP` on the read side. It then emits one `(SKIP, SKIP)`
//! marker column for everything between them.
//!
//! There are two reasons. An intron is routinely tens of kilobases, so one
//! column per base would cost a hundred thousand reference reads and automaton
//! steps for one spliced read. The marker also stops a pattern from crossing
//! the junction by accident. `SKIP` is not a subset of any base, so a trace
//! that reaches the marker dies. Two exonic bases on either side of an intron
//! are not next to each other in the genome, whatever they look like in the
//! read.
//!
//! An intron short enough to fit inside the two context windows is emitted in
//! full, with no marker. A pattern that requires a marker therefore cannot
//! match a short intron.
//!
//! # The two cursors and the offset
//!
//! The walk keeps three positions.
//!
//! - `qpos` indexes `read_seq`, which holds `SEQ` in the order the BAM stores
//!   it. It starts at the end of the read for a backward walk.
//! - `rpos` is the reference coordinate. It starts at the end of the alignment
//!   for a backward walk.
//! - `read_off` counts the query bases that the walk has passed, in emitted
//!   order. It is the `SEQ` index of the column on both strands.
//!
//! A column that consumes no query base carries the nearest 5' offset instead
//! of its own. A column that consumes no reference base carries the nearest 5'
//! coordinate. These columns are a deletion, an intron and an insertion.
//! [`PAD_READ_OFF`] marks a flank column that describes no read base at all,
//! and the outputs write no offset for it.
//!
//! `read_off` is the offset into `SEQ`. It is not the offset into the read as
//! the sequencer read it, because `SEQ` does not hold the hard-clipped bases.
//! [`crate::hits`] adds the hard clip back on when it builds a row. Use
//! [`hard_clips_as_sequenced`] and [`soft_clips_as_sequenced`] for that, and
//! read `docs/design/off-by-one-safeguards.md`.
//!
//! Three `debug_assert` calls at the end of the walk check that both cursors
//! arrived where they should, and that `read_off` reached the length of the
//! read. They catch an operation that consumed the wrong amount.
//!
//! # Reading the reference
//!
//! The walk reads the reference from an [`Aref`], which is a memory-mapped
//! store of `Seq` bytes. It reads through a [`ContigMap`], because a BAM tid
//! and an `Aref` tid index two lists that are ordered independently. Matching
//! them by name is the only thing that keeps the walk on the right chromosome.
//! A record whose contig is not in the reference is an error here. The caller
//! removes those records before the walk.
//!
//! A coordinate past the end of a contig is legal. The `Aref` returns `PAD` for
//! it. This is what lets the flank run off the end of a chromosome without a
//! special case.

use anyhow::{anyhow, Result};
use itertools::{izip, Either};
use rust_htslib::bam::{ext::BamRecordExtensions, record::Cigar, record::CigarStringView, Record};
use std::iter::repeat_n;
use clap::ValueEnum;
use crate::contig_map::ContigMap;
use crate::aref::Aref;
pub use crate::column::{Column, PAD_READ_OFF};
use crate::seq::Seq;
use crate::strand_rule::StrandCall;

/// What to do with inserted bases.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Insertions {
    /// Do not emit them. An inserted base has no reference coordinate. A column
    /// that pairs it against a gap (`Seq::GAP`) therefore sits between two
    /// reference bases that are still adjacent, and it breaks every pattern that
    /// spans it. A reference-anchored caller never sees these columns.
    /// MethylDackel is one example, and any other caller that works from a
    /// pileup is another. The read offsets (`read_5p`, `read_3p`) still count the
    /// inserted bases, so they stay true offsets into SEQ.
    Skip,
    /// Emit them against a gap (`Seq::GAP`) on the reference side. You need this
    /// to write a pattern *about* an insertion, such as "is there a base
    /// inserted next to this cytosine". It also breaks every pattern that spans
    /// an insertion.
    Emit,
}

#[derive(Clone, Copy, Debug)]
pub struct WalkOpts {
    /// The number of extra columns that the walk emits at each end. On the read
    /// side these columns carry [`Seq::CLIP`] next to a soft clip, and
    /// [`Seq::PAD`] past the read. On the reference side they carry the real
    /// reference base. The number is at least `max(k) - 1`, so that a pattern
    /// whose reference part reaches past the read can still match. It is
    /// `max(k)`, unless `--end-context` sets it.
    pub flank: usize,
    /// The number of reference bases that the walk emits at each end of an
    /// intron. These columns carry [`Seq::SKIP`] on the read side. The walk
    /// replaces the rest of the intron with one marker.
    ///
    /// This is the same number as `flank`, for the same reason. An intron edge
    /// and a read end are both places where a pattern runs off the far side of
    /// what the read covers. A pattern that can reach past the end of a read
    /// must reach exactly as far into an intron. The code derives this number
    /// instead of asking for it, so that the two boundaries always agree.
    pub splice_context: usize,
    pub insertions: Insertions,
}

impl Default for WalkOpts {
    fn default() -> Self {
        Self {
            flank: 0,
            splice_context: 0,
            insertions: Insertions::Skip,
        }
    }
}

/// The parts of a record that [`walk_alignment`] reads, as htslib `SAM_*` bits.
///
/// The walk reads FLAG for the strand, and RNAME to find the contig. It reads
/// POS and CIGAR to place the read. It reads SEQ and QUAL to fill the columns. A
/// CRAM reader that serves the walk must decode at least these fields. See
/// [`crate::bam_io::limit_cram_decoding`]. This constant sits beside the
/// function, so that a change to what the walk reads is close to the list that
/// must follow it.
pub const WALK_SAM_FIELDS: u32 = rust_htslib::htslib::sam_fields_SAM_FLAG
    | rust_htslib::htslib::sam_fields_SAM_RNAME
    | rust_htslib::htslib::sam_fields_SAM_POS
    | rust_htslib::htslib::sam_fields_SAM_CIGAR
    | rust_htslib::htslib::sam_fields_SAM_SEQ
    | rust_htslib::htslib::sam_fields_SAM_QUAL;

/// Walk one record's alignment, and call `visit` once for each alignment column.
///
/// `contigs` maps the record's BAM tid onto a reference contig by name. The
/// caller must supply it, because a BAM tid and an [`Aref`] tid are indexes into
/// two lists that are ordered independently. See [`ContigMap`].
///
/// `strand` says which way to walk. The walk emits columns 5'->3' along the
/// conversion strand. A minus-strand call therefore runs backwards along the
/// reference. Both sides are reverse complemented, and `refr_pos` decreases
/// across the record. The caller gives the direction, and the walk does not read
/// it from the FLAG, because each aligner chooses which bits carry it. See
/// [`crate::strand_rule`].
///
/// - `refr_pos` is the reference coordinate of the column. An insertion column
///   has no coordinate of its own, so it carries the nearest 5' coordinate. A
///   coordinate past the end of the contig is emitted as it is, and the
///   reference base there comes back as [`Seq::PAD`].
/// - `read_off` is the 0-based query offset as sequenced, which indexes SEQ. A
///   deletion column and a ref-skip column carry the nearest 5' offset.
/// - The flank at each end is `opts.flank` columns of reference. Where that end
///   of the read is soft-clipped, the columns nearest the aligned part stand for
///   the clipped bases. Each of those columns holds [`Seq::CLIP`] on the read
///   side, the clipped base's `read_off`, and no quality. There are as many of
///   them as the clip is long, up to the width of the flank, so a long clip adds
///   no columns. The rest of the flank holds [`Seq::PAD`], which is past the
///   read, with [`PAD_READ_OFF`]. No offset describes those columns, and the
///   outputs write none.
/// - A `CIGAR N` run does not emit one column for each skipped base. It emits
///   `opts.splice_context` reference bases at each end against [`Seq::SKIP`]. It
///   then emits one `(SKIP, SKIP)` marker for everything between them. An intron
///   that fits inside the two context windows is emitted whole, with no
///   marker.
pub fn walk_alignment<F>(
    record: &Record,
    refr: &Aref,
    contigs: &ContigMap,
    read_seq: &mut Vec<Seq>,
    opts: WalkOpts,
    strand: StrandCall,
    visit: &mut F,
) -> Result<()>
where
    F: FnMut(Column),
{
    let is_ob = strand.walk_reversed();
    // The BAM's tid indexes the BAM's @SQ list; the reference has its own
    // order. Resolving by name is the only thing that keeps the sequence
    // fetched below on the chromosome the record is actually aligned to.
    // Callers filter these records out beforehand, so reaching here is a bug
    // rather than a data condition.
    let tid = contigs.refr_tid(record.tid()).ok_or_else(|| {
        anyhow!(
            "record aligned to BAM tid {} which is not in the reference; \
             it should have been filtered before the walk",
            record.tid()
        )
    })?;
    let cigar = record.cigar();
    let start = record.reference_start();
    let end = cigar.end_pos();
    let skip_ins = opts.insertions == Insertions::Skip;

    // htslib's 4-bit nibble encoding is bit-identical to Seq, so no decode step.
    read_seq.clear();
    let seq = record.seq();
    read_seq.extend((0..seq.len()).map(|i| Seq(seq.encoded_base(i))));
    let read_len = read_seq.len() as i64;

    // Phred scores are stored unencoded and in BAM order, exactly parallel to
    // SEQ, so they slice and orient the same way. A record with '*' for QUAL
    // has no array (or a mismatched one); fall back to -1 throughout.
    let quals = record.qual();
    // htslib stores a '*' QUAL as 0xff in every position, not as an empty array.
    let has_qual = quals.len() == read_seq.len() && quals.first() != Some(&0xff);

    // Query cursor indexes read_seq in BAM (reference) order; the reference
    // cursor walks backwards for is_ob. Both start at the far end in that case.
    let mut qpos: usize = if is_ob { read_seq.len() } else { 0 };
    let mut rpos: i64 = if is_ob { end } else { start };

    // Counts query bases consumed in emitted order — including skipped clipped
    // bases — so this is the SEQ index directly, on both strands.
    let mut read_off: i64 = 0;

    // Soft clips at the walk's first and last ends.
    let (lead_clip, trail_clip) = soft_clip_ends(&cigar, is_ob);

    // Leading flank: the `flank` reference bases immediately 5' of the aligned
    // region in emitted order. For OT that is below `start`; for OB, at `end`.
    // The clipped bases are the walk's first, offsets 0 .. lead_clip.
    if opts.flank > 0 {
        let n = opts.flank as i64;
        let (rs, re) = if is_ob { (end, end + n) } else { (start - n, start) };
        let shown = lead_clip.min(n);
        emit_flank(refr, tid, rs, re, is_ob, Flank::Leading { clip: shown, first_off: lead_clip - shown }, visit)?;
    }

    let ops = if is_ob {
        Either::Left(cigar.iter().rev())
    } else {
        Either::Right(cigar.iter())
    };

    for op in ops {
        match op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let n = *n as usize;
                let q = take_q(&mut qpos, n, is_ob);
                let (rs, re) = take_r(&mut rpos, n as i64, is_ob);
                let read = orient(read_seq[q.clone()].iter().copied(), is_ob);
                let refr = orient(fetch(refr, tid, rs, re)?, is_ob);
                let qual = orient_gen(phred(quals, has_qual, q), is_ob);
                for (r, f, p, qv) in izip!(read, refr, orient_pos(rs..re, is_ob), qual) {
                    visit(Column { read: r, refr: f, refr_pos: p, read_off, qual: qv });
                    read_off += 1;
                }
            }

            Cigar::Del(n) => {
                let n = *n as usize;
                let (rs, re) = take_r(&mut rpos, n as i64, is_ob);
                let read = orient(repeat_n(Seq::GAP, n), is_ob);
                let refr = orient(fetch(refr, tid, rs, re)?, is_ob);
                // Consumes no query: carry the nearest 5' offset.
                let anchor_off = read_off - 1;
                for (r, f, p) in izip!(read, refr, orient_pos(rs..re, is_ob)) {
                    // No read base, so no quality.
                    visit(Column { read: r, refr: f, refr_pos: p, read_off: anchor_off, qual: -1 });
                }
            }

            // An intron. Emitted as `context` reference bases at each end
            // against [`Seq::SKIP`] on the read side, with everything between
            // them collapsed to a single `(SKIP, SKIP)` column.
            //
            // Collapsing is not an optimisation bolted onto the same treatment
            // deletions get -- it is the only tractable one. A `D` is a handful
            // of bases; an `N` is an intron, routinely tens of kilobases, and
            // emitting it per base would mean a hundred thousand reference
            // fetches and automaton steps for one spliced read, plus a hundred
            // thousand chances for a pattern to fire at a position the read
            // never covered.
            //
            // The marker is what keeps a pattern from spanning the junction by
            // accident: `SKIP` is not a subset of any base, so a trace running
            // through it dies unless the pattern asked for a junction there.
            // That is the correct default -- two exonic bases either side of an
            // intron are not adjacent in the genome, whatever they look like in
            // the read.
            Cigar::RefSkip(n) => {
                let n = *n as usize;
                let (rs, re) = take_r(&mut rpos, n as i64, is_ob);
                let anchor_off = read_off - 1;
                let emit_range = |a: i64, z: i64, visit: &mut F| -> Result<()> {
                    let refr = orient(fetch(refr, tid, a, z)?, is_ob);
                    for (f, p) in refr.zip(orient_pos(a..z, is_ob)) {
                        visit(Column {
                            read: Seq::SKIP,
                            refr: f,
                            refr_pos: p,
                            read_off: anchor_off,
                            qual: -1,
                        });
                    }
                    Ok(())
                };

                let ctx = opts.splice_context as i64;
                if n as i64 <= 2 * ctx {
                    // Short enough to show in full, so there is nothing to
                    // elide and no marker. A pattern demanding a marker
                    // therefore cannot match a short intron, which is what
                    // `{N,}`-style alternation is for.
                    emit_range(rs, re, visit)?;
                } else {
                    // In emitted order the 5' context is the low end of the
                    // range for OT and the high end for OB, since the walk runs
                    // backwards along the reference there.
                    let (lead, trail) = if is_ob {
                        ((re - ctx, re), (rs, rs + ctx))
                    } else {
                        ((rs, rs + ctx), (re - ctx, re))
                    };
                    emit_range(lead.0, lead.1, visit)?;
                    // The first elided base in emitted order, so the marker
                    // says where the hole starts rather than where the intron
                    // does; with any context at all the two neighbouring
                    // columns bracket the hole exactly.
                    let hole = if is_ob { re - ctx - 1 } else { rs + ctx };
                    visit(Column {
                        read: Seq::SKIP,
                        refr: Seq::SKIP,
                        refr_pos: hole,
                        read_off: anchor_off,
                        qual: -1,
                    });
                    emit_range(trail.0, trail.1, visit)?;
                }
            }

            Cigar::Ins(n) => {
                let n = *n as usize;
                let q = take_q(&mut qpos, n, is_ob);
                if skip_ins {
                    read_off += n as i64;
                } else {
                    // Consumes no reference: carry the nearest 5' coordinate.
                    // That is the last coordinate emitted, which for is_ob is
                    // `rpos` itself (the cursor moves down), else `rpos - 1`.
                    let anchor_pos = if is_ob { rpos } else { rpos - 1 };
                    let read = orient(read_seq[q.clone()].iter().copied(), is_ob);
                    let refr = orient(repeat_n(Seq::GAP, n), is_ob);
                    let qual = orient_gen(phred(quals, has_qual, q), is_ob);
                    for ((r, f), qv) in read.zip(refr).zip(qual) {
                        visit(Column { read: r, refr: f, refr_pos: anchor_pos, read_off, qual: qv });
                        read_off += 1;
                    }
                }
            }

            // Soft-clipped bases were sequenced but never aligned, so they
            // produce no columns at all; the offset still steps over them.
            Cigar::SoftClip(n) => {
                let n = *n as usize;
                take_q(&mut qpos, n, is_ob);
                read_off += n as i64;
            }

            // Hard clips aren't in SEQ at all, so they advance nothing.
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    // Trailing flank: mirror of the leading one. The clipped bases are the
    // walk's last, offsets read_len - trail_clip .. read_len.
    if opts.flank > 0 {
        let n = opts.flank as i64;
        let (rs, re) = if is_ob { (start - n, start) } else { (end, end + n) };
        let shown = trail_clip.min(n);
        emit_flank(refr, tid, rs, re, is_ob, Flank::Trailing { clip: shown, first_off: read_len - trail_clip }, visit)?;
    }

    if is_ob {
        debug_assert_eq!(qpos, 0);
        debug_assert_eq!(rpos, start);
    } else {
        debug_assert_eq!(qpos, read_seq.len());
        debug_assert_eq!(rpos, end);
    }
    debug_assert_eq!(read_off, read_len);
    Ok(())
}

/// Count the soft-clipped bases at the read's sequenced 5' and 3' ends.
///
/// SEQ holds the soft-clipped bases, and `read_5p` and `read_3p` already count
/// them. `read_5p - soft_clip_5p - hard_clip_5p` is therefore the offset from the
/// first aligned base. A soft clip sits inside any hard clip at the same end. If
/// a soft clip is the only operation in the CIGAR that is not a hard clip, the
/// function counts it once, at the 5' end.
pub fn soft_clips_as_sequenced(record: &Record) -> (i64, i64) {
    const SOFT_CLIP: u32 = 4; // BAM_CSOFT_CLIP
    const HARD_CLIP: u32 = 5; // BAM_CHARD_CLIP
    let raw = record.raw_cigar();
    let first = raw.iter().position(|&v| v & 0xf != HARD_CLIP);
    let last = raw.iter().rposition(|&v| v & 0xf != HARD_CLIP);
    let soft = |i: Option<usize>| match i.map(|i| raw[i]) {
        Some(v) if v & 0xf == SOFT_CLIP => (v >> 4) as i64,
        _ => 0,
    };
    let (first_n, last_n) = (soft(first), if last == first { 0 } else { soft(last) });
    if record.is_reverse() { (last_n, first_n) } else { (first_n, last_n) }
}

/// Count the hard-clipped bases at the read's sequenced 5' and 3' ends.
///
/// SEQ does not hold the hard-clipped bases, but the sequencer read them. An
/// offset that is relative to the read, and not to the record, must therefore
/// add them. The 5' end of a forward record is the first operation in the CIGAR.
/// The 5' end of a reverse record is the last operation, because the BAM stores
/// SEQ and CIGAR in reference orientation. This function reads the raw CIGAR, so
/// it allocates nothing.
pub fn hard_clips_as_sequenced(record: &Record) -> (i64, i64) {
    const HARD_CLIP: u32 = 5; // BAM_CHARD_CLIP
    let raw = record.raw_cigar();
    let hard = |op: Option<&u32>| match op {
        Some(&v) if v & 0xf == HARD_CLIP => (v >> 4) as i64,
        _ => 0,
    };
    let (first, last) = (hard(raw.first()), if raw.len() > 1 { hard(raw.last()) } else { 0 });
    if record.is_reverse() { (last, first) } else { (first, last) }
}

/// Which end a flank is at, and which clip columns it holds. The flank holds
/// `clip` clip columns, nearest the aligned part. The first of them in emitted
/// order is at `first_off`.
#[derive(Clone, Copy)]
enum Flank {
    Leading { clip: i64, first_off: i64 },
    Trailing { clip: i64, first_off: i64 },
}

/// Count the soft-clipped bases at the walk's first end and last end.
fn soft_clip_ends(cigar: &CigarStringView, rc: bool) -> (i64, i64) {
    let ops = || cigar.iter().filter(|op| !matches!(op, Cigar::HardClip(_)));
    let clip = |op: Option<&Cigar>| match op {
        Some(Cigar::SoftClip(n)) => *n as i64,
        _ => 0,
    };
    let (first, last) = (clip(ops().next()), if ops().count() > 1 { clip(ops().last()) } else { 0 });
    if rc { (last, first) } else { (first, last) }
}

/// Emit the flank `[rs, re)` in emitted order. The clipped bases get a real
/// reference base against CLIP. The columns past the read get PAD. A trailing
/// flank runs in the reverse order.
fn emit_flank<F>(
    refr: &Aref,
    tid: usize,
    rs: i64,
    re: i64,
    rc: bool,
    flank: Flank,
    visit: &mut F,
) -> Result<()>
where
    F: FnMut(Column),
{
    let n = re - rs;
    let refr = orient(fetch(refr, tid, rs, re)?, rc);
    for (i, (f, p)) in refr.zip(orient_pos(rs..re, rc)).enumerate() {
        let i = i as i64;
        // Index among the clip columns, if this is one.
        let clip_index = match flank {
            Flank::Leading { clip, .. } if i >= n - clip => Some(i - (n - clip)),
            Flank::Trailing { clip, .. } if i < clip => Some(i),
            _ => None,
        };
        let column = match (clip_index, flank) {
            (Some(k), Flank::Leading { first_off, .. } | Flank::Trailing { first_off, .. }) => {
                // A clipped base is sequenced but not aligned: its offset is
                // real, and it has no observed base or quality here.
                Column { read: Seq::CLIP, refr: f, refr_pos: p, read_off: first_off + k, qual: -1 }
            }
            // Past the read: no read base, no quality, no read offset.
            (None, _) => Column { read: Seq::PAD, refr: f, refr_pos: p, read_off: PAD_READ_OFF, qual: -1 },
        };
        visit(column);
    }
    Ok(())
}

/// Return the phred scores over a query range. Each score is -1 when the record
/// has no QUAL.
#[inline]
fn phred(
    quals: &[u8],
    has_qual: bool,
    q: std::ops::Range<usize>,
) -> impl DoubleEndedIterator<Item = i32> + '_ {
    let n = q.len();
    if has_qual {
        Either::Left(quals[q].iter().map(|&v| v as i32))
    } else {
        Either::Right(repeat_n(-1i32, n))
    }
}

#[inline]
fn fetch(
    refr: &Aref,
    tid: usize,
    rs: i64,
    re: i64,
) -> Result<impl DoubleEndedIterator<Item = Seq> + Clone + '_> {
    refr
        .fetch(tid, rs, re)
        .ok_or_else(|| anyhow!("no contig with reference tid {tid}"))
}

/// Reverse-complement an iterator when `rc` is true. Both arms have one type.
#[inline]
fn orient<I>(it: I, rc: bool) -> impl Iterator<Item = Seq>
where
    I: DoubleEndedIterator<Item = Seq>,
{
    if rc {
        Either::Left(it.rev().map(Seq::complement))
    } else {
        Either::Right(it)
    }
}

/// Reverse an iterator when `rc` is true, and do not complement it. This is for
/// a value that travels beside the bases and does not depend on the strand.
#[inline]
fn orient_gen<T, I>(it: I, rc: bool) -> impl Iterator<Item = T>
where
    I: DoubleEndedIterator<Item = T>,
{
    if rc {
        Either::Left(it.rev())
    } else {
        Either::Right(it)
    }
}

/// Reverse a coordinate range when `rc` is true, so that the positions come out
/// in emitted order.
#[inline]
fn orient_pos<I>(it: I, rc: bool) -> impl Iterator<Item = i64>
where
    I: DoubleEndedIterator<Item = i64>,
{
    orient_gen(it, rc)
}

/// Consume `n` query bases, and return the slice range. In rc mode the function
/// advances to the left.
#[inline]
fn take_q(qpos: &mut usize, n: usize, rc: bool) -> std::ops::Range<usize> {
    if rc { *qpos -= n; }
    let s = *qpos;
    if !rc { *qpos += n; }
    s..s + n
}

/// Consume `n` reference bases, and return `[start, end)`. In rc mode the
/// function advances to the left.
#[inline]
fn take_r(rpos: &mut i64, n: i64, rc: bool) -> (i64, i64) {
    if rc { *rpos -= n; }
    let s = *rpos;
    if !rc { *rpos += n; }
    (s, s + n)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::contig_map::ContigMap;
    use crate::test_support::{
        built, contig_map, header_many, reference, reference_many, round_tripped,
        LAST_IN_TEMPLATE, REVERSE,
    };
    use rust_htslib::bam::record::Cigar;

    const CONTIG: &str = "chr1";
    const GENOME: &str = "ACGTACGTACGTACGTACGT";

    /// The strand these tests walk records on: the directional rule, which is
    /// what the FLAG-reading code the walk used to hold did.
    fn call(rec: &Record) -> StrandCall {
        crate::tags::Library::Directional.call(rec)
    }

    fn walk(rec: &Record, refr: &Aref, opts: WalkOpts) -> Vec<Column> {
        let map = contig_map(refr, CONTIG, GENOME.len());
        let mut seq = Vec::new();
        let mut cols = Vec::new();
        walk_alignment(rec, refr, &map, &mut seq, opts, call(rec), &mut |c| cols.push(c)).unwrap();
        cols
    }

    fn rec(seq: &[u8], cigar: &[Cigar], pos: i64, flags: u16) -> Record {
        built(b"read1", seq, cigar, pos, flags)
    }

    /// The walk direction comes from the call, not from the FLAG.
    ///
    /// This is the whole point of the strand rule: an aligner that writes the
    /// converted reference strand into 0x10 rather than the read's orientation
    /// (BSBolt does) would otherwise be walked the wrong way round, silently,
    /// and every coordinate and offset in the output would be wrong. Here the
    /// same flagless record is walked both ways, and nothing but the call
    /// differs.
    #[test]
    fn the_call_decides_which_way_the_walk_runs() {
        let refr = reference("align_call_direction", CONTIG, GENOME);
        let map = contig_map(&refr, CONTIG, GENOME.len());
        let r = rec(b"ACGT", &[Cigar::Match(4)], 0, 0);

        let go = |strand: StrandCall| {
            let (mut seq, mut cols) = (Vec::new(), Vec::new());
            walk_alignment(&r, &refr, &map, &mut seq, WalkOpts::default(), strand, &mut |c| {
                cols.push(c)
            })
            .unwrap();
            cols.iter().map(|c| c.refr_pos).collect::<Vec<_>>()
        };

        assert_eq!(go(StrandCall::from_origin(crate::tags::Strand::Ot)), [0, 1, 2, 3]);
        assert_eq!(
            go(StrandCall::from_origin(crate::tags::Strand::Ob)),
            [3, 2, 1, 0],
            "a minus-strand call walks backwards though FLAG 0x10 is clear"
        );
    }

    /// An intron is one marker column, not one column per skipped base.
    ///
    /// The 20-base genome makes the point in miniature: a 10-base skip with no
    /// context emits a single column, so the cost of a junction is independent
    /// of its length. On real data that is the difference between one column
    /// and a hundred thousand.
    #[test]
    fn an_intron_collapses_to_one_column() {
        let refr = reference("align_intron", CONTIG, GENOME);
        // 4M 10N 4M: four exonic bases, a ten-base intron, four more.
        let cig = vec![Cigar::Match(4), Cigar::RefSkip(10), Cigar::Match(4)];
        let cols = walk(&rec(b"ACGTACGT", &cig, 0, 0), &refr, WalkOpts::default());

        assert_eq!(cols.len(), 9, "four exonic, one marker, four exonic");
        let m = cols[4];
        assert_eq!(m.read, Seq::SKIP);
        assert_eq!(m.refr, Seq::SKIP, "a marker has no reference base either");
        assert_eq!(m.refr_pos, 4, "the first elided base");
        assert_eq!(m.qual, -1, "no read base, so no quality");
        assert_eq!(m.read_off, 3, "consumes no query: the nearest 5' offset");
        // The exonic columns either side carry their real coordinates, so the
        // junction's extent is readable from the two neighbours.
        assert_eq!(cols[3].refr_pos, 3);
        assert_eq!(cols[5].refr_pos, 14);
    }

    /// With context, the marker sits between the intron's own first and last
    /// reference bases -- which is what a splice-motif pattern matches on.
    #[test]
    fn context_shows_the_intron_edges() {
        let refr = reference("align_intron_ctx", CONTIG, GENOME);
        let cig = vec![Cigar::Match(4), Cigar::RefSkip(10), Cigar::Match(4)];
        let opts = WalkOpts { splice_context: 2, ..WalkOpts::default() };
        let cols = walk(&rec(b"ACGTACGT", &cig, 0, 0), &refr, opts);

        assert_eq!(cols.len(), 4 + 2 + 1 + 2 + 4);
        // Intronic bases are real reference against SKIP on the read side: the
        // read does not cover them, but the reference is still there to match.
        for i in [4, 5, 7, 8] {
            assert_eq!(cols[i].read, Seq::SKIP, "column {i}");
            assert!(cols[i].refr.is_subset_of(Seq::N), "column {i} has a real base");
        }
        assert_eq!(cols[4].refr_pos, 4, "the intron's first base");
        assert_eq!(cols[8].refr_pos, 13, "the intron's last base");
        assert_eq!(cols[6].refr, Seq::SKIP, "the marker, between the two windows");
        assert_eq!(cols[6].refr_pos, 6, "the first elided base");
    }

    /// An intron that fits inside the two context windows is emitted whole,
    /// with no marker at all -- so `,` means "a junction with something elided"
    /// and a pattern demanding one cannot match a short intron.
    #[test]
    fn a_short_intron_has_no_marker() {
        let refr = reference("align_intron_short", CONTIG, GENOME);
        let cig = vec![Cigar::Match(4), Cigar::RefSkip(4), Cigar::Match(4)];
        let opts = WalkOpts { splice_context: 2, ..WalkOpts::default() };
        let cols = walk(&rec(b"ACGTACGT", &cig, 0, 0), &refr, opts);

        assert_eq!(cols.len(), 12, "every intronic base is shown");
        assert!(
            !cols.iter().any(|c| c.refr == Seq::SKIP),
            "nothing was elided, so nothing marks an elision"
        );
        for i in 4..8 {
            assert_eq!(cols[i].read, Seq::SKIP);
            assert_eq!(cols[i].refr_pos, i as i64);
        }
    }

    /// A deletion and an intron are different values now, so a pattern can tell
    /// them apart. They were one arm until the walk learned to.
    #[test]
    fn a_deletion_and_an_intron_are_different_columns() {
        let refr = reference("align_del_vs_n", CONTIG, GENOME);
        let del = walk(
            &rec(b"ACGTACGT", &[Cigar::Match(4), Cigar::Del(2), Cigar::Match(4)], 0, 0),
            &refr,
            WalkOpts::default(),
        );
        let skip = walk(
            &rec(b"ACGTACGT", &[Cigar::Match(4), Cigar::RefSkip(2), Cigar::Match(4)], 0, 0),
            &refr,
            WalkOpts { splice_context: 2, ..WalkOpts::default() },
        );
        assert_eq!(del.len(), skip.len(), "same shape at this length");
        for i in 4..6 {
            assert_eq!(del[i].read, Seq::GAP);
            assert_eq!(skip[i].read, Seq::SKIP);
            assert_eq!(del[i].refr, skip[i].refr, "the reference base is the same either way");
        }
    }

    /// On the bottom strand the walk runs backwards along the reference, so the
    /// 5' context window is the high end of the intron and the marker sits at
    /// the first base elided *in emitted order*.
    #[test]
    fn an_intron_orients_with_the_walk() {
        let refr = reference("align_intron_ob", CONTIG, GENOME);
        let cig = vec![Cigar::Match(4), Cigar::RefSkip(10), Cigar::Match(4)];
        let opts = WalkOpts { splice_context: 2, ..WalkOpts::default() };
        let cols = walk(&rec(b"ACGTACGT", &cig, 0, REVERSE), &refr, opts);

        // Reference positions descend throughout.
        let exonic: Vec<i64> =
            cols.iter().filter(|c| c.refr != Seq::SKIP).map(|c| c.refr_pos).collect();
        assert!(exonic.windows(2).all(|w| w[0] > w[1]), "descending: {exonic:?}");
        assert_eq!(cols[0].refr_pos, 17, "starts at the far end");
        assert_eq!(cols[4].refr_pos, 13, "the 5' context is the high end of the intron");
        assert_eq!(cols[6].refr, Seq::SKIP);
        assert_eq!(cols[6].refr_pos, 11, "the first elided base, counting down");
        assert_eq!(cols[8].refr_pos, 4, "the 3' context is the low end");
    }

    /// A record built in memory walks identically to the same record parsed
    /// from a BAM, so a synthetic fixture is evidence about the real scan path.
    ///
    /// Everything else in the crate's tests builds records with
    /// `test_support::built`; this is what licenses that.
    #[test]
    fn built_records_walk_like_parsed_ones() {
        let refr = reference("align_equiv", CONTIG, GENOME);
        let m8 = vec![Cigar::Match(8)];
        let ins = vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)];
        let del = vec![Cigar::Match(4), Cigar::Del(2), Cigar::Match(4)];
        let skip = vec![Cigar::Match(4), Cigar::RefSkip(6), Cigar::Match(4)];
        let clip = vec![Cigar::SoftClip(2), Cigar::Match(8)];
        let cases: [(&str, &[u8], &[Cigar], i64, u16); 7] = [
            ("plain", b"ACGTACGT", &m8, 0, 0),
            ("reverse", b"ACGTACGT", &m8, 4, REVERSE),
            ("insertion", b"ACGTTTACGT", &ins, 2, 0),
            ("deletion", b"ACGTACGT", &del, 0, 0),
            ("intron", b"ACGTACGT", &skip, 0, 0),
            ("softclip", b"GGACGTACGT", &clip, 6, 0),
            ("mateflags", b"ACGTACGT", &m8, 8, 0x1 | LAST_IN_TEMPLATE),
        ];
        for (name, seq, cigar, pos, flags) in cases {
            let r = rec(seq, cigar, pos, flags);
            let parsed = round_tripped(&r, &format!("align_{name}"), CONTIG, GENOME.len());
            for opts in [
                WalkOpts::default(),
                WalkOpts { flank: 3, ..WalkOpts::default() },
                WalkOpts { insertions: Insertions::Emit, ..WalkOpts::default() },
                WalkOpts { splice_context: 2, ..WalkOpts::default() },
            ] {
                assert_eq!(
                    walk(&r, &refr, opts),
                    walk(&parsed, &refr, opts),
                    "{name} diverged between a built record and a parsed one, opts {opts:?}"
                );
            }
        }
    }

    /// A '*' QUAL (0xff in every position) gives no quality, like a deletion.
    #[test]
    fn a_missing_quality_is_minus_one_not_255() {
        let refr = reference("align_noqual", CONTIG, GENOME);
        let mut r = rec(b"ACGT", &[Cigar::Match(4)], 0, 0);
        r.set(b"r", Some(&rust_htslib::bam::record::CigarString(vec![Cigar::Match(4)])), b"ACGT", &[0xff; 4]);
        assert!(walk(&r, &refr, WalkOpts::default()).iter().all(|c| c.qual == -1));
    }

    /// The bottom strand walks the reference backwards with both sides
    /// complemented, so positions descend across the record.
    #[test]
    fn the_bottom_strand_walks_backwards() {
        let refr = reference("align_ob", CONTIG, GENOME);
        // is_last_in_template() == is_reverse() selects the top strand, so a
        // lone REVERSE is the bottom strand.
        let cols = walk(&rec(b"ACGTACGT", &[Cigar::Match(8)], 0, REVERSE), &refr, WalkOpts::default());
        assert_eq!(cols.len(), 8);
        assert_eq!(cols[0].refr_pos, 7, "OB starts at the 3' end of the span");
        assert_eq!(cols[7].refr_pos, 0);
        // read_off indexes SEQ in emitted order on both strands.
        assert_eq!(cols[0].read_off, 0);
        assert_eq!(cols[7].read_off, 7);

        // Both bits set is the top strand again, and it ascends.
        let top = walk(
            &rec(b"ACGTACGT", &[Cigar::Match(8)], 0, REVERSE | LAST_IN_TEMPLATE),
            &refr,
            WalkOpts::default(),
        );
        assert_eq!(top[0].refr_pos, 0);
        assert_eq!(top[7].refr_pos, 7);
    }

    /// Deletions carry a reference base and no read base; insertions the
    /// reverse. A column with no read base carries `qual = -1`.
    #[test]
    fn indels_have_one_side_only() {
        let refr = reference("align_indel", CONTIG, GENOME);
        let del = vec![Cigar::Match(4), Cigar::Del(2), Cigar::Match(4)];
        let cols = walk(&rec(b"ACGTACGT", &del, 0, 0), &refr, WalkOpts::default());
        assert_eq!(cols.len(), 10, "the two deleted reference bases get columns");
        assert_eq!(cols[4].read, Seq::GAP);
        assert_eq!(cols[4].qual, -1, "no read base means no quality");

        let ins = vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)];
        let opts = WalkOpts { insertions: Insertions::Emit, ..WalkOpts::default() };
        let cols = walk(&rec(b"ACGTTTACGT", &ins, 0, 0), &refr, opts);
        assert_eq!(cols.len(), 10);
        assert_eq!(cols[4].refr, Seq::GAP, "an inserted base has no reference");
    }

    /// Skipping insertions keeps the reference bases flanking one adjacent, so
    /// a pattern can span an insertion. That is the whole reason it is default.
    #[test]
    fn skipped_insertions_leave_the_reference_contiguous() {
        let refr = reference("align_insskip", CONTIG, GENOME);
        let ins = vec![Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)];
        let cols = walk(&rec(b"ACGTTTACGT", &ins, 0, 0), &refr, WalkOpts::default());
        assert_eq!(cols.len(), 8, "the two inserted bases get no column");
        assert_eq!(cols[3].refr_pos, 3);
        assert_eq!(cols[4].refr_pos, 4, "no gap in reference coordinates");
        // read_off still indexes SEQ, so it steps over the skipped bases.
        assert_eq!(cols[3].read_off, 3);
        assert_eq!(cols[4].read_off, 6);
    }

    /// Insertions and soft clips are emitted by two blocks written out
    /// separately, so nothing but this holds them equal.
    ///
    /// Skipping an insertion or a soft clip advances the offset without
    /// emitting a column, so `read_off` stays a true SEQ index.
    #[test]
    fn skipped_insertions_and_soft_clips_advance_the_offset() {
        let refr = reference("align_alike", CONTIG, GENOME);
        let skipped_ins = walk(
            &rec(b"ACGTTTACGT", &[Cigar::Match(4), Cigar::Ins(2), Cigar::Match(4)], 0, 0),
            &refr,
            WalkOpts::default(),
        );
        let skipped_clip = walk(
            &rec(b"TTACGTACGT", &[Cigar::SoftClip(2), Cigar::Match(8)], 0, 0),
            &refr,
            WalkOpts::default(),
        );
        assert_eq!(skipped_ins.len(), 8);
        assert_eq!(skipped_clip.len(), 8);
        assert_eq!(skipped_ins[4].read_off, 6, "the two skipped bases are stepped over");
        assert_eq!(skipped_clip[0].read_off, 2);
    }

    /// Beside a soft clip the flank's inner columns stand for the clipped bases:
    /// CLIP, the real reference base, the clipped base's offset and no quality.
    /// Past the clip, or past a clip longer than the flank, it is PAD with no
    /// offset. Top and bottom walks run the flank from opposite ends.
    #[test]
    fn a_soft_clip_fills_the_flank_beside_the_aligned_part() {
        let refr = reference("align_clip_flank", CONTIG, GENOME);
        let opts = WalkOpts { flank: 3, ..WalkOpts::default() };
        let clip = [Cigar::SoftClip(2), Cigar::Match(6), Cigar::SoftClip(3)];
        let shape = |cols: &[Column]| {
            cols.iter()
                .map(|c| (c.read.to_char(), c.refr_pos, (c.read != Seq::PAD).then_some(c.read_off), c.qual))
                .collect::<Vec<_>>()
        };
        let l = 'L';

        // Top: 2 clipped bases before POS 5, 3 after the aligned end at 11.
        let top = walk(&rec(b"TTACGTACGGG", &clip, 5, 0), &refr, opts);
        assert_eq!(top.len(), 12);
        let (lead, trail) = (shape(&top[..3]), shape(&top[9..]));
        assert_eq!(lead, vec![('X', 2, None, -1), (l, 3, Some(0), -1), (l, 4, Some(1), -1)]);
        assert_eq!(trail, vec![(l, 11, Some(8), -1), (l, 12, Some(9), -1), (l, 13, Some(10), -1)]);
        assert_eq!(top[3].read_off, 2, "the first aligned base follows the clip");
        assert_eq!((top[3].refr, top[1].refr), (Seq::C, Seq::T), "real reference in the clip columns");

        // Bottom: walked from the aligned end, so the 3-base clip comes first.
        let bottom = walk(&rec(b"TTACGTACGGG", &clip, 5, crate::test_support::REVERSE), &refr, opts);
        let (lead, trail) = (shape(&bottom[..3]), shape(&bottom[9..]));
        assert_eq!(lead, vec![(l, 13, Some(0), -1), (l, 12, Some(1), -1), (l, 11, Some(2), -1)]);
        assert_eq!(trail, vec![(l, 4, Some(9), -1), (l, 3, Some(10), -1), ('X', 2, None, -1)]);

        // A clip longer than the flank shows only its innermost bases.
        let long = walk(&rec(b"AAAAAAAAAACGTA", &[Cigar::SoftClip(10), Cigar::Match(4)], 10, 0), &refr, opts);
        assert_eq!(long.len(), 10);
        assert_eq!(shape(&long[..3]).iter().map(|c| c.2).collect::<Vec<_>>(), vec![Some(7), Some(8), Some(9)]);
        assert!(long[7..].iter().all(|c| c.read == Seq::PAD), "no clip at the other end");
    }

    /// The flank supplies real reference context past the aligned region, which
    /// is what lets a query match a cytosine at the last aligned base.
    #[test]
    fn the_flank_extends_past_the_alignment() {
        let refr = reference("align_flank", CONTIG, GENOME);
        let opts = WalkOpts { flank: 2, ..WalkOpts::default() };
        let cols = walk(&rec(b"ACGTACGT", &[Cigar::Match(8)], 4, 0), &refr, opts);
        assert_eq!(cols.len(), 12, "two columns added at each end");
        assert_eq!(cols[0].refr_pos, 2, "starts two before POS");
        assert_eq!(cols[0].read, Seq::PAD, "a flank column has no read base");
        assert_eq!(cols[0].qual, -1);
        assert_eq!(cols[11].refr_pos, 13);
        // Off the start of the contig the reference is a pad rather than an error.
        let cols = walk(&rec(b"ACGTACGT", &[Cigar::Match(8)], 0, 0), &refr, opts);
        assert_eq!(cols[0].refr, Seq::PAD, "before the contig there is no base");
    }

    /// The reference base for a column comes from the contig whose *name*
    /// matches the record's, not from whatever sits at the same ordinal in the
    /// reference.
    ///
    /// This is the regression test for a real bug: a karyotype-ordered BAM
    /// (chr1, chr2, chr10) against a lexicographically ordered FASTA (chr1,
    /// chr10, chr2) sent every chr2 read to chr10's sequence. Nothing crashed
    /// and every base emitted was a real base, so the output looked ordinary
    /// while being wrong everywhere except the one contig whose index happened
    /// to agree. A single-contig fixture cannot see this, which is why this
    /// test needs three.
    #[test]
    fn the_reference_side_comes_from_the_contig_with_the_matching_name() {
        // Distinct sequences so the columns say which contig was read.
        let refr = reference_many(
            "walk_by_name",
            &[("chr1", "AAAAAAAA"), ("chr10", "CCCCCCCC"), ("chr2", "GGGGGGGG")],
        );
        // Karyotype order: chr2 is BAM tid 1, but reference tid 2.
        let head = header_many(&[("chr1", 8), ("chr2", 8), ("chr10", 8)]);
        let map = ContigMap::build(&head, &refr).unwrap();
        assert_eq!(map.refr_tid(1), Some(2), "fixture does not reproduce the reordering");

        let mut rec = built(b"r", b"TTTT", &[Cigar::Match(4)], 0, 0);
        rec.set_tid(1); // chr2 in the BAM

        let mut seq = Vec::new();
        let mut cols = Vec::new();
        walk_alignment(&rec, &refr, &map, &mut seq, WalkOpts::default(), call(&rec), &mut |c| {
            cols.push(c)
        })
        .unwrap();

        let refr_bases: Vec<Seq> = cols.iter().map(|c| c.refr).collect();
        assert_eq!(refr_bases, vec![Seq::G; 4], "read chr10 (C) instead of chr2 (G)");
    }

    /// A record on a contig the reference does not have is a caller error, not
    /// a silent fallback to some other contig.
    #[test]
    fn a_contig_missing_from_the_reference_is_an_error_rather_than_a_guess() {
        let refr = reference_many("walk_missing", &[("chr1", "AAAAAAAA")]);
        let head = header_many(&[("chr1", 8), ("chrUn_scaffold", 8)]);
        let map = ContigMap::build(&head, &refr).unwrap();

        let mut rec = built(b"r", b"TTTT", &[Cigar::Match(4)], 0, 0);
        rec.set_tid(1);

        let mut seq = Vec::new();
        let err = walk_alignment(&rec, &refr, &map, &mut seq, WalkOpts::default(), call(&rec), &mut |_| {})
            .unwrap_err()
            .to_string();
        assert!(err.contains("not in the reference"), "{err}");
    }
}