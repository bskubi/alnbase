//! Define one alignment column. The code that makes a column lives elsewhere.
//!
//! The consumers are the ring, the scanner, the tracer and the overlap resolver.
//! They depend on this module and not on `alignment`, which needs a BAM library
//! to walk a CIGAR. Most of the crate therefore builds and tests without htslib.

use crate::seq::Seq;

/// The `read_off` value of a flank column, which has no read offset. The value
/// lies outside SEQ. A consumer that forgets to check for [`Seq::PAD`]
/// therefore cannot use it to index a real base.
pub const PAD_READ_OFF: i64 = i64::MIN;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Column {
    /// The read base. It is [`Seq::GAP`] for a deletion. In the flank it is
    /// [`Seq::CLIP`] for a soft-clipped base, and [`Seq::PAD`] past the read.
    pub read: Seq,
    /// The reference base. It is [`Seq::GAP`] for an insertion. It is
    /// [`Seq::PAD`] where the coordinate has run off the end of the contig. It
    /// is [`Seq::SKIP`] in the one marker column that represents the interior of
    /// a long intron.
    pub refr: Seq,
    /// The reference coordinate. Where the column has none, this is the
    /// coordinate of the nearest 5' base.
    pub refr_pos: i64,
    /// The offset into the SEQ of the record. For a deletion column or an
    /// intron column it is the offset of the nearest 5' base. For a clip it is
    /// the offset of the clipped base itself. A pad has no offset and carries
    /// [`PAD_READ_OFF`].
    ///
    /// This is not the `off_5p` or `off_3p` value that the outputs report.
    /// Those values count from the end of the read *as the sequencer read it*,
    /// and SEQ does not include the hard-clipped bases. `crate::hits` therefore
    /// adds the hard clip back on when it builds a row.
    pub read_off: i64,
    /// The Phred quality of the read base. It is -1 where the column has no
    /// read base, which is the case for a deletion and for the flank. It is also
    /// -1 when the record carries no quality string.
    pub qual: i32,
}

impl Column {
    /// Whether the read observed this column. A matched read base, an inserted
    /// read base and a deletion are all observations. A flank pad, a clip and an
    /// intron-context column are not observations. They are reference that the
    /// walk shows around the alignment. The reference concordance check uses
    /// this flag, because it compares only what the read saw.
    #[inline]
    pub fn is_observed(&self) -> bool {
        self.read != Seq::PAD && self.read != Seq::CLIP && self.read != Seq::SKIP
    }
}