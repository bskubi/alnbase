//! Read one record's aux block, and parse only the tags that a caller asks for.
//!
//! Two callers need this and they need different answers from it. The parquet
//! columns need the values of a few named tags; the tag writer needs to know
//! whether a tag is already present *and* whether the block can be walked at
//! all. One pass answers both, where `bam_aux_get` walks from the start once
//! per tag and cannot report damage.
//!
//! # Damage
//!
//! An aux block is a flat run of `TAG TYPE VALUE` with no framing, so a field
//! whose type is unknown or whose value is truncated makes everything after it
//! unreachable -- there is no boundary to resynchronise on. The walk stops
//! there and says so. That is also where htslib's `bam_aux_get` gives up, so
//! the two agree about which tags exist; what this adds is that the caller
//! finds out.

use rust_htslib::bam::Record;

/// One aux field, borrowed from the record.
///
/// `value` is everything after the type byte: the fixed-width bytes of a
/// number, the text of a `Z` or `H` without its NUL, or a `B` array's subtype,
/// count and elements.
#[derive(Clone, Copy, Debug)]
pub(crate) struct RawAux<'a> {
    pub(crate) tag: [u8; 2],
    pub(crate) ty: u8,
    pub(crate) value: &'a [u8],
}

/// The aux fields of a record, in file order, skipped over by size.
///
/// Nothing is parsed or validated beyond what finding the next tag needs, so
/// a tag nobody asked for costs a type lookup and, for text, a scan for its
/// NUL. The walk ends at the first field it cannot size -- an unknown type, a
/// truncated value -- because there is no way to find the tag after it. That
/// is also where `bam_aux_get` gives up, so a tag past a damaged field reads
/// as null either way.
pub(crate) struct AuxWalk<'a> {
    rest: &'a [u8],
    /// Set when the walk gave up on a field rather than running out of
    /// fields: the difference between a block that ended and one that broke.
    damaged: bool,
}

impl<'a> AuxWalk<'a> {
    /// The aux block sits after the name, CIGAR, packed sequence and
    /// qualities. The arithmetic is `Record::aux_iter`'s, which cannot be used
    /// directly because it parses and validates every field it passes.
    pub(crate) fn new(r: &'a Record) -> Self {
        let b = r.inner();
        if b.data.is_null() || b.l_data <= 0 {
            return Self { rest: &[], damaged: false };
        }
        // SAFETY: htslib keeps `data` valid for `l_data` bytes for as long as
        // the record is borrowed; `Record::data` builds the same slice.
        let data = unsafe { std::slice::from_raw_parts(b.data, b.l_data as usize) };
        let l_qseq = b.core.l_qseq.max(0) as usize;
        let start = b.core.l_qname as usize
            + b.core.n_cigar as usize * 4
            + l_qseq.div_ceil(2)
            + l_qseq;
        Self { rest: data.get(start..).unwrap_or(&[]), damaged: false }
    }
}

impl<'a> Iterator for AuxWalk<'a> {
    type Item = RawAux<'a>;

    #[inline]
    fn next(&mut self) -> Option<RawAux<'a>> {
        let buf = self.rest;
        // Ends the walk on anything unsizable, so the next call yields None.
        self.rest = &[];
        let give_up = |w: &mut Self| {
            // Bytes remained but could not be read as a field: damage. No
            // bytes at all is simply the end of the block.
            w.damaged = !buf.is_empty();
            None
        };
        if buf.len() < 3 {
            return give_up(self);
        }
        let tag = [buf[0], buf[1]];
        let ty = buf[2];
        let body = &buf[3..];
        let sized = (|| {
            Some(match ty {
                b'A' | b'c' | b'C' => (body.get(..1)?, 1),
                b's' | b'S' => (body.get(..2)?, 2),
                b'i' | b'I' | b'f' => (body.get(..4)?, 4),
                b'd' => (body.get(..8)?, 8),
                b'Z' | b'H' => {
                    let text = std::ffi::CStr::from_bytes_until_nul(body).ok()?.to_bytes();
                    (text, text.len() + 1)
                }
                b'B' => {
                    let count = u32::from_le_bytes(le(body.get(1..5)?)) as usize;
                    let width = match *body.first()? {
                        b'c' | b'C' => 1,
                        b's' | b'S' => 2,
                        b'i' | b'I' | b'f' => 4,
                        _ => return None,
                    };
                    let len = count.checked_mul(width)?.checked_add(5)?;
                    (body.get(..len)?, len)
                }
                _ => return None,
            })
        })();
        let Some((value, used)) = sized else {
            return give_up(self);
        };
        self.rest = &body[used..];
        Some(RawAux { tag, ty, value })
    }

}

/// A fixed-width little-endian value's bytes, as an array. Only ever called
/// with a slice `AuxWalk` sized for the type, so the length always matches.
#[inline]
pub(crate) fn le<const N: usize>(v: &[u8]) -> [u8; N] {
    v[..N].try_into().expect("AuxWalk sized this value for its type")
}

/// An integer-typed aux value, widened. `None` for any other type.
#[inline]
pub(crate) fn int_value(ty: u8, v: &[u8]) -> Option<i64> {
    Some(match ty {
        b'c' => v[0] as i8 as i64,
        b'C' => v[0] as i64,
        b's' => i16::from_le_bytes(le(v)) as i64,
        b'S' => u16::from_le_bytes(le(v)) as i64,
        b'i' => i32::from_le_bytes(le(v)) as i64,
        b'I' => u32::from_le_bytes(le(v)) as i64,
        _ => return None,
    })
}

impl AuxWalk<'_> {
    /// Whether the walk stopped at a field it could not read, rather than at
    /// the end of the block. Only meaningful once the walk is exhausted.
    pub(crate) fn stopped_early(&self) -> bool {
        self.damaged
    }
}

/// What one pass over a record's aux block found.
pub(crate) struct AuxScan {
    /// The walk stopped at a field it could not size, so tags after that point
    /// are unreachable to every reader.
    pub(crate) damaged: bool,
    /// Whether each wanted tag was present, in the order asked for.
    pub(crate) present: Vec<bool>,
}

/// Walk `record`'s aux block once, reporting damage and which of `wanted` is
/// present.
///
/// Replaces one `bam_aux_get` per wanted tag, which walks from the start each
/// time, and adds the answer that call cannot give.
pub(crate) fn scan(record: &Record, wanted: &[[u8; 2]]) -> AuxScan {
    let mut present = vec![false; wanted.len()];
    let mut walk = AuxWalk::new(record);
    let mut seen = 0usize;
    for field in walk.by_ref() {
        for (i, w) in wanted.iter().enumerate() {
            if !present[i] && *w == field.tag {
                present[i] = true;
                seen += 1;
            }
        }
        if seen == wanted.len() && !wanted.is_empty() {
            // Every wanted tag found, but the rest of the block still has to
            // be walked: damage after them is what decides whether a new tag
            // appended at the end would be reachable.
            continue;
        }
    }
    AuxScan { damaged: walk.stopped_early(), present }
}
