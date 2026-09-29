//! Store a base as a set of possible values, and not as a letter. A `Seq` is one
//! byte. It holds four nucleotide bits and four flag bits. All 256 values are
//! legal, and each one has a canonical name. You can therefore write a pattern
//! for any of them, and the code needs no special cases.
//!
//! ```text
//! bit  0  1  2  3   4    5    6    7
//!      A  C  G  T   GAP  PAD  CLIP SKIP
//! ```
//!
//! # The four flag bits are four different absences
//!
//! Each flag says that something is missing. They are easy to confuse, and they
//! do not mean the same thing.
//!
//! - `GAP` — the read is aligned across this reference base, and it does not
//!   carry the base. This is a deletion. The read covers the position and
//!   reports nothing there.
//! - `PAD` — the position is outside the alignment. It is in the flank, or it is
//!   reference past the end of a contig.
//! - `SKIP` — the position is inside the alignment, and a `CIGAR N` skips it.
//!   For RNA this is an intron. The reference has the base, and the read does
//!   not span it.
//! - `CLIP` — the read continues here, and the aligner soft-clipped it. The
//!   sequencer read the base, and the aligner did not align it. The column
//!   therefore shows the reference base and no read base. `CLIP` appears only on
//!   the read side, in the flank next to a soft clip.
//!
//! Keep `GAP` and `SKIP` apart. A `D` and an `N` have the same shape in a CIGAR,
//! and they make different claims. A deletion says that the base is absent from
//! the sequenced fragment. An intron says that the base was spliced out before
//! sequencing. The two were one value until a pattern could tell them apart. See
//! `alignment::walk_alignment`.
//!
//! # Naming
//!
//! A value has two spellings, because two jobs need different things.
//!
//! [`name`](Seq::name) is the canonical spelling, and you write patterns in it.
//! The nucleotide bits use IUPAC (`A C G T R Y S W K M B D H V N`). The flag
//! bits use punctuation: `.` is GAP, `_` is PAD, `:` is CLIP, and `,` is SKIP. A
//! value that holds more than one of these categories is braced, so `{C.}` is a
//! `C` or a deletion. Two values have a spelling of their own. `~` is `{N._:J}`,
//! which is everything a column can hold. `0` is the empty set.
//! [`from_name`](Seq::from_name) parses a name back, so a pattern is text that
//! makes the round trip. The read side and the reference side join with `@`, as
//! in `C~@CG`.
//!
//! Inside a `{}` group, SKIP is spelled `J` and not `,`. The parser rejects a
//! comma in a group. The reason is that `{A,C,G}` is the natural way to guess at
//! a union, and it would otherwise parse as A, C, G *and a junction*. That set
//! is close to the set the user wanted. It behaves the same way on DNA, and it
//! differs only on spliced data. The parser therefore refuses the comma, and the
//! user sees the mistake instead of a set that is quietly wrong.
//!
//! [`to_ascii`](Seq::to_ascii) is the other spelling. It is exactly one byte for
//! every value, which is what a column-aligned grid or a trace dump needs. It
//! uses `Z` for GAP, `X` for PAD, `L` for CLIP and `J` for SKIP. IUPAC does not
//! use these four letters. A value that mixes bases with flags becomes `*`,
//! because it has no form of one character.
//!
//! `to_ascii` uses the letters and not the punctuation because `.` already has
//! another meaning. In a `mark` row, `.` means "do not capture this column", and
//! a `mark` row sits directly above the `read` and `refr` rows that it applies
//! to.
//!
//! [`from_name`](Seq::from_name) accepts both spellings. `Z`, `.` and `-` all
//! parse as GAP. `X` and `_` parse as PAD. `L` and `:` parse as CLIP. `J` and
//! `,` parse as SKIP. The two spellings therefore differ only in the output.

use std::fmt::Write;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, bytemuck::Pod, bytemuck::Zeroable)]
#[repr(transparent)]
pub struct Seq(pub u8);

/// IUPAC letter for each 4-bit nucleotide subset, indexed by the subset itself.
const IUPAC: [u8; 16] = [
    b'0', b'A', b'C', b'M', b'G', b'R', b'S', b'V', b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N',
];

impl Seq {
    pub const EMPTY: Seq = Seq(0);
    pub const A: Seq = Seq(0b0000_0001);
    pub const C: Seq = Seq(0b0000_0010);
    pub const G: Seq = Seq(0b0000_0100);
    pub const T: Seq = Seq(0b0000_1000);
    pub const GAP: Seq = Seq(0b0001_0000);
    pub const PAD: Seq = Seq(0b0010_0000);
    /// A soft-clipped read position in the flank: sequenced, not aligned.
    pub const CLIP: Seq = Seq(0b0100_0000);
    /// Reference skipped by a `CIGAR N`: an intron, in RNA.
    pub const SKIP: Seq = Seq(0b1000_0000);

    // IUPAC ambiguity codes — unions of the base bits.
    pub const R: Seq = Seq(0b0000_0101); // A G
    pub const Y: Seq = Seq(0b0000_1010); // C T
    pub const S: Seq = Seq(0b0000_0110); // C G
    pub const W: Seq = Seq(0b0000_1001); // A T
    pub const K: Seq = Seq(0b0000_1100); // G T
    pub const M: Seq = Seq(0b0000_0011); // A C
    pub const B: Seq = Seq(0b0000_1110); // C G T
    pub const D: Seq = Seq(0b0000_1101); // A G T
    pub const H: Seq = Seq(0b0000_1011); // A C T
    pub const V: Seq = Seq(0b0000_0111); // A C G
    pub const N: Seq = Seq(0b0000_1111); // A C G T

    // Unions with the flag bits. `GAP_PAD` is `{._}`, matching only where the
    // read has no base.
    pub const GAP_PAD: Seq = Seq(0b0011_0000);
    /// Everything that a column the read spans can hold. This is not `~`, which
    /// also admits SKIP. The relational codes `=` and `/` constrain this set,
    /// because two bases that are not there have no relation.
    pub const N_GAP_PAD: Seq = Seq(0b0011_1111);
    /// `~`: everything that a column can hold, including junctions and clips.
    pub const ANY: Seq = Seq(0b1111_1111);

    /// Mask of the four nucleotide bits.
    pub const BASES: u8 = 0b0000_1111;
    /// Mask of the GAP/PAD/CLIP/SKIP bits.
    pub const FLAGS: u8 = 0b1111_0000;

    /// Every `u8` is a legal value now that the fourth flag bit is in use.
    pub const COUNT: usize = 256;

    /// Every legal value, in ascending bit order.
    ///
    /// Ranges over `u8` rather than `0..COUNT as u8`, which would be `0..0` now
    /// that `COUNT` is 256.
    pub fn all() -> impl Iterator<Item = Seq> {
        (0..=u8::MAX).map(Seq)
    }

    /// All eight bits are now assigned, so every bit pattern is a legal set.
    /// This function therefore always returns `true`. It has a name because
    /// callers assert on it, and because the rule is worth stating.
    #[inline]
    pub const fn is_valid(self) -> bool {
        (self.0 as usize) < Self::COUNT
    }

    // ---------------------------------------------------------------- naming

    /// Canonical name, e.g. `C`, `{C.}`, `_`, `~`, `0`.
    ///
    /// One or more characters. For the fixed-width form, see
    /// [`to_ascii`](Self::to_ascii).
    pub fn name(self) -> String {
        let mut s = String::with_capacity(6);
        self.write_name(&mut s);
        s
    }

    /// Append the canonical name without allocating.
    pub fn write_name(self, out: &mut String) {
        
        if self.0 == 0 {
            out.push('0');
            return;
        }
        if self.0 == Self::ANY.0 {
            out.push('~');
            return;
        }
        let bases = self.0 & Self::BASES;
        let has_bases = bases != 0;
        let has_gap = self.0 & Self::GAP.0 != 0;
        let has_pad = self.0 & Self::PAD.0 != 0;
        let has_clip = self.0 & Self::CLIP.0 != 0;
        let has_skip = self.0 & Self::SKIP.0 != 0;
        let has_brackets = (has_bases as u8
            + has_gap as u8
            + has_pad as u8
            + has_clip as u8
            + has_skip as u8)
            >= 2;
        if has_brackets {
            out.push('{')
        }
        if has_bases {
            out.push(IUPAC[bases as usize] as char);
        }
        if has_gap {
            out.push('.');
        }
        if has_pad {
            out.push('_');
        }
        if has_clip {
            out.push(':');
        }
        // `J` inside a group, `,` alone: a comma in a group is rejected, so the
        // braced form has to use the letter or it would not parse back.
        if has_skip {
            out.push(if has_brackets { 'J' } else { ',' });
        }
        if has_brackets {
            out.push('}')
        }
    }

    /// Parse a canonical name. Letters may appear in any order; `0` is empty.
    pub fn from_name(s: &str) -> Option<Seq> {
        if s.is_empty() {
            return None;
        }
        let mut v: u8 = 0;
        for c in s.bytes() {
            v |= match c.to_ascii_uppercase() {
                b'0'|b'{'|b'}' => Seq::EMPTY.0,
                b'Z'|b'.'|b'-' => Seq::GAP.0,
                b'X'|b'_' => Seq::PAD.0,
                b'L' | b':' => Seq::CLIP.0,
                b'J' | b',' => Seq::SKIP.0,
                b'~' => Seq::ANY.0,
                b'U' => Seq::T.0, // RNA
                other => {
                    let mut found = 0u8;
                    let mut i = 1usize;
                    while i < 16 {
                        if IUPAC[i] == other {
                            found = i as u8;
                            break;
                        }
                        i += 1;
                    }
                    if found == 0 {
                        return None;
                    }
                    found
                }
            };
        }
        Some(Seq(v))
    }

    /// One character, for a compact alignment dump. Grids, traces and anything
    /// else that needs a column of exactly one byte use this form.
    ///
    /// `Z` is GAP, `X` is PAD, and `L` is CLIP. [`name`](Self::name) uses `.`,
    /// `_` and `:` for the same three values. This function uses the letters
    /// because `.` is already the "do not capture" marker in a `mark` row. A
    /// value that mixes bases with flags becomes `*`, because it has no form of
    /// one character. This form therefore loses information, and it is not the
    /// inverse of [`from_name`](Self::from_name).
    #[inline]
    pub const fn to_ascii(self) -> u8 {
        Self::TO_ASCII[self.0 as usize]
    }

    #[inline]
    pub const fn to_char(self) -> char {
        self.to_ascii() as char
    }

    const TO_ASCII: [u8; 256] = {
        let mut t = [b'?'; 256];
        let mut v: usize = 0;
        while v < 256 {
            let bases = (v as u8) & Seq::BASES;
            let flags = (v as u8) & Seq::FLAGS;
            t[v] = if v == 0 {
                b'0'
            } else if flags == 0 {
                IUPAC[bases as usize]
            } else if bases == 0 && flags == Seq::GAP.0 {
                b'Z'
            } else if bases == 0 && flags == Seq::PAD.0 {
                b'X'
            } else if bases == 0 && flags == Seq::CLIP.0 {
                b'L'
            } else if bases == 0 && flags == Seq::SKIP.0 {
                b'J'
            } else {
                b'*'
            };
            v += 1;
        }
        t
    };

    /// Decode one character of reference sequence. The character is an IUPAC
    /// letter or `U`, in upper case or lower case. Any other character gives
    /// `None`, and `alnbase index` then refuses the FASTA. A reference that
    /// holds `*`, `-` or a digit is the wrong file, or it is damaged, and no
    /// `Seq` value states that correctly. This table is not
    /// [`from_name`](Self::from_name), which parses pattern syntax.
    pub const FROM_FASTA: [Option<Seq>; 256] = {
        let mut t = [None; 256];
        let mut i: usize = 1;
        while i < 16 {
            let c = IUPAC[i];
            t[c as usize] = Some(Seq(i as u8));
            t[(c | 0x20) as usize] = Some(Seq(i as u8));
            i += 1;
        }
        t[b'U' as usize] = Some(Seq::T);
        t[b'u' as usize] = Some(Seq::T);
        t
    };

    // ------------------------------------------------------------- set logic

    #[inline]
    pub fn complement(self) -> Seq {
        let n = self.0 & Self::BASES;
        let r = ((n & 0b0001) << 3)
            | ((n & 0b0010) << 1)
            | ((n & 0b0100) >> 1)
            | ((n & 0b1000) >> 3);
        Seq(r | (self.0 & !Self::BASES)) // preserve GAP/PAD/CLIP/SKIP
    }

    /// All combinations of set bits
    pub fn bit_combinations(self) -> Vec<Seq> {
        let mask = self.0;
        let mut out = Vec::with_capacity(1 << mask.count_ones());
        let mut s = mask;
        loop {
            out.push(Seq(s));
            if s == 0 {
                break;
            }
            s = (s - 1) & mask;
        }
        out.pop(); // drop the empty set
        out
    }

    /// Returns one sequence per set bit in the mask.
    pub fn individual_bits(self) -> Vec<Seq> {
        let mut mask = self.0;
        let mut out = Vec::with_capacity(mask.count_ones() as usize);
        
        while mask != 0 {
            // Find the index of the lowest set bit
            let bit_idx = mask.trailing_zeros();
            
            // Push that single bit as a new Seq
            out.push(Seq(1 << bit_idx));
            
            // Clear the lowest set bit to advance the loop
            mask &= mask - 1; 
        }
        
        out
    }

    #[inline] pub fn is_subset_of(self, pat: Seq) -> bool { self.0 & !pat.0 == 0 }
    #[inline] pub fn intersects(self, o: Seq) -> bool { self.0 & o.0 != 0 }
    #[inline] pub fn is_empty(self) -> bool { self.0 == 0 }
    #[inline] pub fn ambiguity(self) -> u32 { self.0.count_ones() }
    #[inline] pub fn is_ambiguous(self) -> bool { self.ambiguity() > 1 }
    #[inline] pub fn is_gappy(self) -> bool { self.intersects(Seq::GAP) }
}

macro_rules! bitops {
    ($($trait:ident::$method:ident = $op:tt),* $(,)?) => {
        $(impl std::ops::$trait for Seq {
            type Output = Seq;
            #[inline]
            fn $method(self, rhs: Seq) -> Seq { Seq(self.0 $op rhs.0) }
        })*
    };
}
bitops!(BitAnd::bitand = &, BitOr::bitor = |, BitXor::bitxor = ^);

macro_rules! bitops_assign {
    ($($trait:ident::$method:ident = $op:tt),* $(,)?) => {
        $(impl std::ops::$trait for Seq {
            #[inline]
            fn $method(&mut self, rhs: Seq) { self.0 $op rhs.0; }
        })*
    };
}
bitops_assign!(BitAndAssign::bitand_assign = &=, BitOrAssign::bitor_assign = |=);

impl std::ops::Not for Seq {
    type Output = Seq;
    #[inline]
    fn not(self) -> Seq { Seq(!self.0 & (Seq::BASES | Seq::FLAGS)) }
}

/// Write the one-character form, for an alignment dump.
impl std::fmt::Display for Seq {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_char(self.to_char())
    }
}

impl Default for Seq {
    /// The empty set. It matches nothing, and it differs from every real
    /// observation. A ring slot that nothing has written yet holds this value.
    fn default() -> Self {
        Seq::EMPTY
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_round_trip() {
        for s in Seq::all() {
            let n = s.name();
            assert_eq!(Seq::from_name(&n), Some(s), "failed on {n} ({:#010b})", s.0);
        }
    }

    #[test]
    fn known_names() {
        assert_eq!(Seq::C.name(), "C");
        assert_eq!(Seq::ANY.name(), "~");
        assert_eq!(Seq::GAP.name(), ".");
        assert_eq!(Seq::GAP_PAD.name(), "{._}");
        assert_eq!(Seq::EMPTY.name(), "0");
        assert_eq!((Seq::C | Seq::GAP).name(), "{C.}");
        assert_eq!(Seq::CLIP.name(), ":");
        assert_eq!((Seq::CLIP | Seq::PAD).name(), "{_:}");
        assert_eq!(Seq::SKIP.name(), ",");
        // `,` alone but `J` in a group, because a comma in a group is rejected
        // and the braced form has to parse back.
        assert_eq!((Seq::N | Seq::SKIP).name(), "{NJ}");
        assert_eq!(Seq::N_GAP_PAD.name(), "{N._}", "no longer `~`: that admits SKIP");
    }

    /// Both spellings of every flag value parse, and the two forms of SKIP name
    /// the same set. `~` is the one code whose meaning widened.
    #[test]
    fn junction_names_round_trip() {
        assert_eq!(Seq::from_name(":"), Some(Seq::CLIP));
        assert_eq!(Seq::from_name("L"), Some(Seq::CLIP));
        assert_eq!(Seq::from_name("%"), None, "ERR is gone");
        assert!(Seq::CLIP.is_subset_of(Seq::ANY), "`~` admits a clip");
        assert_eq!(Seq::from_name(","), Some(Seq::SKIP));
        assert_eq!(Seq::from_name("J"), Some(Seq::SKIP));
        assert_eq!(Seq::from_name("j"), Some(Seq::SKIP));
        assert_eq!(Seq::from_name("~"), Some(Seq::ANY));
        assert!(Seq::SKIP.is_subset_of(Seq::ANY), "`~` admits a junction");
        assert!(
            !Seq::SKIP.is_subset_of(Seq::N_GAP_PAD),
            "the relational codes must not match a junction marker"
        );
        for s in Seq::all() {
            assert_eq!(Seq::from_name(&s.name()), Some(s), "{} did not round trip", s.name());
        }
    }

    /// Every value has a one-byte form, and the four flag values have distinct
    /// ones -- a grid column is one character and cannot be ambiguous.
    #[test]
    fn the_flag_values_have_distinct_ascii() {
        assert_eq!(Seq::GAP.to_char(), 'Z');
        assert_eq!(Seq::PAD.to_char(), 'X');
        assert_eq!(Seq::CLIP.to_char(), 'L');
        assert_eq!(Seq::SKIP.to_char(), 'J');
        assert_eq!((Seq::C | Seq::SKIP).to_char(), '*', "mixed values have no one-byte form");
    }

    #[test]
    fn ascii_decode() {
        assert_eq!(Seq::FROM_FASTA[b'a' as usize], Some(Seq::A));
        assert_eq!(Seq::FROM_FASTA[b'N' as usize], Some(Seq::N));
        assert_eq!(Seq::FROM_FASTA[b'h' as usize], Some(Seq::H));
        assert_eq!(Seq::FROM_FASTA[b'u' as usize], Some(Seq::T));
        for bad in [b'-', b'.', b'0', b'*', b'X', b'E', b'J', b'Z', b'%', b'1'] {
            assert_eq!(Seq::FROM_FASTA[bad as usize], None, "{}", bad as char);
        }
    }

    #[test]
    fn complement_preserves_flags() {
        assert_eq!(Seq::C.complement(), Seq::G);
        assert_eq!(Seq::N_GAP_PAD.complement(), Seq::N_GAP_PAD);
        assert_eq!((Seq::C | Seq::PAD).complement(), Seq::G | Seq::PAD);
        assert_eq!(Seq::GAP.complement(), Seq::GAP);
    }
}