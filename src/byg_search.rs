//! Match many patterns at the same time, with one bit-parallel automaton.
//!
//! This is the shift-and algorithm of Baeza-Yates and Gonnet, from "A New
//! Approach to Text Searching" (CACM 1992). The module extends it to hold every
//! pattern in one bitmap.
//!
//! There is one bit per *trace*. A trace is a prefix of one pattern, and it is
//! active when that prefix matches the values that the walk has just read. A
//! pattern of length `n` owns `n` adjacent bits. The higher bits are the longer
//! prefixes, and the patterns lie end to end in the bitmap.
//!
//! To advance by one value the automaton does a small number of bitwise
//! operations over the bitmap, and it does the same operations for any number of
//! patterns. It shifts every trace up by one. It sets the start bit of every
//! pattern. It clears the traces that the new value contradicts. It then reads
//! the traces that have reached the last column of a pattern.
//!
//! One result of this design is important downstream. The *total* length of all
//! of the patterns sets the cost, and the number of patterns does not. When two
//! queries share a pattern, [`crate::query`] keeps one copy of it, and the
//! second query then adds no cost.
//!
//! A pattern that is longer than one word crosses a word boundary. The shift
//! carries the top bit of each word into the bottom of the next word. That carry
//! is the only part of this file where the word layout is visible.

use std::fmt;

#[derive(Debug)]
pub struct Plan {
    pub pattern_lengths: Box<[usize]>,
    pub end_bits: Box<[usize]>,
    /// The total number of trace bits over every pattern. Only the layout tests
    /// read this field. Those tests are what pin the end-to-end layout, so the
    /// field stays.
    #[cfg_attr(not(test), allow(dead_code))]
    pub bit_count: usize,
    pub word_count: usize,
    pub start: Bitmap,
    pub end: Bitmap,
    /// `sym_count * word_count` words, in symbol-major order. Symbol `s` owns
    /// the range `[s * word_count, (s + 1) * word_count)`. Inside that range,
    /// bit `b` says whether the observed value satisfies pattern column `b`.
    ///
    /// This is one allocation, and not one allocation per symbol. With the
    /// fourth flag bit the table holds 256x256 symbols. A `Box<[Bitmap]>` would
    /// then be 65 536 separate heap blocks, and for any run with 64 patterns or
    /// fewer their pointers weigh more than their contents. In the flat table
    /// the words of a frequent symbol sit in two or three cache lines, and
    /// `step` reaches them through one load and not through two dependent
    /// loads.
    compat: Box<[u64]>,
}

impl Plan {
    /// Lay `pattern_lengths` out end to end, with one bit per column.
    ///
    /// A pattern of length zero owns no bits and has no end bit, so it could
    /// never match. Such a pattern is a bug in the caller and not a case to
    /// handle here. The parser rejects an empty pattern long before this point.
    pub fn from_lengths(pattern_lengths: &[usize], sym_count: usize) -> Self {
        debug_assert!(
            pattern_lengths.iter().all(|&n| n > 0),
            "a zero-length pattern has no end bit and can never match"
        );
        let pattern_lengths: Box<[usize]> = Box::from(pattern_lengths);
        let bit_count = pattern_lengths.iter().sum::<usize>();
        let word_count = bit_count / 64 + 1;
        let mut start = Bitmap::new(word_count);
        let mut end = Bitmap::new(word_count);
        let compat = vec![0u64; sym_count * word_count].into_boxed_slice();
        let mut end_bits: Box<[usize]> = Box::from(vec![0; pattern_lengths.len()]);
        let mut cur = 0;
        for (pattern_idx, length) in pattern_lengths.iter().enumerate() {
            start.set(cur);
            cur += length - 1;
            end.set(cur);
            end_bits[pattern_idx] = cur;
            cur += 1;
        }
        Self { pattern_lengths, end_bits, bit_count, word_count, start, end, compat }
    }

    /// The compat words of one observed symbol.
    #[inline(always)]
    pub fn compat_words(&self, compat_idx: usize) -> &[u64] {
        let base = compat_idx * self.word_count;
        &self.compat[base..base + self.word_count]
    }

    /// Record that symbol `compat_idx` satisfies pattern column `bit`.
    ///
    /// This method is the only writer of the table. The symbol-major layout is
    /// therefore written down in one place, and `compat_words` cannot disagree
    /// with it.
    pub fn set_compat(&mut self, compat_idx: usize, bit: usize) {
        self.compat[compat_idx * self.word_count + bit / 64] |= 1u64 << (bit % 64);
    }

    /// Make a hit bitmap of the correct size for this plan. Use this method and
    /// do not build a bitmap by hand, so that one place derives the width.
    pub fn new_hits(&self) -> Hits {
        Hits { inner: Bitmap::new(self.word_count) }
    }

    /// Make a trace bitmap of the correct size for this plan, with every trace
    /// inactive.
    pub fn new_state(&self) -> State<'_> {
        State { plan: self, inner: Bitmap::new(self.word_count) }
    }
}

#[derive(Debug)]
pub struct Hits {
    pub inner: Bitmap
}


pub struct State<'p> {
    pub plan: &'p Plan,
    pub inner: Bitmap
}

impl<'p> State<'p> {
    /// Advance every trace by one observed value, and write the completed
    /// patterns to `hits`.
    ///
    /// `compat_idx` selects the words of a symbol through
    /// [`Plan::compat_words`]. Bit `b` of those words says whether the value
    /// that the walk has just read satisfies pattern column `b`.
    #[inline(always)]
    pub fn step(&mut self, compat_idx: usize, hits: &mut Hits) {
        // Top bit of the word below, waiting to become the bottom bit of this
        // one. Words ascend so the carry always flows into a word not yet read.
        let mut carry = 0u64;
        let compat_words = self.plan.compat_words(compat_idx);

        for w in 0..self.plan.word_count {
            let compat = compat_words[w];
            let cur = self.inner.words[w];
            let start = self.plan.start.words[w];
            let end = self.plan.end.words[w];

            // The top trace of this word shifts into the bottom of the next.
            // Taken before the shift: `cur`'s bit 63 is what belongs at bit 0
            // of the next word, whereas `active`'s bit 63 is a different trace
            // that has already moved this step.
            let carry_out = cur >> 63;

            // Shift extends every live prefix by one column; `| start` opens a
            // fresh trace at each pattern's first column; `& compat` kills the
            // traces this value contradicts.
            let active = ((cur << 1 | carry) | start) & compat;
            self.inner.words[w] = active;

            // `end` holds the last-column bit of every pattern, so this keeps
            // only the traces that just completed.
            hits.inner.words[w] = active & end;

            carry = carry_out;
        }
    }

    pub fn reset(&mut self) {
        self.inner.fill(0);
    }
}

struct ConcatenatedBits<'a>(&'a [u64]);

impl<'a> fmt::Debug for ConcatenatedBits<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut first_word = true;

        // Iterate words from highest index to lowest
        for word in self.0.iter().rev() {
            // Print the word separator if this isn't the very first word we're printing
            if !first_word {
                write!(f, " | ")?;
            }

            let mut first_byte = true;
            
            // Break the 64-bit word into 8 bytes
            for byte in word.to_be_bytes() {
                // Print the byte separator if this isn't the first byte of the current word
                if !first_byte {
                    write!(f, " ")?;
                }
                
                write!(f, "{:08b}", byte)?;
                first_byte = false;
            }
            
            first_word = false;
        }
        Ok(())
    }
}

#[derive(Clone)]
pub struct Bitmap {
    words: Box<[u64]>
}

impl Bitmap {
    pub fn new(nwords: usize) -> Self {
        let words = vec![0; nwords].into_boxed_slice();
        Self {words}
    }

    fn word_idx(bit: usize) -> usize {
        bit / 64
    }

    fn offset_idx(bit: usize) -> usize {
        bit % 64
    }

    fn word(&self, bit: usize) -> u64 {
        let w = Bitmap::word_idx(bit);
        self.words[w]
    }

    fn word_mut(&mut self, bit: usize) -> &mut u64 {
        let w = Bitmap::word_idx(bit);
        &mut self.words[w]
    }

    pub fn set(&mut self, bit: usize) {
        *self.word_mut(bit) |= 1 << Bitmap::offset_idx(bit);
    }

    fn fill(&mut self, value: u64) {
        self.words.fill(value);
    }

    pub fn get(&self, bit: usize) -> u64 {
        self.word(bit) & (1 << Bitmap::offset_idx(bit))
    }

    /// The raw words, for the predicate layer.
    ///
    /// These words and `Plan::end_bits` are the whole interface between the
    /// automaton and the boolean layer. `Plan::end_bits` crosses it at compile
    /// time, and these words cross it at evaluation time.
    #[inline]
    pub fn words(&self) -> &[u64] {
        &self.words
    }

    pub fn any(&self) -> bool {
        for &w in self.words.iter() {
            if w != 0 {
                return true;
            }
        }
        false
    }
}

impl fmt::Debug for Bitmap {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Bitmap")
         .field("words", &ConcatenatedBits(&self.words))
         .finish()
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    /// Build a plan for literal ASCII patterns. It fills `compat` in the way
    /// that `query::fill_compat` fills it: bit `b` of symbol `s` is set when
    /// pattern column `b` accepts `s`.
    fn plan_for(patterns: &[&str]) -> Plan {
        let lengths: Vec<usize> = patterns.iter().map(|p| p.len()).collect();
        let mut plan = Plan::from_lengths(&lengths, 128);
        let mut bit = 0;
        for p in patterns {
            for c in p.bytes() {
                plan.set_compat(c as usize, bit);
                bit += 1;
            }
        }
        plan
    }

    /// Report which patterns complete at each position of `text`.
    fn run(plan: &Plan, text: &str) -> Vec<Vec<usize>> {
        let mut state = plan.new_state();
        let mut hits = plan.new_hits();
        text.bytes()
            .map(|c| {
                state.step(c as usize, &mut hits);
                (0..plan.end_bits.len())
                    .filter(|&i| hits.inner.get(plan.end_bits[i]) != 0)
                    .collect()
            })
            .collect()
    }

    #[test]
    fn patterns_are_laid_out_end_to_end() {
        let plan = Plan::from_lengths(&[2, 4, 1], 4);
        assert_eq!(plan.bit_count, 7);
        // Last column of each pattern: 0+1, 2+3, 6.
        assert_eq!(&*plan.end_bits, &[1, 5, 6]);
        for (i, &start) in [0usize, 2, 6].iter().enumerate() {
            assert_ne!(plan.start.get(start), 0, "pattern {i} has no start bit");
            assert_ne!(plan.end.get(plan.end_bits[i]), 0, "pattern {i} has no end bit");
        }
        // A start bit is not an end bit unless the pattern is one column wide.
        assert_eq!(plan.end.get(0), 0);
        assert_ne!(plan.end.get(6), 0);
    }

    #[test]
    fn a_pattern_completes_at_its_last_column() {
        let plan = plan_for(&["CG"]);
        let fired = run(&plan, "TTCG");
        assert_eq!(fired[0], Vec::<usize>::new());
        assert_eq!(fired[1], Vec::<usize>::new());
        assert_eq!(fired[2], Vec::<usize>::new(), "the C alone is a partial trace");
        assert_eq!(fired[3], vec![0]);
    }

    #[test]
    fn overlapping_occurrences_each_report() {
        // A trace opens at every position, so runs report once per occurrence
        // rather than once per non-overlapping match.
        let plan = plan_for(&["AA"]);
        let fired = run(&plan, "AAAA");
        assert_eq!(fired[0], Vec::<usize>::new());
        assert_eq!(fired[1], vec![0]);
        assert_eq!(fired[2], vec![0]);
        assert_eq!(fired[3], vec![0]);
    }

    #[test]
    fn patterns_of_different_lengths_share_one_automaton() {
        let plan = plan_for(&["CG", "ACGT", "T"]);
        let fired = run(&plan, "ACGT");
        assert_eq!(fired[0], Vec::<usize>::new());
        assert_eq!(fired[1], Vec::<usize>::new());
        assert_eq!(fired[2], vec![0], "CG completes on the G");
        assert_eq!(fired[3], vec![1, 2], "ACGT and T both complete on the T");
    }

    /// The carry is the only place where the word layout is visible, and only a
    /// pattern that crosses a word boundary exercises it.
    #[test]
    fn a_trace_carries_across_a_word_boundary() {
        // 60 bits of filler, then a pattern occupying bits 60..69 -- so its
        // traces must cross from word 0 into word 1 to ever complete.
        let filler = "A".repeat(60);
        let plan = plan_for(&[&filler, "CGCGCGCGCG"]);
        assert_eq!(plan.bit_count, 70);
        assert!(plan.word_count >= 2, "the case needs two words to be a test");
        assert!(plan.end_bits[1] >= 64, "the straddling pattern must end in word 1");

        let fired = run(&plan, "TTCGCGCGCGCG");
        assert_eq!(fired[11], vec![1], "completes despite crossing the boundary");
        // And it does not fire one position early or late.
        assert_eq!(fired[10], Vec::<usize>::new());

        // The long pattern still works in the same automaton.
        let fired = run(&plan, &"A".repeat(60));
        assert_eq!(fired[59], vec![0]);
    }

    #[test]
    fn reset_discards_partial_traces() {
        let plan = plan_for(&["CG"]);
        let mut state = plan.new_state();
        let mut hits = plan.new_hits();

        // Half a match, then reset: the G must not complete the abandoned C.
        state.step(b'C' as usize, &mut hits);
        state.reset();
        state.step(b'G' as usize, &mut hits);
        assert_eq!(hits.inner.get(plan.end_bits[0]), 0);
        assert!(!hits.inner.any());

        // Without the reset it would have.
        let mut state = plan.new_state();
        state.step(b'C' as usize, &mut hits);
        state.step(b'G' as usize, &mut hits);
        assert_ne!(hits.inner.get(plan.end_bits[0]), 0);
    }

    #[test]
    fn an_unsatisfiable_symbol_kills_every_trace() {
        // A symbol no column accepts -- what an unmatched observation looks
        // like -- must leave the bitmap empty rather than merely unadvanced.
        let plan = plan_for(&["CG"]);
        let mut state = plan.new_state();
        let mut hits = plan.new_hits();
        state.step(b'C' as usize, &mut hits);
        state.step(b'X' as usize, &mut hits);
        assert!(!state.inner.any(), "no trace survives an incompatible value");
    }
}