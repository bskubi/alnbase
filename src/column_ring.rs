//! Keep the columns that the walk emitted recently, so that a query can look
//! back at them.
//!
//! A query reports when its last column is the current one, but the coordinates
//! it wants — the anchor, and any captures — belong to earlier columns. The ring
//! holds the last `max_span` of them.
//!
//! # The back-index convention
//!
//! `store` post-increments, so immediately after storing a column `get(1)`
//! returns it. Counting further back, column `i` of a span of `n` is
//! `get(n - i)`; the first column of the span is `get(n)`.
//!
//! Getting this off by one does not fail, it shifts every reported coordinate
//! by one position — output that looks exactly like data. `Query` resolves the
//! back-indices once at compile time for that reason, and the tests below pin
//! the convention.

use crate::column::Column;

pub struct ColumnRing {
    columns: Box<[Column]>,
    /// Total columns stored since the last `clear`, not an index into
    /// `columns`. Masked on access.
    idx: usize,
    capacity: usize,
    mask: usize,
}

impl ColumnRing {
    /// Holds at least `n` columns; rounded up to a power of two so the index
    /// wrap is a mask rather than a division.
    pub fn new(n: usize) -> Self {
        let capacity = n.max(1).next_power_of_two();
        Self {
            columns: vec![Column::default(); capacity].into_boxed_slice(),
            idx: 0,
            capacity,
            mask: capacity - 1,
        }
    }

    /// The column `back` positions ago. `get(1)` is the one just stored.
    #[inline]
    pub fn get(&self, back: usize) -> &Column {
        debug_assert!(back >= 1, "get(0) is the slot about to be written");
        debug_assert!(back <= self.capacity, "looking back past the ring");
        debug_assert!(back <= self.idx, "looking back before the start of the record");
        &self.columns[(self.idx - back) & self.mask]
    }

    #[inline]
    pub fn store(&mut self, col: Column) {
        self.columns[self.idx & self.mask] = col;
        self.idx += 1;
    }

    /// Start a new record. Slots are not zeroed: nothing can read them until
    /// enough columns have been stored, which the `get` assertions enforce.
    #[inline]
    pub fn clear(&mut self) {
        self.idx = 0;
    }

    #[inline]
    #[cfg_attr(not(test), allow(dead_code))] // read by the ring-size test, which pins the rounding
    pub fn capacity(&self) -> usize {
        self.capacity
    }

    /// Columns stored since the last `clear`.
    #[inline]
    pub fn stored(&self) -> usize {
        self.idx
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::Seq;

    fn col(pos: i64) -> Column {
        Column { read: Seq::A, refr: Seq::A, refr_pos: pos, read_off: pos, qual: 30 }
    }

    #[test]
    fn get_one_is_the_most_recent() {
        let mut r = ColumnRing::new(4);
        r.store(col(10));
        assert_eq!(r.get(1).refr_pos, 10);
        r.store(col(11));
        assert_eq!(r.get(1).refr_pos, 11);
        assert_eq!(r.get(2).refr_pos, 10);
    }

    #[test]
    fn span_column_zero_is_get_span() {
        // A span of 4 firing on its last column: column 0 is get(4).
        let mut r = ColumnRing::new(4);
        for p in 0..4 {
            r.store(col(p));
        }
        for i in 0..4i64 {
            assert_eq!(r.get(4 - i as usize).refr_pos, i, "column {i}");
        }
    }

    #[test]
    fn wraps_without_losing_the_window() {
        let mut r = ColumnRing::new(4);
        for p in 0..100 {
            r.store(col(p));
        }
        assert_eq!(r.get(1).refr_pos, 99);
        assert_eq!(r.get(4).refr_pos, 96);
    }

    #[test]
    fn capacity_rounds_up() {
        assert_eq!(ColumnRing::new(5).capacity(), 8);
        assert_eq!(ColumnRing::new(8).capacity(), 8);
        assert_eq!(ColumnRing::new(1).capacity(), 1);
        // A zero-length span cannot happen, but must not panic here.
        assert_eq!(ColumnRing::new(0).capacity(), 1);
    }

    #[test]
    fn clear_resets_the_window() {
        let mut r = ColumnRing::new(4);
        r.store(col(1));
        r.clear();
        assert_eq!(r.stored(), 0);
        r.store(col(9));
        assert_eq!(r.get(1).refr_pos, 9);
    }
}