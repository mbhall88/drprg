use std::cmp::{max, min};
use std::ops::Range;

/// A trait for defining interval operations (intended for use with Range)
pub trait IntervalOp {
    /// Returns the intersection between two Ranges, or None is they do not overlap
    fn intersect(&self, other: &Self) -> Option<Self>
    where
        Self: Sized;
}

impl<Idx: PartialOrd<Idx> + Ord + ToOwned<Owned = Idx>> IntervalOp for Range<Idx> {
    fn intersect(&self, other: &Self) -> Option<Self> {
        if other.start >= self.end || self.start >= other.end {
            None
        } else {
            let start = max(&self.start, &other.start).to_owned();
            let end = min(&self.end, &other.end).to_owned();
            Some(Range { start, end })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interval_op_intersect_no_intersection() {
        let i = 1..3;
        let j = 5..8;

        let overlap = i.intersect(&j);

        assert!(overlap.is_none())
    }

    #[test]
    fn interval_op_intersect_no_intersection_end_off_by_one_check() {
        let i = 1..3;
        let j = 3..8;

        let overlap = i.intersect(&j);

        assert!(overlap.is_none())
    }

    #[test]
    fn interval_op_intersect_no_intersection_start_off_by_one_check() {
        let i = 1..3;
        let j = 0..1;

        let overlap = i.intersect(&j);

        assert!(overlap.is_none())
    }

    #[test]
    fn interval_op_intersect_no_intersection_start_edge() {
        let i = 1..3;
        let j = 0..2;

        let overlap = i.intersect(&j);
        let expected = Some(1..2);

        assert_eq!(overlap, expected)
    }

    #[test]
    fn interval_op_intersect_no_intersection_end_edge() {
        let i = 1..3;
        let j = 2..6;

        let overlap = i.intersect(&j);
        let expected = Some(2..3);

        assert_eq!(overlap, expected)
    }

    #[test]
    fn interval_op_intersect_no_intersection_self() {
        let i = 1..3;
        let j = 1..3;

        let overlap = i.intersect(&j);
        let expected = Some(1..3);

        assert_eq!(overlap, expected)
    }

    #[test]
    fn interval_op_intersect_no_intersection_subset() {
        let i = 1..4;
        let j = 2..3;

        let overlap = i.intersect(&j);
        let expected = Some(2..3);

        assert_eq!(overlap, expected)
    }
}
