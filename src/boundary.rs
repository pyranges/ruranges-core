use crate::{
    ruranges_structs::{GroupType, PositionType},
    sorts,
};

/// Collapse each group's intervals into its outer boundary.
///
/// Returns one entry per group: the index of an input interval belonging to
/// that group, the group's minimum start, its maximum end, and how many input
/// intervals the group contains. The index identifies the group, not a
/// particular interval within it; callers must not depend on *which* member is
/// reported, only that it belongs to the group.
pub fn sweep_line_boundary<G: GroupType, T: PositionType>(
    chrs: &[G],
    starts: &[T],
    ends: &[T],
) -> (Vec<u32>, Vec<T>, Vec<T>, Vec<u32>) {
    let mut out_indices: Vec<u32> = Vec::with_capacity(chrs.len());
    let mut out_starts = Vec::with_capacity(chrs.len());
    let mut out_ends = Vec::with_capacity(chrs.len());
    let mut counts: Vec<u32> = Vec::with_capacity(chrs.len());

    if chrs.is_empty() {
        return (out_indices, out_starts, out_ends, counts);
    };

    let events = sorts::build_sorted_events_single_collection(chrs, starts, ends, T::zero());

    let mut current_chr = events.first().unwrap().chr;
    let mut current_start = events.first().unwrap().pos;
    let final_idx = events.last().unwrap().idx;
    let final_end = events.last().unwrap().pos;
    let mut prev_pos = T::zero();
    let mut prev_idx = 0;
    let mut current_count = 0_u32;

    for e in events {
        if e.chr != current_chr {
            current_chr = e.chr;
            out_indices.push(prev_idx);
            out_starts.push(current_start);
            out_ends.push(prev_pos);
            // Record the finished group before starting the next one's tally.
            counts.push(current_count);
            current_count = 0;
            current_start = e.pos;
        }

        // Each interval contributes a start and an end event, so counting only
        // the start events yields intervals rather than events.
        if e.is_start {
            current_count += 1;
        }

        prev_pos = e.pos;
        prev_idx = e.idx;
    }

    out_indices.push(final_idx);
    out_starts.push(current_start);
    out_ends.push(final_end);
    counts.push(current_count);

    (out_indices, out_starts, out_ends, counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reports_each_group_span_with_a_member_index() {
        // Groups listed in descending id order, each spanning two intervals.
        let groups = [1_u32, 1, 0, 0];
        let starts = [0_i64, 20, 0, 10];
        let ends = [1_i64, 21, 1, 11];

        let (indices, out_starts, out_ends, counts) = sweep_line_boundary(&groups, &starts, &ends);

        assert_eq!(out_starts, vec![0, 0]);
        assert_eq!(out_ends, vec![11, 21]);
        // The index only has to identify the group, never a specific member.
        assert_eq!(groups[indices[0] as usize], 0);
        assert_eq!(groups[indices[1] as usize], 1);
        // Every group is reported, including the ones before the last.
        assert_eq!(counts, vec![2, 2]);
    }

    #[test]
    fn counts_report_intervals_rather_than_sweep_events() {
        // One group of three intervals: three, not six.
        let (_indices, _starts, _ends, counts) =
            sweep_line_boundary(&[0_u32, 0, 0], &[0_i64, 10, 20], &[1_i64, 11, 21]);
        assert_eq!(counts, vec![3]);

        // Overlapping intervals still count once each.
        let (_indices, _starts, _ends, counts) =
            sweep_line_boundary(&[0_u32, 0], &[0_i64, 5], &[10_i64, 15]);
        assert_eq!(counts, vec![2]);
    }

    #[test]
    fn counts_track_groups_of_differing_size() {
        // Group 0 has one interval, group 1 has three.
        let groups = [0_u32, 1, 1, 1];
        let starts = [0_i64, 0, 10, 20];
        let ends = [5_i64, 1, 11, 21];

        let (_indices, _starts, _ends, counts) = sweep_line_boundary(&groups, &starts, &ends);

        assert_eq!(counts, vec![1, 3]);
    }

    #[test]
    fn empty_input_returns_empty_output() {
        let (indices, starts, ends, counts) =
            sweep_line_boundary(&[] as &[u32], &[] as &[i64], &[] as &[i64]);

        assert!(indices.is_empty());
        assert!(starts.is_empty());
        assert!(ends.is_empty());
        assert!(counts.is_empty());
    }
}
