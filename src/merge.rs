use crate::{
    ruranges_structs::{GroupType, PositionType},
    sorts,
};

pub fn sweep_line_merge<G: GroupType, T: PositionType>(
    chrs: &[G],
    starts: &[T],
    ends: &[T],
    slack: T,
) -> (Vec<u32>, Vec<T>, Vec<T>, Vec<u32>) {
    assert_eq!(chrs.len(), starts.len());
    assert_eq!(starts.len(), ends.len());

    let mut out_indices = Vec::with_capacity(chrs.len());
    let mut out_starts = Vec::with_capacity(chrs.len());
    let mut out_ends = Vec::with_capacity(chrs.len());
    let mut counts = Vec::with_capacity(chrs.len());

    if chrs.is_empty() {
        return (out_indices, out_starts, out_ends, counts);
    }

    let intervals = sorts::build_sorted_intervals(chrs, starts, ends, None, T::zero(), false);

    let mut current_group = intervals[0].group;
    let mut current_idx = intervals[0].idx;
    let mut current_start = intervals[0].start;
    let mut current_end = intervals[0].end;
    let mut current_count = 1u32;

    for interval in intervals.iter().skip(1) {
        let within_merge_distance =
            interval.group == current_group && interval.start < current_end + slack;

        if within_merge_distance {
            if interval.end > current_end {
                current_end = interval.end;
            }
            current_count += 1;
            continue;
        }

        out_indices.push(current_idx);
        out_starts.push(current_start);
        out_ends.push(current_end);
        counts.push(current_count);

        current_group = interval.group;
        current_idx = interval.idx;
        current_start = interval.start;
        current_end = interval.end;
        current_count = 1;
    }

    out_indices.push(current_idx);
    out_starts.push(current_start);
    out_ends.push(current_end);
    counts.push(current_count);

    (out_indices, out_starts, out_ends, counts)
}

#[cfg(test)]
mod tests {
    use super::sweep_line_merge;

    #[test]
    fn merge_scans_sorted_intervals_and_keeps_first_index() {
        let groups = [1u32, 1, 1, 2];
        let starts = [10i32, 1, 4, 0];
        let ends = [14i32, 5, 7, 3];

        let (idx, merged_starts, merged_ends, counts) =
            sweep_line_merge(&groups, &starts, &ends, 0);

        assert_eq!(idx, vec![1, 0, 3]);
        assert_eq!(merged_starts, vec![1, 10, 0]);
        assert_eq!(merged_ends, vec![7, 14, 3]);
        assert_eq!(counts, vec![2, 1, 1]);
    }

    #[test]
    fn merge_requires_slack_for_bookended_intervals() {
        let groups = [1u32, 1];
        let starts = [1i32, 6];
        let ends = [6i32, 8];

        let without_slack = sweep_line_merge(&groups, &starts, &ends, 0);
        assert_eq!(without_slack.0, vec![0, 1]);
        assert_eq!(without_slack.1, vec![1, 6]);
        assert_eq!(without_slack.2, vec![6, 8]);
        assert_eq!(without_slack.3, vec![1, 1]);

        let with_slack = sweep_line_merge(&groups, &starts, &ends, 1);
        assert_eq!(with_slack.0, vec![0]);
        assert_eq!(with_slack.1, vec![1]);
        assert_eq!(with_slack.2, vec![8]);
        assert_eq!(with_slack.3, vec![2]);
    }
}
