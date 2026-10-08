use radsort::sort;
use rayon::prelude::*;

use crate::{
    partition,
    ruranges_structs::{GroupType, Interval, PositionType},
    sorts::build_sorted_intervals_par_with_runs,
};

pub fn max_disjoint<G, T>(
    groups: &[G],
    starts: &[T],
    ends: &[T],
    slack: T,
    sort_output: bool,
) -> Vec<u32>
where
    G: GroupType + Send + Sync,
    T: PositionType + Send + Sync,
{
    // Ensure the input slices all have the same length.
    assert_eq!(groups.len(), starts.len());
    assert_eq!(starts.len(), ends.len());

    // Build and sort intervals (group ➜ start ➜ end).
    let (intervals, runs) =
        build_sorted_intervals_par_with_runs(groups, starts, ends, None, slack, true);

    // Different groups are automatically disjoint, so each group is chosen from
    // on its own task and the results concatenated in group order — the order the
    // single serial scan produced them in.
    let batches = partition::batch_runs(&runs, partition::BATCH_ROWS, |&(s, e)| e - s);
    let chosen: Vec<Vec<u32>> = batches
        .into_par_iter()
        .map(|batch| {
            let mut chosen = Vec::new();
            for &(start, end) in &runs[batch] {
                choose_disjoint(&intervals[start..end], slack, &mut chosen);
            }
            chosen
        })
        .collect();
    let mut output: Vec<u32> = chosen.concat();

    if sort_output {
        sort(&mut output);
    }
    output
}

/// Greedily keep each interval of one start-sorted group that starts after the
/// last kept one ends (plus `slack`).
fn choose_disjoint<G: GroupType, T: PositionType>(
    intervals: &[Interval<G, T>],
    slack: T,
    out: &mut Vec<u32>,
) {
    let Some(first) = intervals.first() else {
        return;
    };
    let mut last_end = first.end;
    out.push(first.idx);
    for interval in intervals.iter().skip(1) {
        if interval.start > last_end + slack {
            last_end = interval.end;
            out.push(interval.idx);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sorts::build_sorted_intervals;

    /// The single serial scan the per-group sweep replaced.
    fn serial_reference(groups: &[u32], starts: &[i64], ends: &[i64], slack: i64) -> Vec<u32> {
        let intervals = build_sorted_intervals(groups, starts, ends, None, slack, true);
        let mut output = Vec::new();
        let Some(first) = intervals.first() else {
            return output;
        };
        let (mut current_group, mut last_end) = (first.group, first.end);
        output.push(first.idx);
        for interval in intervals.iter().skip(1) {
            if interval.group != current_group {
                current_group = interval.group;
                last_end = interval.end;
                output.push(interval.idx);
            } else if interval.start > last_end + slack {
                last_end = interval.end;
                output.push(interval.idx);
            }
        }
        output
    }

    #[test]
    fn per_group_sweep_matches_the_serial_scan() {
        let mut state = 0x4D59_5DF4_D0F3_3173_u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        // Either side of PARALLEL_THRESHOLD, few and many groups, both input orders.
        for &n in &[0_usize, 1, 999, 150_000] {
            for &n_groups in &[1_u64, 24, 50_000] {
                for sorted in [false, true] {
                    let mut rows: Vec<(u32, i64, i64)> = (0..n)
                        .map(|_| {
                            let start = (next() % 1_000_000) as i64;
                            (
                                (next() % n_groups) as u32,
                                start,
                                start + 1 + (next() % 500) as i64,
                            )
                        })
                        .collect();
                    if sorted {
                        rows.sort();
                    }
                    let groups: Vec<u32> = rows.iter().map(|r| r.0).collect();
                    let starts: Vec<i64> = rows.iter().map(|r| r.1).collect();
                    let ends: Vec<i64> = rows.iter().map(|r| r.2).collect();
                    for slack in [0, 5] {
                        assert_eq!(
                            max_disjoint(&groups, &starts, &ends, slack, false),
                            serial_reference(&groups, &starts, &ends, slack),
                            "n={n} groups={n_groups} sorted={sorted} slack={slack}"
                        );
                    }
                }
            }
        }
    }
}
