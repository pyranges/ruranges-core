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

    // Build and sort intervals (group ➜ start ➜ end), unshifted: `slack` is
    // applied once, in the overlap test.
    let (intervals, runs) =
        build_sorted_intervals_par_with_runs(groups, starts, ends, None, T::zero(), true);

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

/// Choose a maximum set of pairwise disjoint intervals from one group, sorted
/// by start and then end.
///
/// Two intervals overlap when the later one starts before the earlier one ends
/// plus `slack`, as in `overlaps`: at slack 0 book-ended intervals are disjoint.
///
/// Walking in start order, an interval disjoint from the last one kept is kept.
/// One that overlaps it but ends earlier takes its place: it starts no earlier,
/// so it is still disjoint from everything kept before, and it leaves more room
/// for what follows. This is the earliest-end-first greedy, which is optimal.
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
        if interval.start >= last_end + slack {
            last_end = interval.end;
            out.push(interval.idx);
        } else if interval.end < last_end {
            last_end = interval.end;
            *out.last_mut().expect("each group keeps its first interval") = interval.idx;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sorts::build_sorted_intervals;

    /// The same choice made in one serial scan over every group.
    fn serial_reference(groups: &[u32], starts: &[i64], ends: &[i64], slack: i64) -> Vec<u32> {
        let intervals = build_sorted_intervals(groups, starts, ends, None, 0, true);
        let mut output = Vec::new();
        let mut from = 0;
        while from < intervals.len() {
            let group = intervals[from].group;
            let to = from
                + intervals[from..]
                    .iter()
                    .take_while(|i| i.group == group)
                    .count();
            choose_disjoint(&intervals[from..to], slack, &mut output);
            from = to;
        }
        output
    }

    /// The largest number of pairwise disjoint rows in one group, by brute force.
    fn largest_disjoint_count(starts: &[i64], ends: &[i64], slack: i64) -> u32 {
        let n = starts.len();
        let disjoint = |a: usize, b: usize| {
            let (first, second) = if starts[a] <= starts[b] {
                (a, b)
            } else {
                (b, a)
            };
            starts[second] >= ends[first] + slack
        };
        (0_u32..1 << n)
            .filter(|set| {
                let rows: Vec<usize> = (0..n).filter(|&r| set >> r & 1 == 1).collect();
                rows.iter()
                    .enumerate()
                    .all(|(i, &a)| rows[i + 1..].iter().all(|&b| disjoint(a, b)))
            })
            .map(u32::count_ones)
            .max()
            .unwrap_or(0)
    }

    #[test]
    fn chooses_a_maximum_disjoint_set() {
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..3_000 {
            let n = (next() % 11) as usize;
            let starts: Vec<i64> = (0..n).map(|_| (next() % 40) as i64).collect();
            let ends: Vec<i64> = starts
                .iter()
                .map(|&s| s + 1 + (next() % 15) as i64)
                .collect();
            for slack in [0, 1, 3] {
                let chosen = max_disjoint(&vec![0_u32; n], &starts, &ends, slack, true);
                let rows: Vec<usize> = chosen.iter().map(|&r| r as usize).collect();
                for (i, &a) in rows.iter().enumerate() {
                    for &b in &rows[i + 1..] {
                        let (first, second) = if starts[a] <= starts[b] {
                            (a, b)
                        } else {
                            (b, a)
                        };
                        assert!(starts[second] >= ends[first] + slack, "{a} and {b} overlap");
                    }
                }
                assert_eq!(
                    chosen.len() as u32,
                    largest_disjoint_count(&starts, &ends, slack),
                    "starts={starts:?} ends={ends:?} slack={slack}"
                );
            }
        }
    }

    #[test]
    fn slack_and_book_ended_intervals_follow_overlaps() {
        // Both kept exactly when the gap between them is at least `slack`.
        for slack in 0..4_i64 {
            for gap in 0..8_i64 {
                let kept = max_disjoint(&[0_u32; 2], &[0, 10 + gap], &[10, 20 + gap], slack, true);
                assert_eq!(
                    kept.len(),
                    if gap >= slack { 2 } else { 1 },
                    "slack={slack} gap={gap}"
                );
            }
        }
        // A long interval does not crowd out the two it contains.
        assert_eq!(
            max_disjoint(&[0_u32; 3], &[0, 10, 30], &[100, 20, 40], 0, true),
            [1, 2]
        );
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
