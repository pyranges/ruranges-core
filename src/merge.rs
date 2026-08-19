use rayon::prelude::*;

use crate::{
    partition,
    ruranges_structs::{GroupType, Interval, PositionType},
    sorts,
};

/// `(indices, starts, ends, counts)`, one entry per merged run: the row that
/// opened it, its span, and how many input intervals it swallowed.
pub type MergedRuns<T> = (Vec<u32>, Vec<T>, Vec<T>, Vec<u32>);

pub fn sweep_line_merge<G: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    chrs: &[G],
    starts: &[T],
    ends: &[T],
    slack: T,
) -> MergedRuns<T> {
    assert_eq!(chrs.len(), starts.len());
    assert_eq!(starts.len(), ends.len());

    if chrs.is_empty() {
        return (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    }

    // The builder hands back the group runs because it already knows them;
    // recovering them here would be a serial pass over the whole sorted array.
    let (intervals, runs) =
        sorts::build_sorted_intervals_par_with_runs(chrs, starts, ends, None, T::zero(), false);

    merge_sorted(&intervals, &runs, slack)
}

/// [`sweep_line_merge`] over `i64` coordinates narrow enough to sweep as `i32`,
/// which sorts a third less data for the same answer.
///
/// `None` when the narrowing does not apply — see
/// [`sorts::build_sorted_narrowed_intervals_par`] — leaving the caller to run
/// [`sweep_line_merge`] as it stands.
pub fn sweep_line_merge_narrowed<G: GroupType + Send + Sync>(
    chrs: &[G],
    starts: &[i64],
    ends: &[i64],
    slack: i64,
) -> Option<MergedRuns<i64>> {
    assert_eq!(chrs.len(), starts.len());
    assert_eq!(starts.len(), ends.len());

    if !sorts::slack_fits_narrow(slack) {
        return None;
    }
    let (intervals, runs) = sorts::build_sorted_narrowed_intervals_par(chrs, starts, ends)?;
    let (indices, merged_starts, merged_ends, counts) =
        merge_sorted(&intervals, &runs, slack as i32);

    // One row per merged run, so widening the coordinates back costs nothing
    // like what the narrow sort saved.
    Some((
        indices,
        merged_starts.into_iter().map(i64::from).collect(),
        merged_ends.into_iter().map(i64::from).collect(),
        counts,
    ))
}

/// Merge each group of an already group-major, start-sorted interval array.
///
/// A merged run never spans two groups, so each group's slice can be swept on
/// its own thread and the results concatenated in group order — the order the
/// single serial scan produced them in.
fn merge_sorted<G: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    intervals: &[Interval<G, T>],
    runs: &[(usize, usize)],
    slack: T,
) -> MergedRuns<T> {
    let batches = partition::batch_runs(runs, partition::BATCH_ROWS, |&(s, e)| e - s);
    let collected: Vec<MergedGroup<T>> = batches
        .into_par_iter()
        .map(|batch| {
            let mut merged = MergedGroup::default();
            for &(start, end) in &runs[batch] {
                merge_group(&intervals[start..end], slack, &mut merged);
            }
            merged
        })
        .collect();

    let total = collected.iter().map(|m| m.indices.len()).sum();
    let mut out_indices = Vec::with_capacity(total);
    let mut out_starts = Vec::with_capacity(total);
    let mut out_ends = Vec::with_capacity(total);
    let mut counts = Vec::with_capacity(total);
    for merged in collected {
        out_indices.extend(merged.indices);
        out_starts.extend(merged.starts);
        out_ends.extend(merged.ends);
        counts.extend(merged.counts);
    }

    (out_indices, out_starts, out_ends, counts)
}

/// One task's merged runs, in group order.
struct MergedGroup<T: PositionType> {
    indices: Vec<u32>,
    starts: Vec<T>,
    ends: Vec<T>,
    counts: Vec<u32>,
}

impl<T: PositionType> Default for MergedGroup<T> {
    fn default() -> Self {
        Self {
            indices: Vec::new(),
            starts: Vec::new(),
            ends: Vec::new(),
            counts: Vec::new(),
        }
    }
}

/// Merge one group's intervals, which are sorted by start. Each output run keeps
/// the index of the first input interval that opened it.
fn merge_group<G: GroupType, T: PositionType>(
    intervals: &[Interval<G, T>],
    slack: T,
    out: &mut MergedGroup<T>,
) {
    let Some(first) = intervals.first() else {
        return;
    };
    let mut current_idx = first.idx;
    let mut current_start = first.start;
    let mut current_end = first.end;
    let mut current_count = 1u32;

    for interval in intervals.iter().skip(1) {
        if interval.start < current_end + slack {
            if interval.end > current_end {
                current_end = interval.end;
            }
            current_count += 1;
            continue;
        }

        out.indices.push(current_idx);
        out.starts.push(current_start);
        out.ends.push(current_end);
        out.counts.push(current_count);

        current_idx = interval.idx;
        current_start = interval.start;
        current_end = interval.end;
        current_count = 1;
    }

    out.indices.push(current_idx);
    out.starts.push(current_start);
    out.ends.push(current_end);
    out.counts.push(current_count);
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

#[cfg(test)]
mod parallel_tests {
    use super::sweep_line_merge;

    /// Deterministic SplitMix64, so this test needs no dev-dependency.
    fn next(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// The serial scan the parallel sweep replaced, kept as the reference.
    fn merge_serially(
        chrs: &[u32],
        starts: &[i64],
        ends: &[i64],
        slack: i64,
    ) -> (Vec<u32>, Vec<i64>, Vec<i64>, Vec<u32>) {
        let intervals =
            crate::sorts::build_sorted_intervals(chrs, starts, ends, None, 0, false);
        let (mut idx, mut out_starts, mut out_ends, mut counts) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        if intervals.is_empty() {
            return (idx, out_starts, out_ends, counts);
        }
        let mut group = intervals[0].group;
        let mut current = (
            intervals[0].idx,
            intervals[0].start,
            intervals[0].end,
            1u32,
        );
        for interval in intervals.iter().skip(1) {
            if interval.group == group && interval.start < current.2 + slack {
                current.2 = current.2.max(interval.end);
                current.3 += 1;
                continue;
            }
            idx.push(current.0);
            out_starts.push(current.1);
            out_ends.push(current.2);
            counts.push(current.3);
            group = interval.group;
            current = (interval.idx, interval.start, interval.end, 1);
        }
        idx.push(current.0);
        out_starts.push(current.1);
        out_ends.push(current.2);
        counts.push(current.3);
        (idx, out_starts, out_ends, counts)
    }

    #[test]
    fn the_narrowed_sweep_matches_the_wide_one() {
        let mut state = 0x2C0E;
        // Above PARALLEL_THRESHOLD, or the narrowed builder declines.
        let n = 150_000_usize;
        for &groups in &[24_u64, 20_000] {
            for &slack in &[0_i64, 50] {
                let mut chrs = Vec::with_capacity(n);
                let mut starts = Vec::with_capacity(n);
                let mut ends = Vec::with_capacity(n);
                for _ in 0..n {
                    chrs.push((next(&mut state) % groups) as u32);
                    let start = (next(&mut state) % 250_000_000) as i64;
                    starts.push(start);
                    ends.push(start + 1 + (next(&mut state) % 400) as i64);
                }
                assert_eq!(
                    super::sweep_line_merge_narrowed(&chrs, &starts, &ends, slack),
                    Some(sweep_line_merge(&chrs, &starts, &ends, slack)),
                    "groups={groups} slack={slack}"
                );
            }
        }
    }

    #[test]
    fn a_coordinate_too_wide_to_narrow_declines() {
        let n = 150_000_usize;
        let chrs = vec![0_u32; n];
        let mut starts: Vec<i64> = (0..n as i64).collect();
        let mut ends: Vec<i64> = starts.iter().map(|start| start + 10).collect();
        starts[n / 2] = crate::sorts::NARROW_COORDINATE_LIMIT + 1;
        ends[n / 2] = crate::sorts::NARROW_COORDINATE_LIMIT + 11;

        assert!(super::sweep_line_merge_narrowed(&chrs, &starts, &ends, 0).is_none());
        // And a slack that could overflow the halved range is declined too.
        let small: Vec<i64> = (0..n as i64).collect();
        let small_ends: Vec<i64> = small.iter().map(|start| start + 10).collect();
        assert!(super::sweep_line_merge_narrowed(
            &chrs,
            &small,
            &small_ends,
            crate::sorts::NARROW_COORDINATE_LIMIT
        )
        .is_none());
    }

    #[test]
    fn parallel_merge_matches_the_serial_scan() {
        let mut state = 0x1D0D;
        // Either side of PARALLEL_THRESHOLD, and both a handful of groups and
        // many, which is what decides the batching.
        for &(n, groups) in &[(500_usize, 3_u64), (150_000, 24), (150_000, 20_000)] {
            for &slack in &[0_i64, 50] {
                let mut chrs = Vec::with_capacity(n);
                let mut starts = Vec::with_capacity(n);
                let mut ends = Vec::with_capacity(n);
                for _ in 0..n {
                    chrs.push((next(&mut state) % groups) as u32);
                    let start = (next(&mut state) % 500_000) as i64;
                    starts.push(start);
                    ends.push(start + 1 + (next(&mut state) % 400) as i64);
                }
                assert_eq!(
                    sweep_line_merge(&chrs, &starts, &ends, slack),
                    merge_serially(&chrs, &starts, &ends, slack),
                    "n={n} groups={groups} slack={slack}"
                );
            }
        }
    }
}
