use std::mem::MaybeUninit;

use rayon::prelude::*;

use crate::{
    partition,
    ruranges_structs::{GroupType, Interval, PositionType},
    sorts,
};

/// Label every interval with the id of the connected run it belongs to.
///
/// Returns `(cluster_ids, indices)`: `indices` names the input rows in
/// `(group, start, input order)` order, and `cluster_ids` gives each of them its
/// cluster.
///
/// This is the merge sweep with a per-row output instead of a per-run one — an
/// interval joins the run in progress exactly when `start < max_end + slack`,
/// which is the same test [`crate::merge::sweep_line_merge`] applies. It used to
/// be phrased as a start/end event sweep, which cost twice the elements to sort
/// (one event per interval endpoint instead of one record per interval) and left
/// the scan itself serial across the whole genome. Sorting intervals instead lets
/// each group's scan run on its own thread, because a cluster never spans two
/// groups.
///
/// Cluster ids are assigned exactly as the event sweep assigned them, gaps
/// included: ids run consecutively inside a group, and the counter skips one at
/// every group boundary. Reproducing that takes a serial prefix sum over the
/// per-group cluster counts, which is why the ids are written twice — once
/// locally, in parallel with the scan, then shifted by their group's base.
pub fn sweep_line_cluster<G: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    chrs: &[G],
    starts: &[T],
    ends: &[T],
    slack: T,
) -> (Vec<u32>, Vec<u32>) {
    if chrs.is_empty() {
        return (Vec::new(), Vec::new());
    }

    let (intervals, runs) =
        sorts::build_sorted_intervals_par_with_runs(chrs, starts, ends, None, T::zero(), false);

    cluster_sorted(&intervals, &runs, slack)
}

/// [`sweep_line_cluster`] over `i64` coordinates narrow enough to sweep as
/// `i32`, which sorts a third less data for the same answer.
///
/// `None` when the narrowing does not apply — see
/// [`sorts::build_sorted_narrowed_intervals_par`] — leaving the caller to run
/// [`sweep_line_cluster`] as it stands. Nothing has to be widened afterwards:
/// the output is ids and row numbers, not coordinates.
pub fn sweep_line_cluster_narrowed<G: GroupType + Send + Sync>(
    chrs: &[G],
    starts: &[i64],
    ends: &[i64],
    slack: i64,
) -> Option<(Vec<u32>, Vec<u32>)> {
    if !sorts::slack_fits_narrow(slack) {
        return None;
    }
    let (intervals, runs) = sorts::build_sorted_narrowed_intervals_par(chrs, starts, ends)?;

    Some(cluster_sorted(&intervals, &runs, slack as i32))
}

/// Label each interval of an already group-major, start-sorted array.
fn cluster_sorted<G: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    intervals: &[Interval<G, T>],
    runs: &[(usize, usize)],
    slack: T,
) -> (Vec<u32>, Vec<u32>) {
    let rows = intervals.len();
    let mut cluster_ids: Vec<u32> = Vec::with_capacity(rows);
    let mut indices: Vec<u32> = Vec::with_capacity(rows);
    // One entry per group: how many clusters that group's scan found.
    let mut cluster_counts: Vec<u32> = vec![0; runs.len()];

    let batches = partition::batch_runs(runs, partition::BATCH_ROWS, |&(start, end)| end - start);
    let row_spans: Vec<(usize, usize)> = batches
        .iter()
        .map(|batch| (runs[batch.start].0, runs[batch.end - 1].1))
        .collect();
    let count_spans: Vec<(usize, usize)> = batches
        .iter()
        .map(|batch| (batch.start, batch.end))
        .collect();

    partition::split_runs_mut(&mut cluster_ids.spare_capacity_mut()[..rows], &row_spans)
        .into_par_iter()
        .zip(partition::split_runs_mut(
            &mut indices.spare_capacity_mut()[..rows],
            &row_spans,
        ))
        .zip(partition::split_runs_mut(
            &mut cluster_counts,
            &count_spans,
        ))
        .zip(batches.par_iter())
        .for_each(|(((ids, idxs), counts), batch)| {
            let base = runs[batch.start].0;
            for (count, &(start, end)) in counts.iter_mut().zip(&runs[batch.clone()]) {
                let (from, to) = (start - base, end - base);
                *count = cluster_group(
                    &intervals[start..end],
                    slack,
                    &mut ids[from..to],
                    &mut idxs[from..to],
                );
            }
        });

    // SAFETY: the runs are contiguous and cover 0..rows, and every row in a run
    // was written by the task that owned it.
    unsafe {
        cluster_ids.set_len(rows);
        indices.set_len(rows);
    }

    // The event sweep bumped its counter once per finished cluster and once more
    // at every group boundary, so a group's first id is `sum(counts + 1)` over
    // the groups before it.
    let mut bases = Vec::with_capacity(runs.len());
    let mut next = 0_u32;
    for &count in &cluster_counts {
        bases.push(next);
        next += count + 1;
    }

    if bases.last().is_some_and(|&base| base != 0) {
        partition::split_runs_mut(&mut cluster_ids, &row_spans)
            .into_par_iter()
            .zip(batches.par_iter())
            .for_each(|(ids, batch)| {
                let row_base = runs[batch.start].0;
                for (&base, &(start, end)) in bases[batch.clone()].iter().zip(&runs[batch.clone()])
                {
                    if base == 0 {
                        continue;
                    }
                    for id in &mut ids[start - row_base..end - row_base] {
                        *id += base;
                    }
                }
            });
    }

    (cluster_ids, indices)
}

/// Scan one group's intervals, which are sorted by start, and label each with a
/// cluster id counting from zero. Returns the number of clusters found.
fn cluster_group<G: GroupType, T: PositionType>(
    intervals: &[Interval<G, T>],
    slack: T,
    ids: &mut [MaybeUninit<u32>],
    idxs: &mut [MaybeUninit<u32>],
) -> u32 {
    let Some((first, rest)) = intervals.split_first() else {
        return 0;
    };
    ids[0].write(0);
    idxs[0].write(first.idx);

    let mut cluster = 0_u32;
    let mut current_end = first.end;
    for ((interval, id), idx) in rest
        .iter()
        .zip(ids[1..].iter_mut())
        .zip(idxs[1..].iter_mut())
    {
        if interval.start < current_end + slack {
            if interval.end > current_end {
                current_end = interval.end;
            }
        } else {
            cluster += 1;
            current_end = interval.end;
        }
        id.write(cluster);
        idx.write(interval.idx);
    }

    cluster + 1
}

#[cfg(test)]
mod tests {
    use super::sweep_line_cluster;

    #[test]
    fn clusters_are_numbered_with_a_gap_at_every_group_boundary() {
        let groups = [1u32, 1, 1, 2, 2];
        let starts = [0i32, 4, 10, 0, 1];
        let ends = [5i32, 7, 12, 2, 3];

        let (ids, idx) = sweep_line_cluster(&groups, &starts, &ends, 0);

        // Group 1 has two clusters (0 and 1), group 2 picks up at 3.
        assert_eq!(idx, vec![0, 1, 2, 3, 4]);
        assert_eq!(ids, vec![0, 0, 1, 3, 3]);
    }

    #[test]
    fn bookended_intervals_need_slack_to_join() {
        let groups = [1u32, 1];
        let starts = [1i32, 6];
        let ends = [6i32, 8];

        assert_eq!(sweep_line_cluster(&groups, &starts, &ends, 0).0, vec![0, 1]);
        assert_eq!(sweep_line_cluster(&groups, &starts, &ends, 1).0, vec![0, 0]);
    }

    #[test]
    fn a_contained_interval_does_not_close_the_cluster() {
        let groups = [1u32, 1, 1];
        let starts = [0i32, 2, 8];
        let ends = [20i32, 4, 9];

        let (ids, idx) = sweep_line_cluster(&groups, &starts, &ends, 0);
        assert_eq!(idx, vec![0, 1, 2]);
        assert_eq!(ids, vec![0, 0, 0]);
    }
}

#[cfg(test)]
mod parity_tests {
    use super::sweep_line_cluster;
    use crate::ruranges_structs::Event;

    /// Deterministic SplitMix64, so this test needs no dev-dependency.
    fn next(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// The event sweep the interval sweep replaced, kept as the reference.
    fn cluster_by_events(
        chrs: &[u32],
        starts: &[i64],
        ends: &[i64],
        slack: i64,
    ) -> (Vec<u32>, Vec<u32>) {
        let mut indices = Vec::with_capacity(chrs.len());
        let mut cluster_ids = Vec::with_capacity(chrs.len());
        if chrs.is_empty() {
            return (cluster_ids, indices);
        }

        let events: Vec<Event<u32, i64>> =
            crate::sorts::build_sorted_events_single_collection(chrs, starts, ends, slack);
        let mut current_chr = events.first().unwrap().chr;
        let mut current_cluster = 0;
        let mut active_intervals = 0;

        for e in events {
            if e.chr != current_chr {
                current_cluster += 1;
                active_intervals = 0;
                current_chr = e.chr;
            }

            if e.is_start {
                indices.push(e.idx);
                cluster_ids.push(current_cluster);
                active_intervals += 1;
            } else {
                active_intervals -= 1;
                if active_intervals == 0 {
                    current_cluster += 1;
                }
            }
        }

        (cluster_ids, indices)
    }

    #[test]
    fn interval_sweep_matches_the_event_sweep() {
        let mut state = 0x51D5;
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
                    sweep_line_cluster(&chrs, &starts, &ends, slack),
                    cluster_by_events(&chrs, &starts, &ends, slack),
                    "n={n} groups={groups} slack={slack}"
                );
            }
        }
    }

    #[test]
    fn the_narrowed_sweep_matches_the_wide_one() {
        let mut state = 0x3F19;
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
                    super::sweep_line_cluster_narrowed(&chrs, &starts, &ends, slack),
                    Some(sweep_line_cluster(&chrs, &starts, &ends, slack)),
                    "groups={groups} slack={slack}"
                );
            }
        }
    }

    #[test]
    fn interval_sweep_matches_the_event_sweep_on_presorted_input() {
        // The presorted fast path in the interval builder skips the sort, so it
        // needs its own comparison.
        let mut state = 0x7AC3;
        let n = 150_000_usize;
        let mut chrs = Vec::with_capacity(n);
        let mut starts = Vec::with_capacity(n);
        let mut ends = Vec::with_capacity(n);
        for row in 0..n {
            chrs.push((row / 30_000) as u32);
            let start = (row % 30_000) as i64 * 3;
            starts.push(start);
            ends.push(start + 1 + (next(&mut state) % 8) as i64);
        }
        assert_eq!(
            sweep_line_cluster(&chrs, &starts, &ends, 0),
            cluster_by_events(&chrs, &starts, &ends, 0)
        );
    }
}
