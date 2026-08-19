use radsort::sort_by_key;

use crate::overlaps::{
    pair_fits_narrow, sorted_records_narrowed_par, sorted_records_par, GroupRun, IntervalRecord,
    SortedRecords,
};
use crate::partition;
use crate::ruranges_structs::{GroupType, MinInterval, PositionType};
use rayon::prelude::*;

/// Subtract set 2 from set 1, returning the surviving pieces of each set-1
/// interval as `(idx into set 1, start, end)`.
///
/// Subtracting a set is the same as subtracting its union, so this walks each
/// group's set-2 intervals once into a disjoint cover and then, per set-1
/// interval, emits the gaps in that cover. Two coordinate sorts over `n1` and
/// `n2` records, and every group independent of every other.
///
/// The previous implementation built `2 * (n1 + n2)` sweep events, radix-sorted
/// them three times, and tracked live set-1 intervals in an `FxHashMap` that it
/// scanned in full on every 0↔1 transition of the set-2 depth — so its cost grew
/// with (set-1 pileup depth × set-2 transitions), which is worst exactly where
/// the data is densest.
///
/// Intervals of non-positive length are ignored on the right (they cover
/// nothing) and produce no output on the left (nothing survives a subtraction
/// from an empty span). `polaranges` rejects them before the kernel is reached;
/// direct callers used to get incoherent results for them, because such an
/// interval's end event sorted before its own start event.
pub fn sweep_line_subtract<G: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    chrs1: &[G],
    starts1: &[T],
    ends1: &[T],
    chrs2: &[G],
    starts2: &[T],
    ends2: &[T],
    sort_output: bool,
) -> (Vec<u32>, Vec<T>, Vec<T>) {
    // If either set is empty, set1 is unchanged (or trivially subtracted).
    if chrs1.is_empty() || chrs2.is_empty() {
        return (
            (0..chrs1.len() as u32).collect(),
            starts1.to_vec(),
            ends1.to_vec(),
        );
    }

    let left = sorted_records_par(chrs1, starts1, ends1, false);
    let right = sorted_records_par(chrs2, starts2, ends2, false);
    subtract_from_sorted(&left, &right, sort_output)
}

/// [`sweep_line_subtract`] over `i64` coordinates narrow enough to sweep as
/// `i32`, which sorts a third less data for the same answer.
///
/// `None` when the narrowing does not apply — see
/// [`crate::overlaps::sorted_records_narrowed_par`] — leaving the caller to run
/// [`sweep_line_subtract`] as it stands.
pub fn sweep_line_subtract_narrowed<G: GroupType + Send + Sync>(
    chrs1: &[G],
    starts1: &[i64],
    ends1: &[i64],
    chrs2: &[G],
    starts2: &[i64],
    ends2: &[i64],
    sort_output: bool,
) -> Option<(Vec<u32>, Vec<i64>, Vec<i64>)> {
    if chrs1.is_empty() || chrs2.is_empty() {
        return None;
    }
    // Subtraction carries no slack, so only the coordinates have to fit.
    if !pair_fits_narrow(starts1, ends1, starts2, ends2, 0) {
        return None;
    }

    let left = sorted_records_narrowed_par(chrs1, starts1, ends1, false)?;
    let right = sorted_records_narrowed_par(chrs2, starts2, ends2, false)?;
    let (idxs, starts, ends) = subtract_from_sorted(&left, &right, sort_output);

    // Every surviving piece lies inside its input interval, so widening back is
    // one pass over the output — far less than the narrow sort saved.
    Some((
        idxs,
        starts.into_iter().map(i64::from).collect(),
        ends.into_iter().map(i64::from).collect(),
    ))
}

/// Subtract right from left, given both sides already group-major and
/// start-sorted.
fn subtract_from_sorted<G: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    left: &SortedRecords<G, T>,
    right: &SortedRecords<G, T>,
    sort_output: bool,
) -> (Vec<u32>, Vec<T>, Vec<T>) {
    let paired = left_runs_with_right(left, right);

    let batches = partition::batch_runs(&paired, partition::BATCH_ROWS, |&(l, r)| {
        (l.1 - l.0) + r.map_or(0, |run| run.1 - run.0)
    });

    let collected: Vec<Vec<MinInterval<T>>> = batches
        .into_par_iter()
        .map(|batch| {
            let mut out = Vec::new();
            // One cover buffer per task, reused across the groups in the batch.
            let mut cover: Vec<(T, T)> = Vec::new();
            for &(left_run, right_run) in &paired[batch] {
                let left_group = &left.records[left_run.0..left_run.1];
                match right_run {
                    Some(run) => subtract_group(
                        left_group,
                        &right.records[run.0..run.1],
                        &mut cover,
                        &mut out,
                    ),
                    // No set-2 interval on this group: every set-1 interval
                    // survives whole.
                    None => out.extend(left_group.iter().map(|q| MinInterval {
                        start: q.start,
                        end: q.end,
                        idx: q.idx,
                    })),
                }
            }
            out
        })
        .collect();

    let mut out_events: Vec<MinInterval<T>> =
        Vec::with_capacity(collected.iter().map(Vec::len).sum());
    for batch in collected {
        out_events.extend(batch);
    }

    if sort_output {
        // Stable, so an input interval's surviving pieces stay in ascending
        // coordinate order — the order they were emitted in.
        sort_by_key(&mut out_events, |i| i.idx);
    }

    let mut out_idxs = Vec::with_capacity(out_events.len());
    let mut out_starts = Vec::with_capacity(out_events.len());
    let mut out_ends = Vec::with_capacity(out_events.len());

    for rec in out_events {
        out_idxs.push(rec.idx);
        out_starts.push(rec.start);
        out_ends.push(rec.end);
    }

    (out_idxs, out_starts, out_ends)
}

/// Each left group run paired with the right run of the same group, if any.
///
/// Unlike the overlap sweeps, subtract has to visit left groups the right side
/// does not have: those intervals survive untouched.
fn left_runs_with_right<C: GroupType, T: PositionType>(
    left: &SortedRecords<C, T>,
    right: &SortedRecords<C, T>,
) -> Vec<(GroupRun, Option<GroupRun>)> {
    let mut out = Vec::with_capacity(left.runs.len());
    let mut j = 0usize;
    for &left_run in &left.runs {
        let group = left.records[left_run.0].group;
        while j < right.runs.len() && right.records[right.runs[j].0].group < group {
            j += 1;
        }
        let matched = right
            .runs
            .get(j)
            .filter(|run| right.records[run.0].group == group)
            .copied();
        out.push((left_run, matched));
    }
    out
}

/// Subtract one group's right intervals from its left intervals.
///
/// `cover` is scratch owned by the caller so the allocation survives across
/// groups; it is cleared on entry.
fn subtract_group<C: GroupType, T: PositionType>(
    left: &[IntervalRecord<C, T>],
    right: &[IntervalRecord<C, T>],
    cover: &mut Vec<(T, T)>,
    out: &mut Vec<MinInterval<T>>,
) {
    // The union of the right intervals, as a disjoint ascending cover. `right`
    // is already sorted by start, so this is one linear pass. Touching intervals
    // are merged too: the zero-length gap between them would never be emitted.
    cover.clear();
    for target in right {
        if target.end <= target.start {
            continue;
        }
        match cover.last_mut() {
            Some(last) if target.start <= last.1 => {
                if target.end > last.1 {
                    last.1 = target.end;
                }
            }
            _ => cover.push((target.start, target.end)),
        }
    }

    // `left` is sorted by start, so the first cover interval that can matter
    // only ever moves forward.
    let mut first = 0usize;
    for query in left {
        while first < cover.len() && cover[first].1 <= query.start {
            first += 1;
        }

        let mut cursor = query.start;
        let mut index = first;
        while index < cover.len() && cover[index].0 < query.end {
            if cover[index].0 > cursor {
                out.push(MinInterval {
                    start: cursor,
                    end: cover[index].0,
                    idx: query.idx,
                });
            }
            cursor = cover[index].1;
            if cursor >= query.end {
                break;
            }
            index += 1;
        }

        if cursor < query.end {
            out.push(MinInterval {
                start: cursor,
                end: query.end,
                idx: query.idx,
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::sweep_line_subtract;

    type Group = u32;
    type Pos = i64;

    /// Deterministic SplitMix64, so these tests need no dev-dependency.
    fn next(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn random_intervals(
        state: &mut u64,
        n: usize,
        groups: u64,
        span: u64,
        max_len: u64,
    ) -> (Vec<Group>, Vec<Pos>, Vec<Pos>) {
        let mut chrs = Vec::with_capacity(n);
        let mut starts = Vec::with_capacity(n);
        let mut ends = Vec::with_capacity(n);
        for _ in 0..n {
            chrs.push((next(state) % groups) as Group);
            let start = (next(state) % span) as Pos;
            let len = 1 + (next(state) % max_len) as Pos;
            starts.push(start);
            ends.push(start + len);
        }
        (chrs, starts, ends)
    }

    /// Subtraction computed one coordinate at a time. Independent of the kernel,
    /// so it pins the result rather than restating the algorithm.
    fn brute_force(
        chrs: &[Group],
        starts: &[Pos],
        ends: &[Pos],
        chrs2: &[Group],
        starts2: &[Pos],
        ends2: &[Pos],
    ) -> Vec<(u32, Pos, Pos)> {
        let mut out = Vec::new();
        for q in 0..chrs.len() {
            let mut open: Option<Pos> = None;
            for pos in starts[q]..ends[q] {
                let covered = (0..chrs2.len()).any(|t| {
                    chrs2[t] == chrs[q] && starts2[t] <= pos && pos < ends2[t]
                });
                match (covered, open) {
                    (true, Some(from)) => {
                        out.push((q as u32, from, pos));
                        open = None;
                    }
                    (false, None) => open = Some(pos),
                    _ => {}
                }
            }
            if let Some(from) = open {
                out.push((q as u32, from, ends[q]));
            }
        }
        out
    }

    fn run(
        chrs: &[Group],
        starts: &[Pos],
        ends: &[Pos],
        chrs2: &[Group],
        starts2: &[Pos],
        ends2: &[Pos],
    ) -> Vec<(u32, Pos, Pos)> {
        let (idx, s, e) = sweep_line_subtract(chrs, starts, ends, chrs2, starts2, ends2, true);
        idx.into_iter().zip(s).zip(e).map(|((i, s), e)| (i, s, e)).collect()
    }

    #[test]
    fn subtract_matches_brute_force_on_small_random_input() {
        let mut state = 0x5EED;
        for _ in 0..40 {
            let (chrs, starts, ends) = random_intervals(&mut state, 40, 3, 120, 30);
            let (chrs2, starts2, ends2) = random_intervals(&mut state, 40, 3, 120, 30);
            assert_eq!(
                run(&chrs, &starts, &ends, &chrs2, &starts2, &ends2),
                brute_force(&chrs, &starts, &ends, &chrs2, &starts2, &ends2)
            );
        }
    }

    #[test]
    fn nested_and_touching_targets_are_covered_once() {
        // [5, 100) fully contains [10, 20); [20, 30) touches it. The union is
        // [5, 100), so only the head and tail of the query survive.
        let chrs = [0_u32, 0];
        let starts = [0_i64, 200];
        let ends = [150_i64, 300];
        let chrs2 = [0_u32, 0, 0];
        let starts2 = [5_i64, 10, 20];
        let ends2 = [100_i64, 20, 30];

        assert_eq!(
            run(&chrs, &starts, &ends, &chrs2, &starts2, &ends2),
            vec![(0, 0, 5), (0, 100, 150), (1, 200, 300)]
        );
    }

    #[test]
    fn a_left_group_with_no_right_group_survives_whole() {
        let chrs = [0_u32, 7];
        let starts = [0_i64, 0];
        let ends = [10_i64, 10];
        let chrs2 = [0_u32];
        let starts2 = [2_i64];
        let ends2 = [4_i64];

        assert_eq!(
            run(&chrs, &starts, &ends, &chrs2, &starts2, &ends2),
            vec![(0, 0, 2), (0, 4, 10), (1, 0, 10)]
        );
    }

    #[test]
    fn unsorted_output_holds_the_same_rows_as_sorted_output() {
        let mut state = 0xB0B;
        let (chrs, starts, ends) = random_intervals(&mut state, 5_000, 4, 20_000, 200);
        let (chrs2, starts2, ends2) = random_intervals(&mut state, 5_000, 4, 20_000, 200);

        let mut sorted = run(&chrs, &starts, &ends, &chrs2, &starts2, &ends2);
        let (idx, s, e) =
            sweep_line_subtract(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, false);
        let mut unsorted: Vec<(u32, i64, i64)> =
            idx.into_iter().zip(s).zip(e).map(|((i, s), e)| (i, s, e)).collect();

        sorted.sort_unstable();
        unsorted.sort_unstable();
        assert_eq!(sorted, unsorted);
    }

    /// Above `PARALLEL_THRESHOLD` the records are partitioned and the groups are
    /// swept on separate threads; the result must not depend on that.
    #[test]
    fn parallel_and_serial_sized_inputs_agree() {
        let mut state = 0xFACE;
        let (chrs, starts, ends) = random_intervals(&mut state, 150_000, 5, 400_000, 500);
        let (chrs2, starts2, ends2) = random_intervals(&mut state, 150_000, 5, 400_000, 500);

        let big = run(&chrs, &starts, &ends, &chrs2, &starts2, &ends2);

        // Same data, one group at a time, below the parallel threshold.
        let mut expected: Vec<(u32, i64, i64)> = Vec::new();
        for group in 0..5_u32 {
            let rows: Vec<usize> = (0..chrs.len()).filter(|&i| chrs[i] == group).collect();
            let rows2: Vec<usize> = (0..chrs2.len()).filter(|&i| chrs2[i] == group).collect();
            let (lc, ls, le): (Vec<u32>, Vec<i64>, Vec<i64>) = (
                rows.iter().map(|&i| chrs[i]).collect(),
                rows.iter().map(|&i| starts[i]).collect(),
                rows.iter().map(|&i| ends[i]).collect(),
            );
            let (rc, rs, re): (Vec<u32>, Vec<i64>, Vec<i64>) = (
                rows2.iter().map(|&i| chrs2[i]).collect(),
                rows2.iter().map(|&i| starts2[i]).collect(),
                rows2.iter().map(|&i| ends2[i]).collect(),
            );
            for (i, s, e) in run(&lc, &ls, &le, &rc, &rs, &re) {
                expected.push((rows[i as usize] as u32, s, e));
            }
        }
        expected.sort_unstable();

        let mut got = big;
        got.sort_unstable();
        assert_eq!(got, expected);
    }
}
