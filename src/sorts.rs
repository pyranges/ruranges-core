use radsort::sort_by_key;
use rayon::prelude::*;

use crate::partition;
use crate::ruranges_structs::Event;
use crate::ruranges_structs::GenericEvent;
use crate::ruranges_structs::GroupStruct;
use crate::ruranges_structs::GroupType;
use crate::ruranges_structs::Interval;
use crate::ruranges_structs::MaxEvent;
use crate::ruranges_structs::MinEvent;
use crate::ruranges_structs::PositionType;
use crate::ruranges_structs::SplicedSubsequenceInterval;
use crate::ruranges_structs::SubsequenceInterval;

pub fn build_intervals<C: GroupType, T: PositionType>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    sort_reverse_direction: Option<&[bool]>,
    slack: T,
) -> Vec<Interval<C, T>> {
    let mut intervals = Vec::with_capacity(chrs.len());
    match sort_reverse_direction {
        Some(reverse) => {
            for i in 0..chrs.len() {
                intervals.push(Interval {
                    group: chrs[i],
                    start: if reverse[i] {
                        -(starts[i] - slack)
                    } else {
                        starts[i] - slack
                    },
                    end: if reverse[i] {
                        -(ends[i] + slack)
                    } else {
                        ends[i] + slack
                    },
                    idx: i as u32,
                });
            }
        }
        None => {
            for i in 0..chrs.len() {
                intervals.push(Interval {
                    group: chrs[i],
                    start: starts[i] - slack,
                    end: ends[i] + slack,
                    idx: i as u32,
                });
            }
        }
    };

    intervals
}

pub fn build_subsequence_intervals<G: GroupType, T: PositionType>(
    chrs: &[G],
    starts: &[T],
    ends: &[T],
    strand_flags: &[bool],
) -> Vec<SplicedSubsequenceInterval<G, T>> {
    let mut intervals = Vec::with_capacity(chrs.len());
    for i in 0..chrs.len() {
        intervals.push(SplicedSubsequenceInterval {
            chr: chrs[i],
            start: if strand_flags[i] {
                starts[i]
            } else {
                -starts[i]
            }, // so that negative strand intervals are sorted in the correct direction
            end: if strand_flags[i] { ends[i] } else { -ends[i] }, // we will find the absolute value when using them
            idx: i as u32,
            forward_strand: strand_flags[i],
            temp_cumsum: T::zero(),
            temp_length: T::zero(),
        });
    }

    intervals
}

pub fn build_sequence_intervals(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    strand_flags: &[bool],
    force_plus_strand: bool,
) -> Vec<SubsequenceInterval> {
    let mut intervals: Vec<SubsequenceInterval> = Vec::with_capacity(chrs.len());
    for i in 0..chrs.len() {
        intervals.push(SubsequenceInterval {
            group_id: chrs[i],
            start: if force_plus_strand || strand_flags[i] {
                starts[i]
            } else {
                -starts[i]
            }, // so that negative strand intervals are sorted in the correct direction
            end: if force_plus_strand || strand_flags[i] {
                ends[i]
            } else {
                -ends[i]
            }, // we will find the absolute value when using them
            idx: idxs[i],
            forward_strand: strand_flags[i],
        });
    }

    intervals
}

pub fn build_sorted_intervals<C: GroupType, T: PositionType>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    sort_reverse_direction: Option<&[bool]>,
    slack: T,
    sort_on_ends_too: bool,
) -> Vec<Interval<C, T>> {
    let mut intervals = build_intervals(chrs, starts, ends, sort_reverse_direction, slack);

    if sort_on_ends_too {
        sort_by_key(&mut intervals, |i| i.end);
    };
    sort_by_key(&mut intervals, |i| i.start);
    sort_by_key(&mut intervals, |i| i.group);

    intervals
}

/// Half-open offsets of each group's run in a group-major array.
pub type GroupRuns = Vec<(usize, usize)>;

/// A sorted interval array together with each group's run in it.
pub type SortedIntervals<C, T> = (Vec<Interval<C, T>>, GroupRuns);

/// Parallel counterpart of [`build_sorted_intervals`] with **byte-identical**
/// output, used by the merge and max-disjoint kernels.
///
/// A parallel counting partition makes each chromosome a contiguous run, then
/// each run's coordinate sort runs on its own rayon task — serial radix within a
/// chromosome, parallel across chromosomes (the same shape as the merge sweep).
/// radsort is a stable LSD radix, so equal keys keep input/`idx` order, matching
/// the serial end ➜ start ➜ group passes. Below `PARALLEL_THRESHOLD` it defers
/// to the serial builder, where the rayon fan-out would not be repaid.
pub fn build_sorted_intervals_par<C: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    sort_reverse_direction: Option<&[bool]>,
    slack: T,
    sort_on_ends_too: bool,
) -> Vec<Interval<C, T>> {
    build_sorted_intervals_par_inner(
        chrs,
        starts,
        ends,
        sort_reverse_direction,
        slack,
        sort_on_ends_too,
    )
    .0
}

/// [`build_sorted_intervals_par`] together with the half-open offsets of each
/// group's run in the result.
///
/// Every per-group sweep needs those offsets, and recovering them from the
/// sorted array is a serial pass over `24 * n` bytes — 2.4 GiB at 10^8 rows, on
/// one thread, after all the parallel work is done. The counting partition
/// already knows where each group starts, so on the paths that use it the runs
/// come free; the remaining paths pay a scan of the *group ids* instead of the
/// intervals, which is six times narrower.
pub fn build_sorted_intervals_par_with_runs<
    C: GroupType + Send + Sync,
    T: PositionType + Send + Sync,
>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    sort_reverse_direction: Option<&[bool]>,
    slack: T,
    sort_on_ends_too: bool,
) -> SortedIntervals<C, T> {
    let (intervals, runs) = build_sorted_intervals_par_inner(
        chrs,
        starts,
        ends,
        sort_reverse_direction,
        slack,
        sort_on_ends_too,
    );
    let runs = runs.unwrap_or_else(|| partition::runs_by(&intervals, |i| i.group));
    (intervals, runs)
}

/// The shared body. The group runs are returned whenever the path that ran knew
/// them anyway, so `build_sorted_intervals_par` never pays to compute them.
fn build_sorted_intervals_par_inner<C: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    sort_reverse_direction: Option<&[bool]>,
    slack: T,
    sort_on_ends_too: bool,
) -> (Vec<Interval<C, T>>, Option<GroupRuns>) {
    if chrs.len() < partition::PARALLEL_THRESHOLD {
        return (
            build_sorted_intervals(
                chrs,
                starts,
                ends,
                sort_reverse_direction,
                slack,
                sort_on_ends_too,
            ),
            None,
        );
    }

    // Nothing to sort when the input is already in target order, which
    // coordinate-sorted genomic files usually are. Checked on the inputs, in
    // parallel, and short-circuiting, so it is ~free when they are not. Building
    // only shifts coordinates by `slack`, which preserves order; a
    // reverse-direction build negates them, so that case skips the check.
    let block_list = partition::blocks(chrs.len());
    if sort_reverse_direction.is_none()
        && is_input_group_coord_sorted(chrs, starts, ends, &block_list, sort_on_ends_too)
    {
        // The rows keep their input positions here, so the group ids carry the
        // runs directly and never have to be read back out of the intervals.
        return (
            build_intervals_par(chrs, starts, ends, slack, &block_list),
            Some(partition::runs_by(chrs, |chr| *chr)),
        );
    }

    // Make each group contiguous with the parallel counting partition rather
    // than a serial radix pass over the whole array, then sort each group's
    // coordinates on its own task. Identical output: the partition is stable, so
    // within a group the rows keep input order exactly as the stable radix left
    // them.
    if sort_reverse_direction.is_none() {
        if let Some((layout, block_list)) = partition::plan_group_major(chrs, 1) {
            let runs = layout.group_runs();
            let mut intervals =
                partition::scatter_group_major(chrs, &block_list, layout, 1, |row, slots| {
                    slots[0].write(Interval {
                        group: chrs[row],
                        start: starts[row] - slack,
                        end: ends[row] + slack,
                        idx: row as u32,
                    });
                });
            sort_groups_by_coordinate(&mut intervals, &runs, sort_on_ends_too);
            return (intervals, Some(runs));
        }
    }

    let mut intervals = build_intervals(chrs, starts, ends, sort_reverse_direction, slack);
    sort_by_key(&mut intervals, |i| i.group);
    let runs = partition::runs_by(&intervals, |i| i.group);
    sort_groups_by_coordinate(&mut intervals, &runs, sort_on_ends_too);

    (intervals, Some(runs))
}

/// Largest coordinate, and largest slack, the narrowed builder will accept.
///
/// Halving the `i32` range leaves the sweeps room to add a slack to an end
/// without overflowing: `2^30 + (2^30 - 1)` is exactly `i32::MAX`, and the
/// negative side lands on `i32::MIN`. Genomic coordinates are two orders of
/// magnitude below the bound; protein ones are five.
pub const NARROW_COORDINATE_LIMIT: i64 = 1 << 30;

/// True when a slack is small enough to add to a narrowed coordinate.
pub fn slack_fits_narrow(slack: i64) -> bool {
    slack.unsigned_abs() < NARROW_COORDINATE_LIMIT as u64
}

/// Build the group-major sorted interval array with `i32` coordinates from
/// `i64` inputs, together with each group's run.
///
/// The sort moves a third less data this way — an `Interval<u32, i32>` is 16
/// bytes against 24 — and the narrowing itself is free, because the partition
/// already copies every coordinate once and can convert as it writes. Genomic
/// readers hand out `Int64` columns for coordinates that never leave `i32`, so
/// this is the usual case rather than a special one.
///
/// Returns `None` when the input is too small to be worth partitioning or when
/// a coordinate is outside [`NARROW_COORDINATE_LIMIT`]; the caller then stays on
/// its `i64` path.
///
/// The counting partition declining is *not* one of those cases. It refuses
/// high-cardinality frames on purpose — its histogram is one counter per
/// (block, group), so 500k groups do not fit the budget — but that is a reason
/// to sort differently, not to give up a third of the bytes. This used to
/// return `None` there, which meant merge and cluster silently ran wide on
/// exactly the many-group shapes the narrowing helps most. The whole-array
/// radix fallback below is what [`build_sorted_intervals_par_inner`] and
/// `sorted_records_narrowed_par` already do.
pub fn build_sorted_narrowed_intervals_par<C: GroupType + Send + Sync>(
    chrs: &[C],
    starts: &[i64],
    ends: &[i64],
) -> Option<SortedIntervals<C, i32>> {
    if chrs.len() < partition::PARALLEL_THRESHOLD || !fits_narrow(starts, ends) {
        return None;
    }

    let narrowed = |row: usize| Interval {
        group: chrs[row],
        start: starts[row] as i32,
        end: ends[row] as i32,
        idx: row as u32,
    };

    let (mut intervals, runs) = match partition::plan_group_major(chrs, 1) {
        Some((layout, block_list)) => {
            let runs = layout.group_runs();
            let intervals =
                partition::scatter_group_major(chrs, &block_list, layout, 1, |row, slots| {
                    slots[0].write(narrowed(row));
                });
            (intervals, runs)
        }
        None => {
            let mut intervals: Vec<Interval<C, i32>> = (0..chrs.len()).map(narrowed).collect();
            sort_by_key(&mut intervals, |i| i.group);
            let runs = partition::runs_by(&intervals, |i| i.group);
            (intervals, runs)
        }
    };
    sort_groups_by_coordinate(&mut intervals, &runs, false);

    Some((intervals, runs))
}

/// Largest coordinate a sweep that returns a *distance* will accept.
///
/// [`NARROW_COORDINATE_LIMIT`] bounds each endpoint, which is all a sweep
/// needs when it only ever compares two coordinates. `nearest` subtracts them:
/// its answer is `target - query + 1`, so what has to fit `i32` is the whole
/// span a group can hold, and two coordinates at opposite ends of the wider
/// bound are `2^31` apart — one past `i32::MAX`. Halving the bound halves the
/// span. Human chr1, the longest, is 2.49e8 against this 5.37e8.
pub const NARROW_DISTANCE_LIMIT: i64 = 1 << 29;

/// True when every coordinate is within [`NARROW_COORDINATE_LIMIT`].
pub(crate) fn fits_narrow(starts: &[i64], ends: &[i64]) -> bool {
    fits_narrow_within(starts, ends, NARROW_COORDINATE_LIMIT)
}

/// True when every coordinate is within `limit` of zero.
pub(crate) fn fits_narrow_within(starts: &[i64], ends: &[i64], limit: i64) -> bool {
    let in_range = |value: &i64| (-limit..=limit).contains(value);
    partition::blocks(starts.len())
        .par_iter()
        .all(|&(from, to)| starts[from..to].iter().all(in_range) && ends[from..to].iter().all(in_range))
}

/// Sort each group's slice by start (and end first, when asked), one task per
/// batch of groups.
fn sort_groups_by_coordinate<C: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    intervals: &mut [Interval<C, T>],
    runs: &[(usize, usize)],
    sort_on_ends_too: bool,
) {
    partition::par_groups_mut(intervals, runs, |group| {
        if sort_on_ends_too {
            partition::sort_run_by_key(&mut *group, |i| i.end);
        }
        partition::sort_run_by_key(group, |i| i.start);
    });
}

/// Build the interval array in input order, in parallel.
fn build_intervals_par<C: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    slack: T,
    block_list: &[(usize, usize)],
) -> Vec<Interval<C, T>> {
    let mut intervals: Vec<Interval<C, T>> = Vec::with_capacity(chrs.len());
    // Split by `block_list` rather than by a fixed chunk width: `partition::blocks`
    // divides the rows evenly among the blocks, so a block is not `BLOCK_ROWS`
    // wide and zipping against `par_chunks_mut(BLOCK_ROWS)` would pair each block
    // with the wrong destination.
    partition::split_runs_mut(intervals.spare_capacity_mut(), block_list)
        .into_par_iter()
        .zip(block_list.par_iter())
        .for_each(|(target, &(start, end))| {
            for (slot, row) in target.iter_mut().zip(start..end) {
                slot.write(Interval {
                    group: chrs[row],
                    start: starts[row] - slack,
                    end: ends[row] + slack,
                    idx: row as u32,
                });
            }
        });
    // SAFETY: every slot in 0..len was written above.
    unsafe { intervals.set_len(chrs.len()) };
    intervals
}

/// True when `(group, start[, end])` is already non-decreasing across the inputs.
fn is_input_group_coord_sorted<C: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    block_list: &[(usize, usize)],
    with_ends: bool,
) -> bool {
    block_list.par_iter().all(|&(start, end)| {
        let from = start.saturating_sub(1);
        for row in from..end.saturating_sub(1) {
            let next = row + 1;
            if (chrs[row], starts[row]) > (chrs[next], starts[next]) {
                return false;
            }
            if with_ends
                && chrs[row] == chrs[next]
                && starts[row] == starts[next]
                && ends[row] > ends[next]
            {
                return false;
            }
        }
        true
    })
}

pub fn build_sorted_subsequence_intervals<G: GroupType, T: PositionType>(
    chrs: &[G],
    starts: &[T],
    ends: &[T],
    strand_flags: &[bool],
) -> Vec<SplicedSubsequenceInterval<G, T>> {
    let mut intervals = build_subsequence_intervals(chrs, starts, ends, strand_flags);

    sort_by_key(&mut intervals, |i| i.end);
    sort_by_key(&mut intervals, |i| i.start);
    sort_by_key(&mut intervals, |i| i.chr);

    intervals
}

pub fn build_sorted_sequence_intervals(
    chrs: &[i64],
    starts: &[i64],
    ends: &[i64],
    idxs: &[i64],
    strand_flags: &[bool],
    force_plus_strand: bool,
) -> Vec<SubsequenceInterval> {
    let mut intervals =
        build_sequence_intervals(chrs, starts, ends, idxs, strand_flags, force_plus_strand);

    sort_by_key(&mut intervals, |i| i.end);
    sort_by_key(&mut intervals, |i| i.start);
    sort_by_key(&mut intervals, |i| i.group_id);

    intervals
}

pub fn sort_order_idx<G: GroupType, T: PositionType>(
    chrs: &[G],
    starts: &[T],
    ends: &[T],
    sort_reverse_direction: Option<&[bool]>,
) -> Vec<u32> {
    build_sorted_intervals(chrs, starts, ends, sort_reverse_direction, T::zero(), true)
        .iter()
        .map(|i| i.idx)
        .collect()
}

pub fn build_sorted_events_single_position<C: GroupType, T: PositionType>(
    chrs: &[C],
    pos: &[T],
    start: bool,
    first_set: bool,
    negative_position: bool,
    slack: T,
) -> Vec<Event<C, T>> {
    let mut events = Vec::with_capacity(2 * (chrs.len()));

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        let pos = if start {
            pos[i] - slack
        } else {
            pos[i] + slack
        };
        events.push(Event {
            chr: chrs[i],
            pos: if negative_position { -pos } else { pos },
            is_start: start,
            first_set: first_set,
            idx: i as u32,
        });
    }

    sort_by_key(&mut events, |e| (e.chr, e.pos, e.is_start));

    events
}

pub fn build_sorted_events_single_collection<C: GroupType, T: PositionType>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    slack: T,
) -> Vec<Event<C, T>> {
    let mut events = Vec::with_capacity(2 * (chrs.len()));

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        events.push(Event {
            chr: chrs[i],
            pos: starts[i],
            is_start: true,
            first_set: true,
            idx: i as u32,
        });
        events.push(Event {
            chr: chrs[i],
            pos: ends[i] + slack,
            is_start: false,
            first_set: true,
            idx: i as u32,
        });
    }

    // Sort events by:
    // 1. pos (ascending)
    // 2. is_start before is_end (if pos ties)
    // (We don't strictly need to tie-break by set_id or idx, but we can.)

    sort_by_key(&mut events, |e| e.is_start);
    sort_by_key(&mut events, |e| e.pos);
    sort_by_key(&mut events, |e| e.chr);

    events
}

/// Parallel counterpart of [`build_sorted_events_single_collection`] with
/// **byte-identical** output, used by the single-input sweep kernels (cluster,
/// split, boundary, complement).
///
/// The chr-major counting partition makes each chromosome a contiguous run, then
/// each run's events are sorted on their own rayon task — serial radix within a
/// chromosome, parallel across chromosomes. radsort is a stable LSD radix, so
/// the result matches the serial `is_start → pos → chr` order. Below
/// `PARALLEL_THRESHOLD` it defers to the serial builder. There is no presorted
/// fast path: the event array interleaves start and end events, so it is
/// essentially never already pos-sorted, even when the input intervals are.
pub fn build_sorted_events_single_collection_par<
    C: GroupType + Send + Sync,
    T: PositionType + Send + Sync,
>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    slack: T,
) -> Vec<Event<C, T>> {
    if 2 * chrs.len() < partition::PARALLEL_THRESHOLD {
        return build_sorted_events_single_collection(chrs, starts, ends, slack);
    }

    // Both the counting partition and the radix pass it replaces are stable, and
    // a row's two events stay adjacent and in build order within their
    // chromosome, so the sorted result is identical either way.
    let (mut events, runs) = match partition::plan_group_major(chrs, 2) {
        Some((layout, block_list)) => {
            let runs = layout.group_runs();
            let events =
                partition::scatter_group_major(chrs, &block_list, layout, 2, |row, slots| {
                    slots[0].write(Event {
                        chr: chrs[row],
                        pos: starts[row],
                        is_start: true,
                        first_set: true,
                        idx: row as u32,
                    });
                    slots[1].write(Event {
                        chr: chrs[row],
                        pos: ends[row] + slack,
                        is_start: false,
                        first_set: true,
                        idx: row as u32,
                    });
                });
            (events, runs)
        }
        None => {
            let mut events = Vec::with_capacity(2 * chrs.len());
            for i in 0..chrs.len() {
                events.push(Event {
                    chr: chrs[i],
                    pos: starts[i],
                    is_start: true,
                    first_set: true,
                    idx: i as u32,
                });
                events.push(Event {
                    chr: chrs[i],
                    pos: ends[i] + slack,
                    is_start: false,
                    first_set: true,
                    idx: i as u32,
                });
            }
            sort_by_key(&mut events, |e| e.chr);
            let runs = partition::runs_by(&events, |e| e.chr);
            (events, runs)
        }
    };

    partition::par_groups_mut(&mut events, &runs, |group| {
        partition::sort_run_by_key(&mut *group, |e| e.is_start);
        partition::sort_run_by_key(group, |e| e.pos);
    });

    events
}

/// One event per row, ordered by `(group, pos - slack)`.
///
/// The counting partition makes the groups contiguous on every core and each
/// group's events are then sorted on their own task, which is the same output
/// the two whole-array radix passes below produce: the partition is stable and
/// radsort is stable, so equal `(group, pos)` events keep row order either way.
pub fn build_sorted_events_single_collection_separate_outputs<
    C: GroupType + Send + Sync,
    T: PositionType + Send + Sync,
>(
    chrs: &[C],
    pos: &[T],
    slack: T,
) -> Vec<MinEvent<C, T>> {
    if chrs.len() >= partition::PARALLEL_THRESHOLD {
        if let Some((layout, block_list)) = partition::plan_group_major(chrs, 1) {
            let runs = layout.group_runs();
            let mut events =
                partition::scatter_group_major(chrs, &block_list, layout, 1, |row, slots| {
                    slots[0].write(MinEvent {
                        chr: chrs[row],
                        pos: pos[row] - slack,
                        idx: row as u32,
                    });
                });
            partition::par_groups_mut(&mut events, &runs, |group| {
                partition::sort_run_by_key(group, |e| e.pos)
            });
            return events;
        }
    }

    let mut out_pos: Vec<MinEvent<C, T>> = Vec::with_capacity(chrs.len());

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        out_pos.push(MinEvent {
            chr: chrs[i],
            pos: pos[i] - slack,
            idx: i as u32,
        });
    }

    sort_by_key(&mut out_pos, |e| e.pos);
    sort_by_key(&mut out_pos, |e| e.chr);

    out_pos
}

/// [`build_sorted_events_single_collection_separate_outputs`] over `i64`
/// positions narrow enough to sweep as `i32`.
///
/// Same shape, same order, one difference: the narrowing happens as the events
/// are written rather than in a pass of its own. `nearest` builds four of these
/// views, so doing it here is the difference between four extra passes over the
/// coordinates and none.
///
/// `None` when the counting partition declines or the input is too small to be
/// worth partitioning — the caller then has to narrow up front anyway, so there
/// is nothing to gain from a serial narrowed path here.
pub fn build_sorted_events_narrowed_separate_outputs<C: GroupType + Send + Sync>(
    chrs: &[C],
    pos: &[i64],
    slack: i64,
) -> Option<Vec<MinEvent<C, i32>>> {
    debug_assert_eq!(chrs.len(), pos.len(), "chrs/pos length mismatch");

    if chrs.len() < partition::PARALLEL_THRESHOLD {
        return None;
    }
    let (layout, block_list) = partition::plan_group_major(chrs, 1)?;

    let runs = layout.group_runs();
    let mut events = partition::scatter_group_major(chrs, &block_list, layout, 1, |row, slots| {
        slots[0].write(MinEvent {
            chr: chrs[row],
            // The caller has checked that both the coordinates and the slack
            // fit, so the subtraction is done in `i64` and narrowed after.
            pos: (pos[row] - slack) as i32,
            idx: row as u32,
        });
    });
    partition::par_groups_mut(&mut events, &runs, |group| {
        partition::sort_run_by_key(group, |e| e.pos)
    });
    Some(events)
}

pub fn build_sorted_groups<C: GroupType>(chrs: &[C]) -> Vec<u32> {
    let mut out: Vec<GroupStruct<C>> = (0..chrs.len())
        .map(|i| GroupStruct {
            chr: chrs[i],
            idx: i as u32,
        })
        .collect();

    sort_by_key(&mut out, |e| e.chr);

    // take the chromosome field, cast to u32, collect -----------------------
    out.into_iter().map(|e| e.idx).collect()
}

pub fn build_sorted_events_with_starts_ends<C: GroupType, T: PositionType>(
    chrs: &[C],
    pos: &[T],
    slack: T,
) -> Vec<MinEvent<C, T>> {
    let mut out_pos = Vec::with_capacity(chrs.len());

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        out_pos.push(MinEvent {
            chr: chrs[i],
            pos: pos[i] - slack,
            idx: i as u32,
        });
    }

    sort_by_key(&mut out_pos, |e| e.pos);
    sort_by_key(&mut out_pos, |e| e.chr);

    out_pos
}

pub fn build_sorted_events<C: GroupType, T: PositionType>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
) -> Vec<GenericEvent<C, T>> {
    let mut events = Vec::with_capacity(2 * (chrs.len() + chrs2.len()));

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        events.push(GenericEvent {
            chr: chrs[i],
            pos: if slack < starts[i] {
                starts[i] - slack
            } else {
                T::zero()
            },
            is_start: true,
            first_set: true,
            idx: i as u32,
        });
        events.push(GenericEvent {
            chr: chrs[i],
            pos: ends[i].saturating_add(slack),
            is_start: false,
            first_set: true,
            idx: i as u32,
        });
    }

    for j in 0..chrs2.len() {
        events.push(GenericEvent {
            chr: chrs2[j],
            pos: starts2[j],
            is_start: true,
            first_set: false,
            idx: j as u32,
        });
        events.push(GenericEvent {
            chr: chrs2[j],
            pos: ends2[j],
            is_start: false,
            first_set: false,
            idx: j as u32,
        });
    }

    sort_by_key(&mut events, |e| e.is_start);
    sort_by_key(&mut events, |e| e.pos);
    sort_by_key(&mut events, |e| e.chr);

    events
}

pub fn build_sorted_maxevents_with_starts_ends<C: GroupType, T: PositionType>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
) -> Vec<MaxEvent<C, T>> {
    let mut events = Vec::with_capacity(2 * (chrs.len() + chrs2.len()));

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        events.push(MaxEvent {
            chr: chrs[i],
            pos: starts[i] - slack,
            start: starts[i] - slack,
            end: ends[i] + slack,
            is_start: true,
            first_set: true,
            idx: i as u32,
        });
        events.push(MaxEvent {
            chr: chrs[i],
            pos: ends[i] + slack,
            end: ends[i] + slack,
            start: starts[i] - slack,
            is_start: false,
            first_set: true,
            idx: i as u32,
        });
    }

    for i in 0..chrs2.len() {
        events.push(MaxEvent {
            chr: chrs2[i],
            pos: starts2[i],
            start: starts2[i],
            end: ends2[i],
            is_start: true,
            first_set: false,
            idx: i as u32,
        });
        events.push(MaxEvent {
            chr: chrs2[i],
            pos: ends2[i],
            start: starts2[i],
            end: ends2[i],
            is_start: false,
            first_set: false,
            idx: i as u32,
        });
    }

    sort_by_key(&mut events, |e| e.is_start);
    sort_by_key(&mut events, |e| e.pos);
    sort_by_key(&mut events, |e| e.chr);

    events
}

pub fn build_sorted_events_idxs<C: GroupType, T: PositionType>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
) -> Vec<Event<C, T>> {
    let mut events = Vec::with_capacity(2 * (chrs.len() + chrs2.len()));

    // Convert set1 intervals into events
    for i in 0..chrs.len() {
        events.push(Event {
            chr: chrs[i],
            pos: starts[i] - slack,
            is_start: true,
            first_set: true,
            idx: i as u32,
        });
        events.push(Event {
            chr: chrs[i],
            pos: ends[i] + slack,
            is_start: false,
            first_set: true,
            idx: i as u32,
        });
    }

    for j in 0..chrs2.len() {
        events.push(Event {
            chr: chrs2[j],
            pos: starts2[j],
            is_start: true,
            first_set: false,
            idx: j as u32,
        });
        events.push(Event {
            chr: chrs2[j],
            pos: ends2[j],
            is_start: false,
            first_set: false,
            idx: j as u32,
        });
    }

    sort_by_key(&mut events, |e| e.is_start);
    sort_by_key(&mut events, |e| e.pos);
    sort_by_key(&mut events, |e| e.chr);

    events
}

pub fn build_sorted_events_from_intervals<C: GroupType, T: PositionType>(
    intervals1: &mut [Interval<C, T>],
    intervals2: &mut [Interval<C, T>],
) -> Vec<Event<C, T>> {
    let mut events = Vec::with_capacity(2 * (intervals1.len() + intervals2.len()));

    // Convert set1 intervals into events
    for interval in intervals1 {
        events.push(Event {
            chr: interval.group,
            pos: interval.start,
            is_start: true,
            first_set: true,
            idx: interval.idx,
        });
        events.push(Event {
            chr: interval.group,
            pos: interval.end,
            is_start: false,
            first_set: true,
            idx: interval.idx,
        });
    }

    for interval in intervals2 {
        events.push(Event {
            chr: interval.group,
            pos: interval.start,
            is_start: true,
            first_set: false,
            idx: interval.idx,
        });
        events.push(Event {
            chr: interval.group,
            pos: interval.end,
            is_start: false,
            first_set: false,
            idx: interval.idx,
        });
    }

    // Sort events by:
    // 1. pos (ascending)
    // 2. is_start before is_end (if pos ties)
    // (We don't strictly need to tie-break by set_id or idx, but we can.)
    sort_by_key(&mut events, |e| e.is_start);
    sort_by_key(&mut events, |e| e.pos);

    events
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic SplitMix64, so the parity tests need no dev-dependency.
    fn next(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// `n` random intervals over `groups` chromosomes. `sorted` produces the
    /// coordinate-sorted input the parallel builder has a fast path for.
    fn random_input(
        state: &mut u64,
        n: usize,
        groups: u32,
        sorted: bool,
    ) -> (Vec<u32>, Vec<i64>, Vec<i64>) {
        let mut chrs = Vec::with_capacity(n);
        let mut starts = Vec::with_capacity(n);
        let mut ends = Vec::with_capacity(n);
        for _ in 0..n {
            let start = (next(state) % 1_000_000) as i64;
            chrs.push((next(state) % groups as u64) as u32);
            starts.push(start);
            ends.push(start + 1 + (next(state) % 500) as i64);
        }
        if sorted {
            let mut rows: Vec<(u32, i64, i64)> = chrs
                .iter()
                .zip(&starts)
                .zip(&ends)
                .map(|((&c, &s), &e)| (c, s, e))
                .collect();
            rows.sort();
            for (i, (c, s, e)) in rows.into_iter().enumerate() {
                chrs[i] = c;
                starts[i] = s;
                ends[i] = e;
            }
        }
        (chrs, starts, ends)
    }

    fn interval_key(iv: &Interval<u32, i64>) -> (u32, i64, i64, u32) {
        (iv.group, iv.start, iv.end, iv.idx)
    }

    fn event_key(e: &Event<u32, i64>) -> (u32, i64, bool, bool, u32) {
        (e.chr, e.pos, e.is_start, e.first_set, e.idx)
    }

    #[test]
    fn parallel_interval_builder_matches_serial() {
        let mut state = 0x1234_5678;
        // Sizes either side of PARALLEL_THRESHOLD, and both input orders, so the
        // presorted fast path and the counting partition are both exercised.
        for &n in &[1_usize, 999, 150_000] {
            for sorted in [false, true] {
                for &groups in &[1_u32, 24] {
                    for ends_too in [false, true] {
                        let (chrs, starts, ends) = random_input(&mut state, n, groups, sorted);
                        let serial =
                            build_sorted_intervals(&chrs, &starts, &ends, None, 0, ends_too);
                        let parallel =
                            build_sorted_intervals_par(&chrs, &starts, &ends, None, 0, ends_too);
                        assert_eq!(
                            serial.iter().map(interval_key).collect::<Vec<_>>(),
                            parallel.iter().map(interval_key).collect::<Vec<_>>(),
                            "n={n} sorted={sorted} groups={groups} ends_too={ends_too}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn parallel_interval_builder_matches_serial_with_slack_and_reverse() {
        let mut state = 0x9999;
        let (chrs, starts, ends) = random_input(&mut state, 150_000, 8, false);
        let reverse: Vec<bool> = (0..chrs.len()).map(|i| i % 3 == 0).collect();

        let serial = build_sorted_intervals(&chrs, &starts, &ends, None, 7, false);
        let parallel = build_sorted_intervals_par(&chrs, &starts, &ends, None, 7, false);
        assert_eq!(
            serial.iter().map(interval_key).collect::<Vec<_>>(),
            parallel.iter().map(interval_key).collect::<Vec<_>>()
        );

        let serial = build_sorted_intervals(&chrs, &starts, &ends, Some(&reverse), 0, true);
        let parallel = build_sorted_intervals_par(&chrs, &starts, &ends, Some(&reverse), 0, true);
        assert_eq!(
            serial.iter().map(interval_key).collect::<Vec<_>>(),
            parallel.iter().map(interval_key).collect::<Vec<_>>()
        );
    }

    #[test]
    fn parallel_event_builder_matches_serial() {
        let mut state = 0xABCD;
        for &n in &[0_usize, 1, 999, 150_000] {
            for &groups in &[1_u32, 24] {
                let (chrs, starts, ends) = random_input(&mut state, n, groups, false);
                let serial = build_sorted_events_single_collection(&chrs, &starts, &ends, 3);
                let parallel = build_sorted_events_single_collection_par(&chrs, &starts, &ends, 3);
                assert_eq!(
                    serial.iter().map(event_key).collect::<Vec<_>>(),
                    parallel.iter().map(event_key).collect::<Vec<_>>(),
                    "n={n} groups={groups}"
                );
            }
        }
    }

    /// Group ids the counting partition declines to index, so the parallel
    /// builders fall back to the radix path.
    #[test]
    fn parallel_builders_match_serial_on_a_span_the_partition_declines() {
        let mut state = 0x5A5A;
        let (mut chrs, starts, ends) = random_input(&mut state, 150_000, 4, false);
        chrs[0] = u32::MAX; // a span far wider than MAX_GROUP_SPAN

        let serial = build_sorted_intervals(&chrs, &starts, &ends, None, 0, true);
        let parallel = build_sorted_intervals_par(&chrs, &starts, &ends, None, 0, true);
        assert_eq!(
            serial.iter().map(interval_key).collect::<Vec<_>>(),
            parallel.iter().map(interval_key).collect::<Vec<_>>()
        );

        let serial = build_sorted_events_single_collection(&chrs, &starts, &ends, 0);
        let parallel = build_sorted_events_single_collection_par(&chrs, &starts, &ends, 0);
        assert_eq!(
            serial.iter().map(event_key).collect::<Vec<_>>(),
            parallel.iter().map(event_key).collect::<Vec<_>>()
        );
    }
}
