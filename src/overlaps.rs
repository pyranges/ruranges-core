use std::str::FromStr;

use radsort::sort_by_key;
use rayon::prelude::*;

use crate::partition;
use crate::ruranges_structs::{GroupType, OverlapPair, OverlapType, PositionType};
use crate::sorts;

#[derive(Copy, Clone, Debug)]
pub(crate) struct IntervalRecord<C: GroupType, T: PositionType> {
    pub(crate) group: C,
    pub(crate) start: T,
    pub(crate) end: T,
    pub(crate) idx: u32,
}

/// The serial reference for the two sort keys below. Every kernel takes
/// [`sorted_records_par`] now; this is what that builder is checked against.
#[cfg(test)]
pub(crate) fn sorted_records<C: GroupType, T: PositionType>(
    groups: &[C],
    starts: &[T],
    ends: &[T],
    overlap_type: OverlapType,
) -> Vec<IntervalRecord<C, T>> {
    debug_assert_eq!(groups.len(), starts.len(), "groups/starts length mismatch");
    debug_assert_eq!(groups.len(), ends.len(), "groups/ends length mismatch");

    let mut records = Vec::with_capacity(groups.len());

    for idx in 0..groups.len() {
        records.push(IntervalRecord {
            group: groups[idx],
            start: starts[idx],
            end: ends[idx],
            idx: idx as u32,
        });
    }

    // Sort key choice depends on overlap_type:
    //
    // - `All`: every overlapping pair is emitted, so sweep order does not affect
    //   *which* pairs come out; downstream `overlaps()` then sorts pairs by
    //   `(idx, idx2)` when `sort_output=true`, fully determining the final
    //   ordering. Sorting by `(group, start)` is therefore sufficient and saves
    //   two radsort passes over fat 32-byte records.
    //
    // - `First` / `Last`: the sweep walks the active list in order and emits
    //   only the first/last match, so the in-sweep ordering of records with the
    //   same `(group, start)` decides which target is returned. We need a
    //   tie-breaker that's stable across runs/platforms — including `end` and
    //   `idx` ensures determinism. This was load-bearing for
    //   pyranges/ruranges_py#23 (sort_order=True returning right-side indices
    //   in scrambled order); do not drop those fields here.
    match overlap_type {
        OverlapType::All => sort_by_key(&mut records, |r| (r.group, r.start)),
        OverlapType::First | OverlapType::Last => {
            sort_by_key(&mut records, |r| (r.group, r.start, r.end, r.idx))
        }
    }

    records
}

/// Records made group-contiguous, with the offsets of each group's run.
///
/// The sweeps below need both: the records to walk, and where one chromosome
/// ends and the next begins, so that each group can be handed to its own thread.
pub(crate) struct SortedRecords<C: GroupType, T: PositionType> {
    pub(crate) records: Vec<IntervalRecord<C, T>>,
    /// Half-open offsets of each group's run, ascending by group id.
    pub(crate) runs: Vec<(usize, usize)>,
}

/// Group-major sorted records, in parallel, with the group runs the sweeps need.
///
/// `wide_key` selects the same two keys [`sorted_records`] documents: the narrow
/// `(group, start)`, or `(group, start, end, idx)` when the caller's output
/// depends on which of several equally-positioned targets the sweep reaches
/// first. Output is identical to `sorted_records` for the matching key, because
/// the counting partition is stable and radsort is stable: rows keep `idx` order
/// inside a group exactly as the whole-array radix left them.
pub(crate) fn sorted_records_par<C: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    groups: &[C],
    starts: &[T],
    ends: &[T],
    wide_key: bool,
) -> SortedRecords<C, T> {
    debug_assert_eq!(groups.len(), starts.len(), "groups/starts length mismatch");
    debug_assert_eq!(groups.len(), ends.len(), "groups/ends length mismatch");

    let (mut records, runs) = if groups.len() >= partition::PARALLEL_THRESHOLD {
        match partition::plan_group_major(groups, 1) {
            Some((layout, block_list)) => {
                let runs = layout.group_runs();
                let records =
                    partition::scatter_group_major(groups, &block_list, layout, 1, |row, slots| {
                        slots[0].write(IntervalRecord {
                            group: groups[row],
                            start: starts[row],
                            end: ends[row],
                            idx: row as u32,
                        });
                    });
                (records, runs)
            }
            None => {
                let mut records = build_records(groups, starts, ends);
                sort_by_key(&mut records, |r| r.group);
                let runs = partition::runs_by(&records, |r| r.group);
                (records, runs)
            }
        }
    } else {
        let mut records = build_records(groups, starts, ends);
        sort_by_key(&mut records, |r| r.group);
        let runs = partition::runs_by(&records, |r| r.group);
        (records, runs)
    };

    // Every record in a run shares its group, so the group field drops out of
    // the sort key here — and `idx` drops out of the wide key too, because the
    // partition is stable and the records were built in `idx` order, so a run is
    // already `idx`-ascending before any of these passes.
    partition::par_groups_mut(&mut records, &runs, |group| {
        if wide_key {
            partition::sort_run_by_key(&mut *group, |r| r.end);
        }
        partition::sort_run_by_key(group, |r| r.start);
    });

    SortedRecords { records, runs }
}

/// [`sorted_records_par`] over `i64` coordinates narrow enough to sweep as
/// `i32`, which sorts a third less data for the same answer.
///
/// An `IntervalRecord<u32, i32>` is 16 bytes against 24, and the narrowing is
/// free: the counting partition already copies every coordinate once and
/// converts as it writes. Genomic readers hand out `Int64` columns for
/// coordinates that never leave `i32`, so this is the usual case rather than a
/// special one.
///
/// `None` when the input is too small to be worth partitioning, when a
/// coordinate is outside [`sorts::NARROW_COORDINATE_LIMIT`], or when the
/// counting partition declines; the caller then stays on its `i64` path.
pub(crate) fn sorted_records_narrowed_par<C: GroupType + Send + Sync>(
    groups: &[C],
    starts: &[i64],
    ends: &[i64],
    wide_key: bool,
) -> Option<SortedRecords<C, i32>> {
    debug_assert_eq!(groups.len(), starts.len(), "groups/starts length mismatch");
    debug_assert_eq!(groups.len(), ends.len(), "groups/ends length mismatch");

    if groups.len() < partition::PARALLEL_THRESHOLD || !sorts::fits_narrow(starts, ends) {
        return None;
    }

    let narrowed = |row: usize| IntervalRecord {
        group: groups[row],
        start: starts[row] as i32,
        end: ends[row] as i32,
        idx: row as u32,
    };

    // The same two paths, and the same fallback, as `sorted_records_par`: the
    // counting partition when it pays, the whole-array radix when it does not.
    // Narrowing is worth having either way — it is the *sort* that moves a third
    // less data, not the partition.
    let (mut records, runs) = match partition::plan_group_major(groups, 1) {
        Some((layout, block_list)) => {
            let runs = layout.group_runs();
            let records =
                partition::scatter_group_major(groups, &block_list, layout, 1, |row, slots| {
                    slots[0].write(narrowed(row));
                });
            (records, runs)
        }
        None => {
            let mut records: Vec<IntervalRecord<C, i32>> =
                (0..groups.len()).map(narrowed).collect();
            sort_by_key(&mut records, |r| r.group);
            let runs = partition::runs_by(&records, |r| r.group);
            (records, runs)
        }
    };

    // Same two keys, and the same reasoning, as `sorted_records_par`.
    partition::par_groups_mut(&mut records, &runs, |group| {
        if wide_key {
            partition::sort_run_by_key(&mut *group, |r| r.end);
        }
        partition::sort_run_by_key(group, |r| r.start);
    });

    Some(SortedRecords { records, runs })
}

/// True when both sides of a binary kernel, and its slack, fit the narrowed
/// sweep. Checked before either side is built, so a wide right-hand side does
/// not cost a narrowed sort of the left.
pub(crate) fn pair_fits_narrow(
    starts: &[i64],
    ends: &[i64],
    starts2: &[i64],
    ends2: &[i64],
    slack: i64,
) -> bool {
    sorts::slack_fits_narrow(slack)
        && sorts::fits_narrow(starts, ends)
        && sorts::fits_narrow(starts2, ends2)
}

fn build_records<C: GroupType, T: PositionType>(
    groups: &[C],
    starts: &[T],
    ends: &[T],
) -> Vec<IntervalRecord<C, T>> {
    (0..groups.len())
        .map(|idx| IntervalRecord {
            group: groups[idx],
            start: starts[idx],
            end: ends[idx],
            idx: idx as u32,
        })
        .collect()
}

/// Half-open offsets of one group's run.
pub(crate) type GroupRun = (usize, usize);

/// One group's runs on each side, for a group both sides have.
pub(crate) type PairedGroupRun = (GroupRun, GroupRun);

/// The group runs the two sides have in common, ascending by group id.
///
/// The sweeps only ever visit groups present on both sides, so this merge join
/// over the two run lists is the parallel equivalent of their two-pointer group
/// walk — and it is what makes the per-group work independent.
pub(crate) fn paired_group_runs<C: GroupType, T: PositionType>(
    left: &SortedRecords<C, T>,
    right: &SortedRecords<C, T>,
) -> Vec<PairedGroupRun> {
    let mut out = Vec::new();
    let (mut i, mut j) = (0usize, 0usize);
    while i < left.runs.len() && j < right.runs.len() {
        let left_group = left.records[left.runs[i].0].group;
        let right_group = right.records[right.runs[j].0].group;
        if left_group < right_group {
            i += 1;
        } else if right_group < left_group {
            j += 1;
        } else {
            out.push((left.runs[i], right.runs[j]));
            i += 1;
            j += 1;
        }
    }
    out
}

/// Concatenate per-group results in group order, which is the order the serial
/// two-pointer sweep emitted them in.
fn concat_groups<E: Send>(batches: Vec<Vec<E>>) -> Vec<E> {
    let total = batches.iter().map(Vec::len).sum();
    let mut out = Vec::with_capacity(total);
    for batch in batches {
        out.extend(batch);
    }
    out
}

/// Run `body` on each paired group, batching consecutive pairs into rayon tasks.
/// The batches come back in group order, which is the order the serial group
/// walk emits them in.
pub(crate) fn par_paired_group_batches<E, F>(pairs: &[PairedGroupRun], body: F) -> Vec<Vec<E>>
where
    E: Send,
    F: Fn(GroupRun, GroupRun, &mut Vec<E>) + Send + Sync,
{
    let batches = partition::batch_runs(pairs, partition::BATCH_ROWS, |&(l, r)| {
        (l.1 - l.0) + (r.1 - r.0)
    });
    batches
        .into_par_iter()
        .map(|batch| {
            let mut out = Vec::new();
            for &(left_run, right_run) in &pairs[batch] {
                body(left_run, right_run, &mut out);
            }
            out
        })
        .collect()
}

/// [`par_paired_group_batches`] for a body that writes its own output rather
/// than returning one — the counting kernels scatter straight into a result
/// vector, so there is nothing to collect and concatenate.
pub(crate) fn par_paired_groups_for_each<F>(pairs: &[PairedGroupRun], body: F)
where
    F: Fn(GroupRun, GroupRun) + Send + Sync,
{
    let batches = partition::batch_runs(pairs, partition::BATCH_ROWS, |&(l, r)| {
        (l.1 - l.0) + (r.1 - r.0)
    });
    batches.into_par_iter().for_each(|batch| {
        for &(left_run, right_run) in &pairs[batch] {
            body(left_run, right_run);
        }
    });
}

/// [`par_paired_group_batches`], concatenated.
pub(crate) fn par_paired_groups<E, F>(pairs: &[PairedGroupRun], body: F) -> Vec<E>
where
    E: Send,
    F: Fn(GroupRun, GroupRun, &mut Vec<E>) + Send + Sync,
{
    concat_groups(par_paired_group_batches(pairs, body))
}

#[inline(always)]
pub(crate) fn overlaps_with_slack<T: PositionType>(
    query_start: T,
    query_end: T,
    target_start: T,
    target_end: T,
    slack: T,
) -> bool {
    query_start < target_end.saturating_add(slack) && target_start < query_end.saturating_add(slack)
}

#[inline(always)]
pub(crate) fn query_contained_in_target_with_slack<T: PositionType>(
    query_start: T,
    query_end: T,
    target_start: T,
    target_end: T,
    slack: T,
) -> bool {
    let query_start_slack = query_start.saturating_sub(slack);
    let query_end_slack = query_end.saturating_add(slack);
    query_start_slack >= target_start && query_end_slack <= target_end
}

pub(crate) fn clear_active(active: &mut Vec<usize>, active_head: &mut usize) {
    active.clear();
    *active_head = 0;
}

fn collect_overlap_pairs<C: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
    overlap_type: OverlapType,
    contained: bool,
) -> Vec<OverlapPair> {
    let wide_key = wide_sort_key(overlap_type);
    let left = sorted_records_par(chrs, starts, ends, wide_key);
    let right = sorted_records_par(chrs2, starts2, ends2, wide_key);
    collect_overlap_pairs_par(&left, &right, slack, overlap_type, contained)
}

// ── closed-form counting path ───────────────────────────────────────────────
//
// The sweep counts each query by scanning its live active window, so its cost
// grows with the number of candidate pairs rather than with the input size —
// pathological on deep-coverage data. Counting instead of enumerating is the
// Binary Interval Search (BITS) method:
//
//     Layer RM, Skadron K, Robins G, Hall IM, Quinlan AR. "Binary Interval
//     Search: a scalable algorithm for counting interval intersections."
//     Bioinformatics 29(1):1-7, 2013. doi:10.1093/bioinformatics/bts652
//
// BITS counts the targets that *cannot* overlap — those ending before the query
// begins and those starting after it ends — and subtracts them from the total.
// Written as a difference of two prefix counts, and carrying `slack`:
//
//     count = #{t.start < q.end + slack} - #{t.end + slack <= q.start}
//
// Both terms are `partition_point` searches over sorted columns, so a query
// costs O(log M) no matter how many targets it actually overlaps, and every
// query is independent of every other.
//
// The identity is exact when `slack >= 0` and every target has positive length.
// Otherwise the two excluded sets can intersect — a target with `end <= start`
// can satisfy both — and the difference under-counts, so those groups fall back
// to the sweep.

/// True if any interval in the group has non-positive length (`end <= start`),
/// the case that breaks the counting identity. O(M) scan, cheap beside the end
/// sort the closed form needs anyway.
pub(crate) fn group_has_nonpositive_length<C: GroupType, T: PositionType>(
    group: &[IntervalRecord<C, T>],
) -> bool {
    group.iter().any(|r| r.end <= r.start)
}

/// `right_group`'s ends, shifted by `slack` and sorted — the second prefix count
/// the identity needs. `right_group` is ascending by start, so the first term
/// reads it directly.
pub(crate) fn sorted_end_bounds<C: GroupType, T: PositionType>(
    right_group: &[IntervalRecord<C, T>],
    slack: T,
) -> Vec<T> {
    let mut end_bounds: Vec<T> = right_group
        .iter()
        .map(|r| r.end.saturating_add(slack))
        .collect();
    radsort::sort(&mut end_bounds);
    end_bounds
}

/// How many targets overlap `query`, by the BITS identity.
#[inline]
pub(crate) fn count_for_query<C: GroupType, T: PositionType>(
    query: &IntervalRecord<C, T>,
    right_group: &[IntervalRecord<C, T>],
    end_bounds: &[T],
    slack: T,
) -> usize {
    let query_end_slack = query.end.saturating_add(slack);
    let started = right_group.partition_point(|r| r.start < query_end_slack);
    let ended = end_bounds.partition_point(|&e| e <= query.start);
    started.saturating_sub(ended)
}

/// Whether the closed form is worth its extra end sort for this group.
///
/// The sweep costs O(candidate pairs), the closed form O((N + M) log M), so the
/// formula wins only once the data is dense enough. On sparse input — the common
/// biological case, reads against reads — the end sort is not repaid, so taking
/// the formula unconditionally would be a regression there.
///
/// The pair count is estimated from a uniform sample of queries, counting the
/// targets whose START falls in `[q.start, q.end + slack)`. That is available
/// from the start-sorted targets alone, needs no end sort, and is a lower bound
/// on the true overlap count, so the test only switches when the win is
/// unambiguous. Sampling rather than a uniform-density model is what makes it
/// hold up on real genomic data, which clusters heavily.
pub(crate) fn formula_beats_sweep<C: GroupType, T: PositionType>(
    left_group: &[IntervalRecord<C, T>],
    right_group: &[IntervalRecord<C, T>],
    slack: T,
) -> bool {
    // Each sample costs two binary searches, so the sample has to stay small
    // relative to the sweep it is deciding about — a group of 40k queries is
    // only ~40k sweep steps, which a 4096-sample probe would dominate. A few
    // hundred is enough to separate the regimes given the 4x margin below.
    const SAMPLE: usize = 256;

    let n = left_group.len();
    let m = right_group.len();
    if n == 0 || m == 0 {
        return false;
    }

    let step = (n / n.min(SAMPLE)).max(1);
    let mut sampled: u64 = 0;
    let mut pairs_in_sample: u64 = 0;

    let mut il = 0usize;
    while il < n {
        let query = left_group[il];
        let query_end_slack = query.end.saturating_add(slack);
        let hi = right_group.partition_point(|r| r.start < query_end_slack);
        let lo = right_group.partition_point(|r| r.start < query.start);
        // `hi < lo` for an inverted query (`end < start`), which is malformed
        // input the sweep still tolerates; saturate rather than underflow.
        pairs_in_sample += hi.saturating_sub(lo) as u64;
        sampled += 1;
        il += step;
    }

    // Integer throughout: `usize::BITS - leading_zeros` is the bit length, i.e.
    // floor(log2(m)) + 1, which is all the precision a threshold needs and
    // avoids pulling libm's `log2` into the crate. u128 because the extrapolated
    // pair count can exceed u64 on large dense groups.
    let est_pairs = (pairs_in_sample as u128) * (n as u128) / (sampled as u128);
    let m_log_m = (m as u128) * (usize::BITS - m.leading_zeros()) as u128;
    // The factor keeps a margin around the sparse regime, where a wrong switch
    // costs the end sort; the dense pathology clears it by an order of magnitude.
    est_pairs > 4 * m_log_m
}

/// Whether this overlap type needs the `(group, start, end, idx)` sort key.
///
/// `All` emits every pair, so the sweep order cannot change *which* pairs come
/// out. `First`/`Last` emit one target per query, so the order of equally
/// positioned records decides which one — see [`sorted_records`].
fn wide_sort_key(overlap_type: OverlapType) -> bool {
    match overlap_type {
        OverlapType::All => false,
        OverlapType::First | OverlapType::Last => true,
    }
}

/// Run the overlap sweep against pre-sorted left/right records. Used when
/// callers (e.g. `nearest`) need to share the sorted views with other sweeps.
pub(crate) fn collect_overlap_pairs_from_sorted<C: GroupType, T: PositionType>(
    left: &[IntervalRecord<C, T>],
    right: &[IntervalRecord<C, T>],
    slack: T,
    overlap_type: OverlapType,
    contained: bool,
) -> Vec<OverlapPair> {
    let n1 = left.len();
    let n2 = right.len();

    let mut pairs: Vec<OverlapPair> = Vec::new();

    if n1 == 0 || n2 == 0 {
        return pairs;
    }

    let mut i = 0usize;
    let mut j = 0usize;
    let mut active: Vec<usize> = Vec::new();

    while i < n1 && j < n2 {
        let g1 = left[i].group;
        let g2 = right[j].group;

        if g1 < g2 {
            let group = g1;
            while i < n1 && left[i].group == group {
                i += 1;
            }
            continue;
        }

        if g2 < g1 {
            let group = g2;
            while j < n2 && right[j].group == group {
                j += 1;
            }
            continue;
        }

        let group = g1;

        let i0 = i;
        while i < n1 && left[i].group == group {
            i += 1;
        }
        let i1 = i;

        let j0 = j;
        while j < n2 && right[j].group == group {
            j += 1;
        }
        let j1 = j;

        sweep_pairs_in_group(
            &left[i0..i1],
            &right[j0..j1],
            slack,
            overlap_type,
            contained,
            &mut active,
            &mut pairs,
        );
    }

    pairs
}

/// The overlap sweep for one group, whose left and right records are each a
/// contiguous slice. Every group is independent of every other, so this is the
/// unit both the serial walk above and the parallel driver hand out.
///
/// `active` is passed in only so a serial caller can reuse its allocation across
/// groups; it is cleared on entry.
#[allow(clippy::too_many_arguments)]
fn sweep_pairs_in_group<C: GroupType, T: PositionType>(
    left: &[IntervalRecord<C, T>],
    right: &[IntervalRecord<C, T>],
    slack: T,
    overlap_type: OverlapType,
    contained: bool,
    active: &mut Vec<usize>,
    pairs: &mut Vec<OverlapPair>,
) {
    let mut active_head = 0usize;
    clear_active(active, &mut active_head);
    let mut jr = 0usize;

    for query in left.iter().copied() {
        if contained {
            let query_start_slack = query.start.saturating_sub(slack);

            // A target that starts exactly at the query boundary must be active
            // for containment checks.
            while jr < right.len() && right[jr].start <= query_start_slack {
                active.push(jr);
                jr += 1;
            }

            // Do not expire containment candidates by their end coordinate.
            // In particular, negative slack can contract a query so far that
            // its effective end is before its effective start. Such a query
            // can still satisfy the established pair predicate for a target
            // ending at (or before) the effective start. Expiring candidates
            // here made that result depend on whether an earlier long target
            // prevented the head of this start-sorted list from advancing.
            // Keep every target admitted by the start bound and let the
            // predicate below decide each pair independently.

            match overlap_type {
                OverlapType::All => {
                    for idx in active_head..active.len() {
                        let target = right[active[idx]];
                        if !query_contained_in_target_with_slack(
                            query.start,
                            query.end,
                            target.start,
                            target.end,
                            slack,
                        ) {
                            continue;
                        }

                        pairs.push(OverlapPair {
                            idx: query.idx,
                            idx2: target.idx,
                        });
                    }
                }
                OverlapType::First => {
                    for idx in active_head..active.len() {
                        let target = right[active[idx]];
                        if !query_contained_in_target_with_slack(
                            query.start,
                            query.end,
                            target.start,
                            target.end,
                            slack,
                        ) {
                            continue;
                        }

                        pairs.push(OverlapPair {
                            idx: query.idx,
                            idx2: target.idx,
                        });
                        break;
                    }
                }
                OverlapType::Last => {
                    let mut last_target: Option<IntervalRecord<C, T>> = None;

                    for idx in active_head..active.len() {
                        let target = right[active[idx]];
                        if !query_contained_in_target_with_slack(
                            query.start,
                            query.end,
                            target.start,
                            target.end,
                            slack,
                        ) {
                            continue;
                        }

                        last_target = Some(target);
                    }

                    if let Some(target) = last_target {
                        pairs.push(OverlapPair {
                            idx: query.idx,
                            idx2: target.idx,
                        });
                    }
                }
            }

            continue;
        }

        let query_end_slack = query.end.saturating_add(slack);

        while jr < right.len() && right[jr].start < query_end_slack {
            active.push(jr);
            jr += 1;
        }

        while active_head < active.len() {
            let r = active[active_head];
            if right[r].end.saturating_add(slack) <= query.start {
                active_head += 1;
            } else {
                break;
            }
        }

        if active_head > 0 && active_head * 2 >= active.len() {
            active.drain(0..active_head);
            active_head = 0;
        }

        match overlap_type {
            OverlapType::All => {
                for idx in active_head..active.len() {
                    let target = right[active[idx]];

                    if !overlaps_with_slack(
                        query.start,
                        query.end,
                        target.start,
                        target.end,
                        slack,
                    ) {
                        continue;
                    }

                    pairs.push(OverlapPair {
                        idx: query.idx,
                        idx2: target.idx,
                    });
                }
            }
            OverlapType::First => {
                for idx in active_head..active.len() {
                    let target = right[active[idx]];

                    if !overlaps_with_slack(
                        query.start,
                        query.end,
                        target.start,
                        target.end,
                        slack,
                    ) {
                        continue;
                    }

                    pairs.push(OverlapPair {
                        idx: query.idx,
                        idx2: target.idx,
                    });
                    break;
                }
            }
            OverlapType::Last => {
                let mut last_target: Option<IntervalRecord<C, T>> = None;

                for idx in active_head..active.len() {
                    let target = right[active[idx]];

                    if !overlaps_with_slack(
                        query.start,
                        query.end,
                        target.start,
                        target.end,
                        slack,
                    ) {
                        continue;
                    }

                    last_target = Some(target);
                }

                if let Some(target) = last_target {
                    pairs.push(OverlapPair {
                        idx: query.idx,
                        idx2: target.idx,
                    });
                }
            }
        }
    }
}

/// Parallel driver for [`collect_overlap_pairs_from_sorted`]: one task per batch
/// of groups, results concatenated in group order, which is the order the serial
/// group walk emits them in.
fn collect_overlap_pairs_par<C: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    left: &SortedRecords<C, T>,
    right: &SortedRecords<C, T>,
    slack: T,
    overlap_type: OverlapType,
    contained: bool,
) -> Vec<OverlapPair> {
    let pairs = paired_group_runs(left, right);
    par_paired_groups(&pairs, |left_run, right_run, out| {
        let mut active = Vec::new();
        sweep_pairs_in_group(
            &left.records[left_run.0..left_run.1],
            &right.records[right_run.0..right_run.1],
            slack,
            overlap_type,
            contained,
            &mut active,
            out,
        );
    })
}

#[allow(clippy::too_many_arguments)]
pub fn overlaps<C: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
    overlap_type: &str,
    sort_output: bool,
    contained: bool,
) -> (Vec<u32>, Vec<u32>) {
    let overlap_type = OverlapType::from_str(overlap_type).expect("invalid overlap_type string");

    let mut pairs = collect_overlap_pairs(
        chrs,
        starts,
        ends,
        chrs2,
        starts2,
        ends2,
        slack,
        overlap_type,
        contained,
    );

    if sort_output {
        sort_by_key(&mut pairs, |p| (p.idx, p.idx2));
    }

    pairs.into_iter().map(|pair| (pair.idx, pair.idx2)).unzip()
}

/// [`overlaps`] over `i64` coordinates narrow enough to sweep as `i32`, which
/// sorts a third less data for the same answer.
///
/// `None` when the narrowing does not apply — see
/// [`sorted_records_narrowed_par`]. The result is a pair of row numbers, so
/// nothing widens back.
#[allow(clippy::too_many_arguments)]
pub fn overlaps_narrowed<C: GroupType + Send + Sync>(
    chrs: &[C],
    starts: &[i64],
    ends: &[i64],
    chrs2: &[C],
    starts2: &[i64],
    ends2: &[i64],
    slack: i64,
    overlap_type: &str,
    sort_output: bool,
    contained: bool,
) -> Option<(Vec<u32>, Vec<u32>)> {
    let overlap_type = OverlapType::from_str(overlap_type).expect("invalid overlap_type string");
    if !pair_fits_narrow(starts, ends, starts2, ends2, slack) {
        return None;
    }

    let wide_key = wide_sort_key(overlap_type);
    let left = sorted_records_narrowed_par(chrs, starts, ends, wide_key)?;
    let right = sorted_records_narrowed_par(chrs2, starts2, ends2, wide_key)?;
    let mut pairs =
        collect_overlap_pairs_par(&left, &right, slack as i32, overlap_type, contained);

    if sort_output {
        sort_by_key(&mut pairs, |p| (p.idx, p.idx2));
    }

    Some(pairs.into_iter().map(|pair| (pair.idx, pair.idx2)).unzip())
}

/// How many right intervals each left interval overlaps, one count per left row,
/// in input row order.
///
/// Each group is counted on its own thread, and a group dense enough for the
/// sweep's cost to matter takes the BITS closed form instead — the same two
/// choices [`crate::overlap_indices::overlap_indices_chunked`] makes, for the
/// same reasons. Before this, `count_overlaps` was the last kernel still doing a
/// serial radix sort followed by a single serial pass whose cost grew with the
/// number of candidate *pairs*, which is quadratic on deep-coverage input.
pub fn count_overlaps<C: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
) -> Vec<u32> {
    let mut counts = vec![0_u32; chrs.len()];
    if chrs.is_empty() || chrs2.is_empty() {
        return counts;
    }

    // A count cannot see *which* target was matched, so both sides take the
    // minimal `(group, start)` key.
    let left = sorted_records_par(chrs, starts, ends, false);
    let right = sorted_records_par(chrs2, starts2, ends2, false);
    count_from_sorted(&left, &right, slack, &mut counts);
    counts
}

/// [`count_overlaps`] over `i64` coordinates narrow enough to sweep as `i32`.
///
/// `None` when the narrowing does not apply — see
/// [`sorted_records_narrowed_par`] — leaving the caller on [`count_overlaps`].
/// The counts themselves are `u32` either way, so nothing widens back.
pub fn count_overlaps_narrowed<C: GroupType + Send + Sync>(
    chrs: &[C],
    starts: &[i64],
    ends: &[i64],
    chrs2: &[C],
    starts2: &[i64],
    ends2: &[i64],
    slack: i64,
) -> Option<Vec<u32>> {
    if chrs.is_empty() || chrs2.is_empty() {
        return None;
    }
    if !pair_fits_narrow(starts, ends, starts2, ends2, slack) {
        return None;
    }

    let left = sorted_records_narrowed_par(chrs, starts, ends, false)?;
    let right = sorted_records_narrowed_par(chrs2, starts2, ends2, false)?;

    let mut counts = vec![0_u32; chrs.len()];
    count_from_sorted(&left, &right, slack as i32, &mut counts);
    Some(counts)
}

/// Count each left row against pre-sorted, group-major records.
///
/// Only groups both sides have can contribute, so the rest of `counts` keeps the
/// zero it was created with.
fn count_from_sorted<C: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    left: &SortedRecords<C, T>,
    right: &SortedRecords<C, T>,
    slack: T,
    counts: &mut [u32],
) {
    let pairs = paired_group_runs(left, right);
    let counts = partition::ScatterSlice::new(counts);
    par_paired_groups_for_each(&pairs, |left_run, right_run| {
        count_group(
            &left.records[left_run.0..left_run.1],
            &right.records[right_run.0..right_run.1],
            slack,
            &counts,
        );
    });
}

/// One group's counts, by whichever of the two methods is cheaper for it.
fn count_group<C: GroupType, T: PositionType>(
    left: &[IntervalRecord<C, T>],
    right: &[IntervalRecord<C, T>],
    slack: T,
    counts: &partition::ScatterSlice<'_, u32>,
) {
    // SAFETY, for every `counts.write` below: a record's `idx` is its input row
    // number, each row belongs to exactly one group, and each group is visited by
    // exactly one task — so no index is written twice.
    if slack >= T::zero()
        && !group_has_nonpositive_length(right)
        && formula_beats_sweep(left, right, slack)
    {
        let end_bounds = sorted_end_bounds(right, slack);
        for query in left {
            let count = count_for_query(query, right, &end_bounds, slack);
            unsafe { counts.write(query.idx as usize, saturate_u32(count)) };
        }
        return;
    }

    let mut active: Vec<usize> = Vec::new();
    let mut active_head = 0usize;
    let mut jr = 0usize;

    for query in left.iter().copied() {
        let query_end_slack = query.end.saturating_add(slack);

        while jr < right.len() && right[jr].start < query_end_slack {
            active.push(jr);
            jr += 1;
        }

        while active_head < active.len() {
            let target = right[active[active_head]];
            if target.end.saturating_add(slack) <= query.start {
                active_head += 1;
            } else {
                break;
            }
        }

        if active_head > 0 && active_head * 2 >= active.len() {
            active.drain(0..active_head);
            active_head = 0;
        }

        let mut count = 0_u32;
        for &index in &active[active_head..] {
            let target = right[index];
            if overlaps_with_slack(query.start, query.end, target.start, target.end, slack) {
                count = count.saturating_add(1);
            }
        }

        unsafe { counts.write(query.idx as usize, count) };
    }
}

#[inline]
fn saturate_u32(count: usize) -> u32 {
    count.min(u32::MAX as usize) as u32
}

#[cfg(test)]
mod tests {
    use super::overlaps;

    type Group = u32;
    type Pos = i64;

    #[test]
    fn overlaps_all_returns_expected_pairs() {
        let groups: [Group; 3] = [1, 1, 1];
        let starts: [Pos; 3] = [1, 10, 30];
        let ends: [Pos; 3] = [5, 20, 40];

        let groups2: [Group; 4] = [1, 1, 1, 1];
        let starts2: [Pos; 4] = [3, 11, 18, 35];
        let ends2: [Pos; 4] = [4, 12, 19, 36];

        let (idx1, idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, false,
        );

        assert_eq!(idx1, vec![0, 1, 1, 2]);
        assert_eq!(idx2, vec![0, 1, 2, 3]);
    }

    #[test]
    fn overlaps_issue_23_preserves_order_for_all_first_last() {
        let groups = vec![1_u32; 10];
        let starts: Vec<Pos> = vec![97, 981, 1000, 1227, 1409, 3889, 4398, 4815, 5004, 5047];
        let ends: Vec<Pos> = vec![1097, 1981, 2000, 2227, 2409, 4889, 5398, 5815, 6004, 6047];

        let groups2 = vec![1_u32; 10];
        let starts2: Vec<Pos> = vec![1476, 2110, 2823, 3547, 3578, 4295, 5184, 5545, 7117, 7190];
        let ends2: Vec<Pos> = vec![2476, 3110, 3823, 4547, 4578, 5295, 6184, 6545, 8117, 8190];

        let (all_idx1, all_idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, false,
        );

        assert_eq!(
            all_idx1,
            vec![1, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9]
        );
        assert_eq!(
            all_idx2,
            vec![0, 0, 0, 1, 0, 1, 3, 4, 5, 3, 4, 5, 6, 5, 6, 7, 5, 6, 7, 5, 6, 7]
        );

        let (first_idx1, first_idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "first", true, false,
        );
        assert_eq!(first_idx1, vec![1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(first_idx2, vec![0, 0, 0, 0, 3, 3, 5, 5, 5]);

        let (last_idx1, last_idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "last", true, false,
        );
        assert_eq!(last_idx1, vec![1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(last_idx2, vec![0, 0, 1, 1, 5, 6, 7, 7, 7]);
    }

    #[test]
    fn overlaps_contained_filters_to_left_contained_in_right() {
        let groups: [Group; 2] = [1, 1];
        let starts: [Pos; 2] = [10, 10];
        let ends: [Pos; 2] = [12, 30];

        let groups2: [Group; 2] = [1, 1];
        let starts2: [Pos; 2] = [9, 11];
        let ends2: [Pos; 2] = [40, 29];

        let (idx1, idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, true,
        );

        assert_eq!(idx1, vec![0, 1]);
        assert_eq!(idx2, vec![0, 0]);
    }

    #[test]
    fn overlaps_contained_includes_targets_with_same_start() {
        let groups: [Group; 2] = [1, 1];
        let starts: [Pos; 2] = [1, 6];
        let ends: [Pos; 2] = [3, 9];

        let groups2: [Group; 1] = [1];
        let starts2: [Pos; 1] = [1];
        let ends2: [Pos; 1] = [9];

        let (idx1, idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, true,
        );

        assert_eq!(idx1, vec![0, 1]);
        assert_eq!(idx2, vec![0, 0]);
    }

    #[test]
    fn overlaps_contained_with_negative_slack_matches_legacy_behavior() {
        let groups: [Group; 2] = [1, 1];
        let starts: [Pos; 2] = [1, 4];
        let ends: [Pos; 2] = [3, 9];

        let groups2: [Group; 3] = [1, 1, 1];
        let starts2: [Pos; 3] = [1, 2, 2];
        let ends2: [Pos; 3] = [10, 3, 9];

        let (idx1, idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, -2, "first", true, true,
        );

        assert_eq!(idx1, vec![0, 1]);
        assert_eq!(idx2, vec![0, 0]);
    }

    #[test]
    fn overlaps_contained_negative_slack_is_independent_of_other_targets() {
        let groups: [Group; 1] = [0];
        let starts: [Pos; 1] = [10];
        let ends: [Pos; 1] = [12];

        let target_groups: [Group; 1] = [0];
        let target_starts: [Pos; 1] = [10];
        let target_ends: [Pos; 1] = [12];
        let (alone_left, alone_right) = overlaps(
            &groups,
            &starts,
            &ends,
            &target_groups,
            &target_starts,
            &target_ends,
            -2,
            "all",
            true,
            true,
        );

        let targets_groups: [Group; 2] = [0, 0];
        let targets_starts: [Pos; 2] = [0, 10];
        let targets_ends: [Pos; 2] = [100, 12];
        let (with_other_left, with_other_right) = overlaps(
            &groups,
            &starts,
            &ends,
            &targets_groups,
            &targets_starts,
            &targets_ends,
            -2,
            "all",
            true,
            true,
        );

        assert_eq!((alone_left, alone_right), (vec![0], vec![0]));
        assert_eq!(
            (with_other_left, with_other_right),
            (vec![0, 0], vec![0, 1])
        );
    }

    #[test]
    fn overlaps_handles_unsorted_input() {
        let groups: [Group; 3] = [1, 1, 1];
        let starts: [Pos; 3] = [10, 0, 20];
        let ends: [Pos; 3] = [12, 5, 30];

        let groups2: [Group; 2] = [1, 1];
        let starts2: [Pos; 2] = [11, 3];
        let ends2: [Pos; 2] = [13, 4];

        let (idx1, idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, false,
        );

        assert_eq!(idx1, vec![0, 1]);
        assert_eq!(idx2, vec![0, 1]);
    }

    #[test]
    fn overlaps_respects_slack_for_bookended_intervals() {
        let groups: [Group; 1] = [1];
        let starts: [Pos; 1] = [10];
        let ends: [Pos; 1] = [20];

        let groups2: [Group; 1] = [1];
        let starts2: [Pos; 1] = [20];
        let ends2: [Pos; 1] = [25];

        let (idx1_no_slack, idx2_no_slack) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, false,
        );
        assert!(idx1_no_slack.is_empty());
        assert!(idx2_no_slack.is_empty());

        let (idx1_with_slack, idx2_with_slack) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 1, "all", true, false,
        );
        assert_eq!(idx1_with_slack, vec![0]);
        assert_eq!(idx2_with_slack, vec![0]);
    }
}

#[cfg(test)]
mod sorted_records_tests {
    use super::{sorted_records, sorted_records_par, OverlapType};

    /// Deterministic SplitMix64, so this test needs no dev-dependency.
    fn next(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// The parallel builder claims byte-identical output for the matching key,
    /// including the `First`/`Last` tie-break that pyranges/ruranges_py#23 made
    /// load-bearing. Ties are what this pins, so the coordinates are drawn from
    /// a small range to force plenty of them.
    #[test]
    fn parallel_records_match_the_serial_ones_including_ties() {
        let mut state = 0x2718;
        for &n in &[7_usize, 999, 150_000] {
            for &groups in &[1_u64, 24] {
                let mut chrs = Vec::with_capacity(n);
                let mut starts = Vec::with_capacity(n);
                let mut ends = Vec::with_capacity(n);
                for _ in 0..n {
                    chrs.push((next(&mut state) % groups) as u32);
                    let start = (next(&mut state) % 50) as i64;
                    starts.push(start);
                    ends.push(start + (next(&mut state) % 20) as i64);
                }

                for overlap_type in [OverlapType::All, OverlapType::First, OverlapType::Last] {
                    let wide = !matches!(overlap_type, OverlapType::All);
                    let serial = sorted_records(&chrs, &starts, &ends, overlap_type);
                    let parallel = sorted_records_par(&chrs, &starts, &ends, wide);
                    let key = |r: &super::IntervalRecord<u32, i64>| (r.group, r.start, r.end, r.idx);
                    assert_eq!(
                        serial.iter().map(key).collect::<Vec<_>>(),
                        parallel.records.iter().map(key).collect::<Vec<_>>(),
                        "n={n} groups={groups} type={overlap_type:?}"
                    );
                    // The runs must cover every record, ascending by group.
                    assert_eq!(
                        parallel.runs.iter().map(|&(s, e)| e - s).sum::<usize>(),
                        n
                    );
                }
            }
        }
    }
}

#[cfg(test)]
mod count_overlaps_tests {
    use super::{count_overlaps, count_overlaps_narrowed};

    /// Deterministic SplitMix64, so these tests need no dev-dependency.
    fn next(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn intervals(
        state: &mut u64,
        n: usize,
        groups: u64,
        span: u64,
        max_len: u64,
    ) -> (Vec<u32>, Vec<i64>, Vec<i64>) {
        let mut chrs = Vec::with_capacity(n);
        let mut starts = Vec::with_capacity(n);
        let mut ends = Vec::with_capacity(n);
        for _ in 0..n {
            chrs.push((next(state) % groups) as u32);
            let start = (next(state) % span) as i64;
            starts.push(start);
            ends.push(start + 1 + (next(state) % max_len) as i64);
        }
        (chrs, starts, ends)
    }

    /// Every pair, tested one at a time. Independent of the kernel, so it pins
    /// the answer rather than restating the algorithm.
    fn brute_force(
        chrs: &[u32],
        starts: &[i64],
        ends: &[i64],
        chrs2: &[u32],
        starts2: &[i64],
        ends2: &[i64],
        slack: i64,
    ) -> Vec<u32> {
        (0..chrs.len())
            .map(|q| {
                (0..chrs2.len())
                    .filter(|&t| {
                        chrs2[t] == chrs[q]
                            && starts[q] < ends2[t] + slack
                            && starts2[t] < ends[q] + slack
                    })
                    .count() as u32
            })
            .collect()
    }

    #[test]
    fn counts_match_brute_force_on_small_random_input() {
        let mut state = 0xC0FFEE;
        for slack in [0_i64, 1, 25] {
            for _ in 0..20 {
                let (chrs, starts, ends) = intervals(&mut state, 60, 3, 200, 40);
                let (chrs2, starts2, ends2) = intervals(&mut state, 60, 4, 200, 40);
                assert_eq!(
                    count_overlaps(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack),
                    brute_force(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack),
                    "slack={slack}"
                );
            }
        }
    }

    /// A group deep enough to take the closed form must agree with a group that
    /// takes the sweep, and both with the brute force.
    #[test]
    fn the_closed_form_and_the_sweep_agree_on_a_deep_pileup() {
        let mut state = 0xD3E9;
        // One group, long targets over a short span: every query overlaps most
        // of the targets, which is the regime `formula_beats_sweep` switches in.
        let (chrs, starts, ends) = intervals(&mut state, 4_000, 1, 500, 5);
        let (chrs2, starts2, ends2) = intervals(&mut state, 4_000, 1, 500, 400);
        assert_eq!(
            count_overlaps(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, 0_i64),
            brute_force(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, 0),
        );
    }

    /// Above `PARALLEL_THRESHOLD` the records are partitioned and the groups are
    /// counted on separate threads; the answer must not depend on that, and the
    /// narrowed path must return the same counts as the wide one.
    #[test]
    fn parallel_and_narrowed_paths_agree_with_the_wide_one() {
        let mut state = 0xBEEF;
        let (chrs, starts, ends) = intervals(&mut state, 150_000, 5, 400_000, 500);
        let (chrs2, starts2, ends2) = intervals(&mut state, 150_000, 5, 400_000, 500);

        for slack in [0_i64, 250] {
            let wide = count_overlaps(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack);
            assert_eq!(
                wide,
                brute_force_by_group(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack),
                "slack={slack}"
            );
            assert_eq!(
                count_overlaps_narrowed(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack),
                Some(wide),
                "slack={slack}"
            );
        }
    }

    /// A coordinate, or a slack, too wide for `i32` declines rather than
    /// truncating — silently wrapping a coordinate is what this pins.
    #[test]
    fn coordinates_and_slacks_too_wide_to_narrow_decline() {
        let mut state = 0x5A5A;
        let (chrs, mut starts, mut ends) = intervals(&mut state, 150_000, 4, 400_000, 500);
        let (chrs2, starts2, ends2) = intervals(&mut state, 150_000, 4, 400_000, 500);

        assert!(
            count_overlaps_narrowed(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, 0).is_some()
        );
        assert!(count_overlaps_narrowed(
            &chrs,
            &starts,
            &ends,
            &chrs2,
            &starts2,
            &ends2,
            crate::sorts::NARROW_COORDINATE_LIMIT
        )
        .is_none());

        starts[7] = crate::sorts::NARROW_COORDINATE_LIMIT + 1;
        ends[7] = crate::sorts::NARROW_COORDINATE_LIMIT + 11;
        assert!(
            count_overlaps_narrowed(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, 0).is_none()
        );
    }

    /// The brute force above is O(n*m); at 150k rows it needs the group split.
    fn brute_force_by_group(
        chrs: &[u32],
        starts: &[i64],
        ends: &[i64],
        chrs2: &[u32],
        starts2: &[i64],
        ends2: &[i64],
        slack: i64,
    ) -> Vec<u32> {
        let mut by_group: std::collections::HashMap<u32, Vec<usize>> = Default::default();
        for t in 0..chrs2.len() {
            by_group.entry(chrs2[t]).or_default().push(t);
        }
        (0..chrs.len())
            .map(|q| {
                by_group
                    .get(&chrs[q])
                    .map(|targets| {
                        targets
                            .iter()
                            .filter(|&&t| {
                                starts[q] < ends2[t] + slack && starts2[t] < ends[q] + slack
                            })
                            .count() as u32
                    })
                    .unwrap_or(0)
            })
            .collect()
    }
}

