use std::str::FromStr;

use radsort::sort_by_key;

use crate::ruranges_structs::{GroupType, OverlapPair, OverlapType, PositionType};

#[derive(Copy, Clone, Debug)]
pub(crate) struct IntervalRecord<C: GroupType, T: PositionType> {
    pub(crate) group: C,
    pub(crate) start: T,
    pub(crate) end: T,
    pub(crate) idx: u32,
}

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

fn collect_overlap_pairs<C: GroupType, T: PositionType>(
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
    let left = sorted_records(chrs, starts, ends, overlap_type);
    let right = sorted_records(chrs2, starts2, ends2, overlap_type);
    collect_overlap_pairs_from_sorted(&left, &right, slack, overlap_type, contained)
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
    let mut active_head: usize = 0;

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

        clear_active(&mut active, &mut active_head);

        let mut jr = j0;

        for il in i0..i1 {
            let query = left[il];
            if contained {
                let query_start_slack = query.start.saturating_sub(slack);

                // A target that starts exactly at the query boundary must be active
                // for containment checks.
                while jr < j1 && right[jr].start <= query_start_slack {
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

            while jr < j1 && right[jr].start < query_end_slack {
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

    pairs
}

#[allow(clippy::too_many_arguments)]
pub fn overlaps<C: GroupType, T: PositionType>(
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

pub fn count_overlaps<C: GroupType, T: PositionType>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
) -> Vec<u32> {
    // count_overlaps only returns a count per left interval, so it never cares
    // which target was matched. Use the minimal `All` sort key.
    let left = sorted_records(chrs, starts, ends, OverlapType::All);
    let right = sorted_records(chrs2, starts2, ends2, OverlapType::All);

    let n1 = left.len();
    let n2 = right.len();

    let mut counts = vec![0_u32; n1];

    if n1 == 0 || n2 == 0 {
        return counts;
    }

    let mut i = 0usize;
    let mut j = 0usize;

    let mut active: Vec<usize> = Vec::new();
    let mut active_head: usize = 0;

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

        clear_active(&mut active, &mut active_head);

        let mut jr = j0;

        for il in i0..i1 {
            let query = left[il];
            let query_end_slack = query.end.saturating_add(slack);

            while jr < j1 && right[jr].start < query_end_slack {
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

            let mut count = 0_u32;

            for idx in active_head..active.len() {
                let target = right[active[idx]];

                if overlaps_with_slack(query.start, query.end, target.start, target.end, slack) {
                    count = count.saturating_add(1);
                }
            }

            counts[query.idx as usize] = count;
        }
    }

    counts
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
