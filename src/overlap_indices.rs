use std::str::FromStr;

use crate::ruranges_structs::{GroupType, OverlapType, PositionType};

#[derive(Copy, Clone, Debug)]
struct IntervalRecord<C: GroupType, T: PositionType> {
    group: C,
    start: T,
    end: T,
    idx: u32,
}

fn sorted_records<C: GroupType, T: PositionType>(
    groups: &[C],
    starts: &[T],
    ends: &[T],
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

    radsort::sort_by_key(&mut records, |r| (r.group, r.start, r.end, r.idx));
    records
}

#[inline(always)]
fn overlaps_with_slack<T: PositionType>(
    query_start: T,
    query_end: T,
    target_start: T,
    target_end: T,
    slack: T,
) -> bool {
    query_start < target_end.saturating_add(slack) && target_start < query_end.saturating_add(slack)
}

#[inline(always)]
fn query_contained_in_target_with_slack<T: PositionType>(
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

fn clear_active(active: &mut Vec<usize>, active_head: &mut usize) {
    active.clear();
    *active_head = 0;
}

#[allow(clippy::too_many_arguments)]
pub fn overlap_indices<C: GroupType, T: PositionType>(
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
) -> Vec<u32> {
    let overlap_type = OverlapType::from_str(overlap_type).expect("invalid overlap_type string");
    let keep_all_matches = matches!(overlap_type, OverlapType::All);
    let left = sorted_records(chrs, starts, ends);
    let right = sorted_records(chrs2, starts2, ends2);

    let n1 = left.len();
    let n2 = right.len();

    let mut indices = Vec::new();
    if n1 == 0 || n2 == 0 {
        return indices;
    }

    let mut i = 0usize;
    let mut j = 0usize;
    let mut active: Vec<usize> = Vec::new();
    let mut active_head = 0usize;

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

        for query in left.iter().take(i1).skip(i0).copied() {
            if contained {
                let query_start_slack = query.start.saturating_sub(slack);

                while jr < j1 && right[jr].start <= query_start_slack {
                    active.push(jr);
                    jr += 1;
                }

                while active_head < active.len() {
                    let r = active[active_head];
                    if right[r].end <= query_start_slack {
                        active_head += 1;
                    } else {
                        break;
                    }
                }

                if active_head > 0 && active_head * 2 >= active.len() {
                    active.drain(0..active_head);
                    active_head = 0;
                }

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

                    indices.push(query.idx);
                    if !keep_all_matches {
                        break;
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

            for idx in active_head..active.len() {
                let target = right[active[idx]];
                if !overlaps_with_slack(query.start, query.end, target.start, target.end, slack) {
                    continue;
                }

                indices.push(query.idx);
                if !keep_all_matches {
                    break;
                }
            }
        }
    }

    if sort_output {
        radsort::sort(&mut indices);
    }

    indices
}

#[cfg(test)]
mod tests {
    use super::overlap_indices;
    use crate::overlaps::overlaps;

    type Group = u32;
    type Pos = i64;

    #[test]
    fn overlap_indices_matches_left_side_of_overlaps() {
        let groups: [Group; 5] = [2, 1, 1, 2, 1];
        let starts: [Pos; 5] = [30, 1, 10, 35, 18];
        let ends: [Pos; 5] = [40, 5, 20, 45, 24];
        let groups2: [Group; 5] = [1, 2, 1, 2, 1];
        let starts2: [Pos; 5] = [3, 33, 11, 43, 19];
        let ends2: [Pos; 5] = [4, 34, 12, 44, 23];

        for overlap_type in ["all", "first", "last"] {
            let (left, _) = overlaps(
                &groups,
                &starts,
                &ends,
                &groups2,
                &starts2,
                &ends2,
                0,
                overlap_type,
                true,
                false,
            );
            let left_only = overlap_indices(
                &groups,
                &starts,
                &ends,
                &groups2,
                &starts2,
                &ends2,
                0,
                overlap_type,
                true,
                false,
            );

            assert_eq!(left_only, left);
        }
    }

    #[test]
    fn overlap_indices_matches_contained_overlaps() {
        let groups: [Group; 3] = [1, 1, 1];
        let starts: [Pos; 3] = [5, 10, 20];
        let ends: [Pos; 3] = [8, 15, 25];
        let groups2: [Group; 2] = [1, 1];
        let starts2: [Pos; 2] = [1, 18];
        let ends2: [Pos; 2] = [16, 30];

        let (left, _) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, true,
        );
        let left_only = overlap_indices(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, true,
        );

        assert_eq!(left_only, left);
    }

    #[test]
    fn overlap_indices_can_skip_output_sort() {
        let groups: [Group; 3] = [1, 1, 1];
        let starts: [Pos; 3] = [20, 1, 10];
        let ends: [Pos; 3] = [25, 5, 15];
        let groups2: [Group; 3] = [1, 1, 1];
        let starts2: [Pos; 3] = [2, 11, 21];
        let ends2: [Pos; 3] = [3, 12, 22];

        let left_only = overlap_indices(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", false, false,
        );

        assert_eq!(left_only, vec![1, 2, 0]);
    }
}
