use std::str::FromStr;

use crate::overlaps::{
    clear_active, count_for_query, formula_beats_sweep, group_has_nonpositive_length,
    overlaps_with_slack, pair_fits_narrow, paired_group_runs, par_paired_group_batches,
    query_contained_in_target_with_slack, sorted_end_bounds, sorted_records_narrowed_par,
    sorted_records_par, IntervalRecord, SortedRecords,
};
use crate::ruranges_structs::{GroupType, OverlapType, PositionType};

/// Indexes of the left intervals that overlap something on the right, once per
/// match when `overlap_type` is `"all"`.
#[allow(clippy::too_many_arguments)]
pub fn overlap_indices<C: GroupType + Send + Sync, T: PositionType + Send + Sync>(
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
    let batches = overlap_indices_chunked(
        chrs, starts, ends, chrs2, starts2, ends2, slack, overlap_type, contained,
    );
    let mut indices = Vec::with_capacity(batches.iter().map(Vec::len).sum());
    for batch in batches {
        indices.extend(batch);
    }
    if sort_output {
        radsort::sort(&mut indices);
    }
    indices
}

/// Sort row indexes ascending, as `overlap_indices(.., sort_output = true)`
/// does. Exposed so a caller that took the chunked path can still ask for the
/// sorted order without pulling in a radix sort of its own.
pub fn sort_indices(mut indices: Vec<u32>) -> Vec<u32> {
    radsort::sort(&mut indices);
    indices
}

/// [`overlap_indices`] without the final concatenation, in group order.
///
/// A caller that is going to gather rows with these indexes can use the batches
/// as gather chunks directly. At 10^8 rows against a whole-genome annotation the
/// concatenated vector is several GiB, and both it and the copy into it are pure
/// overhead for that caller.
#[allow(clippy::too_many_arguments)]
pub fn overlap_indices_chunked<C: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
    overlap_type: &str,
    contained: bool,
) -> Vec<Vec<u32>> {
    let overlap_type = OverlapType::from_str(overlap_type).expect("invalid overlap_type string");

    // Unlike `overlaps`, this kernel returns only `query.idx` — it never emits a
    // target. The `(group, start, end, idx)` key that `sorted_records` uses for
    // First/Last exists to decide *which* target the sweep reaches first, which
    // cannot be observed here, so both sides take the narrow `(group, start)`
    // key and skip two radsort passes over 24-byte records per side.
    let left = sorted_records_par(chrs, starts, ends, false);
    let right = sorted_records_par(chrs2, starts2, ends2, false);
    indices_from_sorted(&left, &right, slack, overlap_type, contained)
}

/// [`overlap_indices_chunked`] over `i64` coordinates narrow enough to sweep as
/// `i32`, which sorts a third less data for the same answer.
///
/// `None` when the narrowing does not apply — see
/// [`crate::overlaps::sorted_records_narrowed_par`]. The result is row numbers,
/// so nothing widens back.
#[allow(clippy::too_many_arguments)]
pub fn overlap_indices_chunked_narrowed<C: GroupType + Send + Sync>(
    chrs: &[C],
    starts: &[i64],
    ends: &[i64],
    chrs2: &[C],
    starts2: &[i64],
    ends2: &[i64],
    slack: i64,
    overlap_type: &str,
    contained: bool,
) -> Option<Vec<Vec<u32>>> {
    let overlap_type = OverlapType::from_str(overlap_type).expect("invalid overlap_type string");
    if !pair_fits_narrow(starts, ends, starts2, ends2, slack) {
        return None;
    }

    let left = sorted_records_narrowed_par(chrs, starts, ends, false)?;
    let right = sorted_records_narrowed_par(chrs2, starts2, ends2, false)?;
    Some(indices_from_sorted(
        &left,
        &right,
        slack as i32,
        overlap_type,
        contained,
    ))
}

/// The per-group driver both entry points share: each group takes the closed
/// form or the sweep on its own thread.
fn indices_from_sorted<C: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    left: &SortedRecords<C, T>,
    right: &SortedRecords<C, T>,
    slack: T,
    overlap_type: OverlapType,
    contained: bool,
) -> Vec<Vec<u32>> {
    let keep_all_matches = matches!(overlap_type, OverlapType::All);
    let group_pairs = paired_group_runs(left, right);
    par_paired_group_batches(&group_pairs, move |left_run, right_run, out| {
        let left_group = &left.records[left_run.0..left_run.1];
        let right_group = &right.records[right_run.0..right_run.1];

        // `contained` is a different predicate, which the counting identity does
        // not cover.
        if !contained
            && slack >= T::zero()
            && !group_has_nonpositive_length(right_group)
            && formula_beats_sweep(left_group, right_group, slack)
        {
            count_indices_in_group(left_group, right_group, slack, keep_all_matches, out);
        } else {
            let mut active = Vec::new();
            sweep_indices_in_group(
                left_group,
                right_group,
                slack,
                keep_all_matches,
                contained,
                &mut active,
                out,
            );
        }
    })
}

/// The closed-form path: BITS gives each query its overlap count directly, and
/// the output this kernel produces is that count of copies of `query.idx` — or
/// one copy when the caller only asked whether anything overlapped.
fn count_indices_in_group<C: GroupType, T: PositionType>(
    left_group: &[IntervalRecord<C, T>],
    right_group: &[IntervalRecord<C, T>],
    slack: T,
    keep_all_matches: bool,
    out: &mut Vec<u32>,
) {
    let end_bounds = sorted_end_bounds(right_group, slack);
    for query in left_group {
        let count = count_for_query(query, right_group, &end_bounds, slack);
        if keep_all_matches {
            out.extend(std::iter::repeat_n(query.idx, count));
        } else if count > 0 {
            out.push(query.idx);
        }
    }
}

/// The overlap-filter sweep for one group, whose left and right records are each
/// a contiguous slice. `active` is passed in so a caller can reuse the
/// allocation across groups; it is cleared on entry.
#[allow(clippy::too_many_arguments)]
fn sweep_indices_in_group<C: GroupType, T: PositionType>(
    left: &[IntervalRecord<C, T>],
    right: &[IntervalRecord<C, T>],
    slack: T,
    keep_all_matches: bool,
    contained: bool,
    active: &mut Vec<usize>,
    indices: &mut Vec<u32>,
) {
    let mut active_head = 0usize;
    clear_active(active, &mut active_head);
    let mut jr = 0usize;

    for query in left.iter().copied() {
        if contained {
            let query_start_slack = query.start.saturating_sub(slack);

            while jr < right.len() && right[jr].start <= query_start_slack {
                active.push(jr);
                jr += 1;
            }

            // See `collect_overlap_pairs_from_sorted`: negative slack may make a
            // pair contained even when the target ends at the effective query
            // start. Keep candidates admitted by their start and evaluate the
            // pair predicate below.

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
    fn overlap_indices_contained_negative_slack_is_independent_of_other_targets() {
        let groups: [Group; 1] = [0];
        let starts: [Pos; 1] = [10];
        let ends: [Pos; 1] = [12];

        let target_groups: [Group; 1] = [0];
        let target_starts: [Pos; 1] = [10];
        let target_ends: [Pos; 1] = [12];
        let alone = overlap_indices(
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
        let with_other = overlap_indices(
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

        assert_eq!(alone, vec![0]);
        assert_eq!(with_other, vec![0, 0]);
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

#[cfg(test)]
mod parity_tests {
    use super::overlap_indices;

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

    /// Random intervals over `groups` chromosomes. A small `span` packs the same
    /// interval count into a short region, which is what pushes a group onto the
    /// closed-form path; a large one keeps it on the sweep.
    pub(super) fn random_intervals_for_tests(
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

    /// Every query index, repeated once per overlapping target, in query order.
    /// Independent of both routes through `overlap_indices`, so it pins whichever
    /// one the density test picks.
    fn brute_force(
        chrs: &[Group],
        starts: &[Pos],
        ends: &[Pos],
        chrs2: &[Group],
        starts2: &[Pos],
        ends2: &[Pos],
        slack: Pos,
        keep_all: bool,
    ) -> Vec<u32> {
        let mut out = Vec::new();
        for q in 0..chrs.len() {
            for t in 0..chrs2.len() {
                if chrs[q] == chrs2[t]
                    && starts[q] < ends2[t].saturating_add(slack)
                    && starts2[t] < ends[q].saturating_add(slack)
                {
                    out.push(q as u32);
                    if !keep_all {
                        break;
                    }
                }
            }
        }
        radsort::sort(&mut out);
        out
    }

    /// Both the dense case (which takes the BITS path) and the sparse one (which
    /// stays on the sweep), on inputs small enough to check against O(N*M).
    #[test]
    fn overlap_indices_matches_brute_force_on_random_input() {
        let mut state = 0xC0FFEE;
        for &(n, m, span) in &[
            (500_usize, 500_usize, 2_000_u64), // dense
            (500, 500, 5_000_000),             // sparse
            (2_000, 2_000, 8_000),             // dense, several groups deep
        ] {
            for &slack in &[0_i64, 5] {
                for keep_all in [false, true] {
                    let (chrs, starts, ends) = random_intervals_for_tests(&mut state, n, 4, span, 300);
                    let (chrs2, starts2, ends2) = random_intervals_for_tests(&mut state, m, 4, span, 300);
                    let expected = brute_force(
                        &chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack, keep_all,
                    );
                    let got = overlap_indices(
                        &chrs,
                        &starts,
                        &ends,
                        &chrs2,
                        &starts2,
                        &ends2,
                        slack,
                        if keep_all { "all" } else { "first" },
                        true,
                        false,
                    );
                    assert_eq!(got, expected, "n={n} m={m} span={span} slack={slack} keep_all={keep_all}");
                }
            }
        }
    }

    /// Targets of non-positive length break the counting identity, so those
    /// groups have to fall back to the sweep.
    #[test]
    fn overlap_indices_handles_degenerate_targets_on_a_dense_group() {
        let mut state = 0xDEAD;
        let (chrs, starts, ends) = random_intervals_for_tests(&mut state, 5_000, 1, 20_000, 200);
        let (chrs2, starts2, mut ends2) = random_intervals_for_tests(&mut state, 5_000, 1, 20_000, 200);
        for i in (0..ends2.len()).step_by(97) {
            ends2[i] = starts2[i]; // zero length
        }

        for keep_all in [false, true] {
            let expected =
                brute_force(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, 0, keep_all);
            let got = overlap_indices(
                &chrs,
                &starts,
                &ends,
                &chrs2,
                &starts2,
                &ends2,
                0,
                if keep_all { "all" } else { "first" },
                true,
                false,
            );
            assert_eq!(got, expected, "keep_all={keep_all}");
        }
    }
}

#[cfg(test)]
mod large_parity_tests {
    use super::super::overlaps::overlaps;
    use super::overlap_indices;
    use super::parity_tests::{random_intervals_for_tests as random_intervals};

    /// Above `PARALLEL_THRESHOLD` the records are partitioned and swept per
    /// group, and dense groups take the closed form. `overlaps` is the
    /// unmodified pair sweep, so its left column is an independent reference at
    /// a size brute force cannot reach.
    #[test]
    fn overlap_indices_matches_the_pair_sweep_above_the_parallel_threshold() {
        let mut state = 0xFEED_BEEF;
        for &(n, m, span) in &[
            (60_000_usize, 60_000_usize, 40_000_u64),   // dense: closed form
            (60_000, 60_000, 20_000_000),               // sparse: sweep
        ] {
            for overlap_type in ["all", "first"] {
                let (chrs, starts, ends) = random_intervals(&mut state, n, 6, span, 300);
                let (chrs2, starts2, ends2) = random_intervals(&mut state, m, 6, span, 300);

                let (mut expected, _) = overlaps(
                    &chrs, &starts, &ends, &chrs2, &starts2, &ends2, 0, overlap_type, true, false,
                );
                radsort::sort(&mut expected);

                let got = overlap_indices(
                    &chrs, &starts, &ends, &chrs2, &starts2, &ends2, 0, overlap_type, true, false,
                );
                assert_eq!(got, expected, "n={n} m={m} span={span} type={overlap_type}");
            }
        }
    }
}
