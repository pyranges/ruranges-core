use std::str::FromStr;

use crate::{
    overlaps::{collect_overlap_pairs_from_sorted, sorted_records, IntervalRecord},
    ruranges_structs::{GroupType, MinEvent, Nearest, OverlapPair, OverlapType, PositionType},
    sorts::build_sorted_events_single_collection_separate_outputs,
};

/// Convert pre-sorted IntervalRecord slice into MinEvent vec using the
/// `start` field (optionally shifted by `slack`). Both sort keys —
/// `(group, start)` and `(chr, start - slack)` — order elements
/// identically, so we get a sorted MinEvent vec for free, no radsort
/// passes needed.
fn min_events_from_sorted_starts<C: GroupType, T: PositionType>(
    records: &[IntervalRecord<C, T>],
    slack: T,
) -> Vec<MinEvent<C, T>> {
    records
        .iter()
        .map(|r| MinEvent {
            chr: r.group,
            pos: r.start - slack,
            idx: r.idx,
        })
        .collect()
}

/// Whether a query reports every neighbour at a winning distance or only one
/// of them.
///
/// Every neighbour a query overlaps sits at distance 0, so under
/// [`Ties::All`] "the nearest neighbours" of a query inside a dense region is
/// every interval covering it. [`Ties::First`] reports one row per query per
/// distance instead, which is what `bedtools closest -t first`,
/// `GenomicRanges` `select="arbitrary"` and BEDOPS `--closest` do.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum Ties {
    /// Report every neighbour at each reported distance. The default, and
    /// what [`nearest`] has always done.
    #[default]
    All,
    /// Report a single neighbour per reported distance. *Which* one is not
    /// specified — only that it is at the winning distance, and that the same
    /// input gives the same answer every time.
    First,
}

impl FromStr for Ties {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "all" => Ok(Ties::All),
            "first" => Ok(Ties::First),
            _ => Err("Invalid ties string"),
        }
    }
}

/// For each MinEvent in `sorted_ends`, find up to `k` *unique positions*
/// in `sorted_starts2` that lie to the right (including equal position on the
/// same chromosome). Entries sharing a position share a distance, so under
/// [`Ties::All`] they are all reported and count as one unique position, while
/// under [`Ties::First`] only the first of each position is reported.
pub fn nearest_intervals_to_the_right<C: GroupType, T: PositionType>(
    sorted_ends: &[MinEvent<C, T>],
    sorted_starts2: &[MinEvent<C, T>],
    k: usize,
    ties: Ties,
) -> Vec<Nearest<T>> {
    // We might need more than `sorted_ends.len()` because each end could
    // contribute up to `k` *unique positions* (potentially multiplied by the
    // number of intervals sharing those positions). So we set capacity
    // accordingly.
    // This is not strictly required, but it helps performance to reserve enough space.
    let mut output = Vec::with_capacity(sorted_ends.len().saturating_mul(k));

    let n_starts = sorted_starts2.len();

    // `j` will track our position in sorted_starts2 as we move through sorted_ends.
    let mut j = 0usize;

    // Iterate over each 'end' event
    for end in sorted_ends {
        let end_chr = end.chr;
        let end_pos = end.pos;

        // Advance `j` so that sorted_starts2[j] is the first start
        // that is >= end_pos on the same chrom (or beyond).
        // Because both arrays are sorted, we never need to move `j` backward.
        while j < n_starts {
            let start = &sorted_starts2[j];
            if start.chr < end_chr {
                // still on a smaller chromosome; move j forward
                j += 1;
            } else if start.chr == end_chr && start.pos < end_pos {
                // same chrom but still to the left; move j forward
                j += 1;
            } else {
                // now start.chr > end_chr (i.e. next chromosome) OR
                // start.chr == end_chr && start.pos >= end_pos
                // -> we've reached a region that is "to the right" or next chrom
                break;
            }
        }

        // Now collect up to k unique positions (on the same chromosome).
        let mut unique_count = 0;
        let mut last_pos: Option<T> = None;

        // We'll scan from `j` onward, but we do NOT move `j` itself
        // because the next 'end' might need a similar or slightly advanced position.
        // Instead, we use `local_idx` to look ahead for this specific end.
        let mut local_idx = j;
        while local_idx < n_starts {
            let start = &sorted_starts2[local_idx];

            // If we've passed beyond the chromosome of this end, we won't find
            // any more right-side intervals for this end.
            if start.chr != end_chr {
                break;
            }

            // Check if we're at a new unique position
            let is_new_position = last_pos.map_or(true, |lp| start.pos != lp);
            if is_new_position {
                unique_count += 1;
                if unique_count > k {
                    // we've reached the limit of k unique positions
                    break;
                }
                last_pos = Some(start.pos);
            }

            // A repeated position is a repeated distance, so under `First` the
            // bucket already has its one row.
            if is_new_position || ties == Ties::All {
                let distance = start.pos - end_pos + T::one(); // can be 0 or positive
                output.push(Nearest {
                    distance,
                    idx: end.idx,
                    idx2: start.idx,
                });
            }

            local_idx += 1;
        }
    }

    output
}

/// For each MinEvent in `sorted_ends`, find up to `k` *unique positions*
/// in `sorted_starts2` that lie to the left (strictly smaller position on
/// the same chromosome). Entries sharing a position share a distance, so
/// under [`Ties::All`] they are all reported and count as one unique position
/// in the limit `k`, while under [`Ties::First`] only the first of each
/// position is reported.
pub fn nearest_intervals_to_the_left<C: GroupType, T: PositionType>(
    sorted_ends: &[MinEvent<C, T>],
    sorted_starts2: &[MinEvent<C, T>],
    k: usize,
    ties: Ties,
) -> Vec<Nearest<T>> {
    // The max possible size is (number of ends) * (k + duplicates at each of those k positions).
    // We reserve a rough upper bound for efficiency.
    let mut output = Vec::with_capacity(sorted_ends.len().saturating_mul(k));

    let n_starts = sorted_starts2.len();
    let mut j = 0_usize; // Points into sorted_starts2

    for end in sorted_ends {
        let end_chr = end.chr;
        let end_pos = end.pos;

        // Move `j` forward so that:
        // - All start events at indices < j have start.chr < end_chr
        //   OR (start.chr == end_chr && start.pos <= end_pos).
        // - Equivalently, sorted_starts2[j] is the *first* event that is NOT
        //   strictly to the left of `end`.
        while j < n_starts {
            let start = &sorted_starts2[j];
            if start.chr < end_chr {
                // still a smaller chromosome => definitely to the left
                j += 1;
            } else if start.chr == end_chr && start.pos <= end_pos {
                // same chrom, smaller/equal position => to the left (touching counts)
                j += 1;
            } else {
                // we've reached a start that is not to the left
                break;
            }
        }

        // Now, everything in [0..j) is strictly to the left of `end`.
        // We'll look backwards from j-1 to gather up to k unique positions
        // on the same chromosome.
        if j == 0 {
            // No intervals to the left; skip
            continue;
        }

        let mut local_idx = j - 1;
        let mut unique_count = 0;
        let mut last_pos: Option<T> = None;

        // Descend from j-1 down to 0 (or until we break).
        loop {
            let start = &sorted_starts2[local_idx];

            // Must match the same chromosome
            if start.chr != end_chr {
                break;
            }

            // Check if we have a new (unique) position
            let is_new_position = last_pos.map_or(true, |lp| start.pos != lp);
            if is_new_position {
                unique_count += 1;
                if unique_count > k {
                    break;
                }
                last_pos = Some(start.pos);
            }

            // A repeated position is a repeated distance, so under `First` the
            // bucket already has its one row.
            if is_new_position || ties == Ties::All {
                // Calculate the distance (end.pos - start.pos)
                // Here, start.pos < end.pos by definition if we get here.
                let distance = end_pos - start.pos + T::one();
                output.push(Nearest {
                    distance,
                    idx: end.idx,    // the 'end' event's idx
                    idx2: start.idx, // the 'start' event's idx
                });
            }

            if local_idx == 0 {
                break;
            }
            local_idx -= 1;
        }
    }

    output
}

/// Merges th
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Direction {
    Forward,
    Backward,
    Any,
}

impl FromStr for Direction {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "forward" => Ok(Direction::Forward),
            "backward" => Ok(Direction::Backward),
            "any" => Ok(Direction::Any),
            _ => Err("Invalid direction string"),
        }
    }
}

/// Report the `k` nearest distances of `chrs2/starts2/ends2` for every row of
/// `chrs/starts/ends`, with every neighbour at each of those distances.
///
/// Use [`nearest_with_ties`] to report one neighbour per distance instead.
#[allow(clippy::too_many_arguments)]
pub fn nearest<C: GroupType, T: PositionType>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
    k: usize,
    include_overlaps: bool,
    direction: &str,
    sort_output: bool,
) -> (Vec<u32>, Vec<u32>, Vec<T>)
where
    C: Send + Sync,
    T: Send + Sync,
{
    nearest_with_ties(
        chrs,
        starts,
        ends,
        chrs2,
        starts2,
        ends2,
        slack,
        k,
        include_overlaps,
        direction,
        sort_output,
        Ties::All,
    )
}

/// [`nearest`], with control over how many neighbours a tied distance reports.
///
/// `Ties::First` is not a filter over the `Ties::All` result: the tied rows are
/// never produced. Under `include_overlaps` a query inside a dense region
/// overlaps — and so ties with — every interval covering it, and materialising
/// those before dropping them costs the same memory as keeping them.
#[allow(clippy::too_many_arguments)]
pub fn nearest_with_ties<C: GroupType, T: PositionType>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
    k: usize,
    include_overlaps: bool,
    direction: &str,
    sort_output: bool,
    ties: Ties,
) -> (Vec<u32>, Vec<u32>, Vec<T>)
where
    C: Send + Sync,
    T: Send + Sync,
{
    let dir = Direction::from_str(direction).unwrap();

    let need_left = dir == Direction::Backward || dir == Direction::Any;
    let need_right = dir == Direction::Forward || dir == Direction::Any;

    // Which sorts are needed by which producer:
    //
    //   producer            left order       right order
    //   ------------------  ---------------  ---------------
    //   overlap             (group, start)   (group, start)
    //   nearest_left        (group, start)   (group, end)
    //   nearest_right       (group, end)     (group, start)
    //
    // The `(group, start)` ordering shows up in both overlap+nearest_left
    // (for left) and overlap+nearest_right (for right). We build each
    // shared `IntervalRecord` view at most once and convert it into the
    // `MinEvent` shape that the nearest sweeps want via a cheap O(n) walk
    // (no extra radsort passes, since subtracting a constant `slack`
    // preserves the sort order). The end-based sorts are unique per
    // producer and have to be built independently.
    let want_left_records = include_overlaps || need_left;
    let want_right_records = include_overlaps || need_right;

    // Phase 1: build all sorted views in parallel. The four leaves can run
    // on up to four worker threads.
    let ((left_records, right_records), (left_by_end, right_by_end)) = rayon::join(
        || {
            rayon::join(
                || {
                    if want_left_records {
                        Some(sorted_records(chrs, starts, ends, OverlapType::All))
                    } else {
                        None
                    }
                },
                || {
                    if want_right_records {
                        Some(sorted_records(chrs2, starts2, ends2, OverlapType::All))
                    } else {
                        None
                    }
                },
            )
        },
        || {
            rayon::join(
                || {
                    if need_right {
                        Some(build_sorted_events_single_collection_separate_outputs(
                            chrs, ends, slack,
                        ))
                    } else {
                        None
                    }
                },
                || {
                    if need_left {
                        Some(build_sorted_events_single_collection_separate_outputs(
                            chrs2, ends2, T::zero(),
                        ))
                    } else {
                        None
                    }
                },
            )
        },
    );

    // Phase 2: run the three sweeps in parallel, each reading from the
    // shared views.
    let compute_overlaps = || -> Vec<OverlapPair> {
        if include_overlaps {
            let l = left_records
                .as_ref()
                .expect("want_left_records implied by include_overlaps");
            let r = right_records
                .as_ref()
                .expect("want_right_records implied by include_overlaps");
            // Every overlap of a query is at distance 0, so they are one
            // bucket, and under `First` the sweep must stop at the first of
            // them rather than hand the merge a bucket to throw away. That is
            // exactly `OverlapType::First`.
            //
            // The records were sorted with the `All` key, `(group, start)`.
            // `sorted_records` documents the richer `First` key as what makes
            // *which* target the sweep returns reproducible; radsort is stable
            // and the records are built in `idx` order, so equal
            // `(group, start)` targets still reach the active list in `idx`
            // order and the pick is reproducible here too. Nearest does not
            // promise which tied neighbour it returns anyway, and the two
            // extra radsort passes over 32-byte records would be paid on
            // exactly the large inputs this option exists for.
            let overlap_type = match ties {
                Ties::All => OverlapType::All,
                Ties::First => OverlapType::First,
            };
            // Reproduce `overlaps(... sort_output=true, false)` but
            // skip the redundant sort_records calls inside.
            let mut pairs = collect_overlap_pairs_from_sorted(l, r, slack, overlap_type, false);
            radsort::sort_by_key(&mut pairs, |p| (p.idx, p.idx2));
            pairs
        } else {
            Vec::new()
        }
    };

    let compute_left = || -> Vec<Nearest<T>> {
        if need_left {
            let l = left_records
                .as_ref()
                .expect("want_left_records implied by need_left");
            // Convert pre-sorted left records to MinEvent with pos = start - slack.
            // Sort order is preserved (subtracting a constant doesn't reorder).
            let sorted_starts = min_events_from_sorted_starts(l, slack);
            let sorted_ends2 = right_by_end
                .as_ref()
                .expect("right_by_end set when need_left is true");
            let mut tmp = nearest_intervals_to_the_left(&sorted_starts, sorted_ends2, k, ties);
            // Each `idx` is unique per row and produces one contiguous block
            // whose distances are already non-decreasing (descending local_idx,
            // growing `end_pos - start.pos + 1`). radsort is a stable LSD radix
            // sort, so sorting by `n.idx` alone preserves the within-block
            // distance order.
            radsort::sort_by_key(&mut tmp, |n| n.idx);
            tmp
        } else {
            Vec::new()
        }
    };

    let compute_right = || -> Vec<Nearest<T>> {
        if need_right {
            let r = right_records
                .as_ref()
                .expect("want_right_records implied by need_right");
            // Right side: pos = start - 0, just reuse the shared records.
            let sorted_starts2 = min_events_from_sorted_starts(r, T::zero());
            let sorted_ends = left_by_end
                .as_ref()
                .expect("left_by_end set when need_right is true");
            let mut tmp = nearest_intervals_to_the_right(sorted_ends, &sorted_starts2, k, ties);
            // See comment above — stable sort by `n.idx` is sufficient.
            radsort::sort_by_key(&mut tmp, |n| n.idx);
            tmp
        } else {
            Vec::new()
        }
    };

    // Run the three independent producers in parallel via nested rayon::join.
    // rayon::join is cheap on small inputs: the calling thread runs one closure
    // inline and only yields the other for stealing if a worker is idle.
    let ((overlaps, nearest_left), nearest_right) =
        rayon::join(|| rayon::join(compute_overlaps, compute_left), compute_right);

    merge_three_way_by_index_distance(
        &overlaps,
        &nearest_left,
        &nearest_right,
        k,
        ties,
        sort_output,
    )
}

/// Merges three sources of intervals, grouped by `idx`.
/// For each unique `idx`, returns up to `k` *distinct* distances — with every
/// interval at those distances under [`Ties::All`], and one interval per
/// distance under [`Ties::First`]. Overlaps are treated as distance=0 and form
/// a single bucket per idx.
///
/// All inputs are sorted by `idx` ascending; within each idx group, distances
/// in `nearest_left` and `nearest_right` are non-decreasing.
#[allow(clippy::too_many_arguments)]
pub fn merge_three_way_by_index_distance<T: PositionType>(
    overlaps: &[OverlapPair],
    nearest_left: &[Nearest<T>],
    nearest_right: &[Nearest<T>],
    k: usize,
    ties: Ties,
    sort_output: bool,
) -> (Vec<u32>, Vec<u32>, Vec<T>) {
    // Cap pre-allocation: at most one entry emitted per input row across all
    // three sources. `saturating_mul(k)` over-allocated by ~k× for large k.
    let cap = overlaps.len() + nearest_left.len() + nearest_right.len();

    if sort_output {
        // On the sort path, collect into a single Vec<Nearest> so we can sort
        // by (idx, distance, idx2) to match the prior output ordering exactly.
        let mut results: Vec<Nearest<T>> = Vec::with_capacity(cap);
        merge_inner(
            overlaps,
            nearest_left,
            nearest_right,
            k,
            ties,
            |idx, idx2, distance| {
                results.push(Nearest {
                    idx,
                    idx2,
                    distance,
                });
            },
        );
        // The merge already emits in idx-ascending order with non-decreasing
        // distance inside each idx bucket; only idx2 can be out of order inside
        // a (idx, distance) bucket. That's near-fully-sorted input, which
        // pdqsort handles in close to O(n) thanks to its pattern detection —
        // typically faster than a 3-pass LSD radix sort over 16-byte structs
        // here. Unstable is fine: full tuple key disambiguates every element.
        results.sort_unstable_by_key(|n| (n.idx, n.distance, n.idx2));
        let mut out_idxs = Vec::with_capacity(results.len());
        let mut out_idxs2 = Vec::with_capacity(results.len());
        let mut out_distances = Vec::with_capacity(results.len());
        for rec in results {
            out_idxs.push(rec.idx);
            out_idxs2.push(rec.idx2);
            out_distances.push(rec.distance);
        }
        return (out_idxs, out_idxs2, out_distances);
    }

    let mut out_idxs: Vec<u32> = Vec::with_capacity(cap);
    let mut out_idxs2: Vec<u32> = Vec::with_capacity(cap);
    let mut out_distances: Vec<T> = Vec::with_capacity(cap);
    merge_inner(
        overlaps,
        nearest_left,
        nearest_right,
        k,
        ties,
        |idx, idx2, distance| {
            out_idxs.push(idx);
            out_idxs2.push(idx2);
            out_distances.push(distance);
        },
    );
    (out_idxs, out_idxs2, out_distances)
}

/// Walks the three sources, calling `emit(idx, idx2, distance)` for each kept
/// entry. Tracks the previous distance with `Option<T>` instead of a per-idx
/// `HashSet`, eliminating the per-idx allocation. Overlaps are emitted first
/// (all at distance 0) as a single bucket; remaining buckets come from a
/// 2-way merge of left/right by distance.
fn merge_inner<T: PositionType, F: FnMut(u32, u32, T)>(
    overlaps: &[OverlapPair],
    nearest_left: &[Nearest<T>],
    nearest_right: &[Nearest<T>],
    k: usize,
    ties: Ties,
    mut emit: F,
) {
    let (mut i, mut j, mut r) = (0_usize, 0_usize, 0_usize);

    while i < overlaps.len() || j < nearest_left.len() || r < nearest_right.len() {
        let idx_o = overlaps.get(i).map(|o| o.idx);
        let idx_l = nearest_left.get(j).map(|n| n.idx);
        let idx_r = nearest_right.get(r).map(|n| n.idx);

        let current_idx = match (idx_o, idx_l, idx_r) {
            (None, None, None) => break,
            (Some(a), Some(b), Some(c)) => a.min(b).min(c),
            (Some(a), Some(b), None) => a.min(b),
            (Some(a), None, Some(c)) => a.min(c),
            (None, Some(b), Some(c)) => b.min(c),
            (Some(a), None, None) => a,
            (None, Some(b), None) => b,
            (None, None, Some(c)) => c,
        };

        let i_start = i;
        while i < overlaps.len() && overlaps[i].idx == current_idx {
            i += 1;
        }
        let overlaps_slice = &overlaps[i_start..i];

        let j_start = j;
        while j < nearest_left.len() && nearest_left[j].idx == current_idx {
            j += 1;
        }
        let left_slice = &nearest_left[j_start..j];

        let r_start = r;
        while r < nearest_right.len() && nearest_right[r].idx == current_idx {
            r += 1;
        }
        let right_slice = &nearest_right[r_start..r];

        let mut distinct_count = 0_usize;

        // Overlaps (all distance 0) form one bucket.
        if !overlaps_slice.is_empty() {
            distinct_count = 1;
            match ties {
                Ties::All => {
                    for op in overlaps_slice {
                        emit(op.idx, op.idx2, T::zero());
                    }
                }
                Ties::First => {
                    let op = &overlaps_slice[0];
                    emit(op.idx, op.idx2, T::zero());
                }
            }
            if distinct_count >= k {
                continue;
            }
        }

        // 2-way merge of left/right by distance. `last_dist` tracks the most
        // recent bucket's distance — bucket boundary = distance change.
        let (mut lj, mut rr) = (0_usize, 0_usize);
        let mut last_dist: Option<T> = None;

        loop {
            let dl = left_slice.get(lj).map(|n| n.distance);
            let dr = right_slice.get(rr).map(|n| n.distance);

            let smallest = match (dl, dr) {
                (None, None) => break,
                (Some(a), Some(b)) => {
                    if a <= b {
                        a
                    } else {
                        b
                    }
                }
                (Some(a), None) => a,
                (None, Some(b)) => b,
            };

            if last_dist != Some(smallest) {
                distinct_count += 1;
                if distinct_count > k {
                    break;
                }
                last_dist = Some(smallest);
            }

            if ties == Ties::First {
                // A query five bases from a neighbour on each side has both in
                // this bucket. Take the left one; the choice is arbitrary, and
                // fixing it here is what keeps the answer reproducible.
                let winner = if left_slice.get(lj).map(|n| n.distance) == Some(smallest) {
                    left_slice[lj]
                } else {
                    right_slice[rr]
                };
                emit(winner.idx, winner.idx2, winner.distance);
            }

            while lj < left_slice.len() && left_slice[lj].distance == smallest {
                if ties == Ties::All {
                    let n = left_slice[lj];
                    emit(n.idx, n.idx2, n.distance);
                }
                lj += 1;
            }
            while rr < right_slice.len() && right_slice[rr].distance == smallest {
                if ties == Ties::All {
                    let n = right_slice[rr];
                    emit(n.idx, n.idx2, n.distance);
                }
                rr += 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::{HashMap, HashSet};

    use super::{nearest, nearest_with_ties, Ties};

    /// Deterministic intervals, so a failure is reproducible without a
    /// property-testing dependency.
    fn pseudo_random_intervals(
        n: usize,
        seed: u64,
        groups: u32,
        span: i64,
        max_len: i64,
    ) -> (Vec<u32>, Vec<i64>, Vec<i64>) {
        let mut state = seed;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 33) as i64
        };

        let mut chrs = Vec::with_capacity(n);
        let mut starts = Vec::with_capacity(n);
        let mut ends = Vec::with_capacity(n);
        for _ in 0..n {
            let start = next() % span;
            let len = 1 + next() % max_len;
            chrs.push((next() % i64::from(groups)) as u32);
            starts.push(start);
            ends.push(start + len);
        }
        (chrs, starts, ends)
    }

    /// The one assertion that catches nearly everything: `First` answers
    /// exactly the same queries as `All`, at the same distance, with exactly
    /// one row each.
    #[test]
    fn nearest_first_answers_the_same_queries_as_all_with_one_row_each() {
        // Dense enough that most queries overlap several targets, which is the
        // case the option exists for.
        let (chrs, starts, ends) = pseudo_random_intervals(500, 11, 3, 4_000, 200);
        let (chrs2, starts2, ends2) = pseudo_random_intervals(500, 29, 3, 4_000, 200);

        for include_overlaps in [true, false] {
            for direction in ["any", "forward", "backward"] {
                let (all_idx, _all_idx2, all_dist) = nearest_with_ties(
                    &chrs,
                    &starts,
                    &ends,
                    &chrs2,
                    &starts2,
                    &ends2,
                    0,
                    1,
                    include_overlaps,
                    direction,
                    true,
                    Ties::All,
                );
                let (first_idx, _first_idx2, first_dist) = nearest_with_ties(
                    &chrs,
                    &starts,
                    &ends,
                    &chrs2,
                    &starts2,
                    &ends2,
                    0,
                    1,
                    include_overlaps,
                    direction,
                    true,
                    Ties::First,
                );

                let case = format!("include_overlaps={include_overlaps}, direction={direction}");

                // k=1 leaves one bucket, so every row of a query shares its
                // distance and the map is well defined.
                let winning: HashMap<u32, i64> = all_idx
                    .iter()
                    .copied()
                    .zip(all_dist.iter().copied())
                    .collect();
                assert!(
                    !winning.is_empty(),
                    "{case}: the fixture answers no queries at all"
                );

                let answered: HashSet<u32> = first_idx.iter().copied().collect();
                assert_eq!(
                    answered.len(),
                    first_idx.len(),
                    "{case}: First reported a query more than once"
                );
                assert_eq!(
                    answered,
                    winning.keys().copied().collect::<HashSet<_>>(),
                    "{case}: First and All answer different queries"
                );
                for (idx, dist) in first_idx.iter().zip(first_dist.iter()) {
                    assert_eq!(
                        winning[idx], *dist,
                        "{case}: query {idx} came back at a losing distance"
                    );
                }
                assert!(
                    first_idx.len() < all_idx.len(),
                    "{case}: the fixture has no ties, so it proves nothing"
                );
            }
        }
    }

    #[test]
    fn nearest_first_reports_one_of_the_overlapping_intervals() {
        // One query covering three targets: three rows at distance 0.
        let chrs = vec![1_u32];
        let starts = vec![100_i64];
        let ends = vec![200_i64];

        let chrs2 = vec![1_u32; 3];
        let starts2 = vec![90_i64, 120, 150];
        let ends2 = vec![110_i64, 130, 160];

        let args = (&chrs, &starts, &ends, &chrs2, &starts2, &ends2);
        let (all_idx, _, all_dist) = nearest_with_ties(
            args.0,
            args.1,
            args.2,
            args.3,
            args.4,
            args.5,
            0,
            1,
            true,
            "any",
            true,
            Ties::All,
        );
        assert_eq!(all_idx, vec![0, 0, 0]);
        assert_eq!(all_dist, vec![0, 0, 0]);

        let (first_idx, first_idx2, first_dist) = nearest_with_ties(
            args.0,
            args.1,
            args.2,
            args.3,
            args.4,
            args.5,
            0,
            1,
            true,
            "any",
            true,
            Ties::First,
        );
        assert_eq!(first_idx, vec![0]);
        assert_eq!(first_dist, vec![0]);
        assert_eq!(first_idx2.len(), 1);
        assert!(first_idx2[0] < 3);
    }

    #[test]
    fn nearest_first_breaks_a_two_sided_tie() {
        // 85-95 and 115-125 are the same distance from 100-110, one on each
        // side, so the tie is between the two sweeps rather than inside one.
        let chrs = vec![1_u32];
        let starts = vec![100_i64];
        let ends = vec![110_i64];

        let chrs2 = vec![1_u32; 2];
        let starts2 = vec![85_i64, 115];
        let ends2 = vec![95_i64, 125];

        let (all_idx, _, all_dist) = nearest_with_ties(
            &chrs,
            &starts,
            &ends,
            &chrs2,
            &starts2,
            &ends2,
            0,
            1,
            true,
            "any",
            true,
            Ties::All,
        );
        assert_eq!(all_idx, vec![0, 0]);
        assert_eq!(all_dist, vec![6, 6]);

        let (first_idx, _, first_dist) = nearest_with_ties(
            &chrs,
            &starts,
            &ends,
            &chrs2,
            &starts2,
            &ends2,
            0,
            1,
            true,
            "any",
            true,
            Ties::First,
        );
        assert_eq!(first_idx, vec![0]);
        assert_eq!(first_dist, vec![6]);
    }

    #[test]
    fn nearest_first_prefers_an_overlap_to_a_touching_neighbour() {
        // 95-105 overlaps 100-110 (distance 0); 110-120 merely touches it
        // (distance 1). Picking one row must not pick the wrong bucket.
        let chrs = vec![1_u32];
        let starts = vec![100_i64];
        let ends = vec![110_i64];

        let chrs2 = vec![1_u32; 2];
        let starts2 = vec![95_i64, 110];
        let ends2 = vec![105_i64, 120];

        let (idx, idx2, dist) = nearest_with_ties(
            &chrs,
            &starts,
            &ends,
            &chrs2,
            &starts2,
            &ends2,
            0,
            1,
            true,
            "any",
            true,
            Ties::First,
        );
        assert_eq!(idx, vec![0]);
        assert_eq!(idx2, vec![0]);
        assert_eq!(dist, vec![0]);
    }

    #[test]
    fn nearest_first_reports_one_row_per_distance_when_k_is_two() {
        // Two targets tied at each of two distances. k=2 asks for two
        // distances, so First reports two rows, not four.
        let chrs = vec![1_u32];
        let starts = vec![100_i64];
        let ends = vec![110_i64];

        let chrs2 = vec![1_u32; 4];
        let starts2 = vec![85_i64, 115, 75, 125];
        let ends2 = vec![95_i64, 125, 85, 135];

        let (all_idx, _, all_dist) = nearest_with_ties(
            &chrs,
            &starts,
            &ends,
            &chrs2,
            &starts2,
            &ends2,
            0,
            2,
            true,
            "any",
            true,
            Ties::All,
        );
        assert_eq!(all_idx.len(), 4);
        assert_eq!(all_dist, vec![6, 6, 16, 16]);

        let (first_idx, _, first_dist) = nearest_with_ties(
            &chrs,
            &starts,
            &ends,
            &chrs2,
            &starts2,
            &ends2,
            0,
            2,
            true,
            "any",
            true,
            Ties::First,
        );
        assert_eq!(first_idx, vec![0, 0]);
        assert_eq!(first_dist, vec![6, 16]);
    }

    #[test]
    fn nearest_defaults_to_reporting_every_tied_interval() {
        let (chrs, starts, ends) = pseudo_random_intervals(200, 3, 2, 2_000, 150);
        let (chrs2, starts2, ends2) = pseudo_random_intervals(200, 7, 2, 2_000, 150);

        let bare = nearest(
            &chrs, &starts, &ends, &chrs2, &starts2, &ends2, 0, 1, true, "any", true,
        );
        let explicit = nearest_with_ties(
            &chrs,
            &starts,
            &ends,
            &chrs2,
            &starts2,
            &ends2,
            0,
            1,
            true,
            "any",
            true,
            Ties::All,
        );
        assert_eq!(bare, explicit);
    }

    #[test]
    fn nearest_backward_includes_touching_left_interval_with_distance_one() {
        let chrs = vec![1_u32];
        let starts = vec![5_i64];
        let ends = vec![10_i64];

        let chrs2 = vec![1_u32];
        let starts2 = vec![1_i64];
        let ends2 = vec![5_i64];

        let (idx, idx2, dist) = nearest(
            &chrs, &starts, &ends, &chrs2, &starts2, &ends2, 0, 1, true, "backward", true,
        );

        assert_eq!(idx, vec![0]);
        assert_eq!(idx2, vec![0]);
        assert_eq!(dist, vec![1]);
    }

    #[test]
    fn nearest_forward_keeps_touching_right_interval_with_distance_one() {
        let chrs = vec![1_u32];
        let starts = vec![1_i64];
        let ends = vec![5_i64];

        let chrs2 = vec![1_u32];
        let starts2 = vec![5_i64];
        let ends2 = vec![10_i64];

        let (idx, idx2, dist) = nearest(
            &chrs, &starts, &ends, &chrs2, &starts2, &ends2, 0, 1, true, "forward", true,
        );

        assert_eq!(idx, vec![0]);
        assert_eq!(idx2, vec![0]);
        assert_eq!(dist, vec![1]);
    }
}
