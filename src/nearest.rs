use std::str::FromStr;

use crate::{
    overlaps::overlaps,
    ruranges_structs::{GroupType, MinEvent, Nearest, OverlapPair, PositionType},
    sorts::build_sorted_events_single_collection_separate_outputs,
};

/// For each MinEvent in `sorted_ends`, find up to `k` *unique positions*
/// in `sorted_starts2` that lie to the right (including equal position on the
/// same chromosome). If multiple entries in `sorted_starts2` share the same
/// position, they all get reported, but they count as one unique position.
pub fn nearest_intervals_to_the_right<C: GroupType, T: PositionType>(
    sorted_ends: Vec<MinEvent<C, T>>,
    sorted_starts2: Vec<MinEvent<C, T>>,
    k: usize,
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
    for end in &sorted_ends {
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
            if last_pos.map_or(true, |lp| start.pos != lp) {
                unique_count += 1;
                if unique_count > k {
                    // we've reached the limit of k unique positions
                    break;
                }
                last_pos = Some(start.pos);
            }

            // This start is included in the results
            let distance = start.pos - end_pos + T::one(); // can be 0 or positive
            output.push(Nearest {
                distance,
                idx: end.idx,
                idx2: start.idx,
            });

            local_idx += 1;
        }
    }

    output
}

/// For each MinEvent in `sorted_ends`, find up to `k` *unique positions*
/// in `sorted_starts2` that lie to the left (strictly smaller position on
/// the same chromosome). If multiple entries in `sorted_starts2` share
/// the same position, they all get reported, but they count as one
/// unique position in the limit `k`.
pub fn nearest_intervals_to_the_left<C: GroupType, T: PositionType>(
    sorted_ends: Vec<MinEvent<C, T>>,
    sorted_starts2: Vec<MinEvent<C, T>>,
    k: usize,
) -> Vec<Nearest<T>> {
    // The max possible size is (number of ends) * (k + duplicates at each of those k positions).
    // We reserve a rough upper bound for efficiency.
    let mut output = Vec::with_capacity(sorted_ends.len().saturating_mul(k));

    let n_starts = sorted_starts2.len();
    let mut j = 0_usize; // Points into sorted_starts2

    for end in &sorted_ends {
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
            if last_pos.map_or(true, |lp| start.pos != lp) {
                unique_count += 1;
                if unique_count > k {
                    break;
                }
                last_pos = Some(start.pos);
            }

            // Calculate the distance (end.pos - start.pos)
            // Here, start.pos < end.pos by definition if we get here.
            let distance = end_pos - start.pos + T::one();
            output.push(Nearest {
                distance,
                idx: end.idx,    // the 'end' event's idx
                idx2: start.idx, // the 'start' event's idx
            });

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
) -> (Vec<u32>, Vec<u32>, Vec<T>) {
    let dir = Direction::from_str(direction).unwrap();

    let need_left = dir == Direction::Backward || dir == Direction::Any;
    let need_right = dir == Direction::Forward || dir == Direction::Any;

    let overlaps = if include_overlaps {
        let (idx, idx2) = overlaps(
            chrs, starts, ends, chrs2, starts2, ends2, slack, "all", true, false,
        );
        idx.into_iter()
            .zip(idx2)
            .map(|(idx, idx2)| OverlapPair { idx, idx2 })
            .collect()
    } else {
        Vec::new()
    };

    let nearest_left = if need_left {
        let sorted_starts =
            build_sorted_events_single_collection_separate_outputs(chrs, starts, slack);
        let sorted_ends2 =
            build_sorted_events_single_collection_separate_outputs(chrs2, ends2, T::zero());
        let mut tmp = nearest_intervals_to_the_left(sorted_starts, sorted_ends2, k);
        // Each `idx` is unique per row and produces one contiguous block whose
        // distances are already non-decreasing (descending local_idx, growing
        // `end_pos - start.pos + 1`). radsort is a stable LSD radix sort, so
        // sorting by `n.idx` alone preserves the within-block distance order.
        radsort::sort_by_key(&mut tmp, |n| n.idx);
        tmp
    } else {
        Vec::new()
    };

    let nearest_right = if need_right {
        let sorted_ends =
            build_sorted_events_single_collection_separate_outputs(chrs, ends, slack);
        let sorted_starts2 =
            build_sorted_events_single_collection_separate_outputs(chrs2, starts2, T::zero());
        let mut tmp = nearest_intervals_to_the_right(sorted_ends, sorted_starts2, k);
        // See comment above — stable sort by `n.idx` is sufficient.
        radsort::sort_by_key(&mut tmp, |n| n.idx);
        tmp
    } else {
        Vec::new()
    };

    merge_three_way_by_index_distance(&overlaps, &nearest_left, &nearest_right, k, sort_output)
}

/// Merges three sources of intervals, grouped by `idx`.
/// For each unique `idx`, returns up to `k` *distinct* distances (including all
/// intervals at those distances). Overlaps are treated as distance=0 and form a
/// single bucket per idx.
///
/// All inputs are sorted by `idx` ascending; within each idx group, distances
/// in `nearest_left` and `nearest_right` are non-decreasing.
pub fn merge_three_way_by_index_distance<T: PositionType>(
    overlaps: &[OverlapPair],
    nearest_left: &[Nearest<T>],
    nearest_right: &[Nearest<T>],
    k: usize,
    sort_output: bool,
) -> (Vec<u32>, Vec<u32>, Vec<T>) {
    // Cap pre-allocation: at most one entry emitted per input row across all
    // three sources. `saturating_mul(k)` over-allocated by ~k× for large k.
    let cap = overlaps.len() + nearest_left.len() + nearest_right.len();

    if sort_output {
        // On the sort path, collect into a single Vec<Nearest> so we can sort
        // by (idx, distance, idx2) to match the prior output ordering exactly.
        let mut results: Vec<Nearest<T>> = Vec::with_capacity(cap);
        merge_inner(overlaps, nearest_left, nearest_right, k, |idx, idx2, distance| {
            results.push(Nearest { idx, idx2, distance });
        });
        // Output is already grouped by idx with non-decreasing distance from the
        // merge; this sort only reorders within (idx, distance) buckets.
        radsort::sort_by_key(&mut results, |n| (n.idx, n.distance, n.idx2));
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
    merge_inner(overlaps, nearest_left, nearest_right, k, |idx, idx2, distance| {
        out_idxs.push(idx);
        out_idxs2.push(idx2);
        out_distances.push(distance);
    });
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
            for op in overlaps_slice {
                emit(op.idx, op.idx2, T::zero());
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

            while lj < left_slice.len() && left_slice[lj].distance == smallest {
                let n = left_slice[lj];
                emit(n.idx, n.idx2, n.distance);
                lj += 1;
            }
            while rr < right_slice.len() && right_slice[rr].distance == smallest {
                let n = right_slice[rr];
                emit(n.idx, n.idx2, n.distance);
                rr += 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::nearest;

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
