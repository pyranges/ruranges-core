use rustc_hash::FxHashMap;

use crate::{
    ruranges_structs::{GroupType, PositionType},
    sorts,
};

pub fn sweep_line_complement<G: GroupType, T: PositionType>(
    chrs: &[G],
    starts: &[T],
    ends: &[T],
    slack: T,
    chrom_lens: &FxHashMap<G, T>,
    include_first_interval: bool, // <-- new parameter
) -> (Vec<G>, Vec<T>, Vec<T>, Vec<u32>) {
    let mut out_chrs = Vec::with_capacity(chrs.len());
    let mut out_starts = Vec::with_capacity(chrs.len());
    let mut out_ends = Vec::with_capacity(chrs.len());
    let mut out_idxs = Vec::with_capacity(chrs.len());

    // Early return if no input
    if chrs.is_empty() {
        return (out_chrs, out_starts, out_ends, out_idxs);
    }

    // Build your events array, sorted by chr and pos
    let events = sorts::build_sorted_events_single_collection(chrs, starts, ends, slack);

    // Initialize
    let mut current_chr = events[0].chr;
    let mut active_count = 0_i64;
    // Whether we start "in a hole" (i.e., complement) depends on `include_first_interval`
    let mut in_complement = include_first_interval;
    // Start the first hole at position 0 of the chromosome (only matters if `in_complement == true`)
    let mut current_start = T::zero();
    // Must start at the first event's index, not at 0: row 0 of the caller's
    // input need not belong to the first group in sorted order, in which case a
    // hardcoded 0 attributes the group's gaps to an unrelated row.
    let mut current_index = events[0].idx;

    for e in events {
        // If we hit a new chromosome
        if e.chr != current_chr {
            // If we ended the previous chromosome still in a hole,
            // optionally close it out at the chromosome’s end
            if let Some(chlen) = chrom_lens.get(&current_chr) {
                if in_complement {
                    out_chrs.push(current_chr);
                    out_starts.push(current_start);
                    out_ends.push(*chlen);
                    out_idxs.push(current_index);
                }
            }

            // Reset for new chromosome
            current_chr = e.chr;
            active_count = 0;
            in_complement = include_first_interval;
            current_start = T::zero();
            current_index = e.idx;
        }

        // Process this event
        if e.is_start {
            // coverage X → X + 1
            active_count += 1;
            // If coverage was zero, we just ended a hole
            if active_count == 1 && in_complement && current_start != e.pos {
                // That hole ends at e.pos
                out_chrs.push(current_chr);
                out_starts.push(current_start);
                out_ends.push(e.pos);
                out_idxs.push(current_index);

                // We're no longer in a hole
                in_complement = false;
            }
        } else {
            // coverage X → X - 1
            active_count -= 1;
            // If coverage has just dropped back to zero,
            // we start a new hole here
            if active_count == 0 {
                in_complement = true;
                current_start = e.pos;
            }
        }
    }

    // End of all events: if we finished in a hole and have chromosome lengths
    if let Some(chlen) = chrom_lens.get(&current_chr) {
        if in_complement {
            out_chrs.push(current_chr);
            out_starts.push(current_start);
            out_ends.push(*chlen);
            out_idxs.push(current_index);
        }
    }

    (out_chrs, out_starts, out_ends, out_idxs)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The caller's row 0 need not belong to the group that sorts first, so a
    /// hardcoded starting index attributed that group's gaps to a foreign row.
    #[test]
    fn reported_index_belongs_to_the_group_that_owns_the_gap() {
        // row 0 -> group 1, rows 1 and 2 -> group 0.
        // Group 0 covers [0,1) and [2,3), so its only gap is [1,2).
        let groups = [1_u32, 0, 0];
        let starts = [0_i64, 2, 0];
        let ends = [1_i64, 3, 1];

        let (out_chrs, out_starts, out_ends, out_idxs) =
            sweep_line_complement(&groups, &starts, &ends, 0, &FxHashMap::default(), false);

        assert_eq!(out_starts, vec![1]);
        assert_eq!(out_ends, vec![2]);
        assert_eq!(out_chrs, vec![0]);
        // The gap belongs to group 0, so the reported row must be one of its own.
        assert_eq!(groups[out_idxs[0] as usize], 0);
    }

    #[test]
    fn reported_index_belongs_to_its_group_for_every_group() {
        // Two groups, each with an internal gap, listed in descending id order.
        let groups = [1_u32, 1, 0, 0];
        let starts = [0_i64, 20, 0, 10];
        let ends = [1_i64, 21, 1, 11];

        let (out_chrs, _out_starts, _out_ends, out_idxs) =
            sweep_line_complement(&groups, &starts, &ends, 0, &FxHashMap::default(), false);

        assert_eq!(out_idxs.len(), out_chrs.len());
        for (group, idx) in out_chrs.iter().zip(out_idxs.iter()) {
            assert_eq!(groups[*idx as usize], *group);
        }
    }
}
