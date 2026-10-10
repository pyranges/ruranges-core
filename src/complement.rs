use crate::{
    ruranges_structs::{GroupType, PositionType},
    sorts,
};

use rustc_hash::FxHashSet;

pub fn sweep_line_non_overlaps<G: GroupType, T: PositionType>(
    chrs: &[G],
    starts: &[T],
    ends: &[T],
    chrs2: &[G],
    starts2: &[T],
    ends2: &[T],
    slack: T,
    sort_output: bool,
) -> Vec<u32> {
    let mut no_overlaps = Vec::new();

    // If either set is empty, none can overlap; return everything as “non-overlapping”.
    if chrs.is_empty() || chrs2.is_empty() {
        // Just return all indices as non-overlapping
        return no_overlaps.to_vec();
    }

    // Build up the event list in ascending order (same as before)
    let events = sorts::build_sorted_events_idxs(chrs, starts, ends, chrs2, starts2, ends2, slack);

    let mut overlapped = FxHashSet::default();

    // Active sets
    let mut active1 = FxHashSet::default();
    let mut active2 = FxHashSet::default();

    // Assume the first event determines the “current” chr
    let mut current_chr = events.first().unwrap().chr;

    for e in events {
        // If chromosome changed, clear active sets
        if e.chr != current_chr {
            active1.clear();
            active2.clear();
            current_chr = e.chr;
        }

        if e.is_start {
            // Interval is starting
            if e.first_set {
                // Overlaps with all currently active intervals in set2
                if !active2.is_empty() {
                    overlapped.insert(e.idx);
                }
                // Insert into active1
                active1.insert(e.idx);
            } else {
                // Overlaps with all currently active intervals in set1
                for &idx1 in active1.iter() {
                    overlapped.insert(idx1);
                }
                // Insert into active2
                active2.insert(e.idx);
            }
        } else {
            // Interval is ending
            if e.first_set {
                active1.remove(&e.idx);
                if !overlapped.contains(&e.idx) {
                    no_overlaps.push(e.idx);
                }
                overlapped.remove(&e.idx);
            } else {
                // `overlapped` holds set-1 indices; this one is from set 2.
                active2.remove(&e.idx);
            }
        }
    }

    if sort_output {
        radsort::sort(&mut no_overlaps);
    }
    no_overlaps
}

#[cfg(test)]
mod tests {
    use super::sweep_line_non_overlaps;

    /// Deterministic SplitMix64, so these tests need no dev-dependency.
    fn next(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    #[test]
    fn a_set2_end_does_not_clear_the_set1_row_sharing_its_index() {
        // [0,10) is row 0 of set 1 and overlaps [2,4), row 0 of set 2.
        let free = sweep_line_non_overlaps(&[0u32], &[0i64], &[10], &[0u32], &[2], &[4], 0, true);
        assert!(free.is_empty());
    }

    #[test]
    fn non_overlapping_rows_match_brute_force() {
        let mut state = 0xC0_4E7_u64;
        for trial in 0..2000 {
            let n1 = 1 + (next(&mut state) % 9) as usize;
            let n2 = 1 + (next(&mut state) % 9) as usize;
            let mut make = |n: usize| -> (Vec<u32>, Vec<i64>, Vec<i64>) {
                let c: Vec<u32> = (0..n).map(|_| (next(&mut state) % 2) as u32).collect();
                let s: Vec<i64> = (0..n).map(|_| (next(&mut state) % 40) as i64).collect();
                let e: Vec<i64> = s
                    .iter()
                    .map(|&x| x + 1 + (next(&mut state) % 9) as i64)
                    .collect();
                (c, s, e)
            };
            let (c1, s1, e1) = make(n1);
            let (c2, s2, e2) = make(n2);
            let want: Vec<u32> = (0..n1)
                .filter(|&i| !(0..n2).any(|j| c1[i] == c2[j] && s2[j] < e1[i] && s1[i] < e2[j]))
                .map(|i| i as u32)
                .collect();
            let got = sweep_line_non_overlaps(&c1, &s1, &e1, &c2, &s2, &e2, 0, true);
            assert_eq!(got, want, "trial={trial}");
        }
    }
}
