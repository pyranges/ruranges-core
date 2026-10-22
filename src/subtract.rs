use radsort::sort_by_key;

use crate::ruranges_structs::{GroupType, Interval, PositionType};

/// A right-side interval; subtraction never needs to know which rows formed it.
struct Mask<G: GroupType, T: PositionType> {
    group: G,
    start: T,
    end: T,
}

fn is_sorted_by_group_start<G: GroupType, T: PositionType>(groups: &[G], starts: &[T]) -> bool {
    (1..groups.len()).all(|i| (groups[i - 1], starts[i - 1]) <= (groups[i], starts[i]))
}

/// Sort right intervals by (group, start) and merge overlapping or touching ones
/// within each group. Afterwards both starts and ends are strictly increasing per group.
fn merge_masks<G: GroupType, T: PositionType>(
    groups: &[G],
    starts: &[T],
    ends: &[T],
) -> Vec<Mask<G, T>> {
    let mut masks: Vec<Mask<G, T>> = (0..groups.len())
        .filter(|&i| starts[i] < ends[i])
        .map(|i| Mask {
            group: groups[i],
            start: starts[i],
            end: ends[i],
        })
        .collect();
    if !is_sorted_by_group_start(groups, starts) {
        sort_by_key(&mut masks, |m| m.start);
        sort_by_key(&mut masks, |m| m.group);
    }

    let mut n_merged = 0;
    for i in 0..masks.len() {
        let (group, start, end) = (masks[i].group, masks[i].start, masks[i].end);
        if n_merged > 0 {
            let last = &mut masks[n_merged - 1];
            if last.group == group && start <= last.end {
                if end > last.end {
                    last.end = end;
                }
                continue;
            }
        }
        masks.swap(n_merged, i);
        n_merged += 1;
    }
    masks.truncate(n_merged);
    masks
}

/// Subtract the right intervals from the left intervals.
///
/// Returns, for every surviving fragment, the index of the left row it came from and
/// its start and end. Intervals are half-open. Zero-length intervals are empty: they
/// remove nothing on the right and produce no fragments on the left, also when the
/// other side is empty. (The previous event sweep could emit a bogus fragment for an
/// empty left interval, and could drop a left interval starting at an empty right one.)
///
/// Reversed intervals (`start > end`) are not validated here, as elsewhere in this
/// crate. They are skipped on both sides, like in `tile`, because merging the right
/// side relies on `start < end`.
///
/// With `sort_output` the fragments are in left-row order, and each row's fragments are
/// in ascending coordinate order. Without it they are in (group, start) traversal order.
///
/// After sorting and merging, traversal is O(n + u + k) for n left rows, u merged right
/// intervals and k output fragments. k itself can be O(n * u) when many long left
/// intervals span many separated right intervals.
pub fn sweep_line_subtract<G: GroupType, T: PositionType>(
    chrs1: &[G],
    starts1: &[T],
    ends1: &[T],
    chrs2: &[G],
    starts2: &[T],
    ends2: &[T],
    sort_output: bool,
) -> (Vec<u32>, Vec<T>, Vec<T>) {
    let masks = merge_masks(chrs2, starts2, ends2);

    let mut left: Vec<Interval<G, T>> = (0..chrs1.len())
        .filter(|&i| starts1[i] < ends1[i])
        .map(|i| Interval {
            group: chrs1[i],
            start: starts1[i],
            end: ends1[i],
            idx: i as u32,
        })
        .collect();
    // Already-sorted input keeps its order, so output order needs no restoring.
    let reordered = !is_sorted_by_group_start(chrs1, starts1);
    if reordered {
        sort_by_key(&mut left, |l| l.start);
        sort_by_key(&mut left, |l| l.group);
    }

    let mut out_idxs = Vec::with_capacity(left.len());
    let mut out_starts = Vec::with_capacity(left.len());
    let mut out_ends = Vec::with_capacity(left.len());

    // Left starts are non-decreasing, so masks ending at or before the current left
    // start can never touch a later left interval. Each left interval scans forward
    // from there with its own pointer, since nested left intervals may still need
    // masks that an earlier, longer left interval has scanned past.
    let mut first = 0;
    for l in &left {
        while first < masks.len()
            && (masks[first].group < l.group
                || (masks[first].group == l.group && masks[first].end <= l.start))
        {
            first += 1;
        }

        let mut cursor = l.start;
        let mut j = first;
        while j < masks.len() && masks[j].group == l.group && masks[j].start < l.end {
            if masks[j].start > cursor {
                out_idxs.push(l.idx);
                out_starts.push(cursor);
                out_ends.push(masks[j].start);
            }
            cursor = masks[j].end;
            j += 1;
        }
        if cursor < l.end {
            out_idxs.push(l.idx);
            out_starts.push(cursor);
            out_ends.push(l.end);
        }
    }

    if sort_output && reordered {
        restore_row_order(chrs1.len(), &mut out_idxs, &mut out_starts, &mut out_ends);
    }

    (out_idxs, out_starts, out_ends)
}

/// Stable counting scatter by left-row index: O(n + k), and keeps each row's fragments
/// in the ascending order they were emitted in.
fn restore_row_order<T: PositionType>(
    n_rows: usize,
    idxs: &mut Vec<u32>,
    starts: &mut Vec<T>,
    ends: &mut Vec<T>,
) {
    let mut offsets = vec![0u32; n_rows + 1];
    for &i in idxs.iter() {
        offsets[i as usize + 1] += 1;
    }
    for r in 0..n_rows {
        offsets[r + 1] += offsets[r];
    }

    let k = idxs.len();
    let mut new_idxs = vec![0u32; k];
    let mut new_starts = vec![T::zero(); k];
    let mut new_ends = vec![T::zero(); k];
    for f in 0..k {
        let row = idxs[f] as usize;
        let to = offsets[row] as usize;
        offsets[row] += 1;
        new_idxs[to] = idxs[f];
        new_starts[to] = starts[f];
        new_ends[to] = ends[f];
    }
    *idxs = new_idxs;
    *starts = new_starts;
    *ends = new_ends;
}

#[cfg(test)]
mod tests {
    use super::sweep_line_subtract;

    fn run(
        left: &[(u32, i32, i32)],
        right: &[(u32, i32, i32)],
        sort_output: bool,
    ) -> Vec<(u32, i32, i32)> {
        let col = |v: &[(u32, i32, i32)]| {
            (
                v.iter().map(|r| r.0).collect::<Vec<_>>(),
                v.iter().map(|r| r.1).collect::<Vec<_>>(),
                v.iter().map(|r| r.2).collect::<Vec<_>>(),
            )
        };
        let (g1, s1, e1) = col(left);
        let (g2, s2, e2) = col(right);
        let (idx, starts, ends) = sweep_line_subtract(&g1, &s1, &e1, &g2, &s2, &e2, sort_output);
        (0..idx.len())
            .map(|i| (idx[i], starts[i], ends[i]))
            .collect()
    }

    #[test]
    fn subtract_keeps_touching_and_merges_touching_masks() {
        assert_eq!(run(&[(0, 0, 10)], &[(0, 10, 20)], true), vec![(0, 0, 10)]);
        assert_eq!(run(&[(0, 10, 20)], &[(0, 0, 10)], true), vec![(0, 10, 20)]);
        assert_eq!(
            run(&[(0, 0, 10)], &[(0, 2, 5), (0, 5, 8)], true),
            vec![(0, 0, 2), (0, 8, 10)]
        );
    }

    #[test]
    fn subtract_nested_left_sees_masks_passed_by_longer_left() {
        let left = [(0, 0, 100), (0, 10, 20), (0, 50, 60)];
        let right = [(0, 55, 56), (0, 12, 14), (0, 90, 95)];
        assert_eq!(
            run(&left, &right, true),
            vec![
                (0, 0, 12),
                (0, 14, 55),
                (0, 56, 90),
                (0, 95, 100),
                (1, 10, 12),
                (1, 14, 20),
                (2, 50, 55),
                (2, 56, 60),
            ]
        );
    }

    #[test]
    fn subtract_output_order() {
        let left = [(1, 0, 10), (0, 20, 30), (0, 0, 10)];
        let right = [(0, 4, 6), (1, 4, 6), (0, 24, 26)];
        // Row order, each row's fragments ascending.
        assert_eq!(
            run(&left, &right, true),
            vec![
                (0, 0, 4),
                (0, 6, 10),
                (1, 20, 24),
                (1, 26, 30),
                (2, 0, 4),
                (2, 6, 10)
            ]
        );
        // (group, start) traversal order.
        assert_eq!(
            run(&left, &right, false),
            vec![
                (2, 0, 4),
                (2, 6, 10),
                (1, 20, 24),
                (1, 26, 30),
                (0, 0, 4),
                (0, 6, 10)
            ]
        );
    }

    #[test]
    fn subtract_treats_zero_length_as_empty() {
        assert_eq!(run(&[(0, 5, 5)], &[(0, 8, 9)], true), vec![]);
        assert_eq!(run(&[(0, 5, 10)], &[(0, 5, 5)], true), vec![(0, 5, 10)]);
        assert_eq!(run(&[(0, 5, 5), (0, 1, 2)], &[], true), vec![(1, 1, 2)]);
        assert_eq!(run(&[], &[(0, 1, 2)], true), vec![]);
    }
}
