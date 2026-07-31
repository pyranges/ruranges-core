//! Order-preserving integer codes for string sort keys.
//!
//! Sorting a genomic frame means turning every key column into an ascending
//! integer and handing the result to [`crate::sorts::sort_order_idx`]. The
//! expensive half of that is ordering a column's *distinct* string values —
//! chromosome names, transcript ids — and the ordering rule has to be identical
//! in every library that offers the operation, or the same data sorts two ways.
//! It therefore lives here, next to the kernel that consumes it, rather than
//! once in `pyranges1` and once in `polaranges`.
//!
//! # Natural order
//!
//! [`natural_cmp`] is the specification. A string is read as an alternating
//! sequence of non-digit runs and digit runs; non-digit runs compare as byte
//! strings, digit runs compare as numbers, and a string that runs out of
//! material sorts before one that has more. So `t2 < t10`, `chr9 < chr10 <
//! chrM`, and `a < a1 < a1a < a1b < ab`.
//!
//! This reproduces Python's `natsort.natsorted` under its default settings,
//! whose key is the same alternating `(text, number, text, number, …)` tuple,
//! with two deliberate differences:
//!
//! * **No Unicode normalisation.** `natsort` NFD-decomposes before comparing,
//!   so it orders `e2 < é1 < z1`; here non-digit runs compare by UTF-8 bytes,
//!   which is code-point order, giving `e2 < z1 < é1`. Only strings differing
//!   in composed non-ASCII characters can tell the two apart.
//! * **Digits are ASCII.** Other Unicode decimal digits are ordinary
//!   characters here and compare as bytes.
//!
//! Both libraries' key columns are chromosome names, contig names and feature
//! ids, so neither difference is reachable from realistic data — but they are
//! differences, and they are tested as such.
//!
//! Values that compare equal — `t1` and `t01` are equal to `natsort` too, since
//! leading zeros do not survive the conversion to a number — keep their input
//! order and still receive distinct ranks. That matches what `pyranges1` has
//! always done by natural-sorting a de-duplicated key frame with a stable sort.

use std::cmp::Ordering;

use rayon::prelude::*;

/// Compare two strings in natural order: digit runs as numbers, everything else
/// as bytes.
///
/// Allocation-free — no padded or split key is materialised, which is what makes
/// it affordable to call `n log n` times over a column's distinct values.
///
/// ```
/// use std::cmp::Ordering;
/// use ruranges_core::ranks::natural_cmp;
///
/// assert_eq!(natural_cmp("chr2", "chr10"), Ordering::Less);
/// assert_eq!(natural_cmp("chr10", "chrM"), Ordering::Less);
/// assert_eq!(natural_cmp("t1", "t01"), Ordering::Equal);
/// ```
pub fn natural_cmp(a: &str, b: &str) -> Ordering {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    let (mut i, mut j) = (0usize, 0usize);

    loop {
        // Non-digit run. UTF-8 byte order is code-point order, so comparing the
        // raw slices is the same as comparing the decoded text.
        let a_text = i;
        while i < a.len() && !a[i].is_ascii_digit() {
            i += 1;
        }
        let b_text = j;
        while j < b.len() && !b[j].is_ascii_digit() {
            j += 1;
        }
        match a[a_text..i].cmp(&b[b_text..j]) {
            Ordering::Equal => {}
            other => return other,
        }

        // One or both strings ended. The one with nothing left is the prefix,
        // and a prefix sorts first — `a < a1`.
        match (i == a.len(), j == b.len()) {
            (true, true) => return Ordering::Equal,
            (true, false) => return Ordering::Less,
            (false, true) => return Ordering::Greater,
            (false, false) => {}
        }

        // Digit run, compared as a number: drop leading zeros, then more digits
        // means a larger number, then compare digit by digit. Length-first is
        // what keeps this correct for numbers too large for any integer type.
        while i < a.len() && a[i] == b'0' {
            i += 1;
        }
        while j < b.len() && b[j] == b'0' {
            j += 1;
        }
        let a_num = i;
        while i < a.len() && a[i].is_ascii_digit() {
            i += 1;
        }
        let b_num = j;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        let (da, db) = (&a[a_num..i], &b[b_num..j]);
        match da.len().cmp(&db.len()).then_with(|| da.cmp(db)) {
            Ordering::Equal => {}
            other => return other,
        }
    }
}

/// Rank `values` in natural order: `ranks[i]` is the position `values[i]` would
/// occupy if the slice were sorted.
///
/// Ranks are `0..values.len()`, one per element, so feeding the distinct values
/// of a column through this yields a code column that sorts exactly as the
/// strings do. Elements that compare equal keep their input order.
///
/// ```
/// use ruranges_core::ranks::natural_rank;
///
/// assert_eq!(natural_rank(&["chr10", "chr2", "chrX", "chr1"]), vec![2, 1, 3, 0]);
/// ```
///
/// # Panics
///
/// If `values.len()` exceeds `u32::MAX`.
pub fn natural_rank(values: &[&str]) -> Vec<u32> {
    rank_by(values, natural_cmp)
}

/// Rank `values` in byte-lexical order — `t10` before `t9`.
///
/// The counterpart to [`natural_rank`] for `natsort=False` and for the generic
/// (non-genomic) namespace, which has no chromosome column and so no reason to
/// order strings naturally.
///
/// # Panics
///
/// If `values.len()` exceeds `u32::MAX`.
pub fn lexical_rank(values: &[&str]) -> Vec<u32> {
    rank_by(values, |a, b| a.as_bytes().cmp(b.as_bytes()))
}

fn rank_by(values: &[&str], compare: impl Fn(&str, &str) -> Ordering + Sync) -> Vec<u32> {
    assert!(
        values.len() <= u32::MAX as usize,
        "cannot rank more than u32::MAX values, got {}",
        values.len()
    );

    let mut order: Vec<u32> = (0..values.len() as u32).collect();
    // Breaking ties on the original index makes an unstable sort stable by
    // construction, so this keeps `par_sort_unstable_by`'s speed without
    // leaving the order of equal values up to the scheduler.
    order.par_sort_unstable_by(|&left, &right| {
        compare(values[left as usize], values[right as usize]).then_with(|| left.cmp(&right))
    });

    let mut ranks = vec![0_u32; values.len()];
    for (rank, &position) in order.iter().enumerate() {
        ranks[position as usize] = rank as u32;
    }
    ranks
}

/// Fold one key column's codes into the running group id, so that a whole key
/// list collapses to the single `groups` argument
/// [`crate::sorts::sort_order_idx`] takes.
///
/// `group` holds the codes of every key so far and has `group_cardinality`
/// distinct values; `codes` is the next key, with `codes_cardinality` distinct
/// values. On return `group` orders rows by the old key first and the new key
/// second, and the new cardinality is returned.
///
/// The straightforward `group * codes_cardinality + codes` overflows `u32` once
/// the keys are wide enough, so when the product would not fit, the combined
/// key is re-densified instead: after densifying, cardinality is at most the row
/// count, which callers already cap at `u32`.
///
/// ```
/// use ruranges_core::ranks::fold_ranks;
///
/// // Chromosome codes, then strand codes.
/// let mut group = vec![0, 0, 1, 1];
/// let cardinality = fold_ranks(&mut group, 2, &[1, 0, 1, 0], 2);
/// assert_eq!(group, vec![1, 0, 3, 2]);
/// assert_eq!(cardinality, 4);
/// ```
///
/// # Panics
///
/// If `group` and `codes` have different lengths.
pub fn fold_ranks(
    group: &mut [u32],
    group_cardinality: u32,
    codes: &[u32],
    codes_cardinality: u32,
) -> u32 {
    assert_eq!(
        group.len(),
        codes.len(),
        "group and codes must describe the same rows"
    );
    if group.is_empty() {
        return 0;
    }

    let product = group_cardinality as u64 * codes_cardinality as u64;
    if product <= u32::MAX as u64 {
        for (slot, &code) in group.iter_mut().zip(codes) {
            *slot = *slot * codes_cardinality + code;
        }
        return product as u32;
    }

    // Too wide to pack. Sort the (old key, new key) pairs and number them, which
    // gives the same order in at most `n` distinct values.
    let combined: Vec<u64> = group
        .iter()
        .zip(codes)
        .map(|(&high, &low)| ((high as u64) << 32) | low as u64)
        .collect();
    let mut order: Vec<u32> = (0..group.len() as u32).collect();
    radsort::sort_by_key(&mut order, |&row| combined[row as usize]);

    let mut dense = 0_u32;
    let mut previous = combined[order[0] as usize];
    for &row in &order {
        let key = combined[row as usize];
        if key != previous {
            dense += 1;
            previous = key;
        }
        group[row as usize] = dense;
    }
    dense + 1
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Order a slice with the comparator, the way `rank_by` does.
    fn natural_sorted(values: &[&str]) -> Vec<String> {
        let ranks = natural_rank(values);
        let mut out = vec![String::new(); values.len()];
        for (value, rank) in values.iter().zip(&ranks) {
            out[*rank as usize] = (*value).to_string();
        }
        out
    }

    /// Every expectation here was produced by `natsort.natsorted` 8.4.0 under its
    /// default settings, so this is the agreement test between the two libraries
    /// and the Python implementation they are replacing.
    #[test]
    fn matches_python_natsort() {
        let cases: &[(&[&str], &[&str])] = &[
            (&["t1", "t10", "t9", "t2"], &["t1", "t2", "t9", "t10"]),
            (
                &["t1", "t01", "t001", "t10", "t010"],
                &["t1", "t01", "t001", "t10", "t010"],
            ),
            (
                &["chr10", "chr2", "chrX", "chrM", "chr1", "chrY"],
                &["chr1", "chr2", "chr10", "chrM", "chrX", "chrY"],
            ),
            (
                &["1", "01", "001", "10", "2"],
                &["1", "01", "001", "2", "10"],
            ),
            (&["B", "a", "A", "b"], &["A", "B", "a", "b"]),
            (
                &["a", "a1", "a1b", "a1a", "ab"],
                &["a", "a1", "a1a", "a1b", "ab"],
            ),
            (&["", "a", "1"], &["", "1", "a"]),
            // `-` is not a sign: `a-1` splits as ("a-", 1), and "a" < "a-".
            (&["a-1", "a-2", "a1", "a2"], &["a1", "a2", "a-1", "a-2"]),
            // Default natsort is integer-only, so `1.10` is (1, ".", 10).
            (
                &["v1.5", "v1.10", "v1.2"],
                &["v1.2", "v1.5", "v1.10"],
            ),
            (
                &["a1z", "a1a9", "a1a10"],
                &["a1a9", "a1a10", "a1z"],
            ),
            (&["a1", "a1 ", "a1-"], &["a1", "a1 ", "a1-"]),
        ];
        for (input, expected) in cases {
            assert_eq!(&natural_sorted(input), expected, "input {input:?}");
        }
    }

    /// Numbers wider than any integer type still compare as numbers, because the
    /// digit run is compared by length before content.
    #[test]
    fn compares_arbitrarily_large_numbers() {
        let big_nines = format!("x{}", "9".repeat(30));
        let bigger_ones = format!("x{}", "1".repeat(31));
        let values = [bigger_ones.as_str(), "x1", big_nines.as_str()];
        assert_eq!(
            natural_sorted(&values),
            vec!["x1".to_string(), big_nines.clone(), bigger_ones.clone()],
        );
    }

    /// The documented divergence from Python `natsort`, which NFD-normalises and
    /// would give `e2 < é1 < z1`.
    #[test]
    fn orders_non_ascii_by_code_point_not_nfd() {
        assert_eq!(natural_sorted(&["é1", "e2", "z1"]), vec!["e2", "z1", "é1"]);
    }

    #[test]
    fn equal_values_keep_input_order_and_distinct_ranks() {
        assert_eq!(natural_cmp("t01", "t1"), Ordering::Equal);
        assert_eq!(natural_rank(&["t01", "t1", "t001"]), vec![0, 1, 2]);
        assert_eq!(natural_rank(&["t1", "t001", "t01"]), vec![0, 1, 2]);
    }

    #[test]
    fn lexical_rank_is_byte_order() {
        assert_eq!(lexical_rank(&["t1", "t10", "t9", "t2"]), vec![0, 1, 3, 2]);
        assert_eq!(lexical_rank(&["chr10", "chr2", "chr1"]), vec![1, 2, 0]);
    }

    #[test]
    fn ranks_handle_empty_and_singleton_input() {
        assert_eq!(natural_rank(&[]), Vec::<u32>::new());
        assert_eq!(natural_rank(&["only"]), vec![0]);
        assert_eq!(lexical_rank(&[]), Vec::<u32>::new());
    }

    #[test]
    fn fold_ranks_packs_when_it_fits() {
        let mut group = vec![0, 1, 2, 0];
        let cardinality = fold_ranks(&mut group, 3, &[2, 0, 1, 1], 3);
        assert_eq!(group, vec![2, 3, 7, 1]);
        assert_eq!(cardinality, 9);
    }

    #[test]
    fn fold_ranks_orders_by_old_key_then_new() {
        // Same rows, shuffled: the fold must order by `group` first.
        let mut group = vec![1, 0, 1, 0];
        fold_ranks(&mut group, 2, &[0, 1, 1, 0], 2);
        assert!(group[1] < group[0] && group[0] < group[2]);
        assert!(group[3] < group[1]);
    }

    /// The densifying branch: cardinalities whose product overflows `u32` must
    /// still produce the same *order*, in at most `n` distinct codes.
    #[test]
    fn fold_ranks_densifies_instead_of_overflowing() {
        let mut group = vec![7, 7, 1, 1, 3];
        let cardinality = fold_ranks(&mut group, 100_000, &[5, 2, 9, 9, 0], 100_000);
        // (1,9) (1,9) (3,0) (7,2) (7,5) -> 0 0 1 2 3
        assert_eq!(group, vec![3, 2, 0, 0, 1]);
        assert_eq!(cardinality, 4);
    }

    #[test]
    fn fold_ranks_accepts_empty_input() {
        let mut group: Vec<u32> = Vec::new();
        assert_eq!(fold_ranks(&mut group, 0, &[], 0), 0);
    }
}
