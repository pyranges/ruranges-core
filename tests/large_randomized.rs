//! Large randomized comparisons for the merge and cluster sweeps.
//!
//! The unit tests next to each kernel run at 150k rows, which is enough to cross
//! `PARALLEL_THRESHOLD` and exercise both the counting partition and the radix
//! fallback. These run the same comparisons at millions of rows and across the
//! coordinate ranges that decide whether the narrowed `i32` sweep applies —
//! including ranges where it must *not* apply, so the `i64` path is what gets
//! checked.
//!
//! They are `#[ignore]`d because a full pass moves several GB; run them with
//!
//! ```text
//! cargo test --release --test large_randomized -- --ignored --nocapture
//! ```
//!
//! Three properties are checked, on every shape:
//!
//! 1. `sweep_line_cluster` agrees with the endpoint-event sweep it replaced.
//! 2. `sweep_line_merge` agrees with the serial scan it replaced.
//! 3. Where the narrowed sweeps apply, they agree with the wide ones; where they
//!    do not, they decline rather than truncating a coordinate.

use ruranges_core::cluster::{sweep_line_cluster, sweep_line_cluster_narrowed};
use ruranges_core::merge::{sweep_line_merge, sweep_line_merge_narrowed};
use ruranges_core::sorts::{
    build_sorted_events_single_collection, build_sorted_intervals, NARROW_COORDINATE_LIMIT,
};

/// Deterministic SplitMix64, so these tests need no dev-dependency.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn below(&mut self, bound: u64) -> u64 {
        self.next() % bound.max(1)
    }
}

/// One randomized input, described by the shape it is meant to exercise.
struct Shape {
    name: &'static str,
    rows: usize,
    groups: u64,
    /// Coordinates are drawn from `origin .. origin + span`.
    origin: i64,
    span: u64,
    /// Longest interval, so pileups get deep on the small spans.
    max_length: u64,
    /// Whether the narrowed sweep is expected to take this input.
    narrows: bool,
}

fn shapes(rows: usize) -> Vec<Shape> {
    vec![
        Shape {
            name: "genomic, few groups",
            rows,
            groups: 24,
            origin: 0,
            span: 250_000_000,
            max_length: 400,
            narrows: true,
        },
        Shape {
            name: "genomic, many groups",
            rows,
            groups: 500_000,
            origin: 0,
            span: 20_000,
            max_length: 60,
            narrows: true,
        },
        Shape {
            name: "one group, deep pileup",
            rows,
            groups: 1,
            origin: 0,
            span: 100_000,
            max_length: 5_000,
            narrows: true,
        },
        Shape {
            name: "negative coordinates",
            rows,
            groups: 24,
            // `upstream`/`extend_ranges` legitimately push intervals past an origin.
            origin: -500_000_000,
            span: 400_000_000,
            max_length: 400,
            narrows: true,
        },
        Shape {
            name: "at the narrowing bound",
            rows,
            groups: 8,
            origin: NARROW_COORDINATE_LIMIT - 1_000_000,
            span: 999_000,
            max_length: 500,
            narrows: true,
        },
        Shape {
            name: "past the narrowing bound",
            rows,
            groups: 24,
            // Beyond what `i32` can hold at all, so the `i64` sweep is the one
            // under test here.
            origin: 4 * NARROW_COORDINATE_LIMIT,
            span: 8_000_000_000,
            max_length: 400,
            narrows: false,
        },
        Shape {
            name: "straddling the narrowing bound",
            rows,
            groups: 24,
            origin: NARROW_COORDINATE_LIMIT - 1000,
            span: 2000,
            max_length: 400,
            narrows: false,
        },
    ]
}

fn generate(shape: &Shape, seed: u64) -> (Vec<u32>, Vec<i64>, Vec<i64>) {
    let mut rng = Rng(seed);
    let mut chrs = Vec::with_capacity(shape.rows);
    let mut starts = Vec::with_capacity(shape.rows);
    let mut ends = Vec::with_capacity(shape.rows);
    for _ in 0..shape.rows {
        chrs.push(rng.below(shape.groups) as u32);
        let start = shape.origin + rng.below(shape.span) as i64;
        starts.push(start);
        // Positive length, which is the invariant polaranges enforces before the
        // kernel and the one the sweeps are specified against.
        ends.push(start + 1 + rng.below(shape.max_length) as i64);
    }
    (chrs, starts, ends)
}

/// The endpoint-event sweep `sweep_line_cluster` replaced, kept as the oracle.
fn cluster_by_events(chrs: &[u32], starts: &[i64], ends: &[i64], slack: i64) -> (Vec<u32>, Vec<u32>) {
    let mut indices = Vec::with_capacity(chrs.len());
    let mut cluster_ids = Vec::with_capacity(chrs.len());
    if chrs.is_empty() {
        return (cluster_ids, indices);
    }

    let events = build_sorted_events_single_collection(chrs, starts, ends, slack);
    let mut current_chr = events.first().unwrap().chr;
    let mut current_cluster = 0_u32;
    let mut active = 0_i64;

    for event in events {
        if event.chr != current_chr {
            current_cluster += 1;
            active = 0;
            current_chr = event.chr;
        }
        if event.is_start {
            indices.push(event.idx);
            cluster_ids.push(current_cluster);
            active += 1;
        } else {
            active -= 1;
            if active == 0 {
                current_cluster += 1;
            }
        }
    }

    (cluster_ids, indices)
}

/// The serial scan `sweep_line_merge` replaced, kept as the oracle.
fn merge_serially(
    chrs: &[u32],
    starts: &[i64],
    ends: &[i64],
    slack: i64,
) -> (Vec<u32>, Vec<i64>, Vec<i64>, Vec<u32>) {
    let intervals = build_sorted_intervals(chrs, starts, ends, None, 0, false);
    let (mut idx, mut out_starts, mut out_ends, mut counts) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    if intervals.is_empty() {
        return (idx, out_starts, out_ends, counts);
    }

    let mut group = intervals[0].group;
    let mut current = (intervals[0].idx, intervals[0].start, intervals[0].end, 1_u32);
    for interval in intervals.iter().skip(1) {
        if interval.group == group && interval.start < current.2 + slack {
            current.2 = current.2.max(interval.end);
            current.3 += 1;
            continue;
        }
        idx.push(current.0);
        out_starts.push(current.1);
        out_ends.push(current.2);
        counts.push(current.3);
        group = interval.group;
        current = (interval.idx, interval.start, interval.end, 1);
    }
    idx.push(current.0);
    out_starts.push(current.1);
    out_ends.push(current.2);
    counts.push(current.3);

    (idx, out_starts, out_ends, counts)
}

/// Rows per shape. Large enough that every batch of group runs holds several
/// tasks' worth of work at 12 threads.
const ROWS: usize = 3_000_000;

#[test]
#[ignore = "moves several GB; run with --ignored"]
fn cluster_matches_the_event_sweep_at_scale() {
    for (index, shape) in shapes(ROWS).iter().enumerate() {
        let (chrs, starts, ends) = generate(shape, 0xC10_0000 + index as u64);
        for &slack in &[0_i64, 1, 250] {
            let observed = sweep_line_cluster(&chrs, &starts, &ends, slack);
            let expected = cluster_by_events(&chrs, &starts, &ends, slack);
            assert_eq!(
                observed, expected,
                "shape={} slack={slack}",
                shape.name
            );
        }
        println!("cluster ok: {} ({} rows)", shape.name, shape.rows);
    }
}

#[test]
#[ignore = "moves several GB; run with --ignored"]
fn merge_matches_the_serial_scan_at_scale() {
    for (index, shape) in shapes(ROWS).iter().enumerate() {
        let (chrs, starts, ends) = generate(shape, 0x11E_0000 + index as u64);
        for &slack in &[0_i64, 1, 250] {
            let observed = sweep_line_merge(&chrs, &starts, &ends, slack);
            let expected = merge_serially(&chrs, &starts, &ends, slack);
            assert_eq!(observed, expected, "shape={} slack={slack}", shape.name);
        }
        println!("merge ok: {} ({} rows)", shape.name, shape.rows);
    }
}

/// The narrowed sweeps must either give the wide answer or refuse the input.
/// Silently truncating a coordinate to `i32` is the failure this guards.
#[test]
#[ignore = "moves several GB; run with --ignored"]
fn the_narrowed_sweeps_either_agree_or_decline() {
    for (index, shape) in shapes(ROWS).iter().enumerate() {
        let (chrs, starts, ends) = generate(shape, 0x4A2_0000 + index as u64);
        for &slack in &[0_i64, 250] {
            let narrowed_merge = sweep_line_merge_narrowed(&chrs, &starts, &ends, slack);
            let narrowed_cluster = sweep_line_cluster_narrowed(&chrs, &starts, &ends, slack);
            assert_eq!(
                narrowed_merge.is_some(),
                shape.narrows,
                "merge narrowing on shape={} slack={slack}",
                shape.name
            );
            assert_eq!(
                narrowed_cluster.is_some(),
                shape.narrows,
                "cluster narrowing on shape={} slack={slack}",
                shape.name
            );

            if let Some(narrowed) = narrowed_merge {
                assert_eq!(
                    narrowed,
                    sweep_line_merge(&chrs, &starts, &ends, slack),
                    "merge shape={} slack={slack}",
                    shape.name
                );
            }
            if let Some(narrowed) = narrowed_cluster {
                assert_eq!(
                    narrowed,
                    sweep_line_cluster(&chrs, &starts, &ends, slack),
                    "cluster shape={} slack={slack}",
                    shape.name
                );
            }
        }
        println!(
            "narrowing ok: {} (narrows={})",
            shape.name, shape.narrows
        );
    }
}

/// A slack large enough to overflow the halved `i32` range must be declined even
/// when every coordinate fits, because the sweeps add it to an end.
#[test]
#[ignore = "moves several GB; run with --ignored"]
fn a_slack_that_could_overflow_is_declined() {
    let shape = &shapes(500_000)[0];
    let (chrs, starts, ends) = generate(shape, 0x5AC_0000);

    for &slack in &[NARROW_COORDINATE_LIMIT, NARROW_COORDINATE_LIMIT + 1, i64::MAX] {
        assert!(sweep_line_merge_narrowed(&chrs, &starts, &ends, slack).is_none());
        assert!(sweep_line_cluster_narrowed(&chrs, &starts, &ends, slack).is_none());
        // And the wide sweeps still answer, which is what the caller falls back to.
        assert_eq!(
            sweep_line_merge(&chrs, &starts, &ends, slack),
            merge_serially(&chrs, &starts, &ends, slack)
        );
    }

    // One below the bound is accepted, and must agree with the wide sweep.
    let slack = NARROW_COORDINATE_LIMIT - 1;
    assert_eq!(
        sweep_line_merge_narrowed(&chrs, &starts, &ends, slack),
        Some(sweep_line_merge(&chrs, &starts, &ends, slack))
    );
    assert_eq!(
        sweep_line_cluster_narrowed(&chrs, &starts, &ends, slack),
        Some(sweep_line_cluster(&chrs, &starts, &ends, slack))
    );
}
