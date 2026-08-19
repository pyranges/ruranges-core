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
//! Four properties are checked, on every shape:
//!
//! 1. `sweep_line_cluster` agrees with the endpoint-event sweep it replaced.
//! 2. `sweep_line_merge` agrees with the serial scan it replaced.
//! 3. Where the narrowed sweeps apply, they agree with the wide ones; where they
//!    do not, they decline rather than truncating a coordinate.
//! 4. The same, for the binary kernels — `subtract`, `overlaps`,
//!    `overlap_indices`, `count_overlaps` and `nearest`. `nearest` returns a
//!    *difference* of two coordinates, so it takes the halved
//!    `NARROW_DISTANCE_LIMIT` and declines shapes the others accept; that
//!    difference is one of the things this pins.

use ruranges_core::cluster::{sweep_line_cluster, sweep_line_cluster_narrowed};
use ruranges_core::merge::{sweep_line_merge, sweep_line_merge_narrowed};
use ruranges_core::nearest::{nearest_with_ties, nearest_with_ties_narrowed, Ties};
use ruranges_core::overlap_indices::{
    overlap_indices_chunked, overlap_indices_chunked_narrowed,
};
use ruranges_core::overlaps::{count_overlaps, count_overlaps_narrowed, overlaps, overlaps_narrowed};
use ruranges_core::sorts::{
    build_sorted_events_single_collection, build_sorted_intervals, NARROW_COORDINATE_LIMIT,
    NARROW_DISTANCE_LIMIT,
};
use ruranges_core::subtract::{sweep_line_subtract, sweep_line_subtract_narrowed};

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
    /// Whether the *distance*-producing narrowing is expected to take it, which
    /// needs the halved bound — see [`NARROW_DISTANCE_LIMIT`].
    narrows_distance: bool,
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
            narrows_distance: true,
        },
        Shape {
            name: "genomic, many groups",
            rows,
            groups: 500_000,
            origin: 0,
            span: 20_000,
            max_length: 60,
            narrows: true,
            narrows_distance: true,
        },
        Shape {
            name: "one group, deep pileup",
            rows,
            groups: 1,
            origin: 0,
            span: 100_000,
            max_length: 5_000,
            narrows: true,
            narrows_distance: true,
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
            narrows_distance: true,
        },
        Shape {
            name: "at the narrowing bound",
            rows,
            groups: 8,
            origin: NARROW_COORDINATE_LIMIT - 1_000_000,
            span: 999_000,
            max_length: 500,
            narrows: true,
            narrows_distance: false,
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
            narrows_distance: false,
        },
        Shape {
            name: "straddling the narrowing bound",
            rows,
            groups: 24,
            origin: NARROW_COORDINATE_LIMIT - 1000,
            span: 2000,
            max_length: 400,
            narrows: false,
            narrows_distance: false,
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

// ── the binary kernels ──────────────────────────────────────────────────────
//
// Same question as `the_narrowed_sweeps_either_agree_or_decline`, for the
// kernels that take two sets. The shapes below keep the *coordinate ranges* of
// the unary ones — that is what decides the narrowing — but spread the
// intervals out, because a binary kernel on a deep pileup emits pairs
// quadratically and would be measuring allocation rather than narrowing.

/// Rows per side. Both sides have to clear `PARALLEL_THRESHOLD` for the
/// narrowed builder to be offered the input at all.
const BINARY_ROWS: usize = 300_000;

fn binary_shapes() -> Vec<Shape> {
    let rows = BINARY_ROWS;
    vec![
        Shape {
            name: "genomic, few groups",
            rows,
            groups: 24,
            origin: 0,
            span: 250_000_000,
            max_length: 400,
            narrows: true,
            narrows_distance: true,
        },
        Shape {
            name: "genomic, many groups",
            rows,
            groups: 500_000,
            origin: 0,
            span: 20_000,
            max_length: 60,
            narrows: true,
            narrows_distance: true,
        },
        Shape {
            name: "one group, moderate pileup",
            rows,
            groups: 1,
            origin: 0,
            span: 5_000_000,
            max_length: 5_000,
            narrows: true,
            narrows_distance: true,
        },
        Shape {
            name: "negative coordinates",
            rows,
            groups: 24,
            origin: -500_000_000,
            span: 400_000_000,
            max_length: 400,
            narrows: true,
            narrows_distance: true,
        },
        Shape {
            // Inside the coordinate bound and outside the distance bound, which
            // is the one shape where the two narrowings must disagree.
            name: "at the narrowing bound",
            rows,
            groups: 8,
            origin: NARROW_COORDINATE_LIMIT - 100_000_000,
            span: 99_000_000,
            max_length: 500,
            narrows: true,
            narrows_distance: false,
        },
        Shape {
            name: "past the narrowing bound",
            rows,
            groups: 24,
            origin: 4 * NARROW_COORDINATE_LIMIT,
            span: 8_000_000_000,
            max_length: 400,
            narrows: false,
            narrows_distance: false,
        },
        Shape {
            name: "straddling the narrowing bound",
            rows,
            groups: 24,
            origin: NARROW_COORDINATE_LIMIT - 100_000_000,
            span: 200_000_000,
            max_length: 400,
            narrows: false,
            narrows_distance: false,
        },
    ]
}

fn flatten(chunks: Vec<Vec<u32>>) -> Vec<u32> {
    chunks.into_iter().flatten().collect()
}

/// Every narrowed binary kernel either gives the wide answer or refuses the
/// input — and `nearest` refuses one shape the others take, because a distance
/// needs the whole span to fit `i32`, not just each endpoint.
#[test]
#[ignore = "moves several GB; run with --ignored"]
fn the_narrowed_binary_kernels_either_agree_or_decline() {
    for (index, shape) in binary_shapes().iter().enumerate() {
        let (chrs, starts, ends) = generate(shape, 0xB1_0000 + index as u64);
        let (chrs2, starts2, ends2) = generate(shape, 0xB2_0000 + index as u64);

        // subtract carries no slack.
        let narrowed = sweep_line_subtract_narrowed(
            &chrs, &starts, &ends, &chrs2, &starts2, &ends2, true,
        );
        assert_eq!(
            narrowed.is_some(),
            shape.narrows,
            "subtract narrowing on shape={}",
            shape.name
        );
        if let Some(narrowed) = narrowed {
            assert_eq!(
                narrowed,
                sweep_line_subtract(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, true),
                "subtract shape={}",
                shape.name
            );
        }

        for &slack in &[0_i64, 250] {
            let narrowed =
                count_overlaps_narrowed(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack);
            assert_eq!(
                narrowed.is_some(),
                shape.narrows,
                "count_overlaps narrowing on shape={} slack={slack}",
                shape.name
            );
            if let Some(narrowed) = narrowed {
                assert_eq!(
                    narrowed,
                    count_overlaps(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack),
                    "count_overlaps shape={} slack={slack}",
                    shape.name
                );
            }

            for overlap_type in ["all", "first", "last"] {
                let narrowed = overlaps_narrowed(
                    &chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack, overlap_type, true,
                    false,
                );
                assert_eq!(
                    narrowed.is_some(),
                    shape.narrows,
                    "overlaps narrowing on shape={} type={overlap_type} slack={slack}",
                    shape.name
                );
                if let Some(narrowed) = narrowed {
                    assert_eq!(
                        narrowed,
                        overlaps(
                            &chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack, overlap_type,
                            true, false,
                        ),
                        "overlaps shape={} type={overlap_type} slack={slack}",
                        shape.name
                    );
                }

                // `contained` takes the sweep rather than the closed form, so it
                // is a different path through the same records.
                let narrowed = overlap_indices_chunked_narrowed(
                    &chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack, overlap_type, true,
                );
                assert_eq!(
                    narrowed.is_some(),
                    shape.narrows,
                    "overlap_indices narrowing on shape={} type={overlap_type} slack={slack}",
                    shape.name
                );
                if let Some(narrowed) = narrowed {
                    assert_eq!(
                        flatten(narrowed),
                        flatten(overlap_indices_chunked(
                            &chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack, overlap_type,
                            true,
                        )),
                        "overlap_indices shape={} type={overlap_type} slack={slack}",
                        shape.name
                    );
                }
            }

            for (include_overlaps, direction) in
                [(true, "any"), (false, "any"), (true, "forward"), (true, "backward")]
            {
                let narrowed = nearest_with_ties_narrowed(
                    &chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack, 1, include_overlaps,
                    direction, true, Ties::All,
                );
                assert_eq!(
                    narrowed.is_some(),
                    shape.narrows_distance,
                    "nearest narrowing on shape={} dir={direction} slack={slack}",
                    shape.name
                );
                if let Some(narrowed) = narrowed {
                    assert_eq!(
                        narrowed,
                        nearest_with_ties(
                            &chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack, 1,
                            include_overlaps, direction, true, Ties::All,
                        ),
                        "nearest shape={} dir={direction} slack={slack}",
                        shape.name
                    );
                }
            }
        }

        println!(
            "binary narrowing ok: {} (narrows={}, distance={})",
            shape.name, shape.narrows, shape.narrows_distance
        );
    }
}

/// A slack too large for the narrowed range must be declined by the binary
/// kernels too — and `nearest`, which is bounded by the halved limit, declines
/// a slack the others still accept.
#[test]
#[ignore = "moves several GB; run with --ignored"]
fn a_binary_slack_that_could_overflow_is_declined() {
    let shape = &binary_shapes()[0];
    let (chrs, starts, ends) = generate(shape, 0x5B1_0000);
    let (chrs2, starts2, ends2) = generate(shape, 0x5B2_0000);

    for &slack in &[NARROW_COORDINATE_LIMIT, NARROW_COORDINATE_LIMIT + 1, i64::MAX] {
        assert!(
            count_overlaps_narrowed(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack)
                .is_none()
        );
        assert!(overlaps_narrowed(
            &chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack, "all", true, false
        )
        .is_none());
        assert!(overlap_indices_chunked_narrowed(
            &chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack, "all", false
        )
        .is_none());
    }

    // Accepted by the coordinate-bounded kernels, refused by the distance-bounded
    // one: `nearest` would have to add this to a coordinate and then subtract
    // another from the result.
    let slack = NARROW_COORDINATE_LIMIT - 1;
    assert_eq!(
        count_overlaps_narrowed(&chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack),
        Some(count_overlaps(
            &chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack
        ))
    );
    assert!(nearest_with_ties_narrowed(
        &chrs, &starts, &ends, &chrs2, &starts2, &ends2, slack, 1, true, "any", true, Ties::All
    )
    .is_none());
    assert!(nearest_with_ties_narrowed(
        &chrs,
        &starts,
        &ends,
        &chrs2,
        &starts2,
        &ends2,
        NARROW_DISTANCE_LIMIT,
        1,
        true,
        "any",
        true,
        Ties::All
    )
    .is_none());
}
