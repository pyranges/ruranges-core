use ruranges_core::complement_single::sweep_line_complement;
use rustc_hash::FxHashMap;

macro_rules! complement_gap_tests {
    ($module:ident, $pos:ty) => {
        mod $module {
            use super::*;

            fn run(
                rows: &[(u32, $pos, $pos)],
                lengths: &[(u32, $pos)],
                slack: $pos,
                include_first: bool,
            ) -> Vec<(u32, $pos, $pos)> {
                let groups: Vec<_> = rows.iter().map(|row| row.0).collect();
                let starts: Vec<_> = rows.iter().map(|row| row.1).collect();
                let ends: Vec<_> = rows.iter().map(|row| row.2).collect();
                let lengths: FxHashMap<_, _> = lengths.iter().copied().collect();
                let (out_groups, out_starts, out_ends, indices) =
                    sweep_line_complement(&groups, &starts, &ends, slack, &lengths, include_first);
                assert_eq!(out_groups.len(), out_starts.len());
                assert_eq!(out_groups.len(), out_ends.len());
                assert_eq!(out_groups.len(), indices.len());
                for (&group, &index) in out_groups.iter().zip(&indices) {
                    assert_eq!(groups[index as usize], group);
                }
                out_groups
                    .into_iter()
                    .zip(out_starts)
                    .zip(out_ends)
                    .map(|((group, start), end)| (group, start, end))
                    .collect()
            }

            #[test]
            fn no_empty_or_reversed_terminal_gap_at_either_group_flush() {
                // Reverse group order exercises sorting as well as both flush paths.
                for end in [100, 101] {
                    for slack in [0, 5] {
                        for include_first in [false, true] {
                            let expected = if include_first {
                                vec![(7, 0, 90), (8, 0, 90)]
                            } else {
                                vec![]
                            };
                            assert_eq!(
                                run(
                                    &[(8, 90, end), (7, 90, end)],
                                    &[(7, 100), (8, 100)],
                                    slack,
                                    include_first
                                ),
                                expected,
                            );
                        }
                    }
                }
            }

            #[test]
            fn positive_slack_preserves_surviving_gap_coordinates() {
                assert_eq!(
                    run(&[(7, 30, 40), (7, 10, 20)], &[(7, 100)], 5, false),
                    vec![(7, 20, 30), (7, 40, 100)],
                );
            }

            #[test]
            fn slack_suppresses_internal_and_terminal_gaps_at_the_threshold() {
                // Gaps of 5, 6, and 5 bases; only the 6-base gap survives slack 5.
                let rows = [(7, 0, 10), (7, 15, 20), (7, 26, 30)];
                assert_eq!(run(&rows, &[(7, 35)], 5, false), vec![(7, 20, 26)]);
                assert_eq!(
                    run(&rows, &[(7, 35)], 4, false),
                    vec![(7, 10, 15), (7, 20, 26), (7, 30, 35)],
                );
                assert_eq!(run(&rows, &[(7, 35)], 6, false), vec![]);
            }

            #[test]
            fn slack_preserves_leading_gap_and_nested_interval_union() {
                // The leading gap is optional and has no preceding interval.
                assert_eq!(
                    run(&[(7, 2, 20), (7, 5, 10), (7, 30, 40)], &[(7, 100)], 5, true),
                    vec![(7, 0, 2), (7, 20, 30), (7, 40, 100)],
                );
            }

            #[test]
            fn positive_slack_handles_coordinate_limits_without_overflow() {
                let max = <$pos>::MAX;
                assert_eq!(
                    run(&[(7, max - 10, max - 2)], &[(7, max)], 5, false),
                    vec![],
                );
                assert_eq!(
                    run(&[(7, max - 10, max - 6)], &[(7, max)], 5, false),
                    vec![(7, max - 6, max)],
                );
            }

            #[test]
            fn negative_slack_retains_existing_contracted_end_coordinates() {
                assert_eq!(
                    run(&[(7, 10, 20), (7, 30, 40)], &[(7, 100)], -2, false),
                    vec![(7, 18, 30), (7, 38, 100)],
                );
            }

            #[test]
            fn randomized_gaps_match_an_independent_interval_union() {
                let mut seed = 44_u64;
                let mut next = || {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    seed >> 32
                };
                for _ in 0..300 {
                    let mut rows = Vec::new();
                    for _ in 0..(next() % 25) {
                        let group = 7 + (next() % 3) as u32;
                        let start = (next() % 99) as $pos;
                        let end = start + 1 + (next() % (100 - start as u64)) as $pos;
                        rows.push((group, start, end));
                    }
                    for slack in [0, 1, 5, 20, 100] {
                        for include_first in [false, true] {
                            let mut expected = Vec::new();
                            for group in 7..=9 {
                                let mut intervals: Vec<_> = rows
                                    .iter()
                                    .filter(|row| row.0 == group)
                                    .map(|row| (row.1, row.2))
                                    .collect();
                                intervals.sort_unstable();
                                let mut union: Vec<($pos, $pos)> = Vec::new();
                                for (start, end) in intervals {
                                    if let Some(last) = union.last_mut() {
                                        if start <= last.1 {
                                            last.1 = last.1.max(end);
                                            continue;
                                        }
                                    }
                                    union.push((start, end));
                                }
                                if let Some(first) = union.first() {
                                    if include_first && first.0 > 0 {
                                        expected.push((group, 0, first.0));
                                    }
                                    for pair in union.windows(2) {
                                        if pair[1].0 - pair[0].1 > slack {
                                            expected.push((group, pair[0].1, pair[1].0));
                                        }
                                    }
                                    let end = union.last().unwrap().1;
                                    if 100 - end > slack {
                                        expected.push((group, end, 100));
                                    }
                                }
                            }
                            assert_eq!(
                                run(&rows, &[(7, 100), (8, 100), (9, 100)], slack, include_first),
                                expected,
                                "rows={rows:?}, slack={slack}, include_first={include_first}",
                            );
                        }
                    }
                }
            }
        }
    };
}

complement_gap_tests!(i32_positions, i32);
complement_gap_tests!(i64_positions, i64);
