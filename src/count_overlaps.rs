use radsort::sort_by_key;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::ruranges_structs::{GroupType, PositionType};

#[derive(Clone, Copy)]
struct Target<C, T> {
    group: C,
    start: T,
    end: T,
}

struct Group {
    start: usize,
    end: usize,
    regular: bool,
}

enum GroupIndex<C> {
    Dense(Vec<usize>),
    Sparse(FxHashMap<C, usize>),
}

struct GroupLookup<C> {
    ranges: Vec<Group>,
    index: GroupIndex<C>,
}

impl<C: GroupType> GroupLookup<C> {
    fn new(entries: Vec<(C, Group)>) -> Self {
        // Factorized groups usually have compact nonnegative labels. A direct
        // index avoids hashing every query and keeps many tiny groups local.
        // Arbitrary/signed labels retain the sparse representation.
        let max = entries.last().and_then(|(label, _)| label.to_usize());
        let dense = entries
            .first()
            .is_some_and(|(label, _)| label.to_usize().is_some())
            && max.is_some_and(|max| max < entries.len().saturating_mul(4).saturating_add(64));
        let mut ranges = Vec::with_capacity(entries.len());
        let index = if dense {
            let mut index = vec![0; max.expect("dense labels fit usize") + 1];
            for (label, group) in entries {
                index[label.to_usize().expect("dense labels fit usize")] = ranges.len() + 1;
                ranges.push(group);
            }
            GroupIndex::Dense(index)
        } else {
            let mut index = FxHashMap::default();
            index.reserve(entries.len());
            for (label, group) in entries {
                index.insert(label, ranges.len());
                ranges.push(group);
            }
            GroupIndex::Sparse(index)
        };
        Self { ranges, index }
    }

    fn get(&self, label: &C) -> Option<&Group> {
        let index = match &self.index {
            GroupIndex::Dense(index) => index.get(label.to_usize()?)?.checked_sub(1)?,
            GroupIndex::Sparse(index) => *index.get(label)?,
        };
        self.ranges.get(index)
    }
}

#[derive(Clone, Copy)]
struct Query<C, T> {
    group: C,
    coordinate: T,
    index: usize,
}

struct RankedQuery<T> {
    coordinate: T,
    index: usize,
    rank: usize,
}

#[derive(Clone, Copy)]
pub(crate) enum Strategy {
    Auto,
    Binary,
    Sorted,
}

fn sort_records<R, K: Ord + radsort::Key>(values: &mut [R], key: impl Fn(&R) -> K) {
    if values.is_sorted_by_key(&key) {
        return;
    }
    if values.len() <= 256 {
        values.sort_by_key(key);
    } else {
        sort_by_key(values, key);
    }
}

fn finish_count(count: usize) -> u32 {
    count.min(u32::MAX as usize) as u32
}

/// Count the two strict overlap inequalities without enumerating matching pairs.
#[allow(clippy::too_many_arguments)]
pub(crate) fn count_overlaps<C: GroupType, T: PositionType>(
    groups: &[C],
    starts: &[T],
    ends: &[T],
    groups2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
    strategy: Strategy,
) -> Vec<u32> {
    debug_assert_eq!(groups.len(), starts.len());
    debug_assert_eq!(groups.len(), ends.len());
    debug_assert_eq!(groups2.len(), starts2.len());
    debug_assert_eq!(groups2.len(), ends2.len());
    let mut counts = vec![0; groups.len()];
    if groups.is_empty() || groups2.is_empty() {
        return counts;
    }

    // Avoid sorting/group-map setup when either side is very small.
    if groups.len() <= 4 || groups2.len() <= 4 || groups.len().saturating_mul(groups2.len()) <= 256
    {
        for i in 0..groups.len() {
            let limit = ends[i].saturating_add(slack);
            counts[i] = finish_count(
                (0..groups2.len())
                    .filter(|&j| {
                        groups[i] == groups2[j]
                            && starts[i] < ends2[j].saturating_add(slack)
                            && starts2[j] < limit
                    })
                    .count(),
            );
        }
        return counts;
    }

    let mut reversed = FxHashSet::default();
    let mut targets: Vec<_> = (0..groups2.len())
        .map(|i| {
            if starts2[i] > ends2[i] {
                reversed.insert(groups2[i]);
            }
            Target {
                group: groups2[i],
                start: starts2[i],
                end: ends2[i].saturating_add(slack),
            }
        })
        .collect();
    sort_records(&mut targets, |t| t.group);

    let mut entries = Vec::new();
    let mut start = 0;
    let mut has_large_group = false;
    while start < targets.len() {
        let group = targets[start].group;
        let mut end = start + 1;
        while end < targets.len() && targets[end].group == group {
            end += 1;
        }
        if end - start > 8 {
            sort_records(&mut targets[start..end], |t| t.start);
            has_large_group = true;
        }
        entries.push((
            group,
            Group {
                start,
                end,
                regular: !reversed.contains(&group),
            },
        ));
        start = end;
    }
    let ranges = GroupLookup::new(entries);

    // Small groups use the original predicate directly, so they need no end sort.
    let mut sorted_ends = Vec::new();
    if has_large_group {
        sorted_ends.extend(targets.iter().map(|t| t.end));
        for group in &ranges.ranges {
            if group.end - group.start > 8 {
                sort_records(&mut sorted_ends[group.start..group.end], |&end| end);
            }
        }
    }

    // Sequential rank scans amortize their query sorts over larger batches.
    // Binary searches avoid that setup for small query sets and small inputs.
    let strategy = match strategy {
        Strategy::Auto => {
            if has_large_group && groups.len() >= 4096 && groups.len() >= groups2.len() / 16 {
                Strategy::Sorted
            } else {
                Strategy::Binary
            }
        }
        other => other,
    };
    let sorted = matches!(strategy, Strategy::Sorted);
    if sorted {
        sorted_counts(
            &targets,
            &sorted_ends,
            &ranges,
            groups,
            starts,
            ends,
            slack,
            &mut counts,
        );
        return counts;
    }
    let mut exceptional = Vec::new();
    for i in 0..groups.len() {
        let Some(group) = ranges.get(&groups[i]) else {
            continue;
        };
        let limit = ends[i].saturating_add(slack);
        let slice = &targets[group.start..group.end];
        if slice.len() <= 8 {
            counts[i] = finish_count(
                slice
                    .iter()
                    .filter(|t| starts[i] < t.end && t.start < limit)
                    .count(),
            );
        } else {
            let query = Query {
                group: groups[i],
                coordinate: limit,
                index: i,
            };
            // For non-reversed targets this proves that every target subtracted
            // by the end rank is also included in the start rank. Use saturated
            // arithmetic twice: a single-slack nonempty-window test is insufficient.
            if group.regular && limit.saturating_add(slack) > starts[i] {
                let a = slice.partition_point(|t| t.start < limit);
                let b =
                    sorted_ends[group.start..group.end].partition_point(|&end| end <= starts[i]);
                counts[i] = finish_count(a - b);
            } else {
                exceptional.push(query);
            }
        }
    }

    if !exceptional.is_empty() {
        exact_counts(
            &targets,
            &sorted_ends,
            &ranges,
            starts,
            &mut exceptional,
            &mut counts,
        );
    }
    counts
}

#[allow(clippy::too_many_arguments)]
fn sorted_counts<C: GroupType, T: PositionType>(
    targets: &[Target<C, T>],
    sorted_ends: &[T],
    ranges: &GroupLookup<C>,
    groups: &[C],
    starts: &[T],
    ends: &[T],
    slack: T,
    counts: &mut [u32],
) {
    // Group queries once, avoiding a hash lookup for every individual query.
    let mut queries: Vec<_> = (0..groups.len())
        .map(|i| Query {
            group: groups[i],
            coordinate: ends[i].saturating_add(slack),
            index: i,
        })
        .collect();
    sort_records(&mut queries, |q| q.group);
    let mut exceptional = Vec::new();
    let mut ranked = Vec::new();
    let mut begin = 0;
    while begin < queries.len() {
        let group_id = queries[begin].group;
        let mut end = begin + 1;
        while end < queries.len() && queries[end].group == group_id {
            end += 1;
        }
        let Some(group) = ranges.get(&group_id) else {
            begin = end;
            continue;
        };
        let queries = &mut queries[begin..end];
        if group.end - group.start <= 8 {
            for query in queries {
                counts[query.index] = finish_count(
                    targets[group.start..group.end]
                        .iter()
                        .filter(|t| starts[query.index] < t.end && t.start < query.coordinate)
                        .count(),
                );
            }
        } else if !group.regular
            || queries
                .iter()
                .all(|q| q.coordinate.saturating_add(slack) <= starts[q.index])
        {
            exact_counts(targets, sorted_ends, ranges, starts, queries, counts);
        } else {
            sort_records(queries, |q| q.coordinate);
            let mut position = group.start;
            if queries.is_sorted_by_key(|q| starts[q.index]) {
                // Fixed-length queries commonly have both boundaries in the
                // same order. Both rank counters can then move forward together.
                let mut end_position = group.start;
                for query in queries {
                    while position < group.end && targets[position].start < query.coordinate {
                        position += 1;
                    }
                    while end_position < group.end
                        && sorted_ends[end_position] <= starts[query.index]
                    {
                        end_position += 1;
                    }
                    if query.coordinate.saturating_add(slack) > starts[query.index] {
                        counts[query.index] = finish_count(position - end_position);
                    } else {
                        exceptional.push(*query);
                    }
                }
            } else {
                // Nested/varying-length queries require independently ordered
                // starts and ends. Only this group needs temporary wide ranks.
                ranked.clear();
                for query in queries {
                    while position < group.end && targets[position].start < query.coordinate {
                        position += 1;
                    }
                    ranked.push(RankedQuery {
                        coordinate: starts[query.index],
                        index: query.index,
                        rank: position,
                    });
                }
                sort_records(&mut ranked, |q| q.coordinate);
                position = group.start;
                for query in &ranked {
                    while position < group.end && sorted_ends[position] <= query.coordinate {
                        position += 1;
                    }
                    let limit = ends[query.index].saturating_add(slack);
                    if limit.saturating_add(slack) > starts[query.index] {
                        counts[query.index] = finish_count(query.rank - position);
                    } else {
                        exceptional.push(Query {
                            group: group_id,
                            coordinate: limit,
                            index: query.index,
                        });
                    }
                }
            }
        }
        begin = end;
    }
    if !exceptional.is_empty() {
        exact_counts(
            targets,
            sorted_ends,
            ranges,
            starts,
            &mut exceptional,
            counts,
        );
    }
}

fn add(tree: &mut [usize], mut index: usize) {
    index += 1;
    while index < tree.len() {
        tree[index] += 1;
        index += index & index.wrapping_neg();
    }
}

fn prefix(tree: &[usize], mut end: usize) -> usize {
    let mut count = 0;
    while end > 0 {
        count += tree[end];
        end &= end - 1;
    }
    count
}

/// Offline dominance sweep for queries where rank subtraction is not valid.
/// Inserting targets with start < limit and counting inserted adjusted ends >
/// query.start implements the original predicate even for reversed intervals.
fn exact_counts<C: GroupType, T: PositionType>(
    targets: &[Target<C, T>],
    sorted_ends: &[T],
    groups: &GroupLookup<C>,
    starts: &[T],
    queries: &mut [Query<C, T>],
    counts: &mut [u32],
) {
    sort_records(queries, |q| q.group);
    let mut tree = Vec::new();
    let mut begin = 0;
    while begin < queries.len() {
        let group_id = queries[begin].group;
        let mut end = begin + 1;
        while end < queries.len() && queries[end].group == group_id {
            end += 1;
        }
        let Some(group) = groups.get(&group_id) else {
            begin = end;
            continue;
        };
        let coordinates = &sorted_ends[group.start..group.end];
        let queries = &mut queries[begin..end];
        sort_records(queries, |q| q.coordinate);
        tree.clear();
        tree.resize(coordinates.len() + 1, 0);
        let mut position = group.start;
        for query in queries {
            while position < group.end && targets[position].start < query.coordinate {
                let rank = coordinates.partition_point(|&end| end < targets[position].end);
                add(&mut tree, rank);
                position += 1;
            }
            let rank = coordinates.partition_point(|&end| end <= starts[query.index]);
            counts[query.index] = finish_count(position - group.start - prefix(&tree, rank));
        }
        begin = end;
    }
}

#[cfg(test)]
mod tests {
    use super::{count_overlaps, finish_count, Strategy};
    use crate::ruranges_structs::PositionType;

    fn check<T: PositionType>(left: &[(u32, T, T)], right: &[(u32, T, T)], slack: T) {
        let groups: Vec<_> = left.iter().map(|r| r.0).collect();
        let starts: Vec<_> = left.iter().map(|r| r.1).collect();
        let ends: Vec<_> = left.iter().map(|r| r.2).collect();
        let groups2: Vec<_> = right.iter().map(|r| r.0).collect();
        let starts2: Vec<_> = right.iter().map(|r| r.1).collect();
        let ends2: Vec<_> = right.iter().map(|r| r.2).collect();
        // Independent Cartesian-product oracle: no sorting, endpoint ranks,
        // interval merging, or Fenwick bookkeeping.
        let expected: Vec<_> = left
            .iter()
            .map(|q| {
                right
                    .iter()
                    .filter(|t| {
                        q.0 == t.0
                            && q.1 < t.2.saturating_add(slack)
                            && t.1 < q.2.saturating_add(slack)
                    })
                    .count() as u32
            })
            .collect();
        for strategy in [Strategy::Auto, Strategy::Binary, Strategy::Sorted] {
            assert_eq!(
                count_overlaps(
                    &groups, &starts, &ends, &groups2, &starts2, &ends2, slack, strategy
                ),
                expected
            );
        }
    }

    #[test]
    fn counts_preserve_duplicates_touching_nested_queries_and_missing_groups() {
        let left = [(0, 0i64, 100), (0, 5, 8), (0, 30, 40), (2, 5, 20)];
        let right = [
            (0, -1, 200),
            (0, 0, 5),
            (0, 5, 7),
            (0, 5, 7),
            (0, 8, 9),
            (1, 0, 100),
        ];
        check(&left, &right, 0);
        // Repeat enough to exercise the rank and dominance paths, not just the
        // tiny-input shortcut; duplicates must be counted independently.
        check(&left.repeat(8), &right.repeat(8), 0);
    }

    #[test]
    fn counts_zero_length_negative_slack_reversed_and_saturated_inputs() {
        let left = [
            (0, 5i64, 5),
            (0, 0, 10),
            (0, 15, 5),
            (0, i64::MAX, i64::MAX),
            (0, i64::MIN, i64::MIN + 1),
            (0, i64::MAX - 1, i64::MAX),
        ];
        let ordinary = [(0, 5, 5), (0, 5, 6), (0, 0, 10), (0, i64::MIN, i64::MAX)];
        for slack in [i64::MIN, -6, -1, 0, 1, i64::MAX] {
            check(&left.repeat(4), &ordinary.repeat(4), slack);
            let mut reversed = ordinary.repeat(4);
            reversed.push((0, 8, 2));
            check(&left.repeat(4), &reversed, slack);
        }
    }

    #[test]
    fn randomized_counts_match_brute_force_for_all_strategies() {
        fn next(state: &mut u64) -> u64 {
            *state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            *state
        }
        fn intervals(state: &mut u64, size: usize, normal: bool) -> Vec<(u32, i64, i64)> {
            const POS: [i64; 14] = [
                i64::MIN,
                i64::MIN + 1,
                -50,
                -20,
                -1,
                0,
                0,
                1,
                5,
                10,
                20,
                50,
                i64::MAX - 1,
                i64::MAX,
            ];
            (0..size)
                .map(|_| {
                    let g = (next(state) % 5) as u32;
                    let a = POS[(next(state) >> 16) as usize % POS.len()];
                    let b = POS[(next(state) >> 16) as usize % POS.len()];
                    (
                        g,
                        if normal { a.min(b) } else { a },
                        if normal { a.max(b) } else { b },
                    )
                })
                .collect()
        }
        let mut state = 781;
        for i in 0..150 {
            let left = intervals(&mut state, i % 67, false);
            let right = intervals(&mut state, (i * 7) % 99, i % 2 == 0);
            for slack in [i64::MIN, -20, -1, 0, 1, 20, i64::MAX] {
                check(&left, &right, slack);
            }
        }
    }

    #[test]
    fn counts_small_groups_and_radix_sized_collections() {
        let left: Vec<_> = (0..600)
            .map(|i| ((i % 75) as u32, i as i64 % 19, i as i64 % 19 + 3))
            .collect();
        let right: Vec<_> = (0..600)
            .map(|i| ((i % 75) as u32, i as i64 % 13, i as i64 % 13 + 2))
            .collect();
        check(&left, &right, 0);
        let left: Vec<_> = (0..600)
            .map(|i| ((i % 2) as u32, i as i64 % 19, i as i64 % 19 + 3))
            .collect();
        let right: Vec<_> = (0..600)
            .map(|i| ((i % 2) as u32, i as i64 % 13, i as i64 % 13 + 2))
            .collect();
        check(&left, &right, 0);
        check(&left.repeat(8), &right, 0);
    }

    #[test]
    fn counts_i32_saturation_and_small_side_shortcuts() {
        let intervals = [
            (0, i32::MIN, i32::MIN + 1),
            (0, 0, 10),
            (0, 5, 5),
            (0, i32::MAX - 1, i32::MAX),
            (0, 8, 2),
        ];
        for slack in [i32::MIN, -6, -1, 0, 1, i32::MAX] {
            check(&intervals.repeat(5), &intervals.repeat(5), slack);
            check(&intervals.repeat(5), &intervals[..4].repeat(5), slack);
            check(&intervals[..1], &intervals.repeat(10), slack);
            check(&intervals.repeat(10), &intervals[..1], slack);
        }
    }

    #[test]
    fn counts_support_sparse_large_and_signed_group_labels() {
        let left: Vec<_> = (0..40)
            .map(|i| (u32::MAX - (i % 2), i as i64, i as i64 + 3))
            .collect();
        let right: Vec<_> = (0..40)
            .map(|i| (u32::MAX - (i % 3), i as i64, i as i64 + 2))
            .collect();
        check(&left, &right, 0);
        let groups = vec![-7i32; 20];
        let starts: Vec<_> = (0..20i32).collect();
        let ends: Vec<_> = starts.iter().map(|&s| s + 3).collect();
        let expected: Vec<_> = starts
            .iter()
            .zip(&ends)
            .map(|(&s, &e)| {
                starts
                    .iter()
                    .zip(&ends)
                    .filter(|&(&t, &u)| s < u && t < e)
                    .count() as u32
            })
            .collect();
        for strategy in [Strategy::Auto, Strategy::Binary, Strategy::Sorted] {
            assert_eq!(
                count_overlaps(&groups, &starts, &ends, &groups, &starts, &ends, 0, strategy),
                expected
            );
        }
    }

    #[test]
    fn final_counts_saturate_without_narrowing_intermediate_ranks() {
        assert_eq!(finish_count(u32::MAX as usize), u32::MAX);
        assert_eq!(
            finish_count((u32::MAX as usize).saturating_add(100)),
            u32::MAX
        );
    }
}
