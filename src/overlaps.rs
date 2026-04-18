use std::collections::HashMap;
use std::str::FromStr;
use std::sync::{Arc, Mutex, OnceLock};

use radsort::sort_by_key;
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};

use crate::ruranges_structs::{GroupType, OverlapPair, OverlapType, PositionType};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum OverlapParallelMode {
    Auto,
    Off,
    On,
}

impl Default for OverlapParallelMode {
    fn default() -> Self {
        Self::Auto
    }
}

impl FromStr for OverlapParallelMode {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "auto" => Ok(Self::Auto),
            "off" => Ok(Self::Off),
            "on" => Ok(Self::On),
            other => Err(format!("invalid overlap parallel mode: {other}")),
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct OverlapParallelConfig {
    pub mode: OverlapParallelMode,
    pub threads: Option<usize>,
    pub min_total_intervals: usize,
}

impl OverlapParallelConfig {
    pub const fn new(
        mode: OverlapParallelMode,
        threads: Option<usize>,
        min_total_intervals: usize,
    ) -> Self {
        Self {
            mode,
            threads,
            min_total_intervals,
        }
    }

    pub const fn serial() -> Self {
        Self::new(OverlapParallelMode::Off, Some(1), 0)
    }

    pub const fn parallel(threads: Option<usize>) -> Self {
        Self::new(OverlapParallelMode::On, threads, 0)
    }

    fn normalized_threads(self) -> Option<usize> {
        self.threads.filter(|threads| *threads > 0)
    }
}

impl Default for OverlapParallelConfig {
    fn default() -> Self {
        Self::new(OverlapParallelMode::Auto, None, 0)
    }
}

fn thread_pool_cache() -> &'static Mutex<HashMap<usize, Arc<ThreadPool>>> {
    static THREAD_POOLS: OnceLock<Mutex<HashMap<usize, Arc<ThreadPool>>>> = OnceLock::new();
    THREAD_POOLS.get_or_init(|| Mutex::new(HashMap::new()))
}

fn cached_thread_pool(threads: usize) -> Arc<ThreadPool> {
    let cache = thread_pool_cache();
    let mut pools = cache.lock().expect("thread pool cache mutex poisoned");

    pools
        .entry(threads)
        .or_insert_with(|| {
            Arc::new(
                ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("failed to build overlap thread pool"),
            )
        })
        .clone()
}

fn run_with_optional_thread_pool<R: Send>(
    threads: Option<usize>,
    task: impl FnOnce() -> R + Send,
) -> R {
    match threads.filter(|threads| *threads > 1) {
        Some(threads) => cached_thread_pool(threads).install(task),
        None => task(),
    }
}

fn should_parallelize(
    config: OverlapParallelConfig,
    total_intervals: usize,
    group_pairs_len: usize,
    n_workers: usize,
) -> bool {
    if n_workers <= 1 || group_pairs_len <= 1 {
        return false;
    }

    match config.mode {
        OverlapParallelMode::Off => false,
        OverlapParallelMode::On => true,
        OverlapParallelMode::Auto => total_intervals >= config.min_total_intervals,
    }
}

#[derive(Copy, Clone, Debug)]
struct IntervalRecord<C: GroupType, T: PositionType> {
    group: C,
    start: T,
    end: T,
    idx: u32,
}

fn sorted_records<C: GroupType, T: PositionType>(
    groups: &[C],
    starts: &[T],
    ends: &[T],
) -> Vec<IntervalRecord<C, T>> {
    debug_assert_eq!(groups.len(), starts.len(), "groups/starts length mismatch");
    debug_assert_eq!(groups.len(), ends.len(), "groups/ends length mismatch");

    let mut records = Vec::with_capacity(groups.len());

    for idx in 0..groups.len() {
        records.push(IntervalRecord {
            group: groups[idx],
            start: starts[idx],
            end: ends[idx],
            idx: idx as u32,
        });
    }

    sort_by_key(&mut records, |r| (r.group, r.start, r.end, r.idx));

    records
}

#[derive(Copy, Clone, Debug)]
struct GroupSpan<C: GroupType> {
    group: C,
    start: usize,
    end: usize,
}

#[derive(Copy, Clone, Debug)]
struct GroupPair {
    index: usize,
    left_start: usize,
    left_end: usize,
    right_start: usize,
    right_end: usize,
}

#[derive(Default)]
struct OverlapIndices {
    left: Vec<u32>,
    right: Vec<u32>,
}

impl OverlapIndices {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            left: Vec::with_capacity(capacity),
            right: Vec::with_capacity(capacity),
        }
    }

    fn push(&mut self, left: u32, right: u32) {
        self.left.push(left);
        self.right.push(right);
    }

    fn len(&self) -> usize {
        self.left.len()
    }

    fn extend(&mut self, other: Self) {
        self.left.extend(other.left);
        self.right.extend(other.right);
    }
}

#[inline(always)]
fn overlaps_with_slack<T: PositionType>(
    query_start: T,
    query_end: T,
    target_start: T,
    target_end: T,
    slack: T,
) -> bool {
    query_start < target_end.saturating_add(slack) && target_start < query_end.saturating_add(slack)
}

#[inline(always)]
fn query_contained_in_target_with_slack<T: PositionType>(
    query_start: T,
    query_end: T,
    target_start: T,
    target_end: T,
    slack: T,
) -> bool {
    let query_start_slack = query_start.saturating_sub(slack);
    let query_end_slack = query_end.saturating_add(slack);
    query_start_slack >= target_start && query_end_slack <= target_end
}

fn group_spans<C: GroupType, T: PositionType>(
    records: &[IntervalRecord<C, T>],
) -> Vec<GroupSpan<C>> {
    let mut spans = Vec::new();
    let mut i = 0usize;

    while i < records.len() {
        let group = records[i].group;
        let start = i;
        i += 1;

        while i < records.len() && records[i].group == group {
            i += 1;
        }

        spans.push(GroupSpan {
            group,
            start,
            end: i,
        });
    }

    spans
}

fn matching_group_pairs<C: GroupType>(
    left: &[GroupSpan<C>],
    right: &[GroupSpan<C>],
) -> Vec<GroupPair> {
    let mut pairs = Vec::new();
    let mut i = 0usize;
    let mut j = 0usize;
    let mut index = 0usize;

    while i < left.len() && j < right.len() {
        let g1 = left[i].group;
        let g2 = right[j].group;

        if g1 < g2 {
            i += 1;
            continue;
        }

        if g2 < g1 {
            j += 1;
            continue;
        }

        pairs.push(GroupPair {
            index,
            left_start: left[i].start,
            left_end: left[i].end,
            right_start: right[j].start,
            right_end: right[j].end,
        });
        index += 1;

        i += 1;
        j += 1;
    }

    pairs
}

fn group_pair_weight(pair: &GroupPair) -> usize {
    (pair.left_end - pair.left_start) + (pair.right_end - pair.right_start)
}

fn build_group_batches(group_pairs: &[GroupPair], n_workers: usize) -> Vec<Vec<GroupPair>> {
    if group_pairs.is_empty() {
        return Vec::new();
    }

    let target_batches = group_pairs.len().min(n_workers.max(1) * 4).max(1);

    let mut weighted = group_pairs.to_vec();
    weighted.sort_by(|a, b| {
        group_pair_weight(b)
            .cmp(&group_pair_weight(a))
            .then_with(|| a.index.cmp(&b.index))
    });

    let mut batches = (0..target_batches)
        .map(|_| (0usize, Vec::<GroupPair>::new()))
        .collect::<Vec<_>>();

    for pair in weighted {
        let batch = batches
            .iter_mut()
            .min_by(|(weight_a, pairs_a), (weight_b, pairs_b)| {
                weight_a
                    .cmp(weight_b)
                    .then_with(|| pairs_a.len().cmp(&pairs_b.len()))
            })
            .expect("at least one batch exists");

        batch.0 += group_pair_weight(&pair);
        batch.1.push(pair);
    }

    let mut batches = batches
        .into_iter()
        .filter_map(|(_, mut pairs)| {
            if pairs.is_empty() {
                return None;
            }

            pairs.sort_by_key(|pair| pair.index);
            Some(pairs)
        })
        .collect::<Vec<_>>();

    batches.sort_by_key(|pairs| pairs[0].index);
    batches
}

trait SweepRecord<T: PositionType>: Copy {
    fn start(self) -> T;
    fn end(self) -> T;
    fn idx(self) -> u32;
}

impl<C: GroupType, T: PositionType> SweepRecord<T> for IntervalRecord<C, T> {
    #[inline(always)]
    fn start(self) -> T {
        self.start
    }

    #[inline(always)]
    fn end(self) -> T {
        self.end
    }

    #[inline(always)]
    fn idx(self) -> u32 {
        self.idx
    }
}

fn visit_overlap_pairs_for_group<R: SweepRecord<T>, T: PositionType>(
    left: &[R],
    right: &[R],
    slack: T,
    overlap_type: OverlapType,
    contained: bool,
    emit: &mut dyn FnMut(u32, u32),
) {
    if left.is_empty() || right.is_empty() {
        return;
    }

    let mut active: Vec<usize> = Vec::new();
    let mut active_head: usize = 0;
    let mut jr = 0usize;

    for query in left.iter().copied() {
        let query_start = query.start();
        let query_end = query.end();
        let query_idx = query.idx();

        if contained {
            let query_start_slack = query_start.saturating_sub(slack);

            while jr < right.len() && right[jr].start() <= query_start_slack {
                active.push(jr);
                jr += 1;
            }

            while active_head < active.len() {
                let r = active[active_head];
                if right[r].end() <= query_start_slack {
                    active_head += 1;
                } else {
                    break;
                }
            }

            if active_head > 0 && active_head * 2 >= active.len() {
                active.drain(0..active_head);
                active_head = 0;
            }

            match overlap_type {
                OverlapType::All => {
                    for idx in active_head..active.len() {
                        let target = right[active[idx]];
                        if !query_contained_in_target_with_slack(
                            query_start,
                            query_end,
                            target.start(),
                            target.end(),
                            slack,
                        ) {
                            continue;
                        }

                        emit(query_idx, target.idx());
                    }
                }
                OverlapType::First => {
                    for idx in active_head..active.len() {
                        let target = right[active[idx]];
                        if !query_contained_in_target_with_slack(
                            query_start,
                            query_end,
                            target.start(),
                            target.end(),
                            slack,
                        ) {
                            continue;
                        }

                        emit(query_idx, target.idx());
                        break;
                    }
                }
                OverlapType::Last => {
                    let mut last_target: Option<R> = None;

                    for idx in active_head..active.len() {
                        let target = right[active[idx]];
                        if !query_contained_in_target_with_slack(
                            query_start,
                            query_end,
                            target.start(),
                            target.end(),
                            slack,
                        ) {
                            continue;
                        }

                        last_target = Some(target);
                    }

                    if let Some(target) = last_target {
                        emit(query_idx, target.idx());
                    }
                }
            }

            continue;
        }

        let query_end_slack = query_end.saturating_add(slack);

        while jr < right.len() && right[jr].start() < query_end_slack {
            active.push(jr);
            jr += 1;
        }

        while active_head < active.len() {
            let r = active[active_head];
            if right[r].end().saturating_add(slack) <= query_start {
                active_head += 1;
            } else {
                break;
            }
        }

        if active_head > 0 && active_head * 2 >= active.len() {
            active.drain(0..active_head);
            active_head = 0;
        }

        match overlap_type {
            OverlapType::All => {
                for idx in active_head..active.len() {
                    let target = right[active[idx]];

                    if !overlaps_with_slack(
                        query_start,
                        query_end,
                        target.start(),
                        target.end(),
                        slack,
                    ) {
                        continue;
                    }

                    emit(query_idx, target.idx());
                }
            }
            OverlapType::First => {
                for idx in active_head..active.len() {
                    let target = right[active[idx]];

                    if !overlaps_with_slack(
                        query_start,
                        query_end,
                        target.start(),
                        target.end(),
                        slack,
                    ) {
                        continue;
                    }

                    emit(query_idx, target.idx());
                    break;
                }
            }
            OverlapType::Last => {
                let mut last_target: Option<R> = None;

                for idx in active_head..active.len() {
                    let target = right[active[idx]];

                    if !overlaps_with_slack(
                        query_start,
                        query_end,
                        target.start(),
                        target.end(),
                        slack,
                    ) {
                        continue;
                    }

                    last_target = Some(target);
                }

                if let Some(target) = last_target {
                    emit(query_idx, target.idx());
                }
            }
        }
    }
}

fn collect_overlap_pairs_for_group<R: SweepRecord<T>, T: PositionType>(
    left: &[R],
    right: &[R],
    slack: T,
    overlap_type: OverlapType,
    contained: bool,
) -> Vec<OverlapPair> {
    let mut pairs = Vec::new();
    let mut emit = |idx: u32, idx2: u32| pairs.push(OverlapPair { idx, idx2 });
    visit_overlap_pairs_for_group(left, right, slack, overlap_type, contained, &mut emit);
    pairs
}

fn collect_overlap_indices_for_group<R: SweepRecord<T>, T: PositionType>(
    left: &[R],
    right: &[R],
    slack: T,
    overlap_type: OverlapType,
    contained: bool,
) -> OverlapIndices {
    let mut pairs = OverlapIndices::with_capacity(left.len().saturating_add(right.len()));
    let mut emit = |idx: u32, idx2: u32| pairs.push(idx, idx2);
    visit_overlap_pairs_for_group(left, right, slack, overlap_type, contained, &mut emit);
    pairs
}

fn collect_overlap_pairs_with_config<C: GroupType + Sync + Send, T: PositionType + Sync + Send>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
    overlap_type: OverlapType,
    contained: bool,
    parallel_config: OverlapParallelConfig,
) -> Vec<OverlapPair> {
    let left = sorted_records(chrs, starts, ends);
    let right = sorted_records(chrs2, starts2, ends2);

    if left.is_empty() || right.is_empty() {
        return Vec::new();
    }

    let left_spans = group_spans(&left);
    let right_spans = group_spans(&right);
    let group_pairs = matching_group_pairs(&left_spans, &right_spans);

    if group_pairs.is_empty() {
        return Vec::new();
    }

    let n_workers = parallel_config
        .normalized_threads()
        .unwrap_or_else(rayon::current_num_threads);
    let group_batches = build_group_batches(&group_pairs, n_workers);
    let use_parallel = should_parallelize(
        parallel_config,
        left.len().saturating_add(right.len()),
        group_pairs.len(),
        n_workers,
    );

    let grouped_pairs = if use_parallel && group_batches.len() > 1 {
        run_with_optional_thread_pool(parallel_config.normalized_threads(), move || {
            group_batches
                .par_iter()
                .map(|batch| {
                    let total_pairs_hint = batch.iter().map(group_pair_weight).sum();
                    let mut batch_pairs = Vec::with_capacity(total_pairs_hint);

                    for pair in batch {
                        batch_pairs.extend(collect_overlap_pairs_for_group(
                            &left[pair.left_start..pair.left_end],
                            &right[pair.right_start..pair.right_end],
                            slack,
                            overlap_type,
                            contained,
                        ));
                    }

                    batch_pairs
                })
                .collect::<Vec<_>>()
        })
    } else {
        group_batches
            .iter()
            .map(|batch| {
                let total_pairs_hint = batch.iter().map(group_pair_weight).sum();
                let mut batch_pairs = Vec::with_capacity(total_pairs_hint);

                for pair in batch {
                    batch_pairs.extend(collect_overlap_pairs_for_group(
                        &left[pair.left_start..pair.left_end],
                        &right[pair.right_start..pair.right_end],
                        slack,
                        overlap_type,
                        contained,
                    ));
                }

                batch_pairs
            })
            .collect::<Vec<_>>()
    };

    let total_pairs = grouped_pairs.iter().map(Vec::len).sum();
    let mut pairs = Vec::with_capacity(total_pairs);
    for group in grouped_pairs {
        pairs.extend(group);
    }

    pairs
}

fn collect_overlap_indices_with_config<
    C: GroupType + Sync + Send,
    T: PositionType + Sync + Send,
>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
    overlap_type: OverlapType,
    contained: bool,
    parallel_config: OverlapParallelConfig,
) -> OverlapIndices {
    let left = sorted_records(chrs, starts, ends);
    let right = sorted_records(chrs2, starts2, ends2);

    if left.is_empty() || right.is_empty() {
        return OverlapIndices::default();
    }

    let left_spans = group_spans(&left);
    let right_spans = group_spans(&right);
    let group_pairs = matching_group_pairs(&left_spans, &right_spans);

    if group_pairs.is_empty() {
        return OverlapIndices::default();
    }

    let n_workers = parallel_config
        .normalized_threads()
        .unwrap_or_else(rayon::current_num_threads);
    let group_batches = build_group_batches(&group_pairs, n_workers);
    let use_parallel = should_parallelize(
        parallel_config,
        left.len().saturating_add(right.len()),
        group_pairs.len(),
        n_workers,
    );

    let grouped_pairs = if use_parallel && group_batches.len() > 1 {
        run_with_optional_thread_pool(parallel_config.normalized_threads(), move || {
            group_batches
                .par_iter()
                .map(|batch| {
                    let total_pairs_hint = batch.iter().map(group_pair_weight).sum();
                    let mut batch_pairs = OverlapIndices::with_capacity(total_pairs_hint);

                    for pair in batch {
                        batch_pairs.extend(collect_overlap_indices_for_group(
                            &left[pair.left_start..pair.left_end],
                            &right[pair.right_start..pair.right_end],
                            slack,
                            overlap_type,
                            contained,
                        ));
                    }

                    batch_pairs
                })
                .collect::<Vec<_>>()
        })
    } else {
        group_batches
            .iter()
            .map(|batch| {
                let total_pairs_hint = batch.iter().map(group_pair_weight).sum();
                let mut batch_pairs = OverlapIndices::with_capacity(total_pairs_hint);

                for pair in batch {
                    batch_pairs.extend(collect_overlap_indices_for_group(
                        &left[pair.left_start..pair.left_end],
                        &right[pair.right_start..pair.right_end],
                        slack,
                        overlap_type,
                        contained,
                    ));
                }

                batch_pairs
            })
            .collect::<Vec<_>>()
    };

    let total_pairs = grouped_pairs.iter().map(OverlapIndices::len).sum();
    let mut pairs = OverlapIndices::with_capacity(total_pairs);
    for group in grouped_pairs {
        pairs.extend(group);
    }

    pairs
}

fn count_overlaps_for_group<C: GroupType, T: PositionType>(
    left: &[IntervalRecord<C, T>],
    right: &[IntervalRecord<C, T>],
    slack: T,
) -> Vec<(u32, u32)> {
    let mut counts = Vec::with_capacity(left.len());

    if left.is_empty() || right.is_empty() {
        return counts;
    }

    let mut active: Vec<usize> = Vec::new();
    let mut active_head: usize = 0;
    let mut jr = 0usize;

    for query in left.iter().copied() {
        let query_end_slack = query.end.saturating_add(slack);

        while jr < right.len() && right[jr].start < query_end_slack {
            active.push(jr);
            jr += 1;
        }

        while active_head < active.len() {
            let r = active[active_head];
            if right[r].end.saturating_add(slack) <= query.start {
                active_head += 1;
            } else {
                break;
            }
        }

        if active_head > 0 && active_head * 2 >= active.len() {
            active.drain(0..active_head);
            active_head = 0;
        }

        let mut count = 0_u32;

        for idx in active_head..active.len() {
            let target = right[active[idx]];

            if overlaps_with_slack(query.start, query.end, target.start, target.end, slack) {
                count = count.saturating_add(1);
            }
        }

        counts.push((query.idx, count));
    }

    counts
}

#[allow(clippy::too_many_arguments)]
pub fn overlaps_with_config<C: GroupType + Sync + Send, T: PositionType + Sync + Send>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
    overlap_type: &str,
    sort_output: bool,
    contained: bool,
    parallel_config: OverlapParallelConfig,
) -> (Vec<u32>, Vec<u32>) {
    let overlap_type = OverlapType::from_str(overlap_type).expect("invalid overlap_type string");

    if sort_output {
        let mut pairs = collect_overlap_pairs_with_config(
            chrs,
            starts,
            ends,
            chrs2,
            starts2,
            ends2,
            slack,
            overlap_type,
            contained,
            parallel_config,
        );

        sort_by_key(&mut pairs, |p| (p.idx, p.idx2));
        return pairs.into_iter().map(|pair| (pair.idx, pair.idx2)).unzip();
    }

    let pairs = collect_overlap_indices_with_config(
        chrs,
        starts,
        ends,
        chrs2,
        starts2,
        ends2,
        slack,
        overlap_type,
        contained,
        parallel_config,
    );

    (pairs.left, pairs.right)
}

#[allow(clippy::too_many_arguments)]
pub fn overlaps<C: GroupType + Sync + Send, T: PositionType + Sync + Send>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
    overlap_type: &str,
    sort_output: bool,
    contained: bool,
) -> (Vec<u32>, Vec<u32>) {
    overlaps_with_config(
        chrs,
        starts,
        ends,
        chrs2,
        starts2,
        ends2,
        slack,
        overlap_type,
        sort_output,
        contained,
        OverlapParallelConfig::default(),
    )
}

pub fn count_overlaps_with_config<C: GroupType + Sync + Send, T: PositionType + Sync + Send>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
    parallel_config: OverlapParallelConfig,
) -> Vec<u32> {
    let left = sorted_records(chrs, starts, ends);
    let right = sorted_records(chrs2, starts2, ends2);

    let mut counts = vec![0_u32; left.len()];

    if left.is_empty() || right.is_empty() {
        return counts;
    }

    let left_spans = group_spans(&left);
    let right_spans = group_spans(&right);
    let group_pairs = matching_group_pairs(&left_spans, &right_spans);

    if group_pairs.is_empty() {
        return counts;
    }

    let n_workers = parallel_config
        .normalized_threads()
        .unwrap_or_else(rayon::current_num_threads);
    let group_batches = build_group_batches(&group_pairs, n_workers);
    let use_parallel = should_parallelize(
        parallel_config,
        left.len().saturating_add(right.len()),
        group_pairs.len(),
        n_workers,
    );

    let grouped_counts = if use_parallel && group_batches.len() > 1 {
        run_with_optional_thread_pool(parallel_config.normalized_threads(), move || {
            group_batches
                .par_iter()
                .map(|batch| {
                    let mut batch_counts = Vec::new();

                    for pair in batch {
                        batch_counts.extend(count_overlaps_for_group(
                            &left[pair.left_start..pair.left_end],
                            &right[pair.right_start..pair.right_end],
                            slack,
                        ));
                    }

                    batch_counts
                })
                .collect::<Vec<_>>()
        })
    } else {
        group_batches
            .iter()
            .map(|batch| {
                let mut batch_counts = Vec::new();

                for pair in batch {
                    batch_counts.extend(count_overlaps_for_group(
                        &left[pair.left_start..pair.left_end],
                        &right[pair.right_start..pair.right_end],
                        slack,
                    ));
                }

                batch_counts
            })
            .collect::<Vec<_>>()
    };

    for group_counts in grouped_counts {
        for (idx, count) in group_counts {
            counts[idx as usize] = count;
        }
    }

    counts
}

pub fn count_overlaps<C: GroupType + Sync + Send, T: PositionType + Sync + Send>(
    chrs: &[C],
    starts: &[T],
    ends: &[T],
    chrs2: &[C],
    starts2: &[T],
    ends2: &[T],
    slack: T,
) -> Vec<u32> {
    count_overlaps_with_config(
        chrs,
        starts,
        ends,
        chrs2,
        starts2,
        ends2,
        slack,
        OverlapParallelConfig::default(),
    )
}

#[cfg(test)]
mod tests {
    use super::{
        count_overlaps_with_config, overlaps, overlaps_with_config, OverlapParallelConfig,
        OverlapParallelMode,
    };

    type Group = u32;
    type Pos = i64;

    #[test]
    fn overlaps_all_returns_expected_pairs() {
        let groups: [Group; 3] = [1, 1, 1];
        let starts: [Pos; 3] = [1, 10, 30];
        let ends: [Pos; 3] = [5, 20, 40];

        let groups2: [Group; 4] = [1, 1, 1, 1];
        let starts2: [Pos; 4] = [3, 11, 18, 35];
        let ends2: [Pos; 4] = [4, 12, 19, 36];

        let (idx1, idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, false,
        );

        assert_eq!(idx1, vec![0, 1, 1, 2]);
        assert_eq!(idx2, vec![0, 1, 2, 3]);
    }

    #[test]
    fn overlaps_issue_23_preserves_order_for_all_first_last() {
        let groups = vec![1_u32; 10];
        let starts: Vec<Pos> = vec![97, 981, 1000, 1227, 1409, 3889, 4398, 4815, 5004, 5047];
        let ends: Vec<Pos> = vec![1097, 1981, 2000, 2227, 2409, 4889, 5398, 5815, 6004, 6047];

        let groups2 = vec![1_u32; 10];
        let starts2: Vec<Pos> = vec![1476, 2110, 2823, 3547, 3578, 4295, 5184, 5545, 7117, 7190];
        let ends2: Vec<Pos> = vec![2476, 3110, 3823, 4547, 4578, 5295, 6184, 6545, 8117, 8190];

        let (all_idx1, all_idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, false,
        );

        assert_eq!(
            all_idx1,
            vec![1, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9]
        );
        assert_eq!(
            all_idx2,
            vec![0, 0, 0, 1, 0, 1, 3, 4, 5, 3, 4, 5, 6, 5, 6, 7, 5, 6, 7, 5, 6, 7]
        );

        let (first_idx1, first_idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "first", true, false,
        );
        assert_eq!(first_idx1, vec![1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(first_idx2, vec![0, 0, 0, 0, 3, 3, 5, 5, 5]);

        let (last_idx1, last_idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "last", true, false,
        );
        assert_eq!(last_idx1, vec![1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(last_idx2, vec![0, 0, 1, 1, 5, 6, 7, 7, 7]);
    }

    #[test]
    fn overlaps_contained_filters_to_left_contained_in_right() {
        let groups: [Group; 2] = [1, 1];
        let starts: [Pos; 2] = [10, 10];
        let ends: [Pos; 2] = [12, 30];

        let groups2: [Group; 2] = [1, 1];
        let starts2: [Pos; 2] = [9, 11];
        let ends2: [Pos; 2] = [40, 29];

        let (idx1, idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, true,
        );

        assert_eq!(idx1, vec![0, 1]);
        assert_eq!(idx2, vec![0, 0]);
    }

    #[test]
    fn overlaps_contained_includes_targets_with_same_start() {
        let groups: [Group; 2] = [1, 1];
        let starts: [Pos; 2] = [1, 6];
        let ends: [Pos; 2] = [3, 9];

        let groups2: [Group; 1] = [1];
        let starts2: [Pos; 1] = [1];
        let ends2: [Pos; 1] = [9];

        let (idx1, idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, true,
        );

        assert_eq!(idx1, vec![0, 1]);
        assert_eq!(idx2, vec![0, 0]);
    }

    #[test]
    fn overlaps_contained_with_negative_slack_matches_legacy_behavior() {
        let groups: [Group; 2] = [1, 1];
        let starts: [Pos; 2] = [1, 4];
        let ends: [Pos; 2] = [3, 9];

        let groups2: [Group; 3] = [1, 1, 1];
        let starts2: [Pos; 3] = [1, 2, 2];
        let ends2: [Pos; 3] = [10, 3, 9];

        let (idx1, idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, -2, "first", true, true,
        );

        assert_eq!(idx1, vec![0, 1]);
        assert_eq!(idx2, vec![0, 0]);
    }

    #[test]
    fn overlaps_handles_unsorted_input() {
        let groups: [Group; 3] = [1, 1, 1];
        let starts: [Pos; 3] = [10, 0, 20];
        let ends: [Pos; 3] = [12, 5, 30];

        let groups2: [Group; 2] = [1, 1];
        let starts2: [Pos; 2] = [11, 3];
        let ends2: [Pos; 2] = [13, 4];

        let (idx1, idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, false,
        );

        assert_eq!(idx1, vec![0, 1]);
        assert_eq!(idx2, vec![0, 1]);
    }

    #[test]
    fn overlaps_respects_slack_for_bookended_intervals() {
        let groups: [Group; 1] = [1];
        let starts: [Pos; 1] = [10];
        let ends: [Pos; 1] = [20];

        let groups2: [Group; 1] = [1];
        let starts2: [Pos; 1] = [20];
        let ends2: [Pos; 1] = [25];

        let (idx1_no_slack, idx2_no_slack) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, false,
        );
        assert!(idx1_no_slack.is_empty());
        assert!(idx2_no_slack.is_empty());

        let (idx1_with_slack, idx2_with_slack) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 1, "all", true, false,
        );
        assert_eq!(idx1_with_slack, vec![0]);
        assert_eq!(idx2_with_slack, vec![0]);
    }

    #[test]
    fn overlaps_handles_multiple_groups() {
        let groups: [Group; 4] = [1, 1, 2, 2];
        let starts: [Pos; 4] = [1, 10, 1, 10];
        let ends: [Pos; 4] = [5, 15, 5, 15];

        let groups2: [Group; 4] = [1, 1, 2, 2];
        let starts2: [Pos; 4] = [3, 30, 2, 14];
        let ends2: [Pos; 4] = [4, 40, 3, 20];

        let (idx1, idx2) = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, false,
        );

        assert_eq!(idx1, vec![0, 2, 3]);
        assert_eq!(idx2, vec![0, 2, 3]);
    }

    #[test]
    fn overlaps_unsorted_output_matches_sorted_output_after_manual_sort() {
        let groups: [Group; 6] = [1, 1, 2, 2, 3, 3];
        let starts: [Pos; 6] = [1, 10, 1, 10, 100, 120];
        let ends: [Pos; 6] = [5, 15, 5, 15, 110, 130];

        let groups2: [Group; 6] = [1, 1, 2, 2, 3, 3];
        let starts2: [Pos; 6] = [3, 30, 2, 14, 105, 150];
        let ends2: [Pos; 6] = [4, 40, 3, 20, 111, 160];

        let sorted = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, false,
        );
        let unsorted = overlaps(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", false, false,
        );

        let mut unsorted_pairs = unsorted.0.into_iter().zip(unsorted.1).collect::<Vec<_>>();
        unsorted_pairs.sort_unstable();

        let sorted_pairs = sorted.0.into_iter().zip(sorted.1).collect::<Vec<_>>();
        assert_eq!(unsorted_pairs, sorted_pairs);
    }

    #[test]
    fn overlap_parallel_config_preserves_results() {
        let groups: [Group; 6] = [1, 1, 2, 2, 3, 3];
        let starts: [Pos; 6] = [1, 10, 1, 10, 100, 120];
        let ends: [Pos; 6] = [5, 15, 5, 15, 110, 130];

        let groups2: [Group; 6] = [1, 1, 2, 2, 3, 3];
        let starts2: [Pos; 6] = [3, 30, 2, 14, 105, 150];
        let ends2: [Pos; 6] = [4, 40, 3, 20, 111, 160];

        let serial = OverlapParallelConfig::serial();
        let forced_parallel =
            OverlapParallelConfig::new(OverlapParallelMode::On, Some(2), usize::MAX);

        let serial_pairs = overlaps_with_config(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", true, false, serial,
        );
        let parallel_pairs = overlaps_with_config(
            &groups,
            &starts,
            &ends,
            &groups2,
            &starts2,
            &ends2,
            0,
            "all",
            true,
            false,
            forced_parallel,
        );

        assert_eq!(serial_pairs, parallel_pairs);

        let serial_counts = count_overlaps_with_config(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, serial,
        );
        let parallel_counts = count_overlaps_with_config(
            &groups,
            &starts,
            &ends,
            &groups2,
            &starts2,
            &ends2,
            0,
            forced_parallel,
        );

        assert_eq!(serial_counts, parallel_counts);
    }

    #[test]
    fn overlap_parallel_unsorted_path_preserves_pair_multiset() {
        let groups: [Group; 6] = [1, 1, 2, 2, 3, 3];
        let starts: [Pos; 6] = [1, 10, 1, 10, 100, 120];
        let ends: [Pos; 6] = [5, 15, 5, 15, 110, 130];

        let groups2: [Group; 6] = [1, 1, 2, 2, 3, 3];
        let starts2: [Pos; 6] = [3, 30, 2, 14, 105, 150];
        let ends2: [Pos; 6] = [4, 40, 3, 20, 111, 160];

        let serial = OverlapParallelConfig::serial();
        let forced_parallel =
            OverlapParallelConfig::new(OverlapParallelMode::On, Some(2), usize::MAX);

        let serial_pairs = overlaps_with_config(
            &groups, &starts, &ends, &groups2, &starts2, &ends2, 0, "all", false, false, serial,
        );
        let parallel_pairs = overlaps_with_config(
            &groups,
            &starts,
            &ends,
            &groups2,
            &starts2,
            &ends2,
            0,
            "all",
            false,
            false,
            forced_parallel,
        );

        let mut serial_pairs = serial_pairs
            .0
            .into_iter()
            .zip(serial_pairs.1)
            .collect::<Vec<_>>();
        let mut parallel_pairs = parallel_pairs
            .0
            .into_iter()
            .zip(parallel_pairs.1)
            .collect::<Vec<_>>();

        serial_pairs.sort_unstable();
        parallel_pairs.sort_unstable();

        assert_eq!(serial_pairs, parallel_pairs);
    }
}
