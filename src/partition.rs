//! Group-major partitioning by parallel counting sort.
//!
//! Every sweep kernel starts by making each group (chromosome) contiguous. The
//! radix pass that used to do that runs on one thread over the whole array, and
//! it is what stops the kernels scaling: the per-group work that follows is
//! independent across groups, but it is preceded by a serial pass over the same
//! data.
//!
//! Callers pass group ids that a frame-level factorisation already made dense
//! (`0..n_groups`), so a row's group can index a histogram directly. Two
//! parallel passes then replace the radix sort: one counts each block's rows per
//! group, one writes each row straight to its final slot. Blocks are visited in
//! order and a block's rows are written in order, so the result is stable — the
//! same order the stable radix produced, element for element.
//!
//! The plan declines (returns `None`) when the group ids are too spread out to
//! index or the histogram would be too large; the caller then keeps the radix
//! path.

use std::mem::MaybeUninit;
use std::ops::Range;

use rayon::prelude::*;

use crate::ruranges_structs::GroupType;

/// Largest number of distinct group slots the counting partition will index.
const MAX_GROUP_SPAN: i64 = 16 << 20;

/// Ceiling on the block-by-group histogram, which is `blocks * groups * 4` bytes.
const MAX_HISTOGRAM_BYTES: usize = 512 << 20;

/// Rows per parallel block: enough to amortise a task, small enough that twelve
/// threads still get several blocks each at 10^6 rows.
pub(crate) const BLOCK_ROWS: usize = 128 * 1024;

/// Below this many rows the serial radix passes already beat the rayon fan-out,
/// so the parallel builders defer to their serial counterparts.
pub(crate) const PARALLEL_THRESHOLD: usize = 100_000;

/// Rows a batch of group runs should cover before it becomes its own rayon task.
///
/// A proteome frame has ~634k groups; one task per group would spend more time
/// in the scheduler than in the sweep. Batching consecutive runs up to this many
/// rows keeps both extremes — 24 chromosomes and 634k proteins — at a sensible
/// task count.
pub(crate) const BATCH_ROWS: usize = 64 * 1024;

/// Half-open row ranges to process in parallel, `count` of them.
fn blocks_of(len: usize, count: usize) -> Vec<(usize, usize)> {
    let count = count.max(1);
    let rows = len.div_ceil(count);
    let mut out = Vec::with_capacity(count);
    let mut start = 0;
    while start < len {
        let end = (start + rows).min(len);
        out.push((start, end));
        start = end;
    }
    if out.is_empty() {
        out.push((0, 0));
    }
    out
}

/// Half-open row ranges to process in parallel.
pub(crate) fn blocks(len: usize) -> Vec<(usize, usize)> {
    blocks_of(len, len.div_ceil(BLOCK_ROWS))
}

/// Group consecutive runs into work units of at least `min_rows` rows.
///
/// `weight` reports a run's cost in rows. Returns index ranges into `runs`; a
/// single run heavier than `min_rows` becomes a batch of its own.
pub(crate) fn batch_runs<R>(runs: &[R], min_rows: usize, weight: impl Fn(&R) -> usize) -> Vec<Range<usize>> {
    let mut out = Vec::new();
    let mut start = 0usize;
    let mut carried = 0usize;
    for (index, run) in runs.iter().enumerate() {
        carried += weight(run);
        if carried >= min_rows {
            out.push(start..index + 1);
            start = index + 1;
            carried = 0;
        }
    }
    if start < runs.len() {
        out.push(start..runs.len());
    }
    out
}

/// A raw pointer that may be shared across rayon tasks.
///
/// Sound only under the invariant the scatter upholds: each destination index is
/// written by exactly one task, exactly once, because the write cursors come from
/// an exclusive prefix sum over disjoint (block, group) counts.
struct SendPtr<T>(*mut T);

// Derived `Copy` would carry a `T: Copy` bound; the pointer itself is always Copy.
impl<T> Clone for SendPtr<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for SendPtr<T> {}

impl<T> SendPtr<T> {
    /// The wrapped pointer. A method rather than a field read: closures capture
    /// disjoint fields in edition 2021, and capturing the bare pointer would make
    /// the closure non-`Sync`.
    fn get(self) -> *mut T {
        self.0
    }
}

// SAFETY: dereferenced only at indices reserved for the writing task, so no two
// tasks ever touch the same element.
unsafe impl<T> Send for SendPtr<T> {}
unsafe impl<T> Sync for SendPtr<T> {}

/// Write cursors for a planned group-major partition.
pub(crate) struct GroupLayout {
    /// `cursors[block * groups + group]`: next output index for that pair.
    cursors: Vec<u32>,
    /// Number of group slots, `max - min + 1` over the group ids.
    groups: usize,
    /// Smallest group id seen, subtracted to index the histogram.
    offset: i64,
    /// Total rows the plan covers.
    len: usize,
}

impl GroupLayout {
    /// Half-open offsets of each non-empty group's run in the partitioned array.
    ///
    /// Valid only before [`scatter_group_major`] consumes the layout, which
    /// advances the cursors.
    pub(crate) fn group_runs(&self) -> Vec<(usize, usize)> {
        let mut runs = Vec::new();
        for group in 0..self.groups {
            let start = self.cursors[group] as usize;
            let end = if group + 1 < self.groups {
                self.cursors[group + 1] as usize
            } else {
                self.len
            };
            if end > start {
                runs.push((start, end));
            }
        }
        runs
    }
}

/// Plan a group-major partition of `chrs`, or `None` if counting will not pay.
///
/// Returns the layout together with the block list the scatter must use: the
/// number of blocks is chosen here, because the histogram is one counter per
/// (block, group) and a frame with many groups can only afford a few blocks.
pub(crate) fn plan_group_major<C: GroupType + Send + Sync>(
    chrs: &[C],
    slots_per_row: usize,
) -> Option<(GroupLayout, Vec<(usize, usize)>)> {
    if chrs.is_empty() {
        return None;
    }

    let scan_blocks = blocks(chrs.len());
    let (min, max) = scan_blocks
        .par_iter()
        .map(|&(start, end)| {
            let mut min = i64::MAX;
            let mut max = i64::MIN;
            for value in &chrs[start..end] {
                let value = value.to_i64().unwrap_or(i64::MIN);
                min = min.min(value);
                max = max.max(value);
            }
            (min, max)
        })
        .reduce(|| (i64::MAX, i64::MIN), |a, b| (a.0.min(b.0), a.1.max(b.1)));

    if min == i64::MIN {
        return None; // a group id that does not fit i64
    }
    let span = max.checked_sub(min)?.checked_add(1)?;
    if span <= 0 || span > MAX_GROUP_SPAN {
        return None;
    }
    let groups = span as usize;
    // Cursors are u32, so the whole output has to be indexable by one.
    if chrs.len().checked_mul(slots_per_row)? > u32::MAX as usize {
        return None;
    }

    // As many blocks as the histogram budget allows, never more than one per
    // BLOCK_ROWS rows. With few groups that is the fine blocking; with hundreds of
    // thousands of groups it is a few hundred larger blocks, which still gives
    // every thread work but keeps the histogram to a few hundred MiB.
    let affordable = MAX_HISTOGRAM_BYTES / groups.max(1) / 4;
    if affordable < rayon::current_num_threads() {
        return None; // not even one block per thread: keep the radix path
    }
    let block_list = blocks_of(chrs.len(), affordable.min(scan_blocks.len()));

    // Per-block histograms, then an exclusive prefix sum in group-major order:
    // every row of group g precedes every row of group g+1, and within a group
    // the blocks keep their order, which is what makes the partition stable.
    let mut cursors = vec![0_u32; block_list.len() * groups];
    cursors
        .par_chunks_mut(groups)
        .zip(block_list.par_iter())
        .for_each(|(histogram, &(start, end))| {
            for value in &chrs[start..end] {
                let slot = value.to_i64().expect("group id fits i64") - min;
                histogram[slot as usize] += slots_per_row as u32;
            }
        });

    let mut running = 0_u32;
    for group in 0..groups {
        for block in 0..block_list.len() {
            let slot = block * groups + group;
            let count = cursors[slot];
            cursors[slot] = running;
            running += count;
        }
    }

    Some((
        GroupLayout {
            cursors,
            groups,
            offset: min,
            len: chrs.len() * slots_per_row,
        },
        block_list,
    ))
}

/// Write `slots_per_row` elements per input row into group-major order.
///
/// `build` receives the row index and that row's destination slots, which it must
/// fill. Rows of the same group land contiguously, in input order.
pub(crate) fn scatter_group_major<C, E, F>(
    chrs: &[C],
    block_list: &[(usize, usize)],
    mut layout: GroupLayout,
    slots_per_row: usize,
    build: F,
) -> Vec<E>
where
    C: GroupType + Send + Sync,
    E: Send,
    F: Fn(usize, &mut [MaybeUninit<E>]) + Send + Sync,
{
    let len = layout.len;
    let groups = layout.groups;
    let offset = layout.offset;
    let mut out: Vec<E> = Vec::with_capacity(len);
    let base = SendPtr(out.spare_capacity_mut().as_mut_ptr());

    layout
        .cursors
        .par_chunks_mut(groups)
        .zip(block_list.par_iter())
        .for_each(|(cursors, &(start, end))| {
            for (row, chr) in (start..end).zip(&chrs[start..end]) {
                let group = (chr.to_i64().expect("group id fits i64") - offset) as usize;
                let position = cursors[group] as usize;
                cursors[group] += slots_per_row as u32;
                // SAFETY: the prefix sum reserved `slots_per_row` slots at
                // `position` for this (block, group) pair, and no other task or
                // row is given the same range.
                let slots = unsafe {
                    std::slice::from_raw_parts_mut(base.get().add(position), slots_per_row)
                };
                build(row, slots);
            }
        });

    // SAFETY: the prefix sum reserved exactly `len` slots and every row filled
    // its own, so 0..len is initialised.
    unsafe { out.set_len(len) };
    out
}

/// Split `values` into the disjoint sub-slices named by `runs`, which must be
/// ascending and non-overlapping. Gaps between runs are dropped.
pub(crate) fn split_runs_mut<'a, E>(
    values: &'a mut [E],
    runs: &[(usize, usize)],
) -> Vec<&'a mut [E]> {
    let mut out = Vec::with_capacity(runs.len());
    let mut rest = values;
    let mut consumed = 0usize;
    for &(start, end) in runs {
        let (_, tail) = rest.split_at_mut(start - consumed);
        let (chunk, tail) = tail.split_at_mut(end - start);
        out.push(chunk);
        rest = tail;
        consumed = end;
    }
    out
}

/// Run `body` on each group's slice, batching consecutive groups into rayon
/// tasks so that a frame with hundreds of thousands of groups does not create
/// one task per group.
///
/// `runs` must be ascending and contiguous, as both [`GroupLayout::group_runs`]
/// and [`runs_by`] produce.
pub(crate) fn par_groups_mut<E, F>(values: &mut [E], runs: &[(usize, usize)], body: F)
where
    E: Send,
    F: Fn(&mut [E]) + Send + Sync,
{
    if runs.is_empty() {
        return;
    }
    let batches = batch_runs(runs, BATCH_ROWS, |&(start, end)| end - start);
    let spans: Vec<(usize, usize)> = batches
        .iter()
        .map(|batch| (runs[batch.start].0, runs[batch.end - 1].1))
        .collect();

    split_runs_mut(values, &spans)
        .into_par_iter()
        .zip(batches.par_iter())
        .for_each(|(chunk, batch)| {
            let base = runs[batch.start].0;
            for &(start, end) in &runs[batch.clone()] {
                body(&mut chunk[start - base..end - base]);
            }
        });
}

/// Half-open offsets of each contiguous run of equal key in a slice already
/// grouped by that key.
pub(crate) fn runs_by<E, K: PartialEq>(values: &[E], key: impl Fn(&E) -> K) -> Vec<(usize, usize)> {
    let mut runs = Vec::new();
    if values.is_empty() {
        return runs;
    }
    let mut start = 0usize;
    let mut previous = key(&values[0]);
    for (index, value) in values.iter().enumerate().skip(1) {
        let current = key(value);
        if current != previous {
            runs.push((start, index));
            start = index;
            previous = current;
        }
    }
    runs.push((start, values.len()));
    runs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counting_partition_is_stable_and_group_major() {
        let chrs: Vec<u32> = (0..10_000u32).map(|i| i % 7).collect();
        let (layout, block_list) = plan_group_major(&chrs, 1).expect("plan");
        let runs = layout.group_runs();
        let out: Vec<(u32, u32)> =
            scatter_group_major(&chrs, &block_list, layout, 1, |row, slots| {
                slots[0].write((chrs[row], row as u32));
            });

        assert_eq!(out.len(), chrs.len());
        // Group-major, and input order kept inside each group.
        for &(start, end) in &runs {
            let group = out[start].0;
            for window in out[start..end].windows(2) {
                assert_eq!(window[0].0, group);
                assert!(window[0].1 < window[1].1);
            }
        }
        assert_eq!(runs.len(), 7);
        assert!(runs.windows(2).all(|w| out[w[0].0].0 < out[w[1].0].0));
    }

    #[test]
    fn two_slots_per_row_keeps_the_pair_adjacent() {
        let chrs: Vec<u32> = vec![2, 0, 1, 0, 2];
        let (layout, block_list) = plan_group_major(&chrs, 2).expect("plan");
        let out: Vec<(u32, u32, bool)> =
            scatter_group_major(&chrs, &block_list, layout, 2, |row, slots| {
                slots[0].write((chrs[row], row as u32, true));
                slots[1].write((chrs[row], row as u32, false));
            });

        assert_eq!(out.len(), 10);
        assert!(out.windows(2).all(|w| w[0].0 <= w[1].0));
        for pair in out.chunks(2) {
            assert_eq!(pair[0].1, pair[1].1);
            assert!(pair[0].2 && !pair[1].2);
        }
    }

    #[test]
    fn batching_keeps_run_order_and_covers_everything() {
        let runs: Vec<(usize, usize)> = (0..10).map(|i| (i * 3, i * 3 + 3)).collect();
        let batches = batch_runs(&runs, 10, |&(s, e)| e - s);
        assert_eq!(batches.first().unwrap().start, 0);
        assert_eq!(batches.last().unwrap().end, runs.len());
        for pair in batches.windows(2) {
            assert_eq!(pair[0].end, pair[1].start);
        }
    }

    #[test]
    fn plan_declines_on_a_span_it_cannot_index() {
        let chrs: Vec<i64> = vec![0, i64::MAX / 2];
        assert!(plan_group_major(&chrs, 1).is_none());
    }
}
