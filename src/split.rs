use std::sync::atomic::{AtomicU64, Ordering};

use rayon::prelude::*;

use crate::{
    partition,
    ruranges_structs::{Event, GroupType, PositionType},
    sorts,
};

pub fn sweep_line_split<G: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    chrs: &[G],
    starts: &[T],
    ends: &[T],
    slack: T,
    between: bool,
) -> (Vec<u32>, Vec<T>, Vec<T>) {
    let events = build_events_of_nonempty_rows(chrs, starts, ends, slack);

    // No atom spans two chromosomes, so each chromosome is swept on its own
    // task, batched as merge batches groups. A first pass counts each batch's
    // atoms, so that the second writes them straight to their place in the
    // output, in chromosome order: the order a single sweep emits them in.
    let runs = partition::runs_by(&events, |e| e.chr);
    let batches: Vec<&[Event<G, T>]> =
        partition::batch_runs(&runs, partition::BATCH_ROWS, |&(s, e)| e - s)
            .into_iter()
            .map(|batch| &events[runs[batch.start].0..runs[batch.end - 1].1])
            .collect();
    let counts: Vec<usize> = batches
        .par_iter()
        .map(|&batch| {
            let mut count = 0;
            sweep(batch, between, None, |_, _, _| count += 1);
            count
        })
        .collect();
    let mut slots = Vec::with_capacity(counts.len());
    let mut total = 0;
    for count in counts {
        slots.push((total, total + count));
        total += count;
    }

    let ended: Vec<AtomicU64> = (0..starts.len().div_ceil(64))
        .map(|_| AtomicU64::new(0))
        .collect();
    let mut idxs_out: Vec<u32> = Vec::with_capacity(total);
    let mut starts_out: Vec<T> = Vec::with_capacity(total);
    let mut ends_out: Vec<T> = Vec::with_capacity(total);
    partition::split_runs_mut(&mut idxs_out.spare_capacity_mut()[..total], &slots)
        .into_par_iter()
        .zip(partition::split_runs_mut(
            &mut starts_out.spare_capacity_mut()[..total],
            &slots,
        ))
        .zip(partition::split_runs_mut(
            &mut ends_out.spare_capacity_mut()[..total],
            &slots,
        ))
        .zip(batches.par_iter())
        .for_each(|(((idxs, atom_starts), atom_ends), &batch)| {
            let mut written = 0;
            sweep(batch, between, Some(&ended), |idx, start, end| {
                idxs[written].write(idx);
                atom_starts[written].write(start);
                atom_ends[written].write(end);
                written += 1;
            });
            assert_eq!(written, idxs.len(), "both passes emit the same atoms");
        });
    // SAFETY: the slots tile 0..total, and every one was written above.
    unsafe {
        idxs_out.set_len(total);
        starts_out.set_len(total);
        ends_out.set_len(total);
    }
    (idxs_out, starts_out, ends_out)
}

/// Split one or more whole chromosomes of sorted events into atoms, passing
/// each to `emit` as (index, start, end).
///
/// With `ended`, a covered atom's index names the innermost interval still open
/// over it. `ended` is shared by every task, one bit per row; a task only ever
/// reads the bits of its own chromosomes' rows, which it wrote itself. Without
/// it the sweep only finds the atoms, which is all the counting pass needs.
fn sweep<G: GroupType, T: PositionType>(
    events: &[Event<G, T>],
    between: bool,
    ended: Option<&[AtomicU64]>,
    mut emit: impl FnMut(u32, T, T),
) {
    // Edge case: no intervals
    if events.is_empty() {
        return;
    }

    // The intervals open at the sweep line, innermost last. An interval that
    // has ended stays on the stack until it reaches the top, and is marked in
    // `ended` so it can be recognised then; each is pushed and popped once.
    let mut open: Vec<u32> = Vec::new();
    let is_ended = |ended: &[AtomicU64], idx: u32| {
        ended[idx as usize / 64].load(Ordering::Relaxed) >> (idx % 64) & 1 == 1
    };

    // State for the sweep line
    let mut current_chr = events[0].chr;
    // We initialize coverage to 0, then we will “process” each event,
    // but we need a “last_pos” to track from where we last emitted.
    let mut active_count: u32 = 0;
    let mut last_pos = events[0].pos;
    let mut last_idx = events[0].idx;

    // Decide whether coverage is “on” at the very first position:
    // If the first event is a start, coverage goes from 0 → 1 at that point.
    if events[0].is_start {
        active_count = 1;
        if ended.is_some() {
            open.push(events[0].idx);
        }
    }

    // We iterate from the *second* event onward.
    // At each new event, we emit from last_pos → e.pos if either coverage was > 0 or `between = true`.
    for e_i in 1..events.len() {
        let e = &events[e_i];

        // If chromosome changes, we “jump” to a new chromosome
        // and do *not* produce an interval bridging old->new.
        if e.chr != current_chr {
            // reset
            current_chr = e.chr;
            active_count = if e.is_start { 1 } else { 0 };
            last_pos = e.pos;
            last_idx = e.idx;
            open.clear();
            if e.is_start && ended.is_some() {
                open.push(e.idx);
            }
            continue;
        }

        // same chromosome => we may emit from last_pos..e.pos if it's > 0 length
        // and either coverage>0 or we want the gap (between = true).
        if e.pos > last_pos {
            // If we were in coverage or want gaps, emit the subinterval.
            if active_count > 0 || between {
                // A covered atom belongs to the innermost interval still open;
                // a gap belongs to none, and keeps the index of the event before it.
                let attribution = match ended {
                    Some(ended) if active_count > 0 => {
                        while is_ended(ended, open[open.len() - 1]) {
                            open.pop();
                        }
                        open[open.len() - 1]
                    }
                    _ => last_idx,
                };
                emit(attribution, last_pos, e.pos);
            }
            last_pos = e.pos;
            last_idx = e.idx;
        }

        // Now handle the event itself (this flips coverage up or down).
        if e.is_start {
            active_count += 1;
            if ended.is_some() {
                open.push(e.idx);
            }
        } else {
            // is an end
            if active_count > 0 {
                active_count -= 1;
            }
            if let Some(ended) = ended {
                ended[e.idx as usize / 64].fetch_or(1 << (e.idx % 64), Ordering::Relaxed);
            }
        }
    }
}

/// The sorted events of every row that covers at least one position.
///
/// A row with `end + slack <= start` covers nothing, as in subtract, so it must
/// neither open nor close coverage. Left in, its end event sorts before its own
/// start and is absorbed while nothing else is open, leaving the sweep one
/// interval deep until the end of the chromosome.
fn build_events_of_nonempty_rows<G: GroupType + Send + Sync, T: PositionType + Send + Sync>(
    chrs: &[G],
    starts: &[T],
    ends: &[T],
    slack: T,
) -> Vec<Event<G, T>> {
    let covers_nothing = |row: usize| ends[row] + slack <= starts[row];
    if !(0..starts.len()).into_par_iter().any(covers_nothing) {
        return sorts::build_sorted_events_single_collection_par(chrs, starts, ends, slack);
    }
    let kept: Vec<u32> = (0..starts.len())
        .filter(|&row| !covers_nothing(row))
        .map(|row| row as u32)
        .collect();
    let gather = |values: &[_]| {
        kept.iter()
            .map(|&row| values[row as usize])
            .collect::<Vec<_>>()
    };
    let mut events = sorts::build_sorted_events_single_collection_par(
        &kept
            .iter()
            .map(|&row| chrs[row as usize])
            .collect::<Vec<_>>(),
        &gather(starts),
        &gather(ends),
        slack,
    );
    // Back to the caller's row numbers.
    for event in &mut events {
        event.idx = kept[event.idx as usize];
    }
    events
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Each atom's covering rows, by brute force. A row spans [start, end + slack).
    fn covering(
        starts: &[i64],
        ends: &[i64],
        chrs: &[u32],
        slack: i64,
        (chr, s, e): (u32, i64, i64),
    ) -> Vec<u32> {
        (0..starts.len())
            .filter(|&i| chrs[i] == chr && starts[i] <= s && e <= ends[i] + slack)
            .map(|i| i as u32)
            .collect()
    }

    #[test]
    fn every_atom_is_attributed_to_an_interval_that_covers_it() {
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for &n in &[1_usize, 2, 7, 50, 400] {
            for _ in 0..200 {
                let chrs: Vec<u32> = (0..n).map(|_| (next() % 3) as u32).collect();
                let starts: Vec<i64> = (0..n).map(|_| (next() % 60) as i64).collect();
                // Some zero-length rows, some long ones that nest the rest.
                let ends: Vec<i64> = starts
                    .iter()
                    .map(|&s| s + [0, 1, 3, 8, 40][(next() % 5) as usize])
                    .collect();
                for slack in [0, 1, 3] {
                    let (idxs, s_out, e_out) =
                        sweep_line_split(&chrs, &starts, &ends, slack, false);
                    for k in 0..idxs.len() {
                        let row = idxs[k] as usize;
                        let chr = chrs[row];
                        assert!(
                            covering(&starts, &ends, &chrs, slack, (chr, s_out[k], e_out[k]))
                                .contains(&idxs[k]),
                            "slack {slack}: atom [{}, {}) attributed to row {row} = [{}, {})",
                            s_out[k],
                            e_out[k],
                            starts[row],
                            ends[row]
                        );
                    }
                    // The atoms are exactly the covered positions, once each.
                    let mut want: Vec<(u32, i64)> = (0..n)
                        .flat_map(|i| (starts[i]..ends[i] + slack).map(move |p| (i, p)))
                        .map(|(i, p)| (chrs[i], p))
                        .collect();
                    want.sort();
                    want.dedup();
                    let mut got: Vec<(u32, i64)> = (0..idxs.len())
                        .flat_map(|k| (s_out[k]..e_out[k]).map(move |p| (k, p)))
                        .map(|(k, p)| (chrs[idxs[k] as usize], p))
                        .collect();
                    got.sort();
                    assert_eq!(got, want, "n={n} slack={slack}");
                }
            }
        }
    }

    #[test]
    fn chromosomes_swept_in_parallel_match_one_serial_sweep() {
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        // Well over BATCH_ROWS events, so several tasks share the bitset.
        let n: usize = 300_000;
        let chrs: Vec<u32> = (0..n).map(|_| (next() % 24) as u32).collect();
        let starts: Vec<i64> = (0..n).map(|_| (next() % 1_000_000) as i64).collect();
        let ends: Vec<i64> = starts.iter().map(|&s| s + (next() % 400) as i64).collect();
        for between in [false, true] {
            let events = build_events_of_nonempty_rows(&chrs, &starts, &ends, 0);
            let ended: Vec<AtomicU64> = (0..n.div_ceil(64)).map(|_| AtomicU64::new(0)).collect();
            let mut serial = (Vec::new(), Vec::new(), Vec::new());
            sweep(&events, between, Some(&ended), |idx, start, end| {
                serial.0.push(idx);
                serial.1.push(start);
                serial.2.push(end);
            });
            assert_eq!(
                sweep_line_split(&chrs, &starts, &ends, 0, between),
                serial,
                "between={between}"
            );
        }
    }

    #[test]
    fn zero_length_rows_cover_nothing() {
        let chrs = [0_u32; 4];
        let starts = [10_i64, 20, 40, 70];
        let ends = [10_i64, 30, 50, 80];
        let (_, s, e) = sweep_line_split(&chrs, &starts, &ends, 0, false);
        assert_eq!(
            s.iter().zip(&e).map(|(&s, &e)| (s, e)).collect::<Vec<_>>(),
            [(20, 30), (40, 50), (70, 80)]
        );
        let (_, s, e) = sweep_line_split(&chrs, &starts, &ends, 0, true);
        assert_eq!(
            s.iter().zip(&e).map(|(&s, &e)| (s, e)).collect::<Vec<_>>(),
            [(20, 30), (30, 40), (40, 50), (50, 70), (70, 80)]
        );
    }

    #[test]
    fn nested_intervals_take_the_innermost() {
        let (idxs, s, e) =
            sweep_line_split(&[0_u32; 3], &[0_i64, 10, 12], &[100_i64, 20, 14], 0, false);
        assert_eq!(s, [0, 10, 12, 14, 20]);
        assert_eq!(e, [10, 12, 14, 20, 100]);
        assert_eq!(idxs, [0, 1, 2, 1, 0]);
    }
}
