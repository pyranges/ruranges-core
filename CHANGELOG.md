# Changelog

## 0.1.14

The sweep kernels run on every core instead of one. Nothing in this crate used
`rayon` outside `nearest` and `ranks`; `overlap`, `overlap_indices`, `subtract`
and `merge` were serial end to end, which is why twelve threads bought them
almost nothing (hg38 10^8, one thread → twelve: `overlap(multiple=true)` 37.3 s →
34.8 s, `subtract` 19.8 s → 14.5 s, and on 634k-group proteome input `subtract`
77.2 s → 74.1 s — no scaling at all).

### Kernels

- **A parallel stable counting partition (`partition`, internal).** Group ids
  arrive dense from a frame-level factorisation, so a row's group can index a
  histogram: per-block counts, an exclusive prefix sum in group-major order, then
  a parallel scatter. It replaces the serial radix pass that made each group
  contiguous, and it is *stable*, so the result is element-for-element what the
  radix produced. It declines, and the caller keeps the radix path, when the
  group ids are too spread out to index or the histogram would be too large.
- **Sorts and sweeps run per group.** Once the groups are contiguous, each is
  independent, so both the coordinate sort and the sweep that follows it run on
  their own task — with consecutive groups batched, so that 24 chromosomes and
  634k proteins both get a sensible number of tasks. `build_sorted_intervals_par`
  and `build_sorted_events_single_collection_par` are the new builders; `merge`,
  `max_disjoint`, `cluster`, `split`, `boundary` and `complement` use them.
  Output is unchanged in every case, and the parity is under test.
- **`sweep_line_cluster` sweeps intervals, not events.** Clustering asks the same
  question merging does — an interval joins the run in progress exactly when
  `start < max_end + slack` — so it now sorts one record per interval instead of
  one event per endpoint, and scans each group on its own thread. The event
  phrasing sorted twice the elements and then walked them in a single serial pass
  over the whole genome, which is why only its sort saw more than one core. Ids
  are numbered exactly as before, gaps at group boundaries included; a serial
  prefix sum over the per-group cluster counts reproduces them. At hg38 10^8 the
  sweep goes 4.38 s → 1.35 s on one thread and 1.61 s → 0.32 s on twelve.
- **Coordinates that fit `i32` are swept as `i32`.** An `Interval<u32, i32>` is
  16 bytes against 24, so the sort moves a third less data, and the narrowing is
  free: the counting partition already copies every coordinate once and converts
  as it writes. Genomic readers hand out `Int64` columns for coordinates that
  never leave `i32`, so this is the usual case. `sweep_line_merge_narrowed` and
  `sweep_line_cluster_narrowed` take `i64` slices and return `None` — leaving the
  caller on its existing path — when a coordinate or the slack is outside
  ±2^30, the bound that leaves the sweeps room to add a slack to an end without
  overflowing.
- **The interval builder hands back the group runs.** Recovering them afterwards
  is a serial pass over the sorted array, which is 2.4 GB at 10^8 rows and runs
  after all the parallel work is done. The counting partition already knows where
  each group starts, so `build_sorted_intervals_par_with_runs` returns them;
  `build_sorted_intervals_par` is unchanged and never pays to compute them.
- **`sweep_line_subtract` no longer builds an event list or a hash map.**
  Subtracting a set is the same as subtracting its union, so it now sorts each
  side once, walks each group's right intervals into a disjoint cover, and emits
  the gaps each left interval has in that cover. The old version built
  `2 * (n1 + n2)` events and radix-sorted them three times — at hg38 10^8 that is
  ~3.4 GB moved three times, against `n1 + n2` records sorted once each now — and
  tracked live left intervals in an `FxHashMap` that it scanned *in full* on every
  0↔1 transition of the right-hand depth, so its cost grew as (left pileup depth ×
  right transitions): worst exactly where the data is densest.
- **`overlap_indices` stops paying for a tie-break it discards.** It returns only
  `query.idx` and never emits a target, so the `(group, start, end, idx)` sort key
  that `First`/`Last` need in `overlaps` — to decide *which* target comes back,
  load-bearing for pyranges/ruranges_py#23 — cannot change its output. Both sides
  now take the narrow `(group, start)` key. `overlaps` is unchanged.
- **A closed form for dense groups.** Where a group's pair density makes the
  sweep's per-query candidate scan the dominant cost, `overlap_indices` counts
  with Binary Interval Search (Layer et al. 2013) instead: two `partition_point`
  searches per query, O(log M) regardless of how many targets it hits. Guarded —
  the identity is exact only for `slack >= 0` with positive-length targets, and a
  256-query sample decides whether the extra end sort is repaid, because on sparse
  input it is not. `contained_intervals_only` is a different predicate and keeps
  the sweep.
- **`count_overlaps` was the last serial kernel, and the slowest thing here.** It
  built its records with the serial radix and then walked them in one pass whose
  per-query step scans the live candidate window — so its cost grew with the
  number of candidate *pairs*, not with the input. On 634k-group proteome input at
  10^8 rows that is 229 s against `merge`'s 4.3 s. It now takes the same two paths
  `overlap_indices` does: groups made contiguous by the counting partition, each
  counted on its own task, and the BITS closed form where a group's density
  repays its end sort. Same counts, in the same input row order — every count is
  written by the one task that owns its row's group. At twelve threads the
  operation goes 229 s → 9.1 s on proteome 10^8 and 25.2 s → 6.9 s on hg38 10^8.
- **The binary kernels sweep narrow coordinates as `i32`**, as `merge` and
  `cluster` already did. `sorted_records_narrowed_par` converts while the
  partition copies, so a 16-byte `IntervalRecord<u32, i32>` is sorted where a
  24-byte one used to be, and the entry points below return `None` — leaving the
  caller on its existing path — when a coordinate or the slack is out of range.
  `subtract` widens its surviving pieces back to `i64`; `overlaps`,
  `overlap_indices` and `count_overlaps` return row numbers and counts, so
  nothing widens.
- **`nearest` builds its six sorted views in parallel.** It was still calling the
  serial record builder and sorting its event arrays with two whole-array radix
  passes; both now take the counting partition, and the narrowed entry point
  converts the coordinates once up front rather than per view. At hg38 10^8 on
  twelve threads the build stage goes 2.26 s → 0.47 s, and the operation 6.3 s →
  4.6 s.
- **The counting partition declines when its bookkeeping would outweigh the
  data.** The histogram is one counter per (block, group) and the prefix sum over
  it is serial, so a frame with a few rows per group — a proteome frame at 10^6
  rows has 340k groups, three rows each — spent more time on the histogram than
  on the rows. It now asks for no more blocks than the pool can run several of at
  once, and hands back `None` when `blocks * groups` would exceed the row count;
  the callers all keep a radix path. The narrowed builder takes that path too, so
  the narrowing survives a frame the partition declines.
- **Per-group sorts use insertion sort below 32 elements.** The radix sort's cost
  per *call* dominates on short runs, and short runs are the normal case on a
  high-cardinality frame: building 10^6 records over 340k groups spent 0.127 s in
  the per-group sorts against 0.016 s now. Both are stable, so the output is
  unchanged.

### API

- `overlap_indices_chunked` returns the per-group index batches without
  concatenating them. A caller that is going to gather rows can use them as
  gather chunks; the flat vector is several GiB at 10^8 rows and exists only to be
  split up again. `overlap_indices` is unchanged and now wraps it, and
  `sort_indices` exposes the sort it applies.
- The kernels listed above gained `Send + Sync` bounds on their group and
  position type parameters. Every primitive that satisfies `GroupType` or
  `PositionType` already satisfies them, so no existing call site needs an edit.
- `sweep_line_merge_narrowed` and `sweep_line_cluster_narrowed` are the `i64`
  entry points to the narrowed sweeps, and `sorts::build_sorted_narrowed_intervals_par`
  the builder behind them. `sorts::NARROW_COORDINATE_LIMIT` and
  `sorts::slack_fits_narrow` name the bound they apply.
- `sweep_line_subtract_narrowed`, `overlaps_narrowed`,
  `overlap_indices_chunked_narrowed`, `count_overlaps_narrowed` and
  `nearest_with_ties_narrowed` are the `i64` entry points to the narrowed binary
  kernels. `nearest` takes `sorts::NARROW_DISTANCE_LIMIT`, half the usual bound:
  its answer is a *difference* of two coordinates, so the whole span a group can
  hold has to fit `i32`, not just each endpoint. It therefore declines inputs the
  other four accept.
- `count_overlaps` and `nearest_with_ties` gained `Send + Sync` bounds on their
  type parameters, as the kernels above already had. Every primitive that
  satisfies `GroupType` or `PositionType` already satisfies them.
- `sorts::build_sorted_intervals_par_with_runs` returns each group's run
  alongside the sorted intervals.

### Behaviour

- **Intervals of non-positive length in `subtract`.** They are now ignored on the
  right (they cover nothing) and produce no output on the left (nothing survives
  subtraction from an empty span). Previously such an interval's end event sorted
  *before* its own start event, so a zero-length left interval was never removed
  from the live set and could emit spurious output spanning far beyond itself,
  and an inverted right interval blocked output over the region between its two
  coordinates. `polaranges` rejects these before the kernel, so only direct
  callers were exposed.
- **Unsorted output order.** With `sort_output = false` the row *order* of
  `overlap_indices` and `sweep_line_subtract` may differ from 0.1.13; the rows
  themselves do not. Subtract's old order came from `FxHashMap` iteration and was
  arbitrary; it is now deterministic. `sort_output = true` is unaffected.
- **Intervals of non-positive length in `cluster`.** An interval whose end does
  not exceed its start now opens a cluster like any other, and closes it at its
  start. Under the event sweep its end event sorted before its own start event,
  which drove the live count negative and could leave the following intervals
  labelled with a cluster that had already been closed. As with `subtract`,
  `polaranges` rejects these before the kernel, so only direct callers were
  exposed.

## 0.1.13

- `nearest_with_ties` takes a `Ties` argument: `Ties::All` reports every
  neighbour at each reported distance, as `nearest` always has, and
  `Ties::First` reports one of them. `nearest` itself is unchanged — same
  arity, same results — and now calls the new function with `Ties::All`.
- The saving is real, not cosmetic. With `include_overlaps`, every interval a
  query overlaps is a neighbour at distance 0, so a query inside a dense region
  ties with everything covering it: 100 million hg38-like queries against
  themselves produce 1.1 billion rows under `All` and 7.7 million under
  `First`, 144x fewer. `First` never builds the rows it does not return — the
  overlap sweep stops at the first hit per query (`OverlapType::First`), the
  two directional sweeps emit one row per distinct position, and the merge
  keeps one row per distance bucket. Filtering afterwards would have cost the
  same memory as `All` and more time.
- Which of the tied neighbours comes back is deliberately unspecified, matching
  `bedtools closest -t first`, GenomicRanges `select="arbitrary"` and BEDOPS
  `--closest`. It is reproducible for a given input; nothing is sorted to make
  it so.
- `k` still counts distinct distances. Under `First` that means at most one row
  per distance and so at most `k` rows per query, rather than the combination
  being rejected.
- The three lower-level helpers `nearest_intervals_to_the_left`,
  `nearest_intervals_to_the_right` and `merge_three_way_by_index_distance` each
  gained the same `Ties` argument. They are the internals of `nearest` and have
  no known caller outside this crate, but they are `pub`, so a direct user of
  them needs one edit per call site.

## 0.1.12

- New `ranks` module: `natural_rank`, `lexical_rank` and `natural_cmp` order a
  key column's *distinct* string values, and `fold_ranks` collapses several
  coded keys into the single `group` argument `sorts::sort_order_idx` takes.
  Together with `sorts` these are the whole of a multi-key interval sort that a
  dataframe library cannot do faster itself.
- `natural_cmp` is the definition of natural ordering (`chr2 < chr10 < chrM`)
  for every library built on ruranges. It lives here so that `pyranges1` and
  `polaranges` cannot order the same chromosome names differently -- previously
  each had its own, and they disagreed on keys with non-digit characters next to
  digits (`a-1` versus `a1`). The unit tests check it against expectations
  produced by Python's `natsort` 8.4.0 and document the two places the two
  deliberately differ: no Unicode normalisation, and ASCII digits only.
- Additive: no existing item changed, so 0.1.11 callers need no edits.

## 0.1.11

- `map_to_global_with_exons` reports the annotation row that produced each
  projected fragment as a fifth output. Callers that carry annotation metadata
  onto the result no longer have to rediscover the emitting interval, which
  otherwise costs a scan of every annotation row per query.
- `map_to_global` is unchanged: same arity, same four outputs. The public
  `StrandInterval` is unchanged too, so this release is additive.

## 0.1.9

- Fix contained overlaps with contracted negative slack. `overlaps` and the
  join kernels previously leaked sweep state when `contained` was combined with
  a negative slack, so a query could report a hit that the contracted interval
  no longer supports.

  This changes results for the `contained` + negative-slack combination. Any
  build that picks up 0.1.9 or later will see the corrected behaviour.
