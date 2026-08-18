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

### API

- `overlap_indices_chunked` returns the per-group index batches without
  concatenating them. A caller that is going to gather rows can use them as
  gather chunks; the flat vector is several GiB at 10^8 rows and exists only to be
  split up again. `overlap_indices` is unchanged and now wraps it, and
  `sort_indices` exposes the sort it applies.
- The kernels listed above gained `Send + Sync` bounds on their group and
  position type parameters. Every primitive that satisfies `GroupType` or
  `PositionType` already satisfies them, so no existing call site needs an edit.

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
