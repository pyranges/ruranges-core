# Changelog

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
