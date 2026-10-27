# Changelog

## Unreleased

- Complement omits empty and reversed terminal gaps. Nonnegative slack filters
  internal and terminal gaps without shifting the coordinates of surviving gaps
  or overflowing when interval ends are near the coordinate type's maximum.
  Optional leading gaps and negative-slack end contraction retain their behavior.

## 0.1.14

- Ordinary overlap reporting removes expired targets throughout the candidate
  list and bounds each scan by the current query end. Long targets and nested
  queries no longer cause repeated scans of irrelevant candidates. First/last
  matches avoid enumerating every overlapping pair.
- `count_overlaps` uses endpoint ranks and sorted scans or binary searches,
  with direct checks for tiny inputs. An exact Fenwick-tree fallback preserves
  zero-length/reversed interval, negative-slack, and saturation behavior.
  Some sparse fallback-heavy workloads are slower.
- For `k=1` with `Ties::All`, nearest retains endpoint-bucket references and
  expands only the winning ties. It assembles results directly in query order,
  avoiding discarded candidate vectors and intermediate result sorts.
  `Ties::First`, the general-k API, and established output ordering are retained.
- Subtraction merges right-side coverage and uses a greedy scan instead of
  revisiting every active left interval at each right-side coverage boundary.
  Zero-length and reversed intervals produce no left fragments and remove no
  right-side coverage; this fixes bogus fragments from empty intervals.
- Preserve all fixes and public APIs from 0.1.13, including contracted
  negative-slack containment and nearest tie policies.

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
