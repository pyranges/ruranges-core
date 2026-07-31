# Changelog

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
