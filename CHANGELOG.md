# Changelog

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
