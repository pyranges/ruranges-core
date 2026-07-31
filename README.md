# ruranges-core

Pure Rust core interval algorithms for the `ruranges` Python extension.

The crate is deliberately dataframe-agnostic: everything here takes and returns
plain slices, so both `pyranges1` (pandas) and `polaranges` (polars) can build on
it without either dataframe library appearing in the dependency tree.

## Modules

Most modules are interval algorithms — overlaps, nearest, merge, cluster,
complement, subtract, tile, split. Two are the sort pipeline, and are worth
calling out because they are shared *specification*, not just shared code:

* [`sorts`](src/sorts.rs) — `sort_order_idx` computes the row order for a sort by
  `(group, start, end)` with optional per-row direction reversal. Three stable
  LSD radix passes.
* [`ranks`](src/ranks.rs) — `natural_rank` / `lexical_rank` turn a key column's
  *distinct string values* into ascending integer codes, and `fold_ranks`
  collapses several coded keys into the single `group` argument `sort_order_idx`
  takes.

`ranks::natural_cmp` is the definition of natural ordering (`chr2 < chr10 <
chrM`) for every library built on ruranges. It lives here rather than in each of
them so the same data cannot sort two ways; the unit tests check it against
expectations produced by Python's `natsort` 8.4.0, and document the two places
the two deliberately differ (Unicode normalisation, non-ASCII digits).

## Development

```bash
cargo test
```
