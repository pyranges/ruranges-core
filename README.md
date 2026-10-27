# ruranges-core

Pure Rust core interval algorithms for the `ruranges` Python extension.

Release **0.1.14** improves overlap reporting, overlap counting, subtraction,
and nearest with `k=1`. See the [changelog](CHANGELOG.md) for behavior details.

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

See [further performance opportunities](docs/performance-opportunities.md) for
the nearest optimization plan and other core kernels worth investigating.

## Nearest performance

For `k=1` with `Ties::All`, nearest records the closest endpoint bucket in each
direction and expands only the winning ties. Overlapping queries emit their distance-zero
matches directly. Results are assembled in input order without sorting expanded
directional candidates or a combined result vector. Required ties and the
existing sorted/unsorted output order are preserved; `Ties::First` and general
`k` keep their existing paths.

## Overlap performance

Ordinary overlap reporting unlinks expired targets while preserving start order
and stops at each query's end. This avoids rescanning expired targets behind a
spanning interval and future targets admitted by a previous long query. First
matches use a forward cursor; last matches use a predecessor structure to skip
expired targets without enumerating all matches. Input/output ordering and slack
semantics are preserved. Containment retains its existing traversal.

Small record collections use stable comparison sorting; larger collections use
radix sorting.

## Counting performance

`overlaps::count_overlaps` counts endpoint ranks instead of visiting every
matching pair. Larger query batches use scans in independently sorted endpoint
orders; smaller batches use binary searches. Tiny right-hand groups use direct
predicate checks. Compact group labels use a direct lookup table, while sparse
or signed labels use a hash map.

An exact Fenwick-tree sweep preserves the existing predicate for queries where
rank subtraction is invalid, including zero-length queries, negative slack,
reversed right intervals, and saturation at integer limits. Duplicates count
individually, results stay in input order, and final counts saturate at `u32::MAX`.

## Transcript projection

`map_to_global` and `map_to_global_with_exons` accept unsorted queries and exon
annotations. The core sorts numeric row indices only when necessary and clips
each query to the exon coverage. Exon local spans must not overlap within a
transcript. With `sort_output=true`, fragments follow original query and exon
row order; returned indices always refer to the original inputs. The exon-aware
form additionally returns the annotation row that produced each fragment.
