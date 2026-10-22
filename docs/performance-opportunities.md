# Further performance opportunities in ruranges-core

Initial assessment of the local source tree on 2026-10-22. This document records
optimization opportunities and their rationale. A first `k=1` nearest patch has
since been implemented and measured; see the [nearest report](../benchmarks/nearest.md).
The nearest baseline discussion below describes the implementation before that
patch. Other proposals remain future work with no claimed measured speedup.

Nearest is the next priority because it is an important operation and its common
`k=1` case currently computes results that are later discarded. After nearest,
the clearest remaining algorithmic opportunities are the separate overlap-index
kernel, non-overlap selection, and clustering. Merge already uses a greedy scan.

Keep this work in the Rust core. Python factorization, pandas grouping, copying,
and result construction belong in a separate patch and separate measurements.

## Current position

The local tree already contains improvements to ordinary overlap reporting and
overlap counting. Their measurements and tradeoffs are documented in
[overlap benchmarks](../benchmarks/overlaps.md) and
[counting benchmarks](../benchmarks/count_overlaps.md).

Nearest calls the shared overlap-reporting kernel, so its overlap phase already
inherits that improvement. The separate `overlap_indices` implementation still
uses the previous active-list traversal. Improving one entry point does not
automatically improve every operation with a similar name.

| Priority | Kernel | Main opportunity | Expected beneficiaries |
|---|---|---|---|
| 1 | Nearest | Avoid discarded searches and intermediate result sorting | `k=1`, queries with overlaps, many endpoint ties |
| 2 | Overlap indices | Remove repeated scans of expired/future candidates | Spanning targets and nonmonotonic query ends |
| 3 | Non-overlap selection | Answer existence without repeatedly marking active queries | Long queries against many target intervals |
| 4 | Cluster | Sort intervals instead of two events per interval | Large ordinary interval collections |
| 5 | Merge and shared sorting | Recognize sorted inputs and reduce setup costs | Sorted inputs and small collections |

These priorities reflect source inspection and likely usefulness, not a measured
ranking of time spent in a representative application.

## 1. Nearest

### Baseline algorithm

[nearest.rs](../src/nearest.rs) builds the sorted views needed by three producers:

| Producer | Query order | Target order |
|---|---|---|
| Overlaps | Group, start | Group, start |
| Left neighbors | Group, start | Group, end |
| Right neighbors | Group, end | Group, start |

Shared start-sorted views are built once. Sorting and producer execution already
use Rayon. Each directional sweep advances a shared cursor monotonically and
then visits the requested endpoint positions and their ties. The cursor traversal
is linear after sorting; candidate emission still costs time proportional to the
number of candidates produced.

The producers currently materialize separate result vectors, sort them by query
index, and merge them. With sorted output, the retained results are then sorted
by `(query index, distance, target index)` and copied into output arrays.

This is already a reasonable sweep algorithm. The first question is how much
searching, materialization, and sorting can be avoided.

### A. Skip directional candidates when overlaps supply the answer

The Python `nearest_ranges` defaults are `k=1`, overlaps included, and all ties.
In the core, all overlaps for a query form the distance-zero bucket. When `k=1`
and that bucket exists, directional candidates cannot contribute to the result.

Currently both directional producers still run for those queries. The final
merge discovers the overlap bucket and discards their results. A candidate
implementation should identify queries answered by overlaps and skip their
directional candidate emission. The sweep cursor must still advance correctly
for later queries.

This can save much more than two candidates per query: many target intervals can
share a nearest endpoint, and the current producer emits all of those ties before
the merge discards them. For example, one overlapping target can make a large
bucket of tied non-overlapping neighbors irrelevant to a `k=1` query.

All overlapping targets remain required when the requested tie policy returns
all ties. Skipping directional candidates does not make a genuinely large final
overlap result cheap.

Computing overlap membership first introduces a dependency between phases that
currently execute concurrently. Benchmark the tradeoff, especially when overlaps
are rare. An unconditional extra preprocessing pass may lose on sparse input.

### B. Specialize the `k=1` result pipeline

For each query, retain references to the nearest left and right endpoint buckets
instead of immediately expanding both into `Nearest` records. Compare the two
distances, retain the winning bucket or both tied buckets, and expand only those
needed for output. An overlap bucket takes precedence when overlaps are included.

Per-query descriptors or counts/offsets could let the implementation assemble
results in original query order without sorting every intermediate candidate by
query index. Target-index ordering within equal-distance buckets must still be
preserved when sorted output requires it. Prefix offsets and an extra traversal
also have costs; compare them with the current vector-and-sort approach.

Keep the current sorted sweeps as the baseline for large collections. For a tiny
query collection against many targets, sorting only the required target endpoint
views and binary-searching each query is another candidate. This would avoid
sorting the query side; it is not a reason to replace large sorted sweeps with
random searches everywhere.

### C. Avoid producing `k` buckets independently on both sides

For larger `k`, each direction currently produces up to `k` distinct endpoint
positions, including all ties. The final merge retains up to `k` distinct
distances across overlaps and both directions combined.

A later design could locate each query's left/right boundaries, compare the next
available distance on each side, and expand only buckets that survive the global
limit. If overlaps are present, they already consume one distance bucket.

This may reduce discarded candidates and memory traffic. Retain a separate
general-`k` path until the simpler `k=1` work is measured and validated.

### D. Make setup proportional to the workload

Investigate these separately so each benefit can be measured:

- Detect already-sorted endpoint views before sorting them.
- Compare compact endpoint records or index arrays with full interval records
  where an overlap phase is unnecessary.
- Compare serial execution with Rayon below a measured size threshold.
- Measure comparison sorting versus radix sorting on each actual record shape.
  A threshold measured for overlap records is not automatically right for
  nearest endpoint records or output records.
- Check memory reserved for candidates later discarded, including the final
  merge's capacity based on the sum of all producer lengths.

Changing a hash container is not the main nearest opportunity: its merge already
tracks adjacent distance buckets without a per-query hash set.

### Compatibility requirements

Preserve the following across specialized and general paths:

- `k` counts distinct distances in the final result, rather than target rows.
- All required ties, including ties across left and right directions, survive.
- Duplicated intervals are distinct target rows; do not merge or deduplicate them.
- Included overlaps have distance zero; bookended intervals have distance one at
  zero slack in the existing tests.
- Forward, backward, and either-direction requests retain their current meaning.
- Queries without eligible targets, empty inputs, and missing groups behave as
  before.
- Sorted output keeps the established query/distance/target-index ordering.
  Validate unsorted-output behavior separately rather than assuming its order is
  irrelevant to callers.
- Exercise zero-length and reversed intervals, negative/positive slack, `k=0`, and
  integer limits before deciding which inputs can use a specialized path.

Nearest currently uses ordinary subtraction/addition in some endpoint and
distance calculations; the counting kernel uses saturating arithmetic. Do not
silently transfer counting's arithmetic semantics into nearest during a
performance patch. Characterize existing behavior and handle any intended
semantic change explicitly.

The core entry point inspected here has no `ties` argument, while the Python
call passes one to the NumPy binding. Inspect that binding before claiming an
optimization also avoids work for `ties="first"`.

## 2. Overlap indices

[overlap_indices.rs](../src/overlap_indices.rs) maintains a start-ordered active
vector and expires entries only from its head. A long target can block removal
of expired targets behind it. A long query can also admit targets beyond the end
of a later, shorter query. Those candidates are repeatedly checked.

For ordinary overlap reporting, reuse or extract the traversal techniques now
used by [overlaps.rs](../src/overlaps.rs): permanently unlink expired entries and
stop scanning at the current query's end bound. Prefer a shared traversal with
different output sinks over maintaining two divergent implementations. The index
sink should avoid allocating full pairs just to discard the right index.

The output contract matters: `All` repeats the left index once per matching
target. It is not an existence query. Output time is therefore at least
proportional to the number of matches. First/last selection does not require
enumerating all matches, although only the left index is returned.

Keep containment separate. Its predicate and candidate bounds differ from
ordinary overlap, so the ordinary-overlap stopping rule cannot be copied blindly.

## 3. Non-overlap selection

[complement.rs](../src/complement.rs), despite its module name, contains the
two-collection non-overlap selector. At each target start, it iterates over all
active queries and records that they overlap. A query that has already been
marked can be visited again for every subsequent target start while it remains
active. Replacing `FxHashSet` alone would retain that repeated work.

An existence query has a simpler solution than reporting or counting every
match. For the ordinary overlap predicate, group targets and sort by start. Store
the prefix maximum of their slack-adjusted ends. For query `[s,e)`, let
`limit = e.saturating_add(slack)` and find the prefix with target start `< limit`.
An overlap exists exactly when that prefix is nonempty and its maximum adjusted
end is `> s`.

This supports independent binary-search queries or a sweep over queries ordered
by `limit`. There is no need to enumerate overlapping pairs or maintain an active
set of query identities. A first-match overlap traversal is another option that
may allow more reuse of the existing optimized code.

The prefix-maximum statement is exact for the stated predicate, including
reversed targets. It does not establish equivalence to every edge behavior of
the current event implementation. Check empty inputs, event ties, slack,
degenerate intervals, and output order against a defined oracle before replacing
that implementation. Keep correctness fixes distinct from performance claims.

## 4. Cluster

[cluster.rs](../src/cluster.rs) currently creates two endpoint events per
interval, sorts them, and tracks the number of active intervals.

For ordinary nonempty intervals with zero or positive slack, a greedy scan can
instead sort once by `(group,start)`, track the maximum end of the current
cluster, and start a new cluster when the next interval cannot join it. Use the
maximum end seen in the cluster, not just the preceding interval's end; nested
intervals otherwise cause incorrect splits.

The existing event sort puts ends before starts at equal coordinates. At zero
slack, touching intervals consequently form separate clusters. Preserve that
strict boundary when deriving the greedy condition. Also preserve output order,
original indices, and numeric cluster IDs; the current event loop can leave gaps
in IDs at group transitions.

This removes event duplication and reduces sorting/memory traffic. It does not
establish a different asymptotic sorting complexity. Negative slack can move an
end event before its own start, so it needs an explicit compatibility assessment
or a retained fallback rather than assuming the ordinary greedy proof covers it.

## 5. Merge and shared sorting

[merge.rs](../src/merge.rs) already sorts intervals by group/start and performs a
greedy scan with the current maximum end. It does not have the unnecessary
two-event representation used by cluster.

Likely improvements are detecting sorted input, choosing an appropriate small
sort, and reducing allocation. Preserve the first contributing index, counts,
strict slack boundary, and tie ordering. Keep expectations below those for a
kernel that repeatedly scans irrelevant candidates.

The helpers in [sorts.rs](../src/sorts.rs) serve multiple operations. A shared
change can have a broad effect, including changing equal-key order. Benchmark
representative callers and retain stable ordering where it affects behavior.

## Benchmark and validation plan

Start with nearest, and measure each proposed change independently before
combining them. Use the existing isolated-crate benchmark approach: a fixed
baseline, current candidate, identical dependencies, release builds, and no
replacement of the installed Python extension.

### Nearest workload matrix

| Dimension | Cases |
|---|---|
| Size | Tiny, 100k, and 1M rows per side |
| Query/target ratio | Equal sizes, few queries, few targets |
| Row order | Shuffled; ordered within interleaved groups; fully grouped/sorted |
| Groups | One group, about 25 groups, many tiny groups, unmatched groups |
| Overlap frequency | None, sparse, mostly overlapping, all overlapping |
| Interval shape | Fixed lengths, nested/long intervals, highly variable ends |
| Endpoint ties | Few ties, many duplicate starts/ends, ties across directions |
| Parameters | `k=1`, modest/larger `k`, all directions, overlaps included/excluded |
| Slack | Zero, positive, negative; distinguish valid fast cases from fallbacks |
| Output | Sorted/unsorted; binding tie policies when testing integration |

Use deliberately large discarded tie buckets to expose wasted candidate
materialization. Include dense final-output cases to distinguish avoidable work
from the cost of returning results the caller requested.

For the other kernels, reuse spanning-target and long-query workloads from the
overlap benchmarks. Add nested intervals, touching boundaries, many tiny groups,
already-sorted input, and missing groups for cluster and merge.

### Measurement and acceptance

- Check complete outputs against the baseline before timing. Use an independent
  brute-force oracle on small randomized inputs; do not rely solely on two
  implementations agreeing.
- Preserve meaningful ordering and tie rules in oracle comparisons, including
  the multiplicity of repeated left indices.
- Include grouping, sorting, allocation, traversal, and result construction in
  the primary Rust timing. Exclude input generation and compilation.
- Rotate baseline/candidate timing order and report medians from repeated runs.
- Report Rayon thread counts and compare serial execution where relevant. The
  current nearest implementation is already parallel.
- Record final output size, intermediate candidate counts, and peak memory when
  practical. Phase timings can explain a result but do not replace total time.
- Publish ordinary cases and regressions alongside the favorable workloads.
  A faster pathology alone does not establish a better default for all users.
- Leave Python copying/grouping out of the core claim. Add an end-to-end result
  only when the binding and Python package are deliberately included and tested.

## Suggested patch sequence

1. Add nearest correctness oracles and a reproducible baseline benchmark.
2. Prototype skipping directional candidates for queries answered by overlaps
   at `k=1`; measure sparse and dense-overlap tradeoffs.
3. Prototype per-query bucket descriptors and direct result assembly for `k=1`.
4. Add only sorting/parallelism thresholds supported by those measurements.
5. Extend candidate selection to general `k` if discarded buckets remain a
   significant cost.
6. Address overlap indices, non-overlap selection, and clustering in separate
   patches with their own compatibility tests and benchmark reports.

The strongest initial nearest hypothesis is that the common `k=1` request can
avoid producing, sorting, and discarding directional candidates. The benchmark
should decide how much that matters relative to sorting and required output.
