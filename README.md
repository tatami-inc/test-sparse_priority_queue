# Performance testing for sparse priority queues

## Motivation

See [tatamic-inc/tatami#49](https://github.com/tatami-inc/tatami/issues/49) for the initial suggestion.

The idea here is, when iterating over the secondary dimension of a sparse matrix,
we use a priority queue to keep track of the current secondary index for each primary dimension element.
Then, when we request the next secondary index, we only have to inspect the primary elements with secondary indices that are not greater than the requested one.
This should reduce the search time when compared to a linear scan over the primary dimension elements... in theory.

For comparison, we implement a "simple" linear search that caches the lower bounds for the last requested index in each primary dimension element.
This is used by **tatami** to avoid a costly binary search on each look-up.
We also implement a more complex linear search that short-circuits the iteration across the primary dimension when an all-zero row is encountered.
The short-circuit is achieved by holding the minimum of the lower bound across the primary dimension elements;
if this lower bound is greater than the requested index, we know that a row is all-zero and can skip it.

## Results

In practice, the performance of the queue is disappointing:

```console
$ ./build/testing
Testing a 10000 x 50000 matrix with a density of 0.1
Expecting a sum of 49995342

|               ns/op |                op/s |    err% |     total | benchmark
|--------------------:|--------------------:|--------:|----------:|:----------
|    1,620,330,293.00 |                0.62 |    0.3% |     18.08 | `linear simple`
|    1,595,880,178.00 |                0.63 |    3.6% |     17.64 | `linear shortcircuit`
|    5,359,585,621.00 |                0.19 |    1.2% |     59.00 | `queue`
```

With hindsight, this is perhaps unsurprising.
The queue takes `O(Nz * log(N))` time to depopulate and fill, given `N` columns and `Nz` non-zero elements.
The linear search just takes `O(N)` time, which becomes faster than the queue when `Nz` is a decent fraction of `N`.
Even at 5% density, we still observe a performance degradation:

```console
$ ./build/testing -d 0.05
Testing a 10000 x 50000 matrix with a density of 0.05
Expecting a sum of 24993148

|               ns/op |                op/s |    err% |     total | benchmark
|--------------------:|--------------------:|--------:|----------:|:----------
|    1,074,593,257.00 |                0.93 |    0.2% |     11.87 | `linear simple`
|      942,579,753.00 |                1.06 |    0.4% |     10.58 | `linear shortcircuit`
|    2,829,276,745.00 |                0.35 |    1.9% |     30.94 | `queue`
```

Eventually, at low densities, we see the desired improvement over the linear method.

```console
$ ./build/testing -d 0.01
Testing a 10000 x 50000 matrix with a density of 0.01
Expecting a sum of 4998047

|               ns/op |                op/s |    err% |     total | benchmark
|--------------------:|--------------------:|--------:|----------:|:----------
|      644,002,568.00 |                1.55 |    2.7% |      7.08 | `linear simple`
|      425,376,361.00 |                2.35 |    0.4% |      4.78 | `linear shortcircuit`
|      580,639,986.00 |                1.72 |    0.1% |      6.39 | `queue`

$ ./build/testing -d 0.001
Testing a 10000 x 50000 matrix with a density of 0.001
Expecting a sum of 499806

|               ns/op |                op/s |    err% |     total | benchmark
|--------------------:|--------------------:|--------:|----------:|:----------
|      483,900,455.00 |                2.07 |    0.7% |      5.34 | `linear simple`
|      253,605,505.00 |                3.94 |    0.7% |      2.81 | `linear shortcircuit`
|       55,406,053.00 |               18.05 |    0.2% |      0.62 | `queue`

$ ./build/testing -d 0.0001
Testing a 10000 x 50000 matrix with a density of 0.0001
Expecting a sum of 50016

|               ns/op |                op/s |    err% |     total | benchmark
|--------------------:|--------------------:|--------:|----------:|:----------
|      477,359,379.00 |                2.09 |    0.6% |      5.27 | `linear simple`
|      239,555,460.00 |                4.17 |    0.8% |      2.64 | `linear shortcircuit`
|        5,448,334.00 |              183.54 |    1.0% |      0.06 | `queue`
```

The short-circuit method also hits optimal performance for all-zero matrices, as intended:

```console
$ ./build/testing -d 0
Testing a 10000 x 50000 matrix with a density of 0
Expecting a sum of 0

|               ns/op |                op/s |    err% |     total | benchmark
|--------------------:|--------------------:|--------:|----------:|:----------
|      468,395,090.00 |                2.13 |    0.1% |      5.16 | `linear simple`
|           99,293.40 |           10,071.16 |    0.1% |      0.01 | `linear shortcircuit`
|          120,154.60 |            8,322.61 |    3.2% |      0.01 | `queue`
```

# Discussion

Turns out that linear iteration is pretty fast if there's no action involved other than to hit `continue`. 
Indeed, **tatami** stores its own cache of lower-bounding indices that can be iterated contiguously in memory for optimal access,
unlike the queue where information for adjacent columns are scattered around via push/pop cycles to the underlying heap.
I suspect that the linear search is also very friendly to the branch predictor for very sparse rows where most columns can be skipped.

Note that this demonstration is already a best-case example of using the queue.
If it's not a strict increment, we have to pop the queue until we get to the row we want, incurring wasteful `log(N)` pushes;
or we can scan the underlying heap directly, but then the indices of the non-zero elements won't come out in order and we'll need to resort them.
Moreover, if there is a large jump between the previous and current requested row, 
the information in the queue is useless and we collapse back to `O(N)` time anyway (with extra overhead to rebuild the queue).

I suppose we _could_ switch between the linear and queue methods depending on the sparsity of the matrix, but that seems pretty complicated.
Given the absence of an unqualified performance benefit, I prefer the linear search due to the simplicity of its implementation.
The suboptimality at low densities is tolerable as everything is already fast enough,
though I'll concede that the short-circuiting is useful to skip all-zero rows that are common in genomics data.

As an aside: I have no idea why the short-circuit method is faster than the simple method at moderately low densities (0.1-1%).
There are no all-zero rows in these situations, meaning that none of the short-circuits can actually be taken.
This means that the short-circuit method should be incurring some overhead for zero performance benefit.
Nonetheless, it is faster than the simple method due to some quirk of the compiled code.
¯\\\_(ツ)\_/¯

# Build instructions

Just use the usual CMake process:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```
