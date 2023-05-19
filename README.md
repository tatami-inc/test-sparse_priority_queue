# Performance testing for sparse priority queues

## Motivation

See [tatamic-inc/tatami#49](https://github.com/tatami-inc/tatami/issues/49) for the initial suggestion.

The idea here is, when iterating over the secondary dimension of a sparse matrix,
we use a priority queue to keep track of the current secondary index for each primary dimension element.
Then, when we request the next secondary index, we only have to inspect the primary elements with secondary indices that fall below the requested one.
This should reduce the search time when compared to a linear scan over the primary dimension elements... in theory.

## Results

In practice, the performance is disappointing:

```console
$ ./build/testing
Testing a 10000 x 100000 matrix with a density of 0.1
Linear access time: 5544 for 99997920 non-zero elements
Priority queue access time: 8713 for 99997920 non-zero elements
```

With hindsight, this is perhaps unsurprising.
The heap takes `O(Nz * log(N))` time to depopulate and fill, given `N` columns and `Nz` non-zero elements.
The linear search just takes `O(N)` time, which becomes faster than the heap when `Nz` is a decent fraction of `N`.
Even at 5% density, we still observe a mild performance degradation:

```console
$ ./build/testing -d 0.05
Testing a 10000 x 100000 matrix with a density of 0.05
Linear access time: 4098 for 49999099 non-zero elements
Priority queue access time: 4598 for 49999099 non-zero elements
```

Eventually, at very low densities, we see the desired improvement:

```console
$ ./build/testing -d 0.01
Testing a 10000 x 100000 matrix with a density of 0.01
Linear access time: 2245 for 9997846 non-zero elements
Priority queue access time: 930 for 9997846 non-zero elements
```

# Discussion

Turns out that linear iteration is pretty fast if there's no action involved other than to hit `continue`. 
Indeed, **tatami** stores the `current_indices` vector that can be iterated contiguously in memory for optimal access, unlike the heap where everything is scattered around via push/pop cycles.
I suspect that the linear search is also very friendly to the branch predictor for very sparse rows where most columns can be skipped.

Note that this demonstration is already a best-case example of using the heap.
If it's not a strict increment, we might need to sort the indices of the non-zero elements extracted from each row, because they won't come out of the heap in order.
Moreover, if there is a large jump between the previous and current requested row, the heap information is useless and we collapse back to `O(N)` time anyway (with extra overhead to rebuild the heap).

I suppose we _could_ switch between the linear and heap methods depending on the sparsity of the matrix, but that seems pretty complicated.
Given the absence of an unqualified performance benefit, I will prefer the simplicity of the linear search.

# Build instructions

Just use the usual CMake process:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```
