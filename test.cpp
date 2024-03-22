#include "CLI/App.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"
#include "nanobench.h"

#include <chrono>
#include <vector>
#include <queue>
#include <random>

struct Linear {
private:
    const std::vector<std::vector<int> >& indices;
    int max_index;

private:
    // The cached position of the pointer at each primary element.
    // Specifically, 'indices[i][cached_indptrs[i]]' is the lower bound for 'last_request' in the primary element 'i'.
    std::vector<size_t> cached_indptrs; 

    // This vector contains the cached index being pointed to by 'cached_indptrs'.
    // We store this here as it is more cache-friendly than doing a look-up to 'indices' every time.
    std::vector<int> cached_indices;

    // Closest value in 'cached_indices' to the 'last_request', to see whether we can short-circuit the iteration.
    // This is the minimum of values in 'cached_indices'.
    int closest_cached_index = 0;

    // This is the first position of 'cached_indices' equal to 'closest_cached_index',
    // i.e., 'cached_indices[closest_cached_index_position] == closest_cached_index'.
    size_t closest_cached_index_position = 0;

    // What was the last requested index on the secondary dimension?
    int last_request = 0;

public:
    Linear(const std::vector<std::vector<int> >& idx, int mi) : 
        indices(idx), 
        max_index(mi), 
        cached_indptrs(indices.size()), 
        cached_indices(indices.size())
    {
        for (size_t p = 0, pend = indices.size(); p < pend; ++p) {
            const auto& curi = indices[p];
            cached_indices[p] = (curi.empty() ? max_index : curi.front());
        }
        if (!cached_indices.empty()) {
            auto closestIt = std::min_element(cached_indices.begin(), cached_indices.end());
            closest_cached_index_position = closestIt - cached_indices.begin();
            closest_cached_index = *closestIt;
        }
    }

private:
    template<class Store_>
    void search_above(int secondary, int primary, Store_ store) {
        // Skipping if the curdex (corresponding to curptr) is already higher
        // than secondary. So, we only need to do more work if the request is
        // greater than the stored index. This also catches cases where we're
        // at the end of the dimension, as curdex is set to max_index.
        auto& curdex = cached_indices[primary];
        if (curdex > secondary) {
            return;
        }

        auto& curptr = cached_indptrs[primary];
        if (curdex == secondary) {
            store(primary, cached_indptrs[primary]);
            return;
        }

        // Having a peek at the index of the next non-zero element; the
        // requested index should be equal to or below this, as would be the
        // case for consecutive accesses. An actual implementation would have
        // to account for non-consecutive jumps but we'll keep things simple
        // here for a comparison to an equally simple queue implementation.
        ++curptr;
        const auto& curi = indices[primary];
        auto endptr = curi.size();
        if (curptr == endptr) {
            curdex = max_index;
            return;
        }

        auto iraw = curi.begin();
        auto inext = iraw + curptr;
        curdex = *inext;
        if (curdex == secondary) {
            store(primary, curptr);
            return;
        }
    }

public:
    template<class Store_>
    void search_simple(int secondary, Store_ store) {
        for (size_t p = 0, pend = indices.size(); p < pend; ++p) {
            search_above(secondary, p, store);
        }
        last_request = secondary;
    }

    template<class Store_>
    void search_shortcircuit(int secondary, Store_ store) {
        if (secondary < closest_cached_index) {
            last_request = secondary;
            return; 
        }

        bool found = false;
        for (size_t p = 0, pend = indices.size(); p < pend; ++p) {
            search_above(secondary, p, [&](int i, size_t s) {
                store(i, s);
                found = true;
            });
        }

        if (found) {
            closest_cached_index = secondary;
        } else {
            closest_cached_index = *(std::min_element(cached_indices.begin(), cached_indices.end()));
        }

        last_request = secondary; 
    }
};

struct Pqueue {
private:
    const std::vector<std::vector<int> >& indices;

    typedef std::pair<int, int> HeapElement;
    std::priority_queue<HeapElement, std::vector<HeapElement>, std::greater<HeapElement> > next_heap;

    std::vector<int> hits, tmp_hits;
    std::vector<size_t> state;

public:
    Pqueue(const std::vector<std::vector<int> >& idx) : indices(idx), hits(indices.size()), state(indices.size()) {
        for (size_t c = 0, cend = indices.size(); c < cend; ++c) { // force everything to be re-searched on initialization.
            hits[c] = c;
            --state[c];
        }
    }

public:
    template<class Store_>
    void search(int secondary, Store_ store) {
        tmp_hits.swap(hits);
        hits.clear();

        // Refilling the indices we popped out in the last round.  This gives
        // us an opportunity to check whether they're equal to the current
        // 'secondary' (and thus elide an insertion into the queue).
        for (auto x : tmp_hits) {
            ++state[x];
            const auto& curx = indices[x];
            if (state[x] < curx.size()) {
                int current = curx[state[x]];
                if (current == secondary) {
                    hits.push_back(x);
                } else {
                    next_heap.emplace(current, x);
                }
            }
        }

        // Finding all the queue elements equal to our current position.
        // No need to do anything fancy when we're just incrementing;
        // it's always going to be '>= secondary'.
        while (!next_heap.empty()) {
            auto current_secondary = next_heap.top().first;
            if (current_secondary > secondary) {
                break;
            }
            auto current_primary_index = next_heap.top().second;
            next_heap.pop();
            hits.push_back(current_primary_index);
        }

        /* 
         * We're going to paint the priority queue in the best possible light
         * here by skipping the sort step, which is technically necessary to
         * have 1:1 feature parity with the linear methods.
         */
        // std::sort(hits.begin(), hits.end()); 
        for (auto x : hits) {
            store(x, state[x]);
        }
    } 
};

int main(int argc, char* argv []) {
    CLI::App app{"Sparse priority queue testing"};
    double density;
    app.add_option("-d,--density", density, "Density of the sparse matrix")->default_val(0.1);
    int nr;
    app.add_option("-r,--nrow", nr, "Number of rows")->default_val(10000);
    int nc;
    app.add_option("-c,--ncol", nc, "Number of columns")->default_val(50000);
    CLI11_PARSE(app, argc, argv);

    std::cout << "Testing a " << nr << " x " << nc << " matrix with a density of " << density << std::endl;

    // Simulating a sparse matrix, albeit not very efficiently, but whatever.
    std::vector<std::vector<int> > i(nc);
    std::mt19937_64 generator(1234567);
    std::uniform_real_distribution<double> distu;

    for (size_t c = 0; c < nc; ++c) {
        auto& curi = i[c];
        for (size_t r = 0; r < nr; ++r) {
            if (distu(generator) <= density) {
                curi.push_back(r);
            }
        }
    }

    size_t expected = 0;
    {
        Linear linear(i, nr);
        for (size_t r = 0; r < nr; ++r) {
            linear.search_simple(r, [&](int, size_t) { ++expected; });
        }
        std::cout << "Expecting a sum of " << expected << std::endl;
    }

    // Doing the linear iteration with simple caching.
    ankerl::nanobench::Bench().run("linear simple", [&](){
        Linear linear(i, nr);
        size_t sum = 0;
        for (size_t r = 0; r < nr; ++r) {
            linear.search_simple(r, [&](int, size_t) { ++sum; });
        }
        if (sum != expected) {
            std::cerr << "WARNING: different result from linear access (" << sum << ")" << std::endl;
        }
    });

    // Doing the linear iteration with more caching.
    ankerl::nanobench::Bench().run("linear shortcircuit", [&](){
        Linear linear(i, nr);
        size_t sum = 0;
        for (size_t r = 0; r < nr; ++r) {
            linear.search_shortcircuit(r, [&](int, size_t) { ++sum; });
        }
        if (sum != expected) {
            std::cerr << "WARNING: different result from linear shortcircuit access (" << sum << ")" << std::endl;
        }
    });

    // Comparing to a priority queue.
    ankerl::nanobench::Bench().run("queue", [&](){
        Pqueue pqueue(i);
        size_t sum = 0;
        for (size_t r = 0; r < nr; ++r) {
            pqueue.search(r, [&](int, size_t) { ++sum; });
        }
        if (sum != expected) {
            std::cerr << "WARNING: different result from queue access (" << sum << ")" << std::endl;
        }
    });

    return 0;
}

