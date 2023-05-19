#include "CLI/App.hpp"
#include "CLI/Formatter.hpp"
#include "CLI/Config.hpp"

#include "tatami/tatami.hpp"

#include <chrono>
#include <vector>
#include <queue>
#include <random>

int main(int argc, char* argv []) {
    CLI::App app{"Sparse priority queue testing"};
    double density;
    app.add_option("-d,--density", density, "Density of the sparse matrix")->default_val(0.1);
    int nr;
    app.add_option("-r,--nrow", nr, "Number of rows")->default_val(10000);
    int nc;
    app.add_option("-c,--ncol", nc, "Number of columns")->default_val(100000);
    CLI11_PARSE(app, argc, argv);

    std::cout << "Testing a " << nr << " x " << nc << " matrix with a density of " << density << std::endl;

    // Simulating a sparse matrix, albeit not very efficiently, but whatever.
    std::vector<int> i, j;
    std::vector<double> x;

    std::mt19937_64 generator(1234567);
    std::uniform_real_distribution<double> distu;
    std::normal_distribution<double> distn;

    for (size_t c = 0; c < nc; ++c) {
        for (size_t r = 0; r < nr; ++r) {
            if (distu(generator) <= density) {
                i.push_back(r);
                j.push_back(c);
                x.push_back(distn(generator));
            }
        }
    }

    auto indptrs = tatami::compress_sparse_triplets<false>(nr, nc, x, i, j);
    tatami::ArrayView<double> x_view (x.data(), x.size());
    tatami::ArrayView<int> i_view (i.data(), i.size());
    tatami::ArrayView<size_t> p_view (indptrs.data(), indptrs.size());
    std::shared_ptr<tatami::NumericMatrix> mat(new tatami::CompressedSparseColumnMatrix<double, int, decltype(x_view), decltype(i_view), decltype(p_view)>(nr, nc, x_view, i_view, p_view));

    // Doing the naive linear iteration that's currently in tatami.
    std::vector<double> xbuffer(mat->ncol());
    std::vector<int> ibuffer(mat->ncol());

    auto start = std::chrono::high_resolution_clock::now();
    auto wrk = mat->sparse_row();
    int sum = 0;
    for (int r = 0; r < mat->nrow(); ++r) {
        auto range = wrk->fetch(r, xbuffer.data(), ibuffer.data());
        sum += range.number;
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Linear access time: " << duration.count() << " for " << sum << " non-zero elements" << std::endl;

    // Comparing to a priority queue.
    start = std::chrono::high_resolution_clock::now();

    typedef std::pair<int, int> HeapElement;
    std::priority_queue<HeapElement, std::vector<HeapElement>, std::greater<HeapElement> > next_heap;
    std::vector<int> hits, tmp_hits;
    auto state = indptrs;

    hits.resize(mat->ncol());
    for (int c = 0; c < mat->ncol(); ++c) { // force everything to be re-searched on initialization.
        hits[c] = c;
        --state[c];
    }

    sum = 0;
    for (int r = 0; r < mat->nrow(); ++r) {
        tmp_hits.swap(hits);
        hits.clear();

        for (auto x : tmp_hits) {
            ++state[x];
            if (state[x] < indptrs[x + 1]) {
                int current = i[state[x]];
                if (current == r) {
                    hits.push_back(x);
                } else {
                    next_heap.emplace(current, x);
                }
            }
        }

        while (!next_heap.empty()) {
            auto current_secondary = next_heap.top().first;
            if (current_secondary > r) {
                break;
            }
            auto current_primary_index = next_heap.top().second;
            next_heap.pop();

            if (current_secondary < r) {
                ++state[current_primary_index];
                if (state[current_primary_index] < indptrs[current_primary_index + 1]) {
                    int next_secondary = i[state[current_primary_index]];
                    if (next_secondary == r) {
                        hits.push_back(current_primary_index);
                    } else {
                        next_heap.emplace(next_secondary, current_primary_index);
                    }
                }
            } else {
                hits.push_back(current_primary_index);
            }
        }

        // Copy indices over for a fair comparison.
        for (size_t h = 0; h < hits.size(); ++h) {
            ibuffer[h] = hits[h];
            xbuffer[h] = x[state[h]];
        }

        sum += hits.size();
    } 

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Priority queue access time: " << duration.count() << " for " << sum << " non-zero elements" << std::endl;

    return 0;
}

