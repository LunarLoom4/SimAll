// =============================================================================
// SimAll Beta - Optimization Subsystem
// File   : src/optimization/Doe.cpp
// =============================================================================
#include "optimization/Doe.hpp"

#include <algorithm>
#include <array>
#include <numeric>
#include <random>

namespace simall::optimization
{

std::vector<std::vector<double>> doe_full_factorial(const std::vector<double>& xLow,
                                                    const std::vector<double>& xUp,
                                                    std::size_t levels)
{
    const std::size_t d = xLow.size();
    if (d == 0 || xUp.size() != d || levels < 1)
        return {};
    std::size_t total = 1;
    for (std::size_t i = 0; i < d; ++i)
        total *= levels;
    std::vector<std::vector<double>> out(total, std::vector<double>(d, 0.0));
    for (std::size_t idx = 0; idx < total; ++idx) {
        std::size_t rem = idx;
        for (std::size_t i = 0; i < d; ++i) {
            const std::size_t lev = rem % levels;
            rem /= levels;
            const double frac = (levels == 1) ? 0.5 : double(lev) / double(levels - 1);
            out[idx][i] = xLow[i] + frac * (xUp[i] - xLow[i]);
        }
    }
    return out;
}

std::vector<std::vector<double>> doe_latin_hypercube(const std::vector<double>& xLow,
                                                     const std::vector<double>& xUp,
                                                     std::size_t n,
                                                     std::uint64_t seed)
{
    const std::size_t d = xLow.size();
    if (d == 0 || xUp.size() != d || n == 0)
        return {};
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    std::vector<std::vector<std::size_t>> perm(d, std::vector<std::size_t>(n));
    for (std::size_t i = 0; i < d; ++i) {
        std::iota(perm[i].begin(), perm[i].end(), std::size_t{0});
        std::shuffle(perm[i].begin(), perm[i].end(), rng);
    }
    std::vector<std::vector<double>> out(n, std::vector<double>(d, 0.0));
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t i = 0; i < d; ++i) {
            const double t = (double(perm[i][k]) + u01(rng)) / double(n);
            out[k][i] = xLow[i] + t * (xUp[i] - xLow[i]);
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// Sobol sequence — uses Joe-Kuo direction numbers for the first 6 dimensions.
// Higher-D falls back to LHS to remain library-free.
// ---------------------------------------------------------------------------
namespace
{

constexpr int kSobolMaxBits = 30;

// Polynomials and m-initial values per dim 0..5 (Joe-Kuo Table 1).
const std::array<std::array<unsigned, kSobolMaxBits>, 6> kSobolM = {{
    /* d=0 */ {
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}},
    /* d=1 */ {{1,        3,        5,        15,        17,        51,       85,      255,
                257,      771,      1285,     3855,      4369,      13107,    21845,   65535,
                65537,    196611,   327685,   983055,    1114129,   3342387,  5570645, 16711935,
                16843009, 50529027, 84215045, 252645135, 286331153, 858993459}},
    /* d=2 */ {{1, 1, 7, 11, 13, 61, 67, 79, 465, 721, 823, 4091, 4125, 4141, 28723, 45311, 53505}},
    /* d=3 */ {{1, 3, 7, 5, 7, 43, 49, 147, 439, 1013, 727, 987, 5889, 6915, 16647}},
    /* d=4 */ {{1, 1, 5, 3, 15, 51, 125, 141, 177, 759, 267, 1839, 6929, 16241, 16565}},
    /* d=5 */ {{1, 3, 1, 1, 9, 59, 25, 89, 321, 835, 833, 4033}},
}};

void compute_sobol(std::size_t n, std::size_t d, std::vector<std::vector<double>>& out)
{
    if (d > 6) {
        out.clear();
        return;
    }
    std::vector<std::vector<std::uint32_t>> V(d, std::vector<std::uint32_t>(kSobolMaxBits + 1, 0));
    for (std::size_t k = 0; k < d; ++k)
        for (int j = 1; j <= kSobolMaxBits; ++j)
            V[k][j] = kSobolM[k][j - 1] << (kSobolMaxBits - j);

    out.assign(n, std::vector<double>(d, 0.0));
    std::vector<std::uint32_t> X(d, 0);
    for (std::size_t i = 1; i <= n; ++i) {
        std::uint32_t c = 1;
        std::uint32_t v = i;
        while (v & 1) {
            v >>= 1;
            ++c;
        }
        for (std::size_t k = 0; k < d; ++k) {
            X[k] ^= V[k][c];
            out[i - 1][k] = double(X[k]) / double(1ULL << kSobolMaxBits);
        }
    }
}

} // namespace

std::vector<std::vector<double>> doe_sobol(const std::vector<double>& xLow,
                                           const std::vector<double>& xUp,
                                           std::size_t n)
{
    const std::size_t d = xLow.size();
    if (d == 0 || xUp.size() != d || n == 0)
        return {};
    if (d > 6) {
        return doe_latin_hypercube(xLow, xUp, n); // graceful fallback
    }
    std::vector<std::vector<double>> raw;
    compute_sobol(n, d, raw);
    for (auto& row : raw)
        for (std::size_t i = 0; i < d; ++i)
            row[i] = xLow[i] + row[i] * (xUp[i] - xLow[i]);
    return raw;
}

} // namespace simall::optimization
