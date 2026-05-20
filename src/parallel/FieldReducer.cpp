// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/FieldReducer.cpp
// =============================================================================
#include "parallel/FieldReducer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace simall::parallel
{

double FieldReducer::sum(const util::aligned_vector<double>& v, std::size_t nOwned) const
{
    const std::size_t n = std::min(nOwned, v.size());
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        s += v[i];
    return ctx_.allreduce_sum(s);
}

double FieldReducer::l2(const util::aligned_vector<double>& v, std::size_t nOwned) const
{
    const std::size_t n = std::min(nOwned, v.size());
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        s += v[i] * v[i];
    return std::sqrt(ctx_.allreduce_sum(s));
}

double FieldReducer::min(const util::aligned_vector<double>& v, std::size_t nOwned) const
{
    const std::size_t n = std::min(nOwned, v.size());
    double m = (n == 0) ? std::numeric_limits<double>::infinity() : v[0];
    for (std::size_t i = 1; i < n; ++i)
        if (v[i] < m)
            m = v[i];
    return ctx_.allreduce_min(m);
}

double FieldReducer::max(const util::aligned_vector<double>& v, std::size_t nOwned) const
{
    const std::size_t n = std::min(nOwned, v.size());
    double m = (n == 0) ? -std::numeric_limits<double>::infinity() : v[0];
    for (std::size_t i = 1; i < n; ++i)
        if (v[i] > m)
            m = v[i];
    return ctx_.allreduce_max(m);
}

double FieldReducer::mean(const util::aligned_vector<double>& v, std::size_t nOwned) const
{
    const std::size_t n = std::min(nOwned, v.size());
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        s += v[i];
    const double sg = ctx_.allreduce_sum(s);
    const double cg = ctx_.allreduce_sum(static_cast<double>(n));
    return cg > 0.0 ? sg / cg : 0.0;
}

double FieldReducer::l2_vector(const util::aligned_vector<double>& vx,
                               const util::aligned_vector<double>& vy,
                               const util::aligned_vector<double>& vz,
                               std::size_t nOwned) const
{
    const std::size_t n = std::min({nOwned, vx.size(), vy.size(), vz.size()});
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        s += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
    return std::sqrt(ctx_.allreduce_sum(s));
}

} // namespace simall::parallel
