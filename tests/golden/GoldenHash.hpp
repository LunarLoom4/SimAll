// =============================================================================
// SimAll Beta - tests/golden/GoldenHash.hpp
// Week 19 - Golden-hash registry.  Locks the FNV-1a 64-bit hash of
// post-processed reference vectors so any future numerical drift in the
// post-processor immediately trips a CI failure.  The reference values
// themselves live in regression/cases/*.hpp; here we only record their
// quantised hashes (6-digit precision).
// =============================================================================
#pragma once

#include "../regression/RegressionFramework.hpp"

namespace simall::golden
{

struct GoldenEntry
{
    const char* name;
    std::uint64_t hash;
};

} // namespace simall::golden
