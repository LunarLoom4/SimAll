// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/FieldRegistry.hpp
// Phase  : 6 — solver-side field storage. Mesh-aware; lives separately from
// meshing::Mesh because field identity is solver-domain (ρ, p, k, ε, Yi…).
// =============================================================================
#pragma once

#include "utilities/AlignedAllocator.hpp"

#include <string>
#include <unordered_map>

namespace simall::solver
{

using ScalarField = util::aligned_vector<double>;

struct VectorField
{
    ScalarField x, y, z;
    void resize(std::size_t n)
    {
        x.assign(n, 0);
        y.assign(n, 0);
        z.assign(n, 0);
    }
    std::size_t size() const { return x.size(); }
};

class FieldRegistry
{
public:
    ScalarField& scalar(const std::string& name, std::size_t n)
    {
        auto& f = scalars_[name];
        if (f.size() != n)
            f.assign(n, 0.0);
        return f;
    }
    VectorField& vector(const std::string& name, std::size_t n)
    {
        auto& f = vectors_[name];
        if (f.size() != n)
            f.resize(n);
        return f;
    }
    ScalarField* find_scalar(const std::string& n)
    {
        auto it = scalars_.find(n);
        return it == scalars_.end() ? nullptr : &it->second;
    }
    VectorField* find_vector(const std::string& n)
    {
        auto it = vectors_.find(n);
        return it == vectors_.end() ? nullptr : &it->second;
    }

    /// Iteration accessors (used by checkpoint I/O and post-processing).
    const std::unordered_map<std::string, ScalarField>& scalars() const noexcept
    {
        return scalars_;
    }
    const std::unordered_map<std::string, VectorField>& vectors() const noexcept
    {
        return vectors_;
    }

private:
    std::unordered_map<std::string, ScalarField> scalars_;
    std::unordered_map<std::string, VectorField> vectors_;
};

} // namespace simall::solver
