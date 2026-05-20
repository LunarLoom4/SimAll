// =============================================================================
// SimAll Beta — Scripting subsystem
// File   : src/scripting/UdfHost.cpp
// =============================================================================
#include "scripting/UdfHost.hpp"

#include <algorithm>

namespace simall::scripting
{

UdfHost& UdfHost::instance()
{
    static UdfHost s;
    return s;
}

void UdfHost::register_udf(std::string name, std::shared_ptr<IUdf> udf)
{
    std::lock_guard lk(mu_);
    udfs_[std::move(name)] = std::move(udf);
}

bool UdfHost::unregister_udf(std::string_view name)
{
    std::lock_guard lk(mu_);
    return udfs_.erase(std::string(name)) > 0;
}

bool UdfHost::has(std::string_view name) const
{
    std::lock_guard lk(mu_);
    return udfs_.find(std::string(name)) != udfs_.end();
}

std::shared_ptr<IUdf> UdfHost::get(std::string_view name) const
{
    std::lock_guard lk(mu_);
    auto it = udfs_.find(std::string(name));
    return it == udfs_.end() ? nullptr : it->second;
}

std::vector<std::string> UdfHost::names() const
{
    std::lock_guard lk(mu_);
    std::vector<std::string> n;
    n.reserve(udfs_.size());
    for (auto& [k, _] : udfs_)
        n.push_back(k);
    std::sort(n.begin(), n.end());
    return n;
}

void UdfHost::clear()
{
    std::lock_guard lk(mu_);
    udfs_.clear();
}

double UdfHost::eval_scalar(std::string_view name, double t) const
{
    auto u = get(name);
    return u ? u->eval_scalar(t) : 0.0;
}

double UdfHost::eval_scalar(std::string_view name, double t, double x, double y, double z) const
{
    auto u = get(name);
    if (!u)
        return 0.0;
    if (u->signature() == UdfSignature::ScalarT)
        return u->eval_scalar(t);
    return u->eval_scalar(t, x, y, z);
}

std::array<double, 3> UdfHost::eval_vector(
    std::string_view name, double t, double x, double y, double z) const
{
    auto u = get(name);
    return u ? u->eval_vector(t, x, y, z) : std::array<double, 3>{0, 0, 0};
}

} // namespace simall::scripting
