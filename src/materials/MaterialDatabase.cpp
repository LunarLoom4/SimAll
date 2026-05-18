#include "materials/Material.hpp"

namespace simall::materials {

Material& MaterialDatabase::add(std::string name) {
    auto [it, _] = store_.emplace(name, Material{name});
    return it->second;
}

Material* MaterialDatabase::find(const std::string& name) {
    auto it = store_.find(name);
    return it == store_.end() ? nullptr : &it->second;
}

}  // namespace simall::materials
