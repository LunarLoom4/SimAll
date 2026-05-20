// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/MeshFormats.cpp
// =============================================================================
#include "io/MeshFormats.hpp"

namespace simall::io
{

std::uint8_t vertices_per_element(ElementType t) noexcept
{
    switch (t) {
    case ElementType::Bar2:
        return 2;
    case ElementType::Tri3:
        return 3;
    case ElementType::Quad4:
        return 4;
    case ElementType::Tetra4:
        return 4;
    case ElementType::Pyra5:
        return 5;
    case ElementType::Penta6:
        return 6;
    case ElementType::Hexa8:
        return 8;
    case ElementType::Poly:
        return 0; // variable
    default:
        return 0;
    }
}

const char* element_type_name(ElementType t) noexcept
{
    switch (t) {
    case ElementType::Bar2:
        return "Bar2";
    case ElementType::Tri3:
        return "Tri3";
    case ElementType::Quad4:
        return "Quad4";
    case ElementType::Tetra4:
        return "Tetra4";
    case ElementType::Pyra5:
        return "Pyra5";
    case ElementType::Penta6:
        return "Penta6";
    case ElementType::Hexa8:
        return "Hexa8";
    case ElementType::Poly:
        return "Poly";
    default:
        return "Unknown";
    }
}

std::size_t ElementSection::element_count() const noexcept
{
    if (type == ElementType::Poly) {
        return polyOffsets.empty() ? 0 : polyOffsets.size() - 1;
    }
    const auto v = vertices_per_element(type);
    return v ? nodes.size() / v : 0;
}

std::size_t ImportedMesh::total_nodes() const noexcept
{
    std::size_t n = 0;
    for (auto& z : zones)
        n += z.x.size();
    return n;
}

std::size_t ImportedMesh::total_elements() const noexcept
{
    std::size_t n = 0;
    for (auto& z : zones)
        for (auto& s : z.sections)
            n += s.element_count();
    return n;
}

std::size_t StructuredBlock::point_count() const noexcept
{
    return std::size_t(ni) * std::size_t(nj) * std::size_t(nk);
}
std::size_t StructuredBlock::cell_count() const noexcept
{
    if (ni == 0 || nj == 0 || nk == 0)
        return 0;
    const std::uint32_t ci = ni > 1 ? ni - 1 : 1;
    const std::uint32_t cj = nj > 1 ? nj - 1 : 1;
    const std::uint32_t ck = nk > 1 ? nk - 1 : 1;
    return std::size_t(ci) * std::size_t(cj) * std::size_t(ck);
}

} // namespace simall::io
