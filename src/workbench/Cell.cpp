// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/Cell.cpp
// Phase  : 22 Pass 22.1
// =============================================================================
#include "workbench/Cell.hpp"

namespace simall::workbench
{

Cell::Cell(CellId id, CellKind kind, std::string label)
    : id_(id), kind_(kind), label_(std::move(label))
{
}

PortId Cell::add_input(std::string name, std::string data_type, bool required)
{
    const PortId pid = next_port_id_++;
    ports_.push_back(CellPort{
        .id = pid,
        .direction = PortDirection::Input,
        .name = std::move(name),
        .data_type = std::move(data_type),
        .required = required,
    });
    return pid;
}

PortId Cell::add_output(std::string name, std::string data_type)
{
    const PortId pid = next_port_id_++;
    ports_.push_back(CellPort{
        .id = pid,
        .direction = PortDirection::Output,
        .name = std::move(name),
        .data_type = std::move(data_type),
        .required = false,
    });
    return pid;
}

const CellPort* Cell::find_port(PortId pid) const noexcept
{
    for (const auto& p : ports_) {
        if (p.id == pid)
            return &p;
    }
    return nullptr;
}

} // namespace simall::workbench
