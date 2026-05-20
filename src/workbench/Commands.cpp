// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/Commands.cpp
// Phase  : 22 Pass 22.4
// =============================================================================
#include "workbench/Commands.hpp"

#include "workbench/ChangeJournal.hpp"
#include "workbench/Schematic.hpp"
#include "workbench/StateMachine.hpp"
#include "workbench/WorkflowEngine.hpp"

#include "core/Command.hpp"

#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace simall::workbench {

namespace {

void log(ChangeJournal* j, ChangeDirection dir, std::string text) {
    if (j) j->record(dir, std::move(text));
}

std::string cell_label_or(const Schematic& s, CellId id, const char* fallback) {
    if (const Cell* c = s.cell(id)) return std::string{c->label()};
    return std::string{fallback};
}

std::string port_name_or(const Schematic& s, CellId cid, PortId pid, const char* fallback) {
    if (const Cell* c = s.cell(cid)) {
        if (const CellPort* p = c->find_port(pid)) return p->name;
    }
    return std::string{fallback};
}

// ---------------------------------------------------------------------------
// AddCellCommand
// ---------------------------------------------------------------------------
class AddCellCommand final : public simall::core::ICommand {
public:
    AddCellCommand(Schematic& s, CellKind kind, std::string label,
                   CellId* out_id, ChangeJournal* journal)
        : s_(&s), kind_(kind), label_(std::move(label)),
          out_id_(out_id), journal_(journal) {}

    void execute() override {
        if (!snapshot_) {
            // First-time execute: mint a fresh cell.
            minted_id_ = s_->add_cell(kind_, label_);
        } else {
            // Redo after undo: restore the exact same cell (same id,
            // same ports, same state).
            if (!s_->restore_cell(*snapshot_)) {
                throw std::runtime_error("AddCellCommand: restore_cell failed");
            }
            snapshot_.reset();
        }
        if (out_id_) *out_id_ = minted_id_;
        log(journal_, ChangeDirection::Do, description());
    }

    void undo() override {
        const Cell* live = s_->cell(minted_id_);
        if (!live) throw std::runtime_error("AddCellCommand: cell vanished before undo");
        snapshot_ = *live;     // copy out the full Cell (incl. ports + state)
        std::vector<CellLink> dropped;
        s_->remove_cell(minted_id_, &dropped);
        // A freshly-added cell shouldn't have any links yet, but if some
        // sibling code wired it before undo we save them too so a
        // subsequent redo can put them back via restore_link.
        dropped_links_ = std::move(dropped);
        log(journal_, ChangeDirection::Undo, description());
    }

    std::string description() const override {
        std::ostringstream os;
        os << "Add " << to_string(kind_) << " cell '" << label_ << "'";
        return os.str();
    }

private:
    Schematic*     s_;
    CellKind       kind_;
    std::string    label_;
    CellId*        out_id_;
    ChangeJournal* journal_;
    CellId         minted_id_{kInvalidCellId};
    std::optional<Cell>   snapshot_;
    std::vector<CellLink> dropped_links_;
};

// ---------------------------------------------------------------------------
// RemoveCellCommand
// ---------------------------------------------------------------------------
class RemoveCellCommand final : public simall::core::ICommand {
public:
    RemoveCellCommand(Schematic& s, CellId id, ChangeJournal* journal)
        : s_(&s), id_(id), journal_(journal) {}

    void execute() override {
        const Cell* live = s_->cell(id_);
        if (!live) throw std::runtime_error("RemoveCellCommand: unknown cell");
        snapshot_     = *live;
        descr_cache_  = std::string("Remove cell '") + std::string(live->label()) + "'";
        std::vector<CellLink> dropped;
        s_->remove_cell(id_, &dropped);
        dropped_links_ = std::move(dropped);
        log(journal_, ChangeDirection::Do, description());
    }

    void undo() override {
        if (!snapshot_) throw std::runtime_error("RemoveCellCommand: nothing to restore");
        if (!s_->restore_cell(*snapshot_)) {
            throw std::runtime_error("RemoveCellCommand: id collision on restore");
        }
        for (const CellLink& l : dropped_links_) {
            (void)s_->restore_link(l);    // best-effort; missing endpoints can't happen
        }
        log(journal_, ChangeDirection::Undo, description());
    }

    std::string description() const override {
        return descr_cache_.empty() ? std::string("Remove cell") : descr_cache_;
    }

private:
    Schematic*     s_;
    CellId         id_;
    ChangeJournal* journal_;
    std::optional<Cell>   snapshot_;
    std::vector<CellLink> dropped_links_;
    std::string           descr_cache_;
};

// ---------------------------------------------------------------------------
// ConnectCommand / DisconnectCommand
// ---------------------------------------------------------------------------
class ConnectCommand final : public simall::core::ICommand {
public:
    ConnectCommand(WorkflowEngine& eng, CellLink link, ChangeJournal* journal)
        : eng_(&eng), link_(link), journal_(journal) {}

    void execute() override {
        if (!eng_->connect(link_)) {
            throw std::runtime_error("ConnectCommand: link rejected by Schematic");
        }
        log(journal_, ChangeDirection::Do, description());
    }
    void undo() override {
        eng_->disconnect(link_);
        log(journal_, ChangeDirection::Undo, description());
    }
    std::string description() const override {
        const Schematic& s = eng_->schematic();
        std::ostringstream os;
        os << "Connect "
           << cell_label_or(s, link_.from_cell, "?")
           << "."
           << port_name_or(s, link_.from_cell, link_.from_port, "?")
           << " -> "
           << cell_label_or(s, link_.to_cell, "?")
           << "."
           << port_name_or(s, link_.to_cell, link_.to_port, "?");
        return os.str();
    }
private:
    WorkflowEngine* eng_;
    CellLink        link_;
    ChangeJournal*  journal_;
};

class DisconnectCommand final : public simall::core::ICommand {
public:
    DisconnectCommand(WorkflowEngine& eng, CellLink link, ChangeJournal* journal)
        : eng_(&eng), link_(link), journal_(journal) {}

    void execute() override {
        if (!eng_->disconnect(link_)) {
            throw std::runtime_error("DisconnectCommand: unknown link");
        }
        log(journal_, ChangeDirection::Do, description());
    }
    void undo() override {
        if (!eng_->connect(link_)) {
            throw std::runtime_error("DisconnectCommand: cannot re-add link");
        }
        log(journal_, ChangeDirection::Undo, description());
    }
    std::string description() const override {
        const Schematic& s = eng_->schematic();
        std::ostringstream os;
        os << "Disconnect "
           << cell_label_or(s, link_.from_cell, "?") << " -> "
           << cell_label_or(s, link_.to_cell,   "?");
        return os.str();
    }
private:
    WorkflowEngine* eng_;
    CellLink        link_;
    ChangeJournal*  journal_;
};

// ---------------------------------------------------------------------------
// SetCellLabelCommand
// ---------------------------------------------------------------------------
class SetCellLabelCommand final : public simall::core::ICommand {
public:
    SetCellLabelCommand(Schematic& s, CellId id, std::string new_label,
                        ChangeJournal* journal)
        : s_(&s), id_(id), new_label_(std::move(new_label)), journal_(journal) {}

    void execute() override {
        Cell* c = s_->cell(id_);
        if (!c) throw std::runtime_error("SetCellLabelCommand: unknown cell");
        if (old_label_.empty()) old_label_ = c->label();
        c->set_label(new_label_);
        log(journal_, ChangeDirection::Do, description());
    }
    void undo() override {
        Cell* c = s_->cell(id_);
        if (!c) throw std::runtime_error("SetCellLabelCommand: cell vanished");
        c->set_label(old_label_);
        log(journal_, ChangeDirection::Undo, description());
    }
    std::string description() const override {
        std::ostringstream os;
        os << "Rename cell '" << old_label_ << "' -> '" << new_label_ << "'";
        return os.str();
    }
private:
    Schematic*     s_;
    CellId         id_;
    std::string    new_label_;
    std::string    old_label_;
    ChangeJournal* journal_;
};

// ---------------------------------------------------------------------------
// SetCellStateCommand
// ---------------------------------------------------------------------------
class SetCellStateCommand final : public simall::core::ICommand {
public:
    SetCellStateCommand(Schematic& s, CellId id, CellState new_state,
                        ChangeJournal* journal)
        : s_(&s), id_(id), new_state_(new_state), journal_(journal) {}

    void execute() override {
        Cell* c = s_->cell(id_);
        if (!c) throw std::runtime_error("SetCellStateCommand: unknown cell");
        old_state_ = c->state();
        c->set_state(new_state_);
        log(journal_, ChangeDirection::Do, description());
    }
    void undo() override {
        Cell* c = s_->cell(id_);
        if (!c) throw std::runtime_error("SetCellStateCommand: cell vanished");
        c->set_state(old_state_);
        log(journal_, ChangeDirection::Undo, description());
    }
    std::string description() const override {
        std::ostringstream os;
        os << "Set cell state " << to_string(old_state_)
           << " -> " << to_string(new_state_);
        return os.str();
    }
private:
    Schematic*     s_;
    CellId         id_;
    CellState      new_state_;
    CellState      old_state_{CellState::Unfulfilled};
    ChangeJournal* journal_;
};

// ---------------------------------------------------------------------------
// SetCellAdapterCommand (Pass 22.5)
// ---------------------------------------------------------------------------
class SetCellAdapterCommand final : public simall::core::ICommand {
public:
    SetCellAdapterCommand(Schematic& s, CellId id, std::string new_adapter,
                          ChangeJournal* journal)
        : s_(&s), id_(id),
          new_adapter_(std::move(new_adapter)), journal_(journal) {}

    void execute() override {
        Cell* c = s_->cell(id_);
        if (!c) throw std::runtime_error("SetCellAdapterCommand: unknown cell");
        old_adapter_ = c->adapter_id();
        c->set_adapter_id(new_adapter_);
        log(journal_, ChangeDirection::Do, description());
    }
    void undo() override {
        Cell* c = s_->cell(id_);
        if (!c) throw std::runtime_error("SetCellAdapterCommand: cell vanished");
        c->set_adapter_id(old_adapter_);
        log(journal_, ChangeDirection::Undo, description());
    }
    std::string description() const override {
        std::ostringstream os;
        os << "Bind cell adapter '"
           << (old_adapter_.empty() ? std::string("<none>") : old_adapter_)
           << "' -> '"
           << (new_adapter_.empty() ? std::string("<none>") : new_adapter_)
           << "'";
        return os.str();
    }
private:
    Schematic*     s_;
    CellId         id_;
    std::string    new_adapter_;
    std::string    old_adapter_;
    ChangeJournal* journal_;
};

}  // namespace


// ---------------------------------------------------------------------------
// Public factory functions
// ---------------------------------------------------------------------------
std::unique_ptr<simall::core::ICommand>
make_add_cell_command(Schematic& s, CellKind kind, std::string label,
                      CellId* out_minted_id, ChangeJournal* journal) {
    return std::make_unique<AddCellCommand>(s, kind, std::move(label),
                                            out_minted_id, journal);
}

std::unique_ptr<simall::core::ICommand>
make_remove_cell_command(Schematic& s, CellId id, ChangeJournal* journal) {
    return std::make_unique<RemoveCellCommand>(s, id, journal);
}

std::unique_ptr<simall::core::ICommand>
make_connect_command(WorkflowEngine& eng, CellLink link, ChangeJournal* journal) {
    return std::make_unique<ConnectCommand>(eng, link, journal);
}

std::unique_ptr<simall::core::ICommand>
make_disconnect_command(WorkflowEngine& eng, CellLink link, ChangeJournal* journal) {
    return std::make_unique<DisconnectCommand>(eng, link, journal);
}

std::unique_ptr<simall::core::ICommand>
make_set_cell_label_command(Schematic& s, CellId id, std::string new_label,
                            ChangeJournal* journal) {
    return std::make_unique<SetCellLabelCommand>(s, id, std::move(new_label), journal);
}

std::unique_ptr<simall::core::ICommand>
make_set_cell_state_command(Schematic& s, CellId id, CellState new_state,
                            ChangeJournal* journal) {
    return std::make_unique<SetCellStateCommand>(s, id, new_state, journal);
}

std::unique_ptr<simall::core::ICommand>
make_set_cell_adapter_command(Schematic& s, CellId id,
                              std::string new_adapter_id,
                              ChangeJournal* journal) {
    return std::make_unique<SetCellAdapterCommand>(
        s, id, std::move(new_adapter_id), journal);
}

}  // namespace simall::workbench
