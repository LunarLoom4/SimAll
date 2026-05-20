// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/ChangeJournal.hpp
// Phase  : 22 Pass 22.4
//
// ChangeJournal -- bounded append-only audit log of every mutation a
// workbench Command makes.  Each entry records the direction (`Do` /
// `Undo`), the command's description text, and a monotonically increasing
// sequence number for stable diffing in tests.
//
// The journal is *passive*: commands push entries into it via record();
// it never mutates Schematic / StateMachine.  Plays the same role as the
// undo stack does in core::CommandHistory but with finer per-step
// granularity -- the journal logs both the do-leg of an execute() and
// the undo-leg of an undo(), so a UI panel can render "session history"
// independently of what's currently in the redo bucket.
// =============================================================================
#pragma once

#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <mutex>
#include <string>
#include <vector>

namespace simall::workbench
{

enum class ChangeDirection : std::uint8_t
{
    Do,
    Undo
};

struct ChangeEntry
{
    std::uint64_t sequence{};
    ChangeDirection direction{ChangeDirection::Do};
    std::string description;
};

class ChangeJournal
{
public:
    explicit ChangeJournal(std::size_t capacity = 1024) noexcept : capacity_(capacity) {}

    // Append a new entry.  Older entries are dropped from the front when
    // the buffer is full -- the journal is FIFO-bounded, never blocks.
    void record(ChangeDirection dir, std::string description);

    // Optional sink invoked synchronously inside record(); useful for
    // wiring a UI listener (e.g. status bar log) without bypassing the
    // bounded ring.  Set to nullptr (default) to disable.
    using Sink = std::function<void(const ChangeEntry&)>;
    void set_sink(Sink sink) noexcept;

    [[nodiscard]] std::vector<ChangeEntry> snapshot() const;
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] std::size_t capacity() const noexcept { return capacity_; }
    void clear() noexcept;

private:
    mutable std::mutex mtx_;
    std::deque<ChangeEntry> ring_;
    std::size_t capacity_;
    std::uint64_t next_seq_{0};
    Sink sink_{};
};

} // namespace simall::workbench
