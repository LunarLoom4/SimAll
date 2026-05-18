# Threading Model

Per spec §1.3 the GUI thread MUST NEVER block.

## Three thread classes

| Class | Created by | Purpose | Lifetime |
|---|---|---|---|
| Qt GUI thread | `QApplication` | event loop, painting, dock layout | application |
| Worker pool   | `core::ThreadPool` | meshing, CAD healing, solver iteration, postprocessing | application |
| MPI rank      | `MPI_Init`         | one per process — owns its FieldRegistry partition | application |

## Cross-thread protocol

- Worker → GUI:  emit a `Qt::QueuedConnection` signal or post a
  `core::Event` on the `EventBus` (auto-marshalled to GUI thread).
- GUI → Worker:  enqueue a `core::ICommand` via the project's
  `CommandStack`; the stack dispatches to the pool.
- Worker → Worker:  use `core::TaskGraph` futures.

## Forbidden patterns

- Mutating any `FieldRegistry` array from the GUI thread.
- Reading a worker's in-progress matrix while it is being assembled.
- Calling `MPI_*` from anywhere but the rank's owning worker.

## Pool sizing

`core::ThreadPool` defaults to `std::thread::hardware_concurrency() - 2`
(reserving one for the GUI and one for the kernel). Override via
`SIMALL_THREADS` environment variable.
