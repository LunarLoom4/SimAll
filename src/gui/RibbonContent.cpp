// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/RibbonContent.cpp
// =============================================================================
#include "gui/RibbonContent.hpp"

#include "gui/RibbonBar.hpp"

namespace simall::gui
{

namespace
{
void addBtn(RibbonActions& acts, RibbonBar* r, const QString& tab, const QString& label)
{
    QAction* a = r->add_button(tab, label);
    if (a)
        acts.byName.emplace(label, a);
}
} // namespace

RibbonActions populate_default_ribbon(RibbonBar* r)
{
    RibbonActions a;
    // 1. File
    addBtn(a, r, "File", "New Project");
    addBtn(a, r, "File", "Open Project");
    addBtn(a, r, "File", "Save Project");
    addBtn(a, r, "File", "Save As…");
    addBtn(a, r, "File", "Import");
    addBtn(a, r, "File", "Export");
    addBtn(a, r, "File", "Preferences");
    addBtn(a, r, "File", "Exit");

    // 2. Geometry
    addBtn(a, r, "Geometry", "Import CAD");
    addBtn(a, r, "Geometry", "Heal");
    addBtn(a, r, "Geometry", "Defeature");
    addBtn(a, r, "Geometry", "Booleans");
    addBtn(a, r, "Geometry", "Feature Edges");
    addBtn(a, r, "Geometry", "Tessellate");
    addBtn(a, r, "Geometry", "Persistent IDs");

    // 3. Mesh
    addBtn(a, r, "Mesh", "Surface Mesh");
    addBtn(a, r, "Mesh", "Volume Mesh");
    addBtn(a, r, "Mesh", "Boundary Layer");
    addBtn(a, r, "Mesh", "Refinement Regions");
    addBtn(a, r, "Mesh", "Quality Report");
    addBtn(a, r, "Mesh", "Partition");
    addBtn(a, r, "Mesh", "AMR");

    // 4. Physics
    addBtn(a, r, "Physics", "Models");
    addBtn(a, r, "Physics", "Turbulence");
    addBtn(a, r, "Physics", "Heat Transfer");
    addBtn(a, r, "Physics", "Multiphase");
    addBtn(a, r, "Physics", "Combustion");
    addBtn(a, r, "Physics", "Radiation");
    addBtn(a, r, "Physics", "Particles");
    addBtn(a, r, "Physics", "Acoustics");

    // 5. Materials
    addBtn(a, r, "Materials", "Browser");
    addBtn(a, r, "Materials", "New Material");
    addBtn(a, r, "Materials", "Assign");
    addBtn(a, r, "Materials", "Import Library");

    // 6. Boundary Conditions
    addBtn(a, r, "Boundary Conditions", "Wall");
    addBtn(a, r, "Boundary Conditions", "Inlet");
    addBtn(a, r, "Boundary Conditions", "Outlet");
    addBtn(a, r, "Boundary Conditions", "Symmetry");
    addBtn(a, r, "Boundary Conditions", "Periodic");
    addBtn(a, r, "Boundary Conditions", "Porous Jump");
    addBtn(a, r, "Boundary Conditions", "Fan");
    addBtn(a, r, "Boundary Conditions", "Interface");

    // 7. Solver
    addBtn(a, r, "Solver", "Algorithm");
    addBtn(a, r, "Solver", "Schemes");
    addBtn(a, r, "Solver", "Relaxation");
    addBtn(a, r, "Solver", "Time Stepping");
    addBtn(a, r, "Solver", "Linear Solvers");
    addBtn(a, r, "Solver", "Convergence");

    // 8. Initialization
    addBtn(a, r, "Initialization", "Uniform");
    addBtn(a, r, "Initialization", "Patch");
    addBtn(a, r, "Initialization", "From File");
    addBtn(a, r, "Initialization", "Reset Fields");

    // 9. Run
    addBtn(a, r, "Run", "Start Solver");
    addBtn(a, r, "Run", "Pause");
    addBtn(a, r, "Run", "Stop");
    addBtn(a, r, "Run", "Steady Iterations");
    addBtn(a, r, "Run", "Transient");
    addBtn(a, r, "Run", "Adaptive Δt");

    // 10. Results
    addBtn(a, r, "Results", "Contours");
    addBtn(a, r, "Results", "Iso Surface");
    addBtn(a, r, "Results", "Streamlines");
    addBtn(a, r, "Results", "Vector Glyphs");
    addBtn(a, r, "Results", "Section Cut");
    addBtn(a, r, "Results", "Volume Render");
    addBtn(a, r, "Results", "Screenshot");
    addBtn(a, r, "Results", "Animation");

    // 11. Automation
    addBtn(a, r, "Automation", "Python Console");
    addBtn(a, r, "Automation", "Record Macro");
    addBtn(a, r, "Automation", "Play Macro");
    addBtn(a, r, "Automation", "Run Script");

    // 12. HPC
    addBtn(a, r, "HPC", "MPI Launcher");
    addBtn(a, r, "HPC", "GPU Devices");
    addBtn(a, r, "HPC", "Thread Affinity");
    addBtn(a, r, "HPC", "Job Submit");

    // 13. Plugins
    addBtn(a, r, "Plugins", "Manage");
    addBtn(a, r, "Plugins", "Reload");
    addBtn(a, r, "Plugins", "SDK Docs");
    return a;
}

} // namespace simall::gui
