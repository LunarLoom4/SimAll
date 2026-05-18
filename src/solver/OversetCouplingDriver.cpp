// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/OversetCouplingDriver.cpp
// =============================================================================
#include "solver/OversetCouplingDriver.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::solver {

OversetCouplingDriver::OversetCouplingDriver(
        std::vector<OversetDomain> domains,
        std::vector<OversetLink>   links,
        OversetDriverOptions       opts)
    : domains_(std::move(domains)),
      links_  (std::move(links)),
      opt_    (opts) {}

bool OversetCouplingDriver::build() {
    for (auto& L : links_) {
        if (L.donor < 0 || L.receptor < 0 ||
            L.donor    >= static_cast<int>(domains_.size()) ||
            L.receptor >= static_cast<int>(domains_.size())) {
            SIMALL_LOG_ERROR("Overset",
                "invalid donor/receptor index in OversetLink");
            return false;
        }
        if (!L.interp) L.interp = std::make_unique<meshing::OversetInterpolation>();
        const bool ok = L.interp->build(*domains_[L.donor].mesh,
                                        *domains_[L.receptor].mesh,
                                        L.holeMaskField,
                                        *domains_[L.donor].fields);
        if (!ok) {
            SIMALL_LOG_ERROR("Overset",
                "OversetInterpolation::build failed for link ",
                L.donor, " → ", L.receptor);
            return false;
        }
    }
    lastSnapshot_.assign(links_.size(), {});
    built_ = true;
    return true;
}

void OversetCouplingDriver::interpolate_all_links() {
    for (auto& L : links_) {
        FieldRegistry& dfld = *domains_[L.donor   ].fields;
        FieldRegistry& rfld = *domains_[L.receptor].fields;
        for (const auto& sn : L.scalarVars) {
            const ScalarField* src = dfld.find_scalar(sn);
            ScalarField*       dst = rfld.find_scalar(sn);
            if (!src || !dst) continue;
            L.interp->interpolate_scalar(*src, *dst);
        }
        for (const auto& vn : L.vectorVars) {
            const VectorField* src = dfld.find_vector(vn);
            VectorField*       dst = rfld.find_vector(vn);
            if (!src || !dst) continue;
            L.interp->interpolate_vector(*src, *dst);
        }
    }
}

void OversetCouplingDriver::apply_donor_values_to_bcs() {
    // For each receptor domain, push the now-interpolated cell values
    // into every OversetBc registered on that domain.  We feed a single
    // scalar channel (whichever the BC is monitoring) — the BC's apply()
    // uses it as a Dirichlet target.
    for (const auto& D : domains_) {
        if (D.oversetBcs.empty()) continue;
        // Use the first scalar field name (e.g., "p") as the channel.
        // Real production code may wire one OversetBc per channel.
        const ScalarField* phi = D.fields ? D.fields->find_scalar("p") : nullptr;
        if (!phi) continue;
        std::vector<std::pair<meshing::CellId,double>> map;
        map.reserve(D.receptorCells.size());
        for (auto c : D.receptorCells)
            map.emplace_back(c, (*phi)[c]);
        for (auto* bc : D.oversetBcs)
            if (bc) bc->setDonorValues(map);
    }
}

double OversetCouplingDriver::measure_max_change() {
    double maxChange = 0.0;
    for (std::size_t li = 0; li < links_.size(); ++li) {
        const auto& L = links_[li];
        FieldRegistry& rfld = *domains_[L.receptor].fields;
        // Flatten current state of monitored vars on receptor side.
        std::vector<double> snap;
        for (const auto& sn : L.scalarVars) {
            const ScalarField* f = rfld.find_scalar(sn);
            if (!f) continue;
            for (auto c : domains_[L.receptor].receptorCells) snap.push_back((*f)[c]);
        }
        for (const auto& vn : L.vectorVars) {
            const VectorField* f = rfld.find_vector(vn);
            if (!f) continue;
            for (auto c : domains_[L.receptor].receptorCells) {
                snap.push_back(f->x[c]);
                snap.push_back(f->y[c]);
                snap.push_back(f->z[c]);
            }
        }
        // Compute relative change vs last snapshot.
        auto& last = lastSnapshot_[li];
        if (last.size() == snap.size() && !last.empty()) {
            double num = 0.0, den = 0.0;
            for (std::size_t i = 0; i < snap.size(); ++i) {
                const double d = snap[i] - last[i];
                num += d*d;
                den += snap[i]*snap[i];
            }
            const double rel = std::sqrt(num) / std::sqrt(std::max(den, 1.0e-30));
            if (rel > maxChange) maxChange = rel;
        } else {
            // First iteration: declare large change to ensure another pass.
            maxChange = std::max(maxChange, 1.0);
        }
        last = std::move(snap);
    }
    return maxChange;
}

double OversetCouplingDriver::iterate_once() {
    if (!built_) (void)build();
    interpolate_all_links();
    apply_donor_values_to_bcs();
    for (auto& D : domains_) {
        if (D.solver) (void)D.solver->iterate();
    }
    return measure_max_change();
}

OversetReport OversetCouplingDriver::iterate_to_convergence() {
    OversetReport rep{};
    if (!built_ && !build()) return rep;
    for (int k = 0; k < opt_.maxOuterIters; ++k) {
        const double ch = iterate_once();
        rep.iters       = k + 1;
        rep.finalChange = ch;
        if (opt_.verbose) {
            SIMALL_LOG_INFO("Overset", "  outer iter ", k,
                            "  max-rel-change = ", ch);
        }
        if (ch < opt_.relTol) { rep.converged = true; break; }
    }
    return rep;
}

}  // namespace simall::solver
