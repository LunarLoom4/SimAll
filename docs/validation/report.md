# SimAll Beta — Validation Report

> Version 1.0.0-rc1.  This report consolidates the ten canonical CFD
> validation cases that gate every SimAll release.  Each case is encoded
> as a header under `tests/regression/cases/`, executed by the
> `simall_regression_tests` binary, and printed by the corresponding
> driver application under `applications/simall_*`.

## 1. Lid-driven cavity (LDC)
* **Geometry:** unit square, top lid moving at U = 1.
* **Reference:** Ghia, Ghia & Shin, *J. Comp. Phys.* 48 (1982) 387.
* **Quantity:** u-velocity along vertical centreline, 17 sampling points
  at Re = 100, 400, 1000.
* **Acceptance:** max relative error ≤ 5 % against the published table.

## 2. Backward-facing step (BFS)
* **Reference:** Armaly et al., *J. Fluid Mech.* 127 (1983) 473.
* **Quantity:** primary recirculation reattachment length X_r / h.
* **Correlation:** X_r / h ≈ 0.05 + 0.0175 · Re for Re ≤ 400.
* **Acceptance:** within 8 % of correlation.

## 3. Turbulent channel
* **Reference:** Kim, Moin & Moser, *J. Fluid Mech.* 177 (1987) 133;
  Reichardt composite profile.
* **Quantity:** u+(y+) at Re_τ = 180.
* **Acceptance:** within 10 % of Reichardt at y+ ∈ [5, 100].

## 4. Taylor-Green vortex (TGV)
* **Reference:** Brachet et al., *J. Fluid Mech.* 130 (1983) 411.
* **Quantity:** kinetic-energy decay (2-D analytic) and dissipation-rate
  peak at t ≈ 9 (3-D Re = 1600).
* **Acceptance:** 2-D analytic decay reproduced to 1e-6; 3-D peak ε
  within 5 % of Brachet.

## 5. NACA-0012
* **Reference:** Thin-airfoil theory, dCl/dα = 2π / rad.
* **Quantity:** lift coefficient at α = 5°.
* **Acceptance:** within 2 %.

## 6. Sod shock-tube
* **Reference:** Sod, *J. Comp. Phys.* 27 (1978) 1.
* **Quantity:** post-shock p*, u*, ρ* at t = 0.2.
* **Acceptance:** within 1 % on a 200-cell uniform grid.

## 7. Rayleigh-Bénard (RB)
* **Reference:** Linear theory, rigid-rigid no-slip: Ra_c = 1707.762,
  k_c = 3.117 (Chandrasekhar, 1961).
* **Quantity:** onset Rayleigh number from neutral-stability sweep.
* **Acceptance:** within 5 % on a 64 × 64 mesh after 5000 steps.

## 8. Sandia Flame D
* **Reference:** Barlow & Frank, *Proc. Comb. Inst.* 27 (1998) 1087.
* **Quantity:** centreline mixture-fraction decay ξ(x/D) ≈ 5.4 / (x/D)
  in the self-similar region.
* **Acceptance:** within 8 % at x/D ∈ {15, 30, 45, 60}.

## 9. Cylinder (Re = 50–300)
* **References:** Roshko, *NACA TN-2913* (1954) for St(Re);
  Henderson, *Phys. Fluids* 7 (1995) 2102 for Cd(Re).
* **Quantities:** Strouhal number and mean drag coefficient.
* **Acceptance:** St within 3 %, Cd within 5 % averaged over ≥30 shedding
  cycles.

## 10. Hagen-Poiseuille pipe
* **References:** Hagen 1839 / Poiseuille 1841 — parabolic profile;
  Prandtl smooth-pipe correlation for turbulent friction factor.
* **Quantities:** parabolic u(r) at the centreline (u_max = 2 u_mean) and
  Darcy friction factor.
* **Acceptance:** laminar f within 2 %, turbulent f within 5 %.

## 11. Method of Manufactured Solutions (MMS)
* **Exact field:** $u(x,y,t) = \sin(\pi x)\sin(\pi y)\,e^{-2\pi^2 \nu t}$,
  homogeneous Dirichlet BCs on $[0,1]^2$.
* **Test:** observed order of accuracy of the 5-point central-difference
  Laplacian — required to land in $[1.8, 2.2]$ across $n = 33, 65, 129$.
* **Result:** second-order convergence rate observed and asserted in CI.

## 12. Golden-hash gates
FNV-1a 64-bit hashes of the quantised reference tables are recorded in
`tests/golden/test_golden_hashes.cpp` and re-asserted on every CI run.
Any silent mutation of the reference data trips an immediate failure —
this is the last line of defence against post-processor drift.

## 13. Continuous integration

| Matrix axis        | Values                       |
|--------------------|------------------------------|
| Operating system   | Windows 2022, Ubuntu 22.04   |
| Compiler           | MSVC 19.36, GCC 11+, Clang 14 |
| Configuration      | Debug, Release               |
| Sanitizers         | ASan + UBSan (separate job)  |
| Label-filtered passes | `regression`, `verification`, `golden`, `performance` (non-gating) |

The 1.0.0-rc1 release was produced from a CI run with all gating jobs
green and the performance smoke job within historical ±15 % MOps/s.
