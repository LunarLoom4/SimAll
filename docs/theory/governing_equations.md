# SimAll Beta — Theory Guide: Governing Equations

> Version 1.0.0-rc1.  This document catalogues the partial differential
> equations and constitutive closures implemented by the SimAll Beta
> solver subsystem.  References are abbreviated; the full bibliography
> appears in `docs/validation/report.md`.

## 1. Conservation laws

### 1.1 Mass
$$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0
$$

### 1.2 Momentum (Navier-Stokes)
$$
\frac{\partial (\rho \mathbf{u})}{\partial t} +
  \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) =
  -\nabla p + \nabla \cdot \boldsymbol{\tau} + \rho \mathbf{g} + \mathbf{f}_\text{src}
$$
with the Newtonian viscous stress
$$
\boldsymbol{\tau} = \mu \left( \nabla \mathbf{u} + (\nabla \mathbf{u})^\mathsf{T} \right)
  - \tfrac{2}{3} \mu (\nabla \cdot \mathbf{u}) \mathbf{I}.
$$

### 1.3 Energy (total enthalpy)
$$
\frac{\partial (\rho h_t)}{\partial t} + \nabla \cdot (\rho h_t \mathbf{u}) =
  \frac{\partial p}{\partial t} + \nabla \cdot (k \nabla T)
  + \nabla \cdot (\boldsymbol{\tau} \cdot \mathbf{u}) + S_h .
$$

### 1.4 Species transport (Y_k for k = 1..N_s − 1)
$$
\frac{\partial (\rho Y_k)}{\partial t} + \nabla \cdot (\rho Y_k \mathbf{u}) =
  \nabla \cdot (\rho D_k \nabla Y_k) + \dot{\omega}_k .
$$

## 2. Turbulence closures

| Family                  | Implemented models |
|-------------------------|--------------------|
| RANS — eddy viscosity   | k-ε (Std/RNG/Realizable), k-ω (Std/SST/BSL), SA |
| RANS — Reynolds stress  | LRR Reynolds-Stress Model (LRR-RSM)               |
| Hybrid RANS-LES         | SAS-SST, IDDES (with SA backbone)                  |
| LES                     | Smagorinsky (static + dynamic), WALE               |
| Transition              | k-kL-ω, γ-Re_θ                                     |

Each model resides in `src/turbulence/` and is selectable at runtime via
the `TurbulenceFactory`.

## 3. Multiphase

* **VOF / IsoAdvector** — sharp interface, geometric reconstruction.
* **Eulerian-Eulerian multifluid** — N_p phases, drag/lift/virtual-mass closures.
* **Mixture model + drift-flux** — fast bubbly/slurry approximations.
* **Cavitation** — Kunz, Zwart-Gerber-Belamri.
* **Wall boiling** — RPI partitioning (heat-flux split).

## 4. Combustion

* **FGM** flamelet-generated manifold tables (Z, χ, h).
* **SLFM** steady-laminar flamelet for diffusion flames.
* **Transported-PDF** with IEM micromixing.
* **G-equation** premixed flame-front tracking.
* **TFC** turbulent flame closure.
* **NOx** thermal / prompt / fuel route bookkeeping.
* **Soot** Moss-Brookes and MOMIC (method of moments).

## 5. Radiation

* **Rosseland** (optically thick).
* **Monte-Carlo RTE** (general).
* **Surface-to-surface** view-factor exchange.
* **Solar load** direct + diffuse.
* **WSGG** weighted-sum-of-gray-gases spectral model.

## 6. Reduced-order modelling (ROM)

* **POD/DMD** baseline.
* **Operator inference** (Peherstorfer & Willcox, 2016) — fits Â, B̂, Ĥ
  on snapshot data, no full-order operator access required.
* **ECSW** energy-conserving sampling and weighting hyper-reduction
  (Farhat et al., 2014) via Lawson-Hanson NNLS.

## 7. Adjoint & optimisation

* **Discrete adjoint** — GMRES on Aᵀ with Jacobi preconditioning.
* **Continuous adjoint** — surface and volume shape-sensitivity assembly.
* **MMA** (Svanberg method of moving asymptotes).
* **CMA-ES** with rank-1 + rank-μ updates.
* **DOE** — full-factorial, Latin hypercube, Sobol.
* **Kriging** anisotropic squared-exponential surrogate.

## 8. Discretisation

* Cell-centred finite-volume on arbitrary polyhedra.
* Convection schemes: 1st-order upwind, 2nd-order linear-upwind, central,
  QUICK, MUSCL, MINMOD, van Leer, van Albada, superbee, Koren,
  flux-corrected transport (Zalesak limiter).
* Pressure-velocity coupling: SIMPLE, SIMPLEC, PISO, coupled (block).
* Time integration: 1st/2nd-order BDF, Crank-Nicolson, dual-time stepping
  for steady acceleration of unsteady solver.

## 9. Verification

The `tests/verification/` directory ships a method-of-manufactured-solutions
suite proving second-order convergence of the 5-point Laplacian on a
smooth $\sin(\pi x)\sin(\pi y)$ field; this is the canonical V&V gate.
