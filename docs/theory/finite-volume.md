# Theory — Finite-Volume Discretisation

## Conservation form

Every transport equation in SimAll Beta is integrated in cell-centred
finite-volume form:

$$
\frac{\partial}{\partial t}\!\!\int_{V_P}\!\!\rho\phi\,dV
+\!\!\oint_{\partial V_P}\!\!\!\rho\phi(\mathbf{u}\!\cdot\!\mathbf{n})\,dS
=\!\!\oint_{\partial V_P}\!\!\!\Gamma_\phi(\nabla\phi\!\cdot\!\mathbf{n})\,dS
+\!\!\int_{V_P}\!\!S_\phi\,dV
$$

The semi-discrete form per cell $P$ is

$$
a_P\phi_P + \sum_{N} a_N \phi_N = b_P
$$

where coefficients are assembled by `solver::Convection`,
`solver::Diffusion`, and `solver::TimeIntegrator`, then the right
hand side aggregates contributions from physics source channels.

## Convection schemes (`src/solver/Convection.cpp`)

| Scheme | Order | Bounded | Notes |
|---|---|---|---|
| Upwind            | 1 | yes | default for scalars without `BoundedHR` |
| Central           | 2 | no  | for DNS / LES only |
| QUICK             | 3 | no  | requires structured-ish topology |
| MUSCL + Van Leer  | 2 | yes | TVD, default for compressible |
| Linear-upwind     | 2 | no  | with `Gradient::limiter` for boundedness |

## Gradient reconstruction (`src/solver/Gradient.cpp`)

- Green-Gauss cell-based.
- Green-Gauss node-based (weighted least-squares for nodes).
- Weighted least-squares (default for skewed meshes).

## Limiters (`src/solver/Limiters.cpp`)

Barth-Jespersen, Venkatakrishnan ($\epsilon = K\bar V$), MLP-u.

## Time integration

| Scheme | Order | Use |
|---|---|---|
| Steady            | — | dual-time inner iterations |
| Backward-Euler    | 1 | always stable |
| Crank-Nicolson    | 2 | smooth transients |
| BDF2 / BDF3       | 2/3 | URANS, LES |
| Explicit RK4      | 4 | DNS, acoustics |
