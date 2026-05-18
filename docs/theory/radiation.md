# Theory — Radiation

Models live under [src/radiation/](../../src/radiation).

| Model | Equation | File |
|---|---|---|
| Surface-to-surface (view-factor) | $J_i = \varepsilon_i E_b + (1\!-\!\varepsilon_i)\sum F_{ij}J_j$ | `S2S` |
| Discrete-Ordinates (DOM)         | RTE along $N_\theta\!\times\! N_\phi$ directions | `DiscreteOrdinates` |
| P-1 (spherical harmonics)        | $-\nabla\!\cdot\!(\Gamma\nabla G) + aG = 4a\sigma T^4$ | `P1Radiation` |
| Rosseland diffusion              | $-\nabla\!\cdot\!\big(\tfrac{16\sigma T^3}{3\beta}\nabla T\big)$ | `Rosseland` |
| Monte-Carlo (forward / reverse)  | photon bundles | `MonteCarloRadiation` |

## Coupling with energy

Each model computes a divergence of radiative flux
$\nabla\!\cdot\!\mathbf{q}_r$ and writes its negative to the energy
source channel `S_T`. Wall-radiative-flux contributions are added at
boundary faces.
