# Validation — Taylor-Green Vortex (Re = 1600)

Inviscid initial condition decaying under viscous action in a
triply-periodic $[0, 2\pi]^3$ box.

| Quantity | Reference | Tolerance |
|---|---|---|
| Kinetic-energy dissipation $\varepsilon(t)$ | van Rees et al. 2011 (DNS) | 3 % peak |
| Enstrophy peak time | $t \approx 8 \tau$ | 5 % |
| Resolved $E(\kappa)$ | $\kappa^{-5/3}$ inertial range | qualitative |

Driven by `applications/simall_taylor_green`. Reference data live in
`data/fields/taylor_green_ref_re1600.h5`.
