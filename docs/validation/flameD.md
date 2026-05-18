# Validation — Sandia Flame D

Piloted CH4/air diffusion flame (Barlow & Frank 1998).

| Quantity | Reference | Tolerance |
|---|---|---|
| Radial $\bar T$ at $x/D = 7.5, 15, 30, 45$ | TNF dataset | 50 K |
| Mixture-fraction variance | TNF | 10 % |
| Major species $Y_{CO_2}, Y_{H_2O}, Y_{O_2}$ | TNF | 15 % |

Run via `applications/simall_flameD` with the FPV combustion model
and GRI-3.0 mechanism (`data/mech/gri30.cti`).
