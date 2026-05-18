# Theory — Turbulence Models

All turbulence models live under [src/turbulence/](../../src/turbulence)
and implement `turbulence::ITurbulenceModel`.

## RANS

| Model | Equations | File |
|---|---|---|
| Spalart–Allmaras           | $\tilde\nu$ | `SpalartAllmaras` |
| Standard $k\!-\!\varepsilon$ | $k,\varepsilon$ | `KEpsilon` |
| Realisable $k\!-\!\varepsilon$ | $k,\varepsilon$ | `RealizableKEpsilon` |
| Standard $k\!-\!\omega$      | $k,\omega$ | `KOmega` |
| SST $k\!-\!\omega$           | $k,\omega$ | `KOmegaSST` |
| Reynolds-stress (LRR / SSG)  | 6 + $\varepsilon$ | `ReynoldsStress` |
| $\gamma\text{-}Re_{\theta t}$ transition | $\gamma, Re_{\theta t}$ | `GammaReThetaTransition` |

## Hybrid RANS-LES

- DES, DDES, IDDES on top of SA or SST (`SaDdes`, `SstIddes`).
- Embedded LES with synthetic turbulence inlet (`SyntheticEddyMethod`).

## LES

- Smagorinsky, dynamic Smagorinsky, WALE, sigma model
  (`SmagorinskySGS`, `DynamicSmagorinskySGS`, `WaleSGS`).
- Wall-modelled LES via algebraic wall-stress (`WallFunctions`).

## Coupling with the solver

Each model owns its transport equations as `solver::ScalarTransport`
instances. After every outer iteration it computes
$\mu_t = f(\text{state})$ and writes it to `FieldRegistry`'s `mu_t`
field — momentum equation reads it on the next iteration.
