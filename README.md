# Hydro-fluks (3-charge EoS)

C++ code for a 1D Bjorken-flow hydrodynamic model with conserved baryon, electric
charge and strangeness densities. It evolves a delta-function perturbation in each of
the conserved quantities and records the response (Green's functions / correlation
functions) in spatial rapidity. Python notebooks build the equation-of-state input.

## Contents

| File | Purpose |
|---|---|
| `bj_flow.h`, `bj_flow.cpp` | `bj_flow` class: Bjorken flow solver (grid `Nl = 501`, `eta_max = 5`) |
| `Hydro-fluks.cpp` | Driver: loops over start times, perturbs each charge in turn, writes `Green-<source>-<response>.dat` |
| `Makefile` | `make` builds `Hydro-fluks` |
| `eos-3mu.ipynb` | Hadron-resonance-gas EoS with three chemical potentials (uses `pdg2020`) |
| `lattice-eos-fits-v1/-v2.ipynb` | Polynomial/Pade fits of lattice EoS data (v2 is the later version) |
| `pdg2020` | PDG hadron list (mass, degeneracy, B, Q, S, ...) |
| `*-n-1n0.dat` | Tabulated EoS quantities (T, e, p, s, mu, chi, c_B^2, c_Q^2, c_S^2) vs. density |

## Build and run

```
make
./Hydro-fluks
```

The solver reads transport/EoS tables named `<name_rho>-{Mst,dMst,etas,dpde,dpdrhoB,dpdrhoQ,dpdrhoS,D}.dat`
(`name_rho` is set in `bj_flow.h`, default `1n0`) from the working directory. Example tables
for `0n0` are in the `eta2tau13` project.

## Known issues (left unchanged on purpose)

- `delta()` in `Hydro-fluks.cpp` uses `PI = 3.14`; use `M_PI` if the exact Gaussian normalisation matters.
- The `Noutput` variable is unused (compiler warnings).
- The notebooks load data from `../datasets/` which is not part of this repository.
