# Hydro-fluks: linearised Bjorken hydrodynamics with B, Q, S

## Physics

Perturbations (`delta B, delta Q, delta S, delta e, delta u`, with `u` the longitudinal flow rapidity-velocity perturbation) around a
boost-invariant Bjorken background are evolved in proper time `t` and spatial rapidity `eta` (grid `Nl = 501` points on `[-5, 5]`,
`deta = 0.02`, periodic). With `h = deta`, `d1 f = (f[i+1] - f[i-1]) / 2h`, `d2 f = (f[i+1] + f[i-1] - 2 f[i]) / h^2`,
and all coefficients evaluated at time `t`:

```
dt B = -B/t + D/t^2 d2 B + (rhoB + D rhoB') / t * d1 u
dt Q = -Q/t + D/t^2 d2 Q + (rhoQ + D rhoQ') / t * d1 u
dt S = -S/t + D/t^2 d2 S
dt e = -(e + dP)/t - (M - 4 etas/(3t)) / t * d1 u,        dP = cs2 e + dBP B + dBP Q + dSP S
dt u = -(2/t + M'/M) u - (cs2 d1 e + dBP d1 B + dQP d1 Q + dSP d1 S) / (M t) + 4 etas/(3 t^2 M) d2 u
```

`rhoB(t) = rhoB0 * t0/t`, `rhoQ = 0.4 rhoB`, `t0 = 1` (set in `bj_flow.h`). `D`, `M` (`Mst`), `etas`, `cs2 = dp/de`, `dBP = dp/drhoB`,
`dQP`, `dSP` are polynomials in `t` read from table files (below).

The Green's function drivers put a Gaussian (`sigma^2 = 0.01`) source `delta(eta)/t1` into one field at start time `t1`
(500 values in `[1, 11]`), evolve to `t = 11`, and store the response of all five fields.

## Code

| File | Role |
|---|---|
| `bj_flow.h/.cpp` | Explicit solver (forward Euler, central differences, `dt = 2e-5`) |
| `bj_flow_cn.h/.cpp` | **Crank-Nicolson** solver, derived from `bj_flow` (same coefficients, same stencil) |
| `Hydro-fluks.cpp` | Driver for the explicit solver, writes `Green-<src>-<resp>.dat` |
| `bj_flow_rk4.h/.cpp` | **Runge-Kutta 4** solver, derived from `bj_flow_cn` (shares the operator `A(t)`) |
| `Hydro-fluks-cn.cpp` | Same output with CN: `./Hydro-fluks-cn [dt] [Ntstart] [table_prefix] [out_dir]`; compiled with `-DUSE_RK4` (`make Hydro-fluks-rk4`) it runs RK4 |
| `tests/test_cn.cpp`, `tests/test_rk4.cpp` | Validation of CN and RK4 (`make test`, `make test-rk4`) |
| `compare/` | Three-solver Green-function comparison and plot script (`make compare`), figures in `docs/figures/` |

### Table files

`<prefix>-{Mst,dMst,etas,dpde,dpdrhoB,dpdrhoQ,dpdrhoS,D}.dat`: 11 numbers each, polynomial coefficients **highest power first**
(`value(t) = sum_k c[10-k] t^k`). `prefix` is `1n0` by default (`bj_flow("name")` changes it). `dMst` is read but the code
differentiates the `Mst` polynomial itself.
### The `1n0` tables (`tables/`, `make tables`)

`tables/make_1n0_tables.py` builds the eight coefficient files from the HRG trajectory tables `T, e, p, mu, cB2, cQ2, cS2 -n-1n0.dat`
(from `eos-3mu.ipynb`; `t = 1 ... 11 fm`, 1000 points, units fm^-1 / fm^-4), following the recipe of `eta2tau13/eos-polynomial-fits.ipynb`
(degree-10 polynomial in `t`, weighted for uniform relative accuracy, fitted on Chebyshev nodes):

* trajectory: ideal Bjorken flow, `rhoB = 0.16/t fm^-3`, `rhoQ = 0.4 rhoB`, `rhoS = 0` (checked: `t de/dt = -(e+p)` to 2e-4, `s t` constant);
* `s = (e + p - muB rhoB)/T` (entropy density); `eta = 0.197 * 2/(4 pi) * s` (eta/s = 2/(4 pi)); `Mst = 0.197 (e+p) - (4/3) eta/t`; `D = 2/(2 pi T)`;
* `dpdrho{B,Q,S} = 0.197 * c{B,Q,S}2` (assuming `cB2, cQ2, cS2` are `dp/drho` in fm^-1); `dpde = cs2 = -(t dp/dt + cB2 rhoB + cQ2 rhoQ)/(e+p)` (0.145 ... 0.149).

Fit errors (relative): `Mst` 2e-3, `etas` 3e-3, `dpde` 3e-4, `D` 4e-5, `dpdrho*` below 1e-3, `dMst` 1.5e-2 (the solver differentiates the `Mst` polynomial instead).
The `Mst` polynomial derivative deviates from the data derivative by up to 0.3 GeV/fm^3/fm, mostly at the late-time end (`t = 11`: -0.040 vs -0.050).

Caveats found in the notebooks (not changed): `s-n-1n0.dat` is not the entropy density (it is `e + p + muB rhoB`, without division by `T`);
`muQS` in `eos-3mu.ipynb` sets `muQ = 0.4 rhoB`, `muS = 0` instead of solving for `rhoQ = 0.4 rhoB`, `rhoS = 0` (small effect);
the `usecols` for `pdg2020` do not match that file's column layout. The tables were therefore derived from the existing outputs and the HRG was not recomputed.
`tests/test-*.dat` (0n0 tables plus synthetic constant `dpdrho*`) remain for the numerical tests.

## Crank-Nicolson scheme

The system is `dy/dt = A(t) y` with `y = (B, Q, S, e, u)`. Each step solves

```
(I - dt/2 A(t + dt/2)) y^{n+1} = (I + dt/2 A(t + dt/2)) y^n
```

* second-order in `dt`, unconditionally stable (the explicit scheme needs `dt ~ 2e-5` because forward Euler with central differences is only marginally stable for the wave-like part);
* `A` has constant coefficients in `eta`, so the left side is a **periodic block-tridiagonal matrix with 5x5 blocks**; it is solved in `O(Nl)`
  by block Thomas elimination with the cyclic closure treated through the last unknown (`bj_flow_cn.cpp`, `step()`);
* the last step of `evolve_to(t_end)` is shortened to land on `t_end`.

### Validation (`tests/test_cn`, `t = 1 -> 11` via `./test_cn 11` (`make test` stops at 6), delta source in each of the five fields, tables in `tests/`)

| CN `dt` | difference to CN(`dt = 0.00125`) | ratio to previous |
|---|---|---|
| 0.04 | 1.1e-4 | - |
| 0.02 | 2.7e-5 | 4.0 |
| 0.01 | 6.8e-6 | 4.1 |
| 0.005 | 1.6e-6 | 4.2 |

Error falls 4x per halved step (second order). The difference between CN and the explicit code saturates at about 2.3e-5 (relative L2 over all fields):
that is the first-order error of the explicit scheme itself. Run time for one start time: explicit 2.2 s, CN with `dt = 0.01` 0.21 s.

## Runge-Kutta 4

`bj_flow_rk4` integrates `y' = A(t) y` with the classical RK4 scheme (stages at `t`, `t + dt/2` twice, `t + dt`), using the same operator blocks as the CN solver.

* 4th order in `dt`, but **explicit**: `dt` must stay below the diffusion stability limit, roughly `2.78 h^2 t^2 / (4 max(D, 4 etas/(3M)))`
  (about 1.2e-3 at `t = 1` for the test tables, growing like `t^2`; the default `dt = 2e-4` has a safety factor of about 6).
* With the test tables RK4 at `dt = 1.1e-3 ... 2e-4` agrees with RK4 at `dt = 5e-5` to about 1e-12 (relative L2, `t = 1 -> 6`), i.e. the time-integration error is at round-off level;
  the remaining 2.4e-5 difference to the explicit code is the forward-Euler error. (`make test-rk4`)
* Cost: about 4 operator applications per step, so RK4 at `dt = 5e-4` runs about 3.6x faster than forward Euler at `dt = 2e-5`; Crank-Nicolson at `dt = 0.01` is faster still (about 11x) but only second order.

## Solver comparison (`make compare`)

Green functions at `t = 11` for a delta source in each field, started at `t1 = 1` and `t1 = 6`; forward Euler (`dt = 2e-5`), CN (`dt = 0.01`), RK4 (`dt = 5e-4`).
Figures (with their parameters) are in `docs/figures/`: `green_comparison_t1_1.png`, `green_comparison_t1_6.png` (5x5 Green-function matrices with the three solvers overlaid)
and `green_comparison_differences.png` (max difference per source/response pair, relative to the largest Euler value).
Differences between solvers are at most 3.4e-5 (`t1 = 1`) and 1.1e-5 (`t1 = 6`), run times for ten runs 17.6 s (Euler), 4.9 s (RK4), 1.6 s (CN).
They use the real `1n0` tables (`PREFIX=tables/1n0`, the default of `make compare`).

## Things to know

* Grid is periodic (point `Nl-1` is a neighbour of point 0) although `initial_conditions` zeroes the end points; harmless for localised sources.
* **Fixed:** `dP = cs2 de + dBP dB + dQP dQ + dSP dS`. The explicit code used to have `dBP` in front of `Q` (typo; `|dQP|` is only about 1% of `|dBP|` in the `1n0` tables); all solvers now use `dQP`. `fix_dQP_in_pressure = false` in `bj_flow_cn`/`bj_flow_rk4` reproduces the old behaviour. Effect on the Green functions (CN, relative to the maximum): 1.7% (`t1 = 1`) and 0.7% (`t1 = 6`) for the `Q` source, at most 1.2e-4 for the other sources. The comparison figures were regenerated with the fix.
* `Hydro-fluks.cpp` uses `PI = 3.14` in the Gaussian; `Hydro-fluks-cn.cpp` does too, to stay comparable.

## Publications

* O. Savchuk and S. Pratt, *Correlations of Conserved Quantities at Finite Baryon Density*, Phys. Rev. C 109, 024910 (2024), arXiv:2311.02046.
* O. Savchuk, *Sensitivity of transverse momentum correlations to early-stage and thermal fluctuations*, arXiv:2402.12504.
