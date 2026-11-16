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
| `Hydro-fluks-cn.cpp` | Same output with CN: `./Hydro-fluks-cn [dt] [Ntstart] [table_prefix] [out_dir]` |
| `tests/test_cn.cpp` | Explicit vs CN comparison and CN convergence study (`make test`) |

### Table files

`<prefix>-{Mst,dMst,etas,dpde,dpdrhoB,dpdrhoQ,dpdrhoS,D}.dat`: 11 numbers each, polynomial coefficients **highest power first**
(`value(t) = sum_k c[10-k] t^k`). `prefix` is `1n0` by default (`bj_flow("name")` changes it). `dMst` is read but the code
differentiates the `Mst` polynomial itself.
The repository does not contain the `1n0-*` tables and the `0n0` set in `eta2tau13` has no `dpdrho*` files.
`tests/test-*.dat` are copies of the `0n0` tables plus **synthetic constant** `dpdrho*` tables, for numerical testing only (not physical).

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

## Things to know

* Grid is periodic (point `Nl-1` is a neighbour of point 0) although `initial_conditions` zeroes the end points; harmless for localised sources.
* `dP` uses `dBP` in front of `Q` (not `dQP`) in the explicit code; `bj_flow_cn` reproduces this by default, `fix_dQP_in_pressure = true` uses `dQP`. Please check which is intended.
* `Hydro-fluks.cpp` uses `PI = 3.14` in the Gaussian; `Hydro-fluks-cn.cpp` does too, to stay comparable.
