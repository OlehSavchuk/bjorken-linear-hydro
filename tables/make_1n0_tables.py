"""Builds the solver input tables  tables/1n0-{Mst,dMst,etas,dpde,dpdrhoB,dpdrhoQ,dpdrhoS,D}.dat
from the HRG trajectory tables of this repository (T,e,p,mu,cB2,cQ2,cS2 -n-1n0.dat, produced by eos-3mu.ipynb).

Trajectory: ideal Bjorken expansion, rhoB = 0.16/t fm^-3, rhoQ = 0.4 rhoB, rhoS = 0, t in [1, 11] fm (1000 points).
Units in the input files: T, mu in fm^-1; e, p in fm^-4 (multiply by hbar c = 0.197 for GeV/fm^3); cB2, cQ2, cS2 in fm^-1.

Recipe (same recipe as eta2tau13/eos-polynomial-fits.ipynb for the 0n0 set: degree-10 polynomials in t, highest power first; here with weights for uniform relative accuracy):
  s      = (e + p - muB rhoB) / T          real entropy density [fm^-3]  (muQ rhoQ is ~1e-4 of e+p and neglected, rhoS = 0)
  etas   = 0.197 * 2/(4 pi) * s            shear viscosity, eta/s = 2/(4 pi)         [GeV/fm^2]
  Mst    = 0.197 (e + p) - 4/3 etas / t                                              [GeV/fm^3]
  dMst   = d Mst / dt
  D      = 2 / (2 pi T)                    charge diffusion, D T = 2/(2 pi)          [fm]
  dpdrho{B,Q,S} = 0.197 * c{B,Q,S}2                                                  [GeV]
  dpde   = cs2 = -(t dp/dt + cB2 rhoB + cQ2 rhoQ) / (e + p)   from  dp = cs2 de + sum_i dp/drho_i drho_i  with
           de/dt = -(e+p)/t, drho_i/dt = -rho_i/t along the (ideal) trajectory.
Run from this folder:  python3 make_1n0_tables.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..")
HBARC = 0.197
DEG = 10

T = np.loadtxt(os.path.join(SRC, "T-n-1n0.dat"))
t, Tf = T[0], T[1]
e, p, mu, cB, cQ, cS = (np.loadtxt(os.path.join(SRC, "%s-n-1n0.dat" % n)) for n in ["e", "p", "mu", "cB2", "cQ2", "cS2"])
rhoB = 0.16 / t
rhoQ = 0.4 * rhoB

s = (e + p - mu * rhoB) / Tf
etas = HBARC * 2 / (4 * np.pi) * s
Mst = HBARC * (e + p) - 4. / 3. * etas / t
dMst = np.gradient(Mst, t, edge_order=2)
D = 2. / (2 * np.pi * Tf)
cs2 = -(t * np.gradient(p, t, edge_order=2) + cB * rhoB + cQ * rhoQ) / (e + p)
k = np.arange(2000)
tc = 0.5 * (t[0] + t[-1]) + 0.5 * (t[-1] - t[0]) * np.cos(np.pi * (k + 0.5) / len(k))   # Chebyshev nodes on [1, 11]
quantities = {"Mst": Mst, "dMst": dMst, "etas": etas, "dpde": cs2, "D": D,
              "dpdrhoB": HBARC * cB, "dpdrhoQ": HBARC * cQ, "dpdrhoS": HBARC * cS}

print("checks on the input tables:")
print("  t de/dt / -(e+p) in [%.4f, %.4f] (1 = ideal Bjorken)" % ((t * np.gradient(e, t, edge_order=2) / -(e + p))[2:-2].min(), (t * np.gradient(e, t, edge_order=2) / -(e + p))[2:-2].max()))
print("  s*t = %.4f ... %.4f (constant for ideal flow)" % ((s * t).min(), (s * t).max()))
for name, y in quantities.items():
    # weighted fit so that the *relative* error is uniform in t (the quantities fall by 10x over [1, 11]);
    # the weight is floored at 5% of max|y| for quantities that cross zero
    scale = np.maximum(np.abs(y), 0.05 * np.max(np.abs(y)))
    # fit on Chebyshev nodes (near-minimax; smaller end-point errors than a uniform grid), data linearly interpolated
    yc, sc = np.interp(tc, t, y), np.interp(tc, t, scale)
    c = np.polyfit(tc, yc, DEG, w=1. / sc)
    fit = np.polyval(c, t)
    print("  %-8s range [% .4g, % .4g]  max relative fit error = %.1e" % (name, y.min(), y.max(), np.max(np.abs(fit - y) / scale)))
    np.savetxt(os.path.join(HERE, "1n0-%s.dat" % name), c)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"xtick.top": True, "ytick.right": True, "xtick.direction": "in", "ytick.direction": "in",
                         "xtick.minor.visible": True, "ytick.minor.visible": True, "font.size": 9})
    fig, ax = plt.subplots(2, 4, figsize=(14, 6.5))
    labels = {"Mst": "Mst [GeV/fm$^3$]", "dMst": "dMst/dt [GeV/fm$^3$/fm]", "etas": "eta [GeV/fm$^2$]", "dpde": "dp/de = cs$^2$",
              "D": "D [fm]", "dpdrhoB": "dp/d$\\rho_B$ [GeV]", "dpdrhoQ": "dp/d$\\rho_Q$ [GeV]", "dpdrhoS": "dp/d$\\rho_S$ [GeV]"}
    for a, (name, y) in zip(ax.flat, quantities.items()):
        a.plot(t, y, "k", lw=2.4, alpha=0.5, label="from HRG tables")
        a.plot(t, np.polyval(np.loadtxt(os.path.join(HERE, "1n0-%s.dat" % name)), t), "tab:red", lw=1.2, ls="--", label="degree-%d fit (solver input)" % DEG)
        a.set_xlabel("t [fm]"); a.set_ylabel(labels[name])
    ax.flat[0].legend(frameon=False, fontsize=8)
    fig.text(0.01, 0.005, "HRG (pdg2020) trajectory: ideal Bjorken, rhoB = 0.16/t fm^-3, rhoQ = 0.4 rhoB, rhoS = 0; T(t=1) = %.0f MeV -> T(t=11) = %.0f MeV; "
             "eta/s = 2/(4 pi), D T = 2/(2 pi); 1000 points, t in [1, 11] fm" % (HBARC * 1000 * Tf[0], HBARC * 1000 * Tf[-1]), fontsize=8, family="monospace")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    os.makedirs(os.path.join(HERE, "..", "docs", "figures"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "..", "docs", "figures", "eos_1n0_tables.png"), dpi=130)
except ImportError:
    pass

# The solver takes dMst from the derivative of the Mst polynomial, not from the dMst file.
dpoly = np.polyval(np.polyder(np.loadtxt(os.path.join(HERE, "1n0-Mst.dat"))), t)
print("  d/dt of the Mst polynomial vs data derivative: max abs diff %.3f GeV/fm^3/fm (at t = 11: %.3f vs %.3f)" % (np.max(np.abs(dpoly - dMst)), dpoly[-1], dMst[-1]))
