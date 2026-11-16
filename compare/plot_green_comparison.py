"""Plots the Green functions from compare_solvers (forward Euler vs Crank-Nicolson vs RK4).

Reads compare/results/*.dat and params.txt, writes docs/figures/*.png.
Run from the repository root:  python3 compare/plot_green_comparison.py
Every figure carries the run parameters read from params.txt.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
OUT = os.path.join(HERE, "..", "docs", "figures")
os.makedirs(OUT, exist_ok=True)

# figure style: ticks on all four sides, pointing inward
plt.rcParams.update({"xtick.top": True, "ytick.right": True, "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.minor.visible": True, "ytick.minor.visible": True, "font.size": 9})

P = dict(l.strip().split("=", 1) for l in open(os.path.join(RES, "params.txt")) if "=" in l)
Nl, eta_max = int(P["Nl"]), float(P["eta_max"])
eta = np.linspace(-eta_max, eta_max, Nl)
SRC = ["B", "Q", "S", "E", "U"]
SRC_LABEL = [r"$\delta B$", r"$\delta Q$", r"$\delta S$", r"$\delta e$", r"$\delta u$"]
RESP_LABEL = [r"$B$", r"$Q$", r"$S$", r"$e$", r"$u$"]
SOLVERS = [("euler", "forward Euler  (dt = %s)" % P["dt_euler"], dict(color="k", lw=2.6, alpha=0.55)),
           ("cn", "Crank-Nicolson  (dt = %s)" % P["dt_cn"], dict(color="tab:red", lw=1.4, ls="--")),
           ("rk4", "RK4  (dt = %s)" % P["dt_rk4"], dict(color="tab:blue", lw=1.4, ls=":"))]
if P["table_prefix"] == "test":
    TABLES = "tables: tests/test-*.dat  (test set: 0n0 EoS + synthetic constant dp/drho = 0.05, 0.02, 0.03; not physical)"
else:
    TABLES = ("tables: %s-*.dat  (HRG EoS, ideal trajectory rhoB = 0.16/t, rhoQ = 0.4 rhoB, rhoS = 0; "
              "eta/s = 2/(4 pi), D T = 2/(2 pi))" % P["table_prefix"])


def param_text(t1):
    return ("Linearised Bjorken hydro, periodic eta grid: Nl = %s points on [-%g, %g], deta = %.4f;  delta source: sigma^2 = %s, PI = %s;\n"
            "evolution t = %g -> %g;  rhoB0 = 0.16, rhoQ0 = 0.4 rhoB0, t0 = 1 (bj_flow.h);  %s;\n"
            "run time over all 10 sources: Euler %.1f s, CN %.1f s, RK4 %.1f s"
            % (P["Nl"], eta_max, eta_max, 2 * eta_max / (Nl - 1), P["source_sigma2"], P["source_PI"], t1, float(P["t_end"]), TABLES,
               float(P["seconds_euler"]), float(P["seconds_cn"]), float(P["seconds_rk4"])))


def load(solver, t1, src):
    return np.loadtxt(os.path.join(RES, "green_%s_t1_%g_src%s.dat" % (solver, t1, src)))   # (5 responses, Nl)


t1_values = [float(x) for x in P["t1_values"].split(",")]
data = {(s, t1, a): load(s, t1, a) for s, _, _ in SOLVERS for t1 in t1_values for a in SRC}

# ---------------------------------------------------------------- figure 1: Green functions, 5 x 5
for t1 in t1_values:
    fig, ax = plt.subplots(5, 5, figsize=(14, 11.5), sharex=True)
    for i, a in enumerate(SRC):
        for j in range(5):
            axx = ax[i, j]
            for s, lab, kw in SOLVERS:
                axx.plot(eta, data[(s, t1, a)][j], label=lab if (i, j) == (0, 0) else None, **kw)
            axx.set_xlim(-eta_max, eta_max)
            axx.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))
            if np.max(np.abs(data[("euler", t1, a)][j])) == 0:
                axx.text(0.5, 0.5, "identically 0\n(decoupled)", transform=axx.transAxes, ha="center", va="center", fontsize=9)
                axx.set_yticks([])
            if i == 0: axx.set_title("response " + RESP_LABEL[j])
            if j == 0: axx.set_ylabel("source " + SRC_LABEL[i] + "\n" + r"$G$")
            if i == 4: axx.set_xlabel(r"$\eta$")
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.965), fontsize=10)
    fig.suptitle("Green functions at t = %g after a delta source at t1 = %g" % (float(P["t_end"]), t1), y=0.995, fontsize=12)
    fig.text(0.01, 0.003, param_text(t1), fontsize=7.5, va="bottom", family="monospace")
    fig.tight_layout(rect=(0, 0.045, 1, 0.965))
    fig.savefig(os.path.join(OUT, "green_comparison_t1_%g.png" % t1), dpi=130)
    plt.close(fig)

# ---------------------------------------------------------------- figure 2: how large are the differences?
pairs = [("cn", "euler", "CN - Euler"), ("rk4", "euler", "RK4 - Euler"), ("cn", "rk4", "CN - RK4")]
fig, ax = plt.subplots(len(t1_values), 3, figsize=(13, 4.6 * len(t1_values)), squeeze=False)
for r, t1 in enumerate(t1_values):
    for c, (s1, s2, title) in enumerate(pairs):
        M = np.full((5, 5), np.nan)
        for i, a in enumerate(SRC):
            for j in range(5):
                ref = data[("euler", t1, a)][j]
                scale = np.max(np.abs(ref))
                if scale > 0:
                    M[i, j] = max(np.max(np.abs(data[(s1, t1, a)][j] - data[(s2, t1, a)][j])) / scale, 1e-16)
        im = ax[r, c].imshow(np.log10(M), cmap="viridis", vmin=-8, vmax=-2)
        for i in range(5):
            for j in range(5):
                if np.isfinite(M[i, j]):
                    ax[r, c].text(j, i, "%.0e" % M[i, j], ha="center", va="center", color="w", fontsize=8)
                else:
                    ax[r, c].text(j, i, "G = 0", ha="center", va="center", color="gray", fontsize=8)
        ax[r, c].set_xticks(range(5)); ax[r, c].set_xticklabels(RESP_LABEL)
        ax[r, c].set_yticks(range(5)); ax[r, c].set_yticklabels(SRC_LABEL)
        ax[r, c].tick_params(which="both", top=True, right=True, direction="in")
        ax[r, c].minorticks_off()
        ax[r, c].set_title("%s   (t1 = %g)" % (title, t1))
        ax[r, c].set_xlabel("response"); ax[r, c].set_ylabel("source")
cax = fig.add_axes([0.90, 0.2, 0.015, 0.6])
fig.colorbar(im, cax=cax, label=r"$\log_{10}\,[\max|\Delta G| / \max|G_{\rm Euler}|]$")
cax.tick_params(which="both", direction="in", right=True)
fig.suptitle("Relative differences between solvers (per source/response pair)", fontsize=12)
fig.text(0.01, 0.003, param_text(t1_values[0]).replace("t = %g -> %g" % (t1_values[0], float(P["t_end"])), "t = t1 -> %g" % float(P["t_end"])),
         fontsize=7.5, va="bottom", family="monospace")
fig.subplots_adjust(left=0.06, right=0.88, bottom=0.12, top=0.9, wspace=0.3, hspace=0.35)
fig.savefig(os.path.join(OUT, "green_comparison_differences.png"), dpi=130)
print("wrote", sorted(f for f in os.listdir(OUT) if f.startswith("green_comparison")))
