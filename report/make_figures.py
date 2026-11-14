"""Draw the figures used in the report, sized for the page and without titles (the LaTeX captions carry those).

Run from the project root, after report.hawser_sweep:
    python3 -m report.make_figures

Long sweeps are not repeated: they are read from data/E05_limit_cycles, data/E08_control and report/data. Everything
else (eigenvalues, Hopf curves, rest positions, the wind sweeps) is recomputed with the project's own functions.
Writes PNG files (300 dpi) to report/figures/.
"""

import json
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc

from experiments.E03_pitchfork.run import draw_pitchfork
from experiments.E10_wind.run import draw_diagram, no_wind_reference, s_curve, signed_beam, sweep, windy
from fishtail import continuation as co
from fishtail import control, hydro
from fishtail import model_rigid as mr
from fishtail.params import SEA_WATER, SIMOS_FIG53, VIDAL_100
from fishtail.plotting import COLORS, INK, INK_2, MUTED, ROOT, draw_mooring, plot_branch, use_style

OUT = Path(__file__).resolve().parent / "figures"
WIDTH = 6.9  # inches: the figure width that maps to \linewidth at (almost) natural font size
DASH = (0, (4, 2.5))


def style():
    use_style()
    plt.rcParams.update({"font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5, "legend.fontsize": 7.5,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "lines.linewidth": 1.6})


def load(experiment):
    return json.loads((ROOT / "data" / experiment / "summary.json").read_text().replace("NaN", "null"))


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


def panel(ax, letter):
    ax.set_title(f"({letter})", loc="left", fontsize=9)


# --- Fig. 1: plan view -----------------------------------------------------------------------------------------------

def fig_schematic():
    s = mr.Rigid(VIDAL_100, x_A_L=0.32, L_H_L=0.62, C_F=0.03)
    gamma, psi = np.radians(32.0), np.radians(12.0)
    state = np.array([gamma, psi, 0.0, 0.0])
    fig, ax = plt.subplots(figsize=(3.45, 1.85))
    ax.set_aspect("equal")
    ax.axis("off")
    draw_mooring(ax, state, s, COLORS[0], lw=1.3)
    e_x = np.array([np.cos(psi), np.sin(psi)])
    e_y = np.array([-np.sin(psi), np.cos(psi)])
    A = s.L_H_L * np.array([-np.cos(gamma), np.sin(gamma)])
    O = A - s.x_A_L * e_x
    ax.plot([0.22, -1.42], [0, 0], color=MUTED, lw=0.8, ls=(0, (5, 3)))  # current line through the buoy
    ax.plot([A[0], A[0] + 0.32], [A[1], A[1]], color=MUTED, lw=0.8, ls=(0, (5, 3)))  # parallel to the current at A
    ax.plot(*A, "o", ms=4.5, color=COLORS[1], zorder=7)
    ax.plot(*O, "o", ms=4, color=INK, zorder=7)
    ax.plot([O[0], A[0] + 0.32 * e_x[0]], [O[1], A[1] + 0.32 * e_x[1]], color=INK_2, lw=0.8)  # ship axis
    ax.add_patch(Arc((0, 0), 0.62, 0.62, theta1=180 - np.degrees(gamma), theta2=180, color=INK, lw=0.9))
    ax.text(-0.39, 0.09, r"$\gamma$", fontsize=9, ha="center", va="center")
    ax.add_patch(Arc(A, 0.5, 0.5, theta1=0, theta2=np.degrees(psi), color=INK, lw=0.9))
    ax.text(A[0] + 0.30, A[1] + 0.03, r"$\psi$", fontsize=9, va="center")
    ax.text(*(0.5 * A + np.array([0.03, 0.05])), r"hawser, $L_H$", fontsize=7.5, ha="left", va="bottom")
    ax.text(0.0, -0.06, "buoy", fontsize=7.5, ha="center", va="top", color=INK_2)
    ax.annotate("attachment\npoint A", A, xytext=(A[0] + 0.03, A[1] + 0.17), fontsize=7.5, ha="left", color=INK_2,
                arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.6))
    ax.annotate("midship O", O, xytext=(O[0] - 0.12, O[1] + 0.33), fontsize=7.5, ha="center", color=INK_2,
                arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.6))
    p1, p2 = O + 0.145 * (-e_y), A + 0.145 * (-e_y)
    ax.annotate("", xy=p1, xytext=p2, arrowprops=dict(arrowstyle="<->", color=INK, lw=0.8, shrinkA=0, shrinkB=0))
    ax.text(*(0.5 * (p1 + p2) - 0.045 * e_y), r"$x_A$", fontsize=9, ha="center", va="top")
    for y in (-0.16, 0.42):
        ax.annotate("", xy=(0.02, y), xytext=(0.28, y), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.2))
    ax.text(0.15, 0.45, r"current $U$", fontsize=7.5, ha="center", va="bottom", color=INK_2)
    ax.set_xlim(-1.45, 0.32)
    ax.set_ylim(-0.25, 0.68)
    save(fig, "fig01_schematic")


# --- Fig. 2: pitchfork -----------------------------------------------------------------------------------------------

def fig_pitchfork():
    fig, ax = plt.subplots(figsize=(3.7, 2.1))
    draw_pitchfork(ax, VIDAL_100, COLORS[0], "rest heading", 0.5, 0.05, show_hopf=False)
    ax.plot([], [], color=INK_2, label="stable")
    ax.plot([], [], color=INK_2, ls=DASH, label="unstable")
    ax.plot([], [], "^", color=INK_2, ms=7, ls="none", label="pitchfork")
    ax.set_xlim(0, 0.5)
    ax.set_ylim(-80, 80)
    ax.set_yticks(range(-75, 76, 25))
    ax.set_xlabel(r"attachment point $x_A/L$ (0 = midship, 0.5 = bow)")
    ax.set_ylabel(r"heading at rest $\psi$ [deg]")
    ax.legend(loc="upper right", ncol=2, fontsize=6.8, columnspacing=0.8, handlelength=1.6, borderaxespad=0.1)
    save(fig, "fig02_pitchfork")


# --- Fig. 3: limit-cycle amplitude against friction --------------------------------------------------------------------

def fig_amplitude():
    data = load("E05_limit_cycles")
    C_crit = data["C_F_hopf"]
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 1.7))
    for key, marker, label, color in [("down", "v", r"$C_F$ decreasing", COLORS[0]),
                                      ("up", "^", r"$C_F$ increasing", COLORS[1])]:
        rows = data["fig01"][key]
        C = np.array([r["C_F"] for r in rows])
        kw = dict(color=color, ms=4.5, mfc="none", mew=1.0, ls="none")
        axes[0].plot(C, [r["gamma_amp_deg"] for r in rows], marker, label=label, **kw)
        axes[1].plot(C, [r["psi_amp_deg"] for r in rows], marker, **kw)
        axes[2].plot(C, np.array([r["F_H_max"] for r in rows]) / C, marker, **kw)
        axes[2].plot(C, np.array([r["F_H_min"] for r in rows]) / C, marker, color=color, ms=3.5, mew=0.6, ls="none")
    for ax in axes[:2]:
        ax.plot(C_crit, 0, "o", ms=6, color=INK, mec="white", mew=1.0, zorder=6)
        ax.set_ylim(bottom=-1)
    axes[2].axhline(0, color=INK_2, lw=0.9)
    axes[2].axhline(1, color=MUTED, lw=0.7)
    axes[0].set_ylabel(r"amplitude of $\gamma$ [deg]")
    axes[1].set_ylabel(r"amplitude of $\psi$ [deg]")
    axes[2].set_ylabel(r"max and min of $F_H/F_{H0}$")
    for ax, letter in zip(axes, "abc"):
        ax.set_xlabel(r"$C_F$")
        panel(ax, letter)
    axes[0].legend(loc="upper right")
    fig.tight_layout(pad=0.4)
    save(fig, "fig03_friction")


# --- Fig. 5: hawser length ---------------------------------------------------------------------------------------------

def fig_hawser():
    L_H = np.linspace(0.02, 2.0, 100)
    sweeps = {r["direction"]: r for r in json.loads((OUT.parent / "data" / "hawser_sweep.json").read_text())}
    L_cr, C_tab = sweeps["up"]["L_cr"], sweeps["up"]["C_F"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(WIDTH, 2.0))
    a1.plot(L_H, co.hopf_curve(mr.Rigid(VIDAL_100, 0.5, 0.5, 0.03), L_H), color=COLORS[0], lw=1.8, label="this model")
    a1.plot(SIMOS_FIG53["L_H_L"], SIMOS_FIG53["100"]["C_F"], "o", color=INK, ms=3.6, mfc="white", mew=0.9,
            label="Simos (2001), Fig. 5.3", ls="none")
    a1.axhline(C_tab, color=MUTED, lw=0.8, ls=":")
    a1.plot(L_cr, C_tab, "o", color=INK, ms=5.5, zorder=6)
    a1.set_xlim(0, 2)
    a1.set_ylim(0, 0.06)
    a1.text(1.5, 0.054, "stable: trails straight", ha="center", fontsize=8, color=INK_2)
    a1.text(1.5, 0.022, "unstable: fishtails", ha="center", fontsize=8, color=INK_2)
    a1.set_xlabel(r"hawser length $L_H/L$")
    a1.set_ylabel(r"Hopf friction $C_F$")
    a1.legend(loc="lower right")
    panel(a1, "a")

    up = sweeps["up"]["rows"]
    L = np.array([r["L_H_L"] for r in up])
    a2.plot(L, [r["gamma_amp_deg"] for r in up], color=COLORS[0], label=r"swing of $\gamma$")
    a2.plot(L, [r["psi_amp_deg"] for r in up], color=COLORS[1], label=r"swing of $\psi$")
    for key, marker, label in [("up", "^", "length increasing"), ("down", "v", "length decreasing")]:
        rows = sweeps[key]["rows"]
        for c, var in [(COLORS[0], "gamma_amp_deg"), (COLORS[1], "psi_amp_deg")]:
            a2.plot([r["L_H_L"] for r in rows], [r[var] for r in rows], marker, color=c, ms=3.2, mfc="none", mew=0.7)
        a2.plot([], [], marker, color=INK_2, ms=3.2, mfc="none", mew=0.7, ls="none", label=label)
    a2.plot(L_cr, 0, "o", color=INK, ms=5.5, mec="white", mew=0.8, zorder=6, label="Hopf point")
    a2.set_xlim(0, 1.5)
    a2.set_ylim(-1, 50)
    a2.set_xlabel(r"hawser length $L_H/L$")
    a2.set_ylabel("swing amplitude [deg]")
    a2.legend(loc="upper right", ncol=3, fontsize=6.4, columnspacing=0.7, handlelength=1.4, handletextpad=0.3,
               borderaxespad=0.1)
    panel(a2, "b")
    fig.tight_layout(pad=0.4)
    save(fig, "fig04_hawser_length")


# --- Fig. 6: wind ------------------------------------------------------------------------------------------------------

def fig_wind():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(WIDTH, 1.95))
    no_wind_reference(a1)
    draw_diagram(a1, windy(0.20, 90.0), lw=1.7)
    a1.plot([], [], color=MUTED, lw=1.0, label="no wind")
    a1.plot([], [], color=INK, label="stable")
    a1.plot([], [], color=INK, ls=DASH, label="unstable")
    a1.plot([], [], "D", color=INK, mfc="white", ms=5, ls="none", label="fold")
    a1.set_xlim(0, 0.5)
    a1.set_ylim(-60, 60)
    a1.set_xlabel(r"attachment point $x_A/L$")
    a1.set_ylabel(r"heading at rest $\psi$ [deg]")
    a1.legend(loc="upper right", ncol=2, columnspacing=0.8, handlelength=1.6, borderaxespad=0.1)
    panel(a1, "a")

    psi, h, stable = s_curve()
    keep = np.abs(h) <= 0.35
    plot_branch(a2, h[keep], np.degrees(psi[keep]), stable[keep], INK, lw=1.7)
    dh = np.diff(h)
    for k in range(dh.size - 1):
        if dh[k] * dh[k + 1] < 0:
            a2.plot(h[k + 1], np.degrees(psi[k + 1]), "D", ms=5, color=INK, mfc="white", mew=1.1, zorder=5)
    grid = np.round(np.arange(-0.30, 0.3001, 0.01), 3)
    up = sweep(grid, co.wind_equilibria(signed_beam(-0.30))[0])
    down = sweep(grid[::-1], co.wind_equilibria(signed_beam(0.30))[-1])
    a2.plot(grid, up, "o", ms=2.6, color=COLORS[0], label=r"simulation, $F_w$ rising", ls="none")
    a2.plot(grid[::-1], down, "s", ms=2.8, mfc="none", mew=0.7, color=COLORS[1], label=r"simulation, $F_w$ falling", ls="none")
    a2.axvline(0, color=MUTED, lw=0.7)
    a2.set_xlabel(r"side wind force $F_w/F_{H0}$ (> 0: from port)")
    a2.set_ylabel(r"heading at rest $\psi$ [deg]")
    a2.set_ylim(-16, 23)
    a2.legend(loc="upper left", handletextpad=0.2, ncol=2, columnspacing=0.8)
    panel(a2, "b")
    fig.tight_layout(pad=0.4)
    save(fig, "fig05_wind")


# --- Figs. 7 and 8: control --------------------------------------------------------------------------------------------

def fig_tug():
    data = load("E08_control")
    base = mr.Rigid(VIDAL_100, x_A_L=0.5, L_H_L=0.5, C_F=0.03)
    ratios = np.array(data["tug"]["ratios"])
    amps = np.array(data["tug"]["amp_deg"])
    crit = data["tug"]["hopf_ratio_at_base"]
    from experiments.E08_control.run import tug
    stable = [co.is_stable(co.eigenvalues(np.zeros(4), tug(base, float(r)))) for r in ratios]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(WIDTH, 1.85))
    plot_branch(a1, ratios, np.zeros(ratios.size), stable, INK, label="trailing straight", lw=1.6)
    has = amps > 0
    xs, ys = np.concatenate([ratios[has], [crit]]), np.concatenate([amps[has], [0.0]])
    a1.plot(xs, ys, color=COLORS[0], label=r"swing amplitude of $\gamma$")
    a1.plot(xs, -ys, color=COLORS[0])
    a1.plot(crit, 0, "o", color=INK, ms=5.5, zorder=5, label=rf"Hopf, $F_{{tug}}$ = {crit:.3f} $F_{{H0}}$")
    a1.set_xlabel(r"tug pull $F_{tug}/F_{H0}$")
    a1.set_ylabel(r"$\gamma$ [deg]")
    a1.legend(loc="upper right")
    panel(a1, "a")

    d = data["fig07"]
    L = np.array(d["L_H_L"])
    names = [("no control", "no control"), ("tug 0.5 F_H0", r"tug, 0.5 $F_{H0}$"),
             ("stern thruster, ★ gains", "thruster, ★ gains")]
    for c, (key, label) in zip(COLORS, names):
        curve = np.array([np.nan if v is None else v for v in d[key]])
        a2.plot(L, curve, color=c, label=label)
        none = np.isnan(curve)
        if none.any():
            a2.plot(L[none], np.full(none.sum(), 0.0009), "x", color=c, ms=3.5, mew=1.0)
    for U in (0.5, 1.0, 2.0):
        a2.axhline(hydro.friction_coefficient(U, VIDAL_100, SEA_WATER), color=MUTED, lw=0.7, ls=":")
    a2.text(1.5, 0.0118, "real ship, 0.5–2 m/s", fontsize=7.2, ha="right", color=INK_2)
    a2.text(0.60, 0.0032, r"×: no fishtailing at any $C_F$", fontsize=7.2, color=INK_2)
    a2.set_xlabel(r"hawser length $L_H/L$")
    a2.set_ylabel(r"Hopf friction $C_F$")
    a2.set_ylim(0, 0.068)
    a2.legend(loc="upper right", ncol=3, fontsize=6.8, columnspacing=0.8, handlelength=1.5, borderaxespad=0.1)
    panel(a2, "b")
    fig.tight_layout(pad=0.4)
    save(fig, "fig06_tug")


def fig_thruster():
    data = load("E08_control")
    base = mr.Rigid(VIDAL_100, x_A_L=0.5, L_H_L=0.5, C_F=0.03)
    star = data["star"]
    P, D = np.linspace(0, 0.04, 81), np.linspace(0, 0.1, 101)
    S = np.array([[control.decay_rate(co.eigenvalues(np.zeros(4), replace(base, k_p=float(p), k_d=float(d))))
                   for d in D] for p in P])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(WIDTH, 1.9), gridspec_kw={"width_ratios": [1.15, 1]})
    m = a1.pcolormesh(D, P, S, cmap="viridis", shading="auto", rasterized=True)
    fig.colorbar(m, ax=a1, pad=0.02, label=r"decay rate $\sigma$ [$U/L$]")
    cs = a1.contour(D, P, S, levels=[0.0, 0.05, 0.1], colors=["white", INK, INK], linewidths=[1.8, 1.3, 0.7])
    a1.clabel(cs, fmt=lambda v: "Hopf" if v == 0 else rf"$\sigma$ = {v:g}", fontsize=6.8)
    a1.plot(star["k_d"], star["k_p"], "*", color="white", mec=INK, ms=11, zorder=6)
    a1.set_xlim(0, 0.06)
    a1.text(0.049, 0.02, r"$\sigma$ capped at 0.22", ha="center", fontsize=7.5, color=INK)
    a1.set_xlabel(r"derivative gain $k_d$")
    a1.set_ylabel(r"proportional gain $k_p$")
    panel(a1, "a")

    d = data["saturation"]
    a2.semilogx(d["limits"], d["amp_from_cycle"], "o-", color=COLORS[1], ms=3, lw=1.3, label="start: full fishtail")
    a2.semilogx(d["limits"], d["amp_from_nudge"], "s--", color=COLORS[0], ms=3, lw=1.3, label=r"start: small nudge")
    bist = d["bistable_limits"]
    a2.axvspan(min(bist), max(bist), color=COLORS[3], alpha=0.15, lw=0, label="bistable")
    a2.set_xlabel(r"thrust limit / $F_{H0}$")
    a2.set_ylabel(r"final swing of $\gamma$ [deg]")
    a2.legend(loc="upper right")
    panel(a2, "b")
    fig.tight_layout(pad=0.4)
    save(fig, "fig07_thruster")


def main():
    style()
    fig_schematic()
    fig_pitchfork()
    fig_amplitude()
    fig_hawser()
    fig_wind()
    fig_tug()
    fig_thruster()


if __name__ == "__main__":
    main()
