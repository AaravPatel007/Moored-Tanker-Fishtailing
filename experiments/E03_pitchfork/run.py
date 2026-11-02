"""E03 - When does the ship stop trailing straight and sit at an angle? (pitchfork bifurcation)

Run from the project root:
    /home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E03_pitchfork.run

Knob: attachment point x_A/L (forward of midship; 0.5 = at the bow). Rigid hawser, Vidal de Negreiros 100 % and 40 %.
Writes figures to figures/E03_pitchfork/ and numbers to data/E03_pitchfork/summary.json.
Sign convention (Simos 2001): psi > 0 = bow turned to port, gamma > 0 = ship swung to port. Time t* = U t / L.
"""

import json
from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from fishtail import continuation as co
from fishtail import hydro
from fishtail import model_rigid as mr
from fishtail.params import TABLE_5_2, TABLE_5_3_KW, VIDAL_100, VIDAL_40
from fishtail.plotting import (COLORS, INK, INK_2, MUTED, ROOT, current_arrow, draw_mooring, plot_branch, save,
                               use_style)

EXP = "E03_pitchfork"
LOADS = [(VIDAL_100, "100 % load"), (VIDAL_40, "40 % load (ballast)")]
L_H_L = 0.5  # hawser length used unless stated
C_HIGH = 0.05  # high friction (like a small tank model): the pitchfork is "clean", nothing fishtails
C_LOW = 0.015  # low friction (full-scale ship, about 1-2 m/s; E01 fig08)
X_GRID = np.linspace(0.0, 0.6, 241)


def branches(hull, L_H_L, C_F, x_grid=X_GRID):
    """Follow both equilibrium branches along x_A/L. Returns a dict of arrays:
    x, psi (angled heading, NaN where it does not exist), stable_straight, stable_angled, gamma (angled), F_H (angled),
    and the leading eigenvalues on each branch."""
    out = {k: [] for k in ["x", "psi", "gamma", "F_H", "stable_straight", "stable_angled", "eig_straight",
                           "eig_angled"]}
    for xa in x_grid:
        s = mr.Rigid(hull, xa, L_H_L, C_F)
        eig0 = co.eigenvalues(np.zeros(4), s)
        heads = co.equilibrium_headings(hull, xa)
        out["x"].append(xa)
        out["stable_straight"].append(co.is_stable(eig0))
        out["eig_straight"].append(eig0)
        if heads:
            state = co.equilibrium_state(heads[0], s)  # one angled pair in this range (checked: never more)
            eig = co.eigenvalues(state, s)
            out["psi"].append(state[1])
            out["gamma"].append(state[0])
            out["F_H"].append(mr.hawser_force_star(state, s))
            out["stable_angled"].append(co.is_stable(eig))
            out["eig_angled"].append(eig)
        else:
            for k in ["psi", "gamma", "F_H"]:
                out[k].append(np.nan)
            out["stable_angled"].append(False)
            out["eig_angled"].append(np.full(4, np.nan + 0j))
    return {k: np.array(v) for k, v in out.items()}


def angled_state(hull, xa, L_H_L, C_F):
    s = mr.Rigid(hull, xa, L_H_L, C_F)
    return co.equilibrium_state(co.equilibrium_headings(hull, xa)[0], s), s


def hopf_points(hull, L_H_L, C_F):
    """x_A/L values where a complex eigenvalue pair crosses the imaginary axis, on each branch (Hopf points)."""
    x_cr = hydro.critical_turret_position(hull)
    straight = co.find_crossings(lambda xa: co.max_oscillatory_growth(mr.Rigid(hull, xa, L_H_L, C_F)),
                                 np.linspace(x_cr + 1e-3, 0.6, 60))
    angled = co.find_crossings(lambda xa: co.oscillatory_growth(*angled_state(hull, xa, L_H_L, C_F)),
                               np.linspace(1e-3, x_cr - 1e-3, 60))
    return straight, angled


def draw_pitchfork(ax, hull, color, label, L_H_L, C_F, summary_key=None, summary=None, show_hopf=True):
    b = branches(hull, L_H_L, C_F)
    x_cr = hydro.critical_turret_position(hull)
    plot_branch(ax, b["x"], np.zeros_like(b["x"]), b["stable_straight"], color)
    ax.plot([], [], color=color, lw=2.0, label=label)  # legend entry: colour = loading condition
    ok = ~np.isnan(b["psi"])
    for sign in (1, -1):
        plot_branch(ax, b["x"][ok], sign * np.degrees(b["psi"][ok]), b["stable_angled"][ok], color)
    ax.plot(x_cr, 0, "^", ms=10, color=color, mec="white", mew=1.4, zorder=6)
    hopf = ([], [])
    if show_hopf:
        hopf = hopf_points(hull, L_H_L, C_F)
        for xh in hopf[0]:
            ax.plot(xh, 0, "o", ms=8, color=color, mec="white", mew=1.4, zorder=6)
        for xh in hopf[1]:
            psi_h = np.degrees(co.equilibrium_headings(hull, xh)[0])
            ax.plot([xh, xh], [psi_h, -psi_h], "o", ms=8, color=color, mec="white", mew=1.4, zorder=6)
    if summary is not None:
        summary[summary_key] = {"x_Tcr": x_cr, "hopf_straight_x_A": hopf[0], "hopf_angled_x_A": hopf[1]}
    return b


def legend_keys(ax):
    ax.plot([], [], color=INK_2, ls="-", label="stable")
    ax.plot([], [], color=INK_2, ls=(0, (4, 2.5)), label="unstable")
    ax.plot([], [], "^", color=INK_2, ms=9, ls="none", label="pitchfork ($x_A = x_{Tcr}$)")


def fig01_pitchfork(summary):
    """The pitchfork diagram at high friction: heading at rest against attachment point, both loads."""
    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    for color, (hull, label) in zip(COLORS, LOADS):
        draw_pitchfork(ax, hull, color, label, L_H_L, C_HIGH, f"fig01 {label}", summary, show_hopf=True)
    legend_keys(ax)
    ax.plot([], [], "o", color=INK_2, ms=8, ls="none", label="Hopf point")
    ax.axvline(0.5, color=MUTED, lw=0.8)
    ax.text(0.505, -80, "bow", color=INK_2, fontsize=9)
    ax.set_xlim(0, 0.6)
    ax.set_ylim(-85, 85)
    ax.set_yticks(range(-75, 76, 25))
    ax.set_xlabel("attachment point  $x_A/L$  (forward of midship; 0.5 = bow)")
    ax.set_ylabel("heading at rest  ψ [deg]   (> 0: bow to port)")
    ax.set_title(f"Pitchfork: attach the hawser too far aft and the ship sits at an angle\n"
                 f"rigid hawser $L_H/L$ = {L_H_L}, high friction $C_F$ = {C_HIGH}", loc="left")
    ax.legend(loc="upper right", fontsize=8.5, ncol=2)
    fig.tight_layout()
    save(fig, EXP, "fig01_pitchfork")


def fig02_eigenvalues(summary):
    """Real parts of the eigenvalues along both branches (100 % load): one real eigenvalue passes through zero."""
    hull = VIDAL_100
    b = branches(hull, L_H_L, C_HIGH)
    x_cr = hydro.critical_turret_position(hull)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6), sharey=True)
    for ax, key, title in [(ax1, "eig_straight", "(a) Trailing straight (ψ = 0)"),
                           (ax2, "eig_angled", "(b) Sitting at an angle (ψ ≠ 0)")]:
        eig = b[key]
        real_only = np.abs(eig.imag) < 1e-7
        for k in range(4):
            ax.plot(b["x"][real_only[:, k]], eig[real_only[:, k], k].real, ".", ms=3.5, color=COLORS[0])
            ax.plot(b["x"][~real_only[:, k]], eig[~real_only[:, k], k].real, ".", ms=3.5, color=COLORS[1])
        ax.axhline(0, color=INK_2, lw=1.0)
        ax.axvline(x_cr, color=MUTED, lw=0.8)
        ax.plot(x_cr, 0, "^", ms=10, color=INK, mec="white", mew=1.4, zorder=6)
        ax.set_xlabel("attachment point  $x_A/L$")
        ax.set_title(title, loc="left")
        ax.set_xlim(0, 0.6)
    ax1.plot([], [], "o", color=COLORS[0], label="real eigenvalue")
    ax1.plot([], [], "o", color=COLORS[1], label="complex pair (real part)")
    ax1.legend(loc="lower left")
    ax1.set_ylabel("real part of eigenvalue  [units of U/L]")
    ax1.set_ylim(-2.8, 0.6)
    ax1.text(x_cr + 0.015, 0.25, "one real eigenvalue\ncrosses zero at $x_{Tcr}$", ha="left", color=INK_2,
             fontsize=9)
    fig.suptitle(f"Eigenvalues along the two branches, 100 % load, $L_H/L$ = {L_H_L}, $C_F$ = {C_HIGH}: "
                 f"> 0 means a small disturbance grows", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig02_eigenvalues")


def fig03_low_friction(summary):
    """The same diagram at full-scale friction: Hopf points appear on both branches (fishtailing)."""
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0), sharey=True)
    for ax, color, (hull, label) in zip(axes, COLORS, LOADS):
        draw_pitchfork(ax, hull, color, label, L_H_L, C_LOW, f"fig03 {label}", summary, show_hopf=True)
        x_cr = hydro.critical_turret_position(hull)
        hs = summary[f"fig03 {label}"]["hopf_straight_x_A"]
        # Past x_Tcr the straight position oscillates unstably from the first Hopf point on (or from x_Tcr itself
        # if its complex pair is already unstable there).
        x_fish = hs[0] if hs else x_cr
        ax.axvspan(x_fish, 0.6, color=COLORS[1], alpha=0.07, lw=0)
        ax.text(0.5 * (x_fish + 0.6), 60, "straight position\nfishtails", ha="center", color=INK_2, fontsize=9)
        ax.set_title(f"({'a' if hull is VIDAL_100 else 'b'}) {label}: $x_{{Tcr}}/L$ = {x_cr:.3f}", loc="left")
        ax.set_xlim(0, 0.6)
        ax.set_xlabel("attachment point  $x_A/L$")
    legend_keys(axes[0])
    axes[0].plot([], [], "o", color=INK_2, ms=8, ls="none", label="Hopf point")
    axes[0].legend(loc="lower right", fontsize=8.5)
    axes[0].set_ylabel("heading at rest  ψ [deg]")
    axes[0].set_ylim(-85, 85)
    fig.suptitle(f"Low friction ($C_F$ = {C_LOW}, real ship): the same pitchfork, but now the branches lose "
                 f"stability to oscillations (Hopf ●);  $L_H/L$ = {L_H_L}", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig03_low_friction")


def fig04_hawser(summary):
    """Hawser angle and force at the angled equilibrium: the heading does not depend on friction, these do."""
    hull = VIDAL_100
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(13.0, 4.4))
    x_cr = hydro.critical_turret_position(hull)
    for color, C_F in zip(COLORS, [0.015, 0.03, 0.05]):
        b = branches(hull, L_H_L, C_F, x_grid=np.linspace(0, x_cr, 121))
        ax1.plot(b["x"], np.degrees(b["psi"]), color=color, lw=[4.0, 2.4, 1.2][[0.015, 0.03, 0.05].index(C_F)],
                 label=f"$C_F$ = {C_F}")
        ax2.plot(b["x"], np.degrees(b["gamma"]), color=color, label=f"$C_F$ = {C_F}")
        ax3.plot(b["x"], b["F_H"] / C_F, color=color, label=f"$C_F$ = {C_F}")
    ax1.set_ylabel("heading ψ [deg]")
    ax1.set_title("(a) Heading: the same for every $C_F$", loc="left")
    ax2.set_ylabel("hawser angle γ [deg]")
    ax2.set_title("(b) Hawser angle: larger at low friction", loc="left")
    ax3.set_ylabel("hawser force ÷ straight-line value  $F_H / F_{H0}$")
    ax3.set_title("(c) Hawser force rises when the ship\nsits at an angle", loc="left")
    for ax in (ax1, ax2, ax3):
        ax.set_xlabel("attachment point  $x_A/L$")
        ax.axvline(x_cr, color=MUTED, lw=0.8)
        ax.set_xlim(0, 0.4)
        ax.legend(loc="upper right", fontsize=8.5)
    ax3.set_yscale("log")
    fig.suptitle(f"Angled equilibrium (port branch), 100 % load, $L_H/L$ = {L_H_L}.  "
                 f"$F_{{H0}}$ = $C_F\\,\\frac{{1}}{{2}}\\rho U^2 L T$ = tension when trailing straight",
                 x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig04_hawser_angle_and_force")


def fig05_time_runs(summary):
    """Time runs either side of x_Tcr: the ship settles at an angle (aft attachment) or straight (forward)."""
    hull = VIDAL_100
    x_cr = hydro.critical_turret_position(hull)
    cases = [(0.25, COLORS[0]), (0.45, COLORS[2])]
    start = [np.radians(5.0), np.radians(2.0), 0.0, 0.0]
    fig = plt.figure(figsize=(12.5, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15])
    ax_psi, ax_gam = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    summary["fig05"] = {}
    for col, (xa, color) in enumerate(cases):
        s = mr.Rigid(hull, xa, L_H_L, C_HIGH)
        t, y = mr.simulate(s, start, t_end_star=120.0, n_out=2401)
        ax_psi.plot(t, np.degrees(y[1]), color=color, label=f"$x_A/L$ = {xa}")
        ax_gam.plot(t, np.degrees(y[0]), color=color, label=f"$x_A/L$ = {xa}")
        summary["fig05"][str(xa)] = {"psi_end_deg": float(np.degrees(y[1, -1])),
                                     "gamma_end_deg": float(np.degrees(y[0, -1]))}
        heads = co.equilibrium_headings(hull, xa)
        if heads:
            ax_psi.axhline(np.degrees(heads[0]), color=color, lw=0.8, ls=":")
        ax = fig.add_subplot(gs[1, col])
        for k, i in enumerate([0, 150, 400, 2400]):
            draw_mooring(ax, y[:, i], s, color, alpha=0.25 + 0.75 * k / 3, fill=(i == 2400))
        mid = np.array([mr.positions(y[:, i], s)["O"] for i in range(0, t.size, 4)])
        ax.plot(mid[:, 0], mid[:, 1], color=color, lw=0.8, alpha=0.6)
        current_arrow(ax, 0.25, 0.45)
        ax.set_aspect("equal")
        ax.set_xlim(-1.6, 0.4)
        ax.set_ylim(-0.6, 0.65)
        ax.set_xlabel("X / L  (upstream →)")
        ax.set_ylabel("Y / L  (port ↑)")
        verdict = "sits at an angle" if xa < x_cr else "returns to trailing straight"
        ax.set_title(f"({'c' if col == 0 else 'd'}) Top view, $x_A/L$ = {xa}: {verdict}\n"
                     f"faint = start, dark = t* = 120; thin line = path of midship", loc="left")
    ax_psi.set_ylabel("heading ψ [deg]")
    ax_psi.set_title(f"(a) Heading  (dotted: predicted angled heading)", loc="left")
    ax_gam.set_ylabel("hawser angle γ [deg]")
    ax_gam.set_title("(b) Hawser angle", loc="left")
    for ax in (ax_psi, ax_gam):
        ax.set_xlabel("time  t* = U t / L")
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.legend(loc="upper right")
    fig.suptitle(f"Released from γ = 5°, ψ = 2° either side of $x_{{Tcr}}/L$ = {x_cr:.3f}  (100 % load, "
                 f"$L_H/L$ = {L_H_L}, $C_F$ = {C_HIGH})", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig05_time_runs")


def normal_form_prediction(hull, x):
    """Small-angle prediction of the angled heading [rad] (our expansion of eqs. 3.21 and 3.29 about psi = 0):

        moment about attachment  ~  c2 (x_Tcr - x) psi  -  q psi|psi|,
        c2 = (pi T/L)(1 + 0.4 C_B B/T),   q = l_p C_Y + x (C_Y - pi T / 2L)
    so psi_e = c2 (x_Tcr - x) / q : a STRAIGHT line in (x_Tcr - x), not a square root.
    """
    TL = hull.T / hull.L
    c2 = np.pi * TL * (1 + 0.4 * hull.C_B * hull.B / hull.T)
    q = hull.lp_L * hull.C_Y + x * (hull.C_Y - np.pi * TL / 2)
    return c2 * (hydro.critical_turret_position(hull) - x) / q


def fig06_normal_form(summary):
    """How fast the angle grows past the pitchfork: proportional to the distance, not its square root."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.6))
    summary["fig06_slopes"] = {}
    for color, (hull, label) in zip(COLORS, LOADS):
        x_cr = hydro.critical_turret_position(hull)
        dist = np.geomspace(1e-4, 0.15, 40)
        psi = np.array([co.equilibrium_headings(hull, x_cr - d)[0] for d in dist])
        ax1.loglog(dist, np.degrees(psi), "o", ms=4.5, color=color, label=f"{label}: computed")
        ax1.loglog(dist, np.degrees(normal_form_prediction(hull, x_cr - dist)), color=color, lw=1.2,
                   label=f"{label}: formula ψ = $c_2 Δ/q$")
        slope = np.polyfit(np.log(dist[:15]), np.log(psi[:15]), 1)[0]
        summary["fig06_slopes"][label] = slope
        x = np.linspace(0, x_cr, 200)
        ax2.plot(x, np.degrees([co.equilibrium_headings(hull, xx)[0] if xx < x_cr else 0 for xx in x]),
                 color=color, label=f"{label}: computed")
        ax2.plot(x, np.degrees(normal_form_prediction(hull, x)), color=color, lw=1.0, ls=":",
                 label=f"{label}: straight-line formula")
    d = np.geomspace(1e-4, 0.15, 10)
    ax1.loglog(d, 30 * np.sqrt(d), color=MUTED, lw=1.0, ls="--", label="slope ½ (square-root law), for comparison")
    ax1.set_xlabel("distance past the pitchfork  Δ = $(x_{Tcr} - x_A)/L$")
    ax1.set_ylabel("heading at rest ψ [deg]")
    ax1.set_title("(a) Log-log: slope 1, so ψ ∝ Δ (not √Δ)", loc="left")
    ax1.legend(loc="lower right", fontsize=8)
    ax2.set_xlabel("attachment point  $x_A/L$")
    ax2.set_ylabel("heading at rest ψ [deg]")
    ax2.set_title("(b) The straight-line formula is good up to about 20° of heading", loc="left")
    ax2.set_ylim(0, 60)
    ax2.legend(loc="upper right", fontsize=8)
    fig.suptitle("Why the branches leave in a straight line: the sideways drag grows like ψ|ψ| (quadratic), "
                 "so the normal form is ẋ = rx − x|x|, not Strogatz's ẋ = rx − x³", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig06_normal_form")


def fig07_simos_static_tests(summary):
    """Compare with Simos's static-bifurcation tests (Table 5-2) and his Froude correction k_w (Table 5-3)."""
    hull = VIDAL_40
    x = np.linspace(0.0, 0.45, 301)
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    summary["fig07"] = {}
    for color, U, kw, psi_meas in zip(COLORS, TABLE_5_3_KW["U_model"], TABLE_5_3_KW["k_w"], TABLE_5_2["psi_deg"]):
        # Simos eq. 5.4: C6C -> C6C (1 + k_w). Equilibrium: (1 + k_w) C6C - (x_A/L) C2C = 0.
        f = lambda psi, xa: (1 + kw) * hydro.C6C(np.pi - psi, hull) - xa * hydro.C2C(np.pi - psi, hull)
        psi_e = []
        for xa in x:
            grid = np.linspace(1e-6, np.pi / 2, 2001)
            vals = f(grid, xa)
            idx = np.where(vals[:-1] * vals[1:] < 0)[0]
            psi_e.append(np.degrees(brentq(f, grid[idx[0]], grid[idx[0] + 1], args=(xa,))) if idx.size else 0.0)
        psi_e = np.array(psi_e)
        ax.plot(x, psi_e, color=color, label=f"model with $k_w$ = {kw}  ({U} m/s)")
        ax.plot(TABLE_5_2["x_A_L"], psi_meas, "o", color=color, ms=9, mec="white", mew=1.5, zorder=6)
        summary["fig07"][str(U)] = {"k_w": kw, "psi_model_at_0.20": float(np.interp(0.20, x, psi_e)),
                                    "psi_measured": psi_meas}
    ax.plot([], [], "o", color=INK_2, ms=9, label="measured, Simos Table 5-2 (1:90 model)")
    ax.axvline(TABLE_5_2["x_A_L"], color=MUTED, lw=0.8)
    ax.set_xlabel("attachment point  $x_A/L$")
    ax.set_ylabel("heading at rest ψ [deg]")
    ax.set_xlim(0, 0.45)
    ax.set_ylim(0, 60)
    ax.set_title("40 % load, static-bifurcation tests at $x_A/L$ = 0.20 ($L_H/L$ = 0.287).\n"
                 "Our model alone ($k_w$ = 0) matches only the slowest test; Simos's wave correction $k_w$ "
                 "(fitted to these tests) moves the pitchfork aft", loc="left", fontsize=10)
    ax.legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    save(fig, EXP, "fig07_simos_static_tests")


def main():
    use_style()
    summary = {}
    fig01_pitchfork(summary)
    fig02_eigenvalues(summary)
    fig03_low_friction(summary)
    fig04_hawser(summary)
    fig05_time_runs(summary)
    fig06_normal_form(summary)
    fig07_simos_static_tests(summary)
    folder = ROOT / "data" / EXP
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
