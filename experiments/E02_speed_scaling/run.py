"""E02 - Does the current speed change what the moored ship does?

Run from the project root:
    /home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E02_speed_scaling.run

Writes figures to figures/E02_speed_scaling/ and numbers to data/E02_speed_scaling/summary.json.
Sign convention (Simos 2001): psi > 0 = bow turned to port, gamma > 0 = ship swung to port. t* = U t / L.

Idea (notes/01_why_current_speed_matters.md): every force grows like U^2, so a faster current replays the same motion
faster - unless the friction coefficient C_F changes, and it does change with the Reynolds number U L / nu.
"""

import json
from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from fishtail import continuation as co
from fishtail import hydro
from fishtail import model_rigid as mr
from fishtail.params import FRESH_WATER, SEA_WATER, VIDAL_100, froude_scaled
from fishtail.plotting import COLORS, INK, INK_2, MUTED, ROOT, save, use_style

EXP = "E02_speed_scaling"
START = [np.radians(5.0), 0.0, 0.0, 0.0]  # released with the hawser 5 deg to port, ship at rest


def fig01_fixed_friction(summary):
    """Same friction coefficient, three current speeds: different in seconds, identical in t* = U t / L."""
    base = mr.Rigid(VIDAL_100, x_A_L=0.5, L_H_L=0.5, C_F=0.035)  # unstable set-up, so the swing grows
    speeds = [0.5, 1.0, 2.0]  # full-scale current speeds [m/s]
    t_star_end = 120.0
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.5, 7.0))
    curves = []
    for color, U in zip(COLORS, speeds):
        s = replace(base, U=U)
        t_end = t_star_end * s.time_unit  # [s]
        t = np.linspace(0, t_end, 3001)
        sol = solve_ivp(mr.rhs, (0, t_end), START, args=(s,), t_eval=t, rtol=1e-9, atol=1e-12)  # in SECONDS
        psi_deg = np.degrees(sol.y[1])
        ax1.plot(t / 3600, psi_deg, color=color, label=f"U = {U} m/s")
        ax2.plot(t / s.time_unit, psi_deg, color=color, lw=[4.5, 2.6, 1.2][speeds.index(U)],
                 label=f"U = {U} m/s")
        curves.append(psi_deg)
    summary["fixed_C_F_max_difference_deg"] = float(max(np.abs(c - curves[0]).max() for c in curves))
    ax1.set_xlabel("time t [hours], full scale")
    ax1.set_ylabel("heading ψ [deg]")
    ax1.set_title("(a) In real time: a faster current makes the same swing happen faster", loc="left")
    ax2.set_xlabel("non-dimensional time  t* = U t / L")
    ax2.set_ylabel("heading ψ [deg]")
    ax2.set_title("(b) In t* = U t / L the three runs are the same curve (drawn thick to thin)", loc="left")
    for ax in (ax1, ax2):
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.legend(loc="upper left")
    fig.suptitle(f"Friction coefficient held fixed ($C_F$ = {base.C_F}): Vidal de Negreiros 100 %, "
                 f"$x_A/L$ = {base.x_A_L}, $L_H/L$ = {base.L_H_L}, released with γ = 5°",
                 x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig01_fixed_friction")


def fig02_friction_from_speed(summary):
    """Friction from the ITTC line (depends on speed and scale): the tank model at 0.2/0.4/0.6 m/s and the real ship.

    Same set-up for all four runs; only C_F differs. Hawser length chosen so that the stability boundary lies inside
    the range of friction values the tests span (L_H/L = 0.15).
    """
    rho_ratio = FRESH_WATER.rho / SEA_WATER.rho
    model = froude_scaled(VIDAL_100, rho_ratio)
    cases = [(f"1:90 model, {U} m/s", froude_scaled(VIDAL_100, rho_ratio), FRESH_WATER, U) for U in (0.2, 0.4, 0.6)]
    cases.append(("full-scale ship, 1.90 m/s", VIDAL_100, SEA_WATER, 1.90))
    x_A_L, L_H_L = 0.5, 0.15
    s_ref = mr.Rigid(VIDAL_100, x_A_L, L_H_L, C_F=0.02)
    C_crit = co.hopf_C_F(s_ref)
    summary["fig02"] = {"L_H_L": L_H_L, "x_A_L": x_A_L, "C_F_hopf": C_crit, "runs": {}}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.6), gridspec_kw={"width_ratios": [1.6, 1]})
    for color, (label, hull, fluid, U) in zip(COLORS, cases):
        C_F = float(hydro.friction_coefficient(U, hull, fluid))
        s = mr.Rigid(VIDAL_100, x_A_L, L_H_L, C_F)  # in t* only C_F matters, so the full-scale hull is used
        t, y = mr.simulate(s, START, t_end_star=150.0, n_out=3001)
        ax1.plot(t, np.degrees(y[1]), color=color, label=f"{label}: $C_F$ = {C_F:.4f}")
        ax2.plot(L_H_L, C_F, "o", color=color, ms=8, mec="white", mew=1.5, zorder=5)
        summary["fig02"]["runs"][label] = {"C_F": C_F, "stable": bool(C_F > C_crit),
                                           "psi_end_deg": float(np.degrees(np.abs(y[1, -600:]).max()))}
    ax1.axhline(0, color=MUTED, lw=0.8)
    ax1.set_xlabel("non-dimensional time  t* = U t / L")
    ax1.set_ylabel("heading ψ [deg]")
    ax1.set_title("(a) Same ship and hawser; only the friction differs", loc="left")
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=2, fontsize=8.5)  # below, clear of curves

    L_H = np.linspace(0.02, 0.6, 59)
    ax2.plot(L_H, co.hopf_curve(s_ref, L_H), color=INK, lw=2.0)
    ax2.axvline(L_H_L, color=MUTED, lw=0.8)
    ax2.text(0.35, 0.040, "stable\n(trails straight)", color=INK_2, fontsize=9, ha="center")
    ax2.text(0.42, 0.012, "unstable\n(fishtails)", color=INK_2, fontsize=9, ha="center")
    ax2.set_xlim(0, 0.6)
    ax2.set_ylim(0, 0.05)
    ax2.set_xlabel("hawser length  $L_H/L$")
    ax2.set_ylabel("friction coefficient $C_F$")
    ax2.set_title("(b) Where the four runs sit on the\nstability boundary (Hopf curve, E06)", loc="left")
    fig.suptitle(f"Friction from the ITTC line, eq. 3.20: model tests have more friction, so they look more stable "
                 f"(100 %, $x_A/L$ = {x_A_L}, $L_H/L$ = {L_H_L})", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig02_friction_from_speed")


def fig03_growth_vs_speed(summary):
    """How fast a small swing grows or dies, against full-scale current speed, for a few hawser lengths."""
    U_full = np.linspace(0.25, 3.0, 56)
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    summary["fig03"] = {}
    for color, L_H_L in zip(COLORS, [0.05, 0.1, 0.2, 0.5]):
        growth = []
        for U in U_full:
            C_F = float(hydro.friction_coefficient(U, VIDAL_100, SEA_WATER))
            growth.append(co.max_oscillatory_growth(mr.Rigid(VIDAL_100, 0.5, L_H_L, C_F)))
        ax.plot(U_full, growth, color=color, label=f"$L_H/L$ = {L_H_L}")
        summary["fig03"][str(L_H_L)] = {"U": list(U_full), "growth": [float(g) for g in growth]}
    ax.axhline(0, color=INK_2, lw=1.0)
    ax.text(2.95, 0.003, "grows: fishtails", ha="right", va="bottom", color=INK_2, fontsize=9)
    ax.text(2.95, -0.003, "dies out: trails straight", ha="right", va="top", color=INK_2, fontsize=9)
    ax.set_xlabel("current speed U [m/s], full-scale ship  (enters only through $C_F$(Re))")
    ax.set_ylabel("growth rate of the swing  Re λ  [units of U/L]")
    ax.set_title("Full-scale Vidal de Negreiros 100 %, $x_A/L$ = 0.5: faster current → less friction → "
                 "less stable", loc="left")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))  # outside, so no curve is hidden
    fig.tight_layout()
    save(fig, EXP, "fig03_growth_vs_speed")


def main():
    use_style()
    summary = {}
    fig01_fixed_friction(summary)
    fig02_friction_from_speed(summary)
    fig03_growth_vs_speed(summary)
    folder = ROOT / "data" / EXP
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps({k: v for k, v in summary.items() if k != "fig03"}, indent=2))


if __name__ == "__main__":
    main()
