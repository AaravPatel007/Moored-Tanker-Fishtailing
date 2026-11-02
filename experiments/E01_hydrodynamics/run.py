"""E01 - Do our water forces match the published model, and is the bare hull unstable in yaw?

Run from the project root:
    /home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E01_hydrodynamics.run

Writes figures to figures/E01_hydrodynamics/ and key numbers to data/E01_hydrodynamics/summary.json.
Sign convention (Simos 2001): x forward, y to port, alpha = 180 deg = current from dead ahead, psi = pi - alpha,
psi > 0 = bow turned to port.
"""

import json
from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from fishtail import hydro, model_pinned
from fishtail.params import (FRESH_WATER, FULL_SCALE_CF_CHECK, SEA_WATER, TABLE_5_1_CF, TABLE_5_1_CT,
                             TABLE_5_1_U_MODEL, VIDAL_100, VIDAL_100_IPT, VIDAL_40, X_TCR_40_CHECK, froude_scaled)
from fishtail.plotting import COLORS, INK, INK_2, MUTED, ROOT, save, use_style

EXP = "E01_hydrodynamics"
RHO_RATIO = FRESH_WATER.rho / SEA_WATER.rho  # model tests in fresh water, ship in sea water
U_FULL = 1.90  # full-scale current speed used by Simos [m/s] (thesis p. 139 and p. 195)
U_MODEL = 0.20  # matching 1:90 model towing speed [m/s] (1.90 / sqrt(90); Table 5-1, p. 139)
LOADS = [(VIDAL_100, "100 % load"), (VIDAL_40, "40 % load (ballast)")]

# Simos (2001) Fig. 4.9, p. 142: his THEORY points (the x marks) for the yaw-rotating test of the loaded Vidal de
# Negreiros, U = 1.90 m/s and r = 0.33 deg/s full scale, forces at model scale. Read by eye by us from the printed
# figure (enlarged), so each value is only good to about +-1.5 N (or N m) and +-2 deg. Pairs: (heading deg, value).
SIMOS_FIG49_SWAY_N = [(0, 6), (10, 8), (21, 10.5), (31, 13), (42, 16), (52, 17.5), (94, 17.5), (104, 16.5),
                      (114, 15), (124, 12), (135, 8.5), (145, 5.5), (155, 1.5), (166, -2), (176, -4), (187, -6),
                      (207, -8.5), (218, -10.5), (228, -13), (239, -15.5), (249, -16.5), (276, -13.5), (290, -12),
                      (301, -10), (311, -7), (322, -4), (332, -1), (342, 2.5), (353, 4.5)]
SIMOS_FIG49_YAW_NM = [(0, -6), (10, -5), (21, -4.5), (31, -4.5), (42, -4.5), (52, -5.5), (63, -6.5), (73, -8),
                      (83, -9.5), (94, -11.5), (104, -14), (114, -17), (124, -18), (135, -18.5), (145, -17.5),
                      (156, -15), (166, -11.5), (177, -8), (187, -4), (198, -1.5), (208, 0), (218, 1), (229, 1),
                      (239, 0.5), (249, -0.5), (260, -2.5), (270, -4.5), (280, -7.5), (291, -9), (301, -10.5),
                      (312, -11), (322, -11), (332, -10), (343, -8.5), (353, -7)]


def fig01_static_coefficients(summary):
    """C1C, C2C, C6C against angle of attack, both loading conditions."""
    deg = np.linspace(0, 180, 721)
    alpha = np.radians(deg)
    fig, axes = plt.subplots(3, 1, figsize=(11.5, 9.0), sharex=True)
    for color, (hull, label) in zip(COLORS, LOADS):
        C0_model = hydro.friction_coefficient(U_MODEL, froude_scaled(hull, RHO_RATIO), FRESH_WATER)
        C0_full = hydro.friction_coefficient(U_FULL, hull, SEA_WATER)
        axes[0].plot(deg, hydro.C1C(alpha, hull, C0_model), color=color,
                     label=f"{label}, 1:90 model at {U_MODEL} m/s ($C_F$ = {C0_model:.4f})")
        axes[0].plot(deg, hydro.C1C(alpha, hull, C0_full), color=color, ls=":",
                     label=f"{label}, full scale at {U_FULL} m/s ($C_F$ = {C0_full:.4f})")
        axes[1].plot(deg, hydro.C2C(alpha, hull), color=color, label=f"{label} ($C_Y$ = {hull.C_Y})")
        axes[2].plot(deg, hydro.C6C(alpha, hull), color=color,
                     label=f"{label} ($l_p/L$ = {hull.lp_L:.3f} aft)")
        axes[1].plot(90, hull.C_Y, "o", color=color, ms=6, mec="white", mew=1.5)
        summary.setdefault("friction_C_F", {})[label] = {"model_0.20": C0_model, "full_1.90": C0_full}

    axes[0].set_ylabel("surge  $C_{1C} = F_X / (\\frac{1}{2}\\rho U^2 LT)$")
    axes[1].set_ylabel("sway  $C_{2C} = F_Y / (\\frac{1}{2}\\rho U^2 LT)$")
    axes[2].set_ylabel("yaw  $C_{6C} = N / (\\frac{1}{2}\\rho U^2 L^2 T)$")
    axes[0].set_title("Static current coefficients, Vidal de Negreiros (Leite et al. 1998 model)\n"
                      "eqs. 3.14, 3.21, 3.29; dots: $C_{2C}$(90°) = $C_Y$ (measured)", loc="left")
    axes[2].set_xlabel("angle of attack α [deg]\n(0° = current from astern, 90° = from starboard, "
                       "180° = from dead ahead)")
    axes[2].set_xticks(range(0, 181, 30))
    for ax in axes:
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))  # outside, so no curve is hidden
    fig.tight_layout()
    save(fig, EXP, "fig01_static_coefficients")


def fig02_yaw_moment_parts(summary):
    """The three parts of C6C (loaded ship), and the moment near head-on current as a function of heading."""
    hull = VIDAL_100
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    deg = np.linspace(0, 180, 721)
    cross_flow, munk, lift = hydro.C6C_terms(np.radians(deg), hull)
    ax1.plot(deg, munk, color=COLORS[1], label="Munk moment (turns the hull broadside)")
    ax1.plot(deg, lift, color=COLORS[0], label="hull lift + rudder (Clarke et al. 1982)")
    ax1.plot(deg, cross_flow, color=COLORS[2], label="cross-flow drag (acting $l_p$ aft)")
    ax1.plot(deg, cross_flow + munk + lift, color=INK, lw=2.4, label="total $C_{6C}$")
    ax1.axhline(0, color=MUTED, lw=0.8)
    ax1.set_xticks(range(0, 181, 30))
    ax1.set_xlabel("angle of attack α [deg]   (180° = current from dead ahead)")
    ax1.set_ylabel("$C_{6C} = N / (\\frac{1}{2}\\rho U^2 L^2 T)$")
    ax1.set_title("(a) Parts of the yaw moment, 100 % load", loc="left")
    ax1.legend(loc="upper left")

    psi_deg = np.linspace(-20, 20, 401)
    cross_flow, munk, lift = hydro.C6C_terms(np.pi - np.radians(psi_deg), hull)
    ax2.plot(psi_deg, munk, color=COLORS[1], label="Munk moment")
    ax2.plot(psi_deg, lift, color=COLORS[0], label="hull lift + rudder")
    ax2.plot(psi_deg, cross_flow, color=COLORS[2], label="cross-flow drag")
    ax2.plot(psi_deg, cross_flow + munk + lift, color=INK, lw=2.4, label="total")
    ax2.axhline(0, color=MUTED, lw=0.8)
    ax2.axvline(0, color=MUTED, lw=0.8)
    ax2.set_xlabel("heading ψ [deg]   (0° = bow straight into the current, > 0 = bow turned to port)")
    ax2.set_ylabel("$C_{6C}$  (> 0 turns the bow to port)")
    ax2.set_title("(b) Near head-on current: the total pushes the bow further away", loc="left")
    ax2.annotate("bow turned to port (ψ > 0)\n→ moment to port (> 0)\n→ turns further: UNSTABLE",
                 xy=(12, float(hydro.C6C(np.pi - np.radians(12), hull))), xytext=(-19, 0.045),
                 color=INK, fontsize=9, arrowprops=dict(arrowstyle="->", color=INK_2, lw=0.8))
    ax2.legend(loc="lower right")

    # Slopes d(C6C)/d(psi) at psi = 0, per radian. Since psi = pi - alpha, d/dpsi = -d/dalpha.
    h = 1e-8
    parts_plus = hydro.C6C_terms(np.pi - h, hull)
    parts_minus = hydro.C6C_terms(np.pi + h, hull)
    names = ["cross_flow", "munk", "lift"]
    slopes = {n: (p - m) / (2 * h) for n, p, m in zip(names, parts_plus, parts_minus)}
    slopes["total"] = sum(slopes.values())
    summary["dC6C_dpsi_at_head_on_100"] = slopes
    fig.tight_layout()
    save(fig, EXP, "fig02_yaw_moment_parts")


def pinned_fixed_points(hull):
    """Headings where the pinned hull is at rest: r = 0 and no yaw moment. Returns a sorted list in radians."""
    f = lambda psi: float(hydro.C6C(np.pi - psi, hull))
    grid = np.linspace(1e-3, np.pi - 1e-3, 3001)  # look strictly between 0 and 180 deg for sign changes
    values = [f(p) for p in grid]
    roots = [brentq(f, a, b) for a, b, fa, fb in zip(grid[:-1], grid[1:], values[:-1], values[1:]) if fa * fb < 0]
    return sorted([0.0, np.pi] + roots + [-x for x in roots])


def fig03_pinned_hull(summary):
    """Hull held by a pin at midship, free to turn: time history and phase portrait."""
    hull, fluid, U = VIDAL_100, SEA_WATER, 1.0  # results in t* = U t / L do not depend on U
    L = hull.L
    to_star = L / U  # time unit L/U [s]

    # Fixed points and their eigenvalues (in units of U/L so they are non-dimensional)
    fixed = []
    for psi in pinned_fixed_points(hull):
        eig = np.linalg.eigvals(model_pinned.jacobian(np.array([psi, 0.0]), hull, fluid, U)) * to_star
        fixed.append({"psi_deg": float(np.degrees(psi)), "eigenvalues_UoverL": [complex(e) for e in eig],
                      "stable": bool(np.all(eig.real < 0))})
    summary["pinned_hull_fixed_points_100"] = [
        {**p, "eigenvalues_UoverL": [[e.real, e.imag] for e in p["eigenvalues_UoverL"]]} for p in fixed]

    # Two runs starting 1 deg either side of head-on
    t_end = 25 * to_star
    runs = []
    for psi0 in [1.0, -1.0]:
        sol = solve_ivp(model_pinned.rhs, (0, t_end), [np.radians(psi0), 0.0], args=(hull, fluid, U),
                        rtol=1e-9, atol=1e-12, max_step=0.05 * to_star, dense_output=True)
        t = np.linspace(0, t_end, 2001)
        psi, r = sol.sol(t)
        runs.append((psi0, t / to_star, np.degrees(psi), r * to_star))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [1, 1.25]})
    for color, (psi0, t_star, psi_deg, r_star) in zip(COLORS, runs):
        ax1.plot(t_star, psi_deg, color=color, label=f"start ψ = {psi0:+.0f}°, ψ' = 0")
    stable = [p for p in fixed if p["stable"]]
    for p in stable:
        ax1.axhline(p["psi_deg"], color=MUTED, lw=0.8)
    ax1.text(25, stable[-1]["psi_deg"] + 4, f"stable heading ±{stable[-1]['psi_deg']:.0f}°", ha="right",
             color=INK_2, fontsize=9)
    ax1.set_xlabel("time  t* = U t / L")
    ax1.set_ylabel("heading ψ [deg]")
    ax1.set_title("(a) Released 1° off head-on, the hull swings away", loc="left")
    ax1.set_yticks(range(-90, 91, 30))
    ax1.legend(loc="center right")

    # Phase portrait: streamlines of (dpsi/dt*, dr*/dt*) on a grid
    r_max = 1.3 * max(np.abs(run[3]).max() for run in runs)
    psi_grid = np.radians(np.linspace(-180, 180, 73))
    r_grid = np.linspace(-r_max, r_max, 41)
    dpsi = np.zeros((r_grid.size, psi_grid.size))
    dr = np.zeros_like(dpsi)
    for i, r_star in enumerate(r_grid):
        for j, psi in enumerate(psi_grid):
            psi_dot, r_dot = model_pinned.rhs(0, [psi, r_star * U / L], hull, fluid, U)
            dpsi[i, j] = np.degrees(psi_dot * to_star)  # deg per unit t*
            dr[i, j] = r_dot * to_star**2  # change of r* per unit t*
    ax2.streamplot(np.degrees(psi_grid), r_grid, dpsi, dr, color=MUTED, density=1.1, linewidth=0.6,
                   arrowsize=0.8)
    for color, (psi0, t_star, psi_deg, r_star) in zip(COLORS, runs):
        ax2.plot(psi_deg, r_star, color=color, lw=2.0)
    for p in fixed:
        # psi = +180 and -180 deg are the same heading (stern into the current): draw it at both edges
        for x in ([-180.0, 180.0] if abs(p["psi_deg"]) == 180 else [p["psi_deg"]]):
            ax2.plot(x, 0, "o", ms=8, mew=1.6, color=INK, mfc=INK if p["stable"] else "white", zorder=5,
                     clip_on=False)
    ax2.plot([], [], "o", ms=8, color=INK, label="stable fixed point")
    ax2.plot([], [], "o", ms=8, mew=1.6, color=INK, mfc="white", label="unstable fixed point (saddle)")
    ax2.set_xlim(-180, 180)
    ax2.set_ylim(-r_max, r_max)
    ax2.set_xticks(range(-180, 181, 45))
    ax2.set_xlabel("heading ψ [deg]   (0° = bow into the current, ±180° = stern into the current)")
    ax2.set_ylabel("yaw rate  r* = r L / U")
    ax2.set_title("(b) Phase portrait: head-on is a saddle", loc="left")
    ax2.legend(loc="upper right", fontsize=8.5)
    fig.suptitle("Bare hull held by a pin at midship (no hawser), 100 % load", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig03_pinned_hull")


def fig04_turret_position(summary):
    """Yaw 'restoring stiffness' near head-on for a hull held by a pivot at x_T: where must the pivot be?"""
    xT = np.linspace(0, 0.5, 201)
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    summary["x_Tcr_over_L"] = {}
    for color, (hull, label) in zip(COLORS, LOADS):
        dC2, dC6 = hydro.head_on_slopes(hull)
        # Moment about the pivot: N - x_T F_Y. Restoring stiffness = -d(moment)/d(psi) = dC6/dalpha - (x_T/L) dC2/dalpha
        stiffness = dC6 - xT * dC2
        x_cr = hydro.critical_turret_position(hull)
        summary["x_Tcr_over_L"][label] = x_cr
        ax.plot(xT, stiffness, color=color, label=label)
        ax.plot(x_cr, 0, "^", ms=9, color=color, mec="white", mew=1.2, zorder=5)
        ax.annotate(f"$x_{{Tcr}}/L$ = {x_cr:.3f}", xy=(x_cr, 0), xytext=(x_cr + 0.012, -0.045 if x_cr > 0.3 else 0.03),
                    fontsize=9, color=INK)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.text(0.01, 0.075, "> 0: pushed back, points into the current (stable)", color=INK_2, fontsize=9)
    ax.text(0.27, -0.125, "< 0: turned away (unstable)", color=INK_2, fontsize=9)
    ax.set_xlabel("pivot (turret) position $x_T/L$, forward of midship")
    ax.set_ylabel("yaw restoring stiffness  $-\\,dC_{6C,pivot}/d\\psi$")
    ax.set_title("Hull held by a pivot at $x_T$, head-on current: static yaw stability\n"
                 f"▲ = critical position (eq. 3.37); thesis gives ≈ {X_TCR_40_CHECK} for 40 % (p. 198)", loc="left")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.92))
    fig.tight_layout()
    save(fig, EXP, "fig04_turret_position")


def fig05_sectional_drag(summary):
    """Our straight-line sideways drag distribution C_D(x) (decision D6)."""
    x = np.linspace(-0.5, 0.5, 201)
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    for color, (hull, label) in zip(COLORS, LOADS):
        ax.plot(x, hydro.sectional_drag(x, hull), color=color,
                label=f"{label}: $C_Y$ = {hull.C_Y}, $l_p/L$ = {hull.lp_L:.3f}")
        ax.axhline(hull.C_Y, color=color, lw=0.9, ls=":")
        ax.plot(-hull.lp_L, hull.C_Y, "o", color=color, ms=7, mec="white", mew=1.5, zorder=5)
    ax.plot([], [], ":", color=MUTED, label="dotted: uniform drag with the same $C_Y$")
    ax.plot([], [], "o", color=MUTED, label="dot: centre of the sideways drag ($l_p$ aft)")
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0, None)
    ax.set_xlabel("position along the hull  x/L   (−0.5 = stern, 0 = midship, +0.5 = bow)")
    ax.set_ylabel("sectional sideways drag coefficient $C_D(x)$")
    ax.set_title("Sideways drag along the hull: straight line through the two measured numbers\n"
                 "area under the line = $C_Y$; its centre is $l_p$ aft of midship", loc="left")
    ax.legend(loc="lower left", fontsize=8.5)
    fig.tight_layout()
    save(fig, EXP, "fig05_sectional_drag")


def fig06_yaw_rate_terms(summary):
    """Extra force and moment caused by turning, ship facing the current (u = U, v = 0), 100 % load."""
    hull, fluid, U = VIDAL_100, SEA_WATER, 1.0
    rho, L, T = fluid.rho, hull.L, hull.T
    r_star = np.linspace(-0.5, 0.5, 201)
    parts = [hydro.yaw_rate_terms(U, 0.0, rs * U / L, hull, fluid) for rs in r_star]
    get = lambda key: np.array([p[key] for p in parts])
    force_unit = 0.5 * rho * U**2 * L * T
    moment_unit = force_unit * L

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))
    ax1.plot(r_star, get("Y_wing") / force_unit, color=COLORS[0], label="wing / lift (Clarke $Y_r$ + $M_{11}ur$)")
    ax1.plot(r_star, get("Y_cross_flow") / force_unit, color=COLORS[2], label="cross-flow drag")
    ax1.plot(r_star, (get("Y_wing") + get("Y_cross_flow")) / force_unit, color=INK, lw=2.4, label="total")
    ax1.set_ylabel("extra sway force  $\\Delta F_Y / (\\frac{1}{2}\\rho U^2 LT)$")
    ax1.set_title("(a) Sway force from turning", loc="left")

    ax2.plot(r_star, get("N_wing") / moment_unit, color=COLORS[0], label="wing / lift (Clarke $N_r$ + $M_{26}ur$)")
    # With v = 0 and our straight-line C_D these two parts are exactly equal (both -(C_Y/32) r|r| in these units),
    # so the cross-flow line is drawn wider underneath the vortex line.
    ax2.plot(r_star, get("N_cross_flow") / moment_unit, color=COLORS[2], lw=4.5, label="cross-flow drag")
    ax2.plot(r_star, get("N_vortex") / moment_unit, color=COLORS[1], lw=1.6,
             label="bow/stern vortices (equal to cross-flow here)")
    total = (get("N_wing") + get("N_cross_flow") + get("N_vortex")) / moment_unit
    ax2.plot(r_star, total, color=INK, lw=2.4, label="total")
    ax2.set_ylabel("extra yaw moment  $\\Delta N / (\\frac{1}{2}\\rho U^2 L^2 T)$")
    ax2.set_title("(b) Yaw moment from turning: always resists the turn", loc="left")
    for ax in (ax1, ax2):
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.axvline(0, color=MUTED, lw=0.8)
        ax.set_xlabel("yaw rate  r* = r L / U")
        ax.legend(loc="best", fontsize=8.5)
    fig.suptitle("Ship facing the current (u = U, v = 0) and turning, 100 % load (Simos 2001, eq. 4.37)",
                 x=0.01, ha="left", fontsize=11)
    i = np.argmin(np.abs(r_star - 0.05))
    summary["yaw_damping_slope_100"] = float(total[i] / r_star[i])  # dN'/dr* for small r*
    fig.tight_layout()
    save(fig, EXP, "fig06_yaw_rate_terms")


def fig07_yaw_rotating_check(summary):
    """Reproduce the theory curve of Simos Fig. 4.9: model towed at constant speed while spinning about midship.

    Kinematics (4.39): u = U cos psi, v = -U sin psi, r constant. The figure shows the water's force on the hull.
    From (4.40) the rig must supply -(M22 - M11) U r cos psi - F_Y,R, so the water's total sideways force on the
    hull is F_Y,R + (M22 - M11) U r cos psi (the second part is the ideal-fluid inertia force, which (3.9) keeps on
    the left-hand side). The moment about midship is just N_R.
    """
    model = froude_scaled(VIDAL_100_IPT, RHO_RATIO)  # Simos used the IPT 1998 C_Y, l_p C_Y here (Table 4-5)
    fluid = FRESH_WATER
    U = U_MODEL
    r = np.radians(0.33) * np.sqrt(VIDAL_100.scale)  # 0.33 deg/s full scale; model time runs sqrt(90) times faster
    psi = np.radians(np.linspace(0, 360, 721))
    FY, N, FY_static, N_static = [], [], [], []
    for p in psi:
        u, v = U * np.cos(p), -U * np.sin(p)
        _, fy, n = hydro.current_forces(u, v, r, model, fluid, C0=0.0)  # C0 only affects F_X, not plotted
        _, fy0, n0 = hydro.static_forces(u, v, model, fluid, C0=0.0)
        inertia = (model.M22 - model.M11) * U * r * np.cos(p)
        FY.append(fy + inertia)
        N.append(n)
        FY_static.append(fy0)
        N_static.append(n0)
    deg = np.degrees(psi)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.4, 7.0), sharex=True)
    for ax, ours, static, simos, unit in [(ax1, FY, FY_static, SIMOS_FIG49_SWAY_N, "N"),
                                          (ax2, N, N_static, SIMOS_FIG49_YAW_NM, "N m")]:
        ax.plot(deg, static, color=MUTED, lw=1.2, ls=":", label="ours, static part only (r = 0)")
        ax.plot(deg, ours, color=COLORS[0], label="ours, full model (static + yaw-rate terms)")
        sx, sy = np.array(simos).T
        ax.plot(sx, sy, "x", color=INK, ms=6, mew=1.4, label="Simos Fig. 4.9 theory (read by eye, ±1.5)")
        ax.axhline(0, color=MUTED, lw=0.8)
        ours_at_points = np.interp(sx, deg, ours)
        summary.setdefault("fig49_rms_difference", {})[unit] = float(np.sqrt(np.mean((ours_at_points - sy) ** 2)))
    ax1.set_ylim(-40, 40)
    ax2.set_ylim(-40, 10)
    ax1.set_ylabel("sway force on hull [N]")
    ax2.set_ylabel("yaw moment about midship [N m]")
    ax2.set_xlabel("heading ψ [deg]  (hull spinning at constant rate while towed)")
    ax2.set_xticks(range(0, 361, 45))
    ax1.legend(loc="upper right", fontsize=8.5)
    ax1.set_title(f"Check against Simos Fig. 4.9: yaw-rotating test, Vidal de Negreiros 100 %, model 1:90\n"
                  f"U = {U} m/s, r = {np.degrees(r):.2f} deg/s (full scale 1.90 m/s, 0.33 deg/s); "
                  f"$C_Y$ = {model.C_Y}, $l_pC_Y$ = 0.044", loc="left")
    fig.tight_layout()
    save(fig, EXP, "fig07_yaw_rotating_check")


def fig08_friction(summary):
    """Friction coefficient C_F against current speed, model and full scale, with Simos's Table 5-1."""
    U_full = np.linspace(0.3, 6.0, 300)  # full-scale current speed [m/s]
    lam = VIDAL_100.scale
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), sharey=True)
    for ax, (hull, label), key in zip(axes, LOADS, ["100", "40"]):
        model = froude_scaled(hull, RHO_RATIO)
        model_k025 = froude_scaled(replace(hull, k_form=0.25), RHO_RATIO)  # Leite et al.'s form factor
        ax.plot(U_full, hydro.friction_coefficient(U_full, hull, SEA_WATER), color=COLORS[0],
                label="full ship (sea water), k = 0.39")
        ax.plot(U_full, hydro.friction_coefficient(U_full / np.sqrt(lam), model, FRESH_WATER), color=COLORS[1],
                label="1:90 model (fresh water), k = 0.39")
        ax.plot(U_full, hydro.friction_coefficient(U_full / np.sqrt(lam), model_k025, FRESH_WATER), color=COLORS[1],
                lw=1.0, ls=":", label="1:90 model, k = 0.25")
        U_tab = np.array(TABLE_5_1_U_MODEL) * np.sqrt(lam)
        lows, highs = np.array(TABLE_5_1_CT[key]).T
        ax.vlines(U_tab, lows, highs, color=INK_2, lw=3, alpha=0.35, label="Table 5-1: measured $C_T$ range (model)")
        ax.plot(U_tab, TABLE_5_1_CF[key], "x", color=INK, ms=7, mew=1.5, label="Table 5-1: Simos's $C_F$ (model)")
        if key == "40":
            ax.plot(FULL_SCALE_CF_CHECK["U"], FULL_SCALE_CF_CHECK["C_F"], "o", color=INK, mfc="white", ms=8,
                    mew=1.5, label="p. 195: full scale, 1.90 m/s → 0.015")
        ax.set_title(f"({'a' if key == '100' else 'b'}) {label}", loc="left")
        ax.set_xlabel("current speed, full scale [m/s]   (model speed = full ÷ √90)")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=8)  # below, clear of curves
    axes[0].set_ylabel("friction coefficient $C_F$ (= $C_0$, eq. 3.20)")
    axes[0].set_ylim(0, 0.06)
    fig.suptitle("Faster current → larger Reynolds number → less friction. Model tests (small Re) have more "
                 "friction than the real ship.", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig08_friction")


def main():
    use_style()
    summary = {}
    fig01_static_coefficients(summary)
    fig02_yaw_moment_parts(summary)
    fig03_pinned_hull(summary)
    fig04_turret_position(summary)
    fig05_sectional_drag(summary)
    fig06_yaw_rate_terms(summary)
    fig07_yaw_rotating_check(summary)
    fig08_friction(summary)

    folder = ROOT / "data" / EXP
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
