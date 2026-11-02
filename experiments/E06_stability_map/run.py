"""E06 - Where in the design space is the moored ship safe? (two-parameter stability maps; main result)

Run from the project root:
    /home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E06_stability_map.run

Maps in the planes (L_H/L, C_F), (x_A/L, L_H/L) and (x_A/L, C_F); validation against Simos (2001) Fig. 5.3.
Writes figures to figures/E06_stability_map/ and numbers to data/E06_stability_map/summary.json.
Sign convention (Simos 2001): psi > 0 = bow turned to port, gamma > 0 = ship swung to port. Time t* = U t / L.
"""

import json
from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from fishtail import continuation as co
from fishtail import hydro
from fishtail import model_rigid as mr
from fishtail.params import (FRESH_WATER, SEA_WATER, SIMOS_FIG53, TABLE_5_1_CT, TABLE_5_1_U_MODEL, VIDAL_100,
                             VIDAL_40, froude_scaled)
from fishtail.plotting import COLORS, INK, INK_2, MUTED, ROOT, save, use_style

EXP = "E06_stability_map"
LOADS = {"100": (VIDAL_100, "100 % load"), "40": (VIDAL_40, "40 % load (ballast)")}
REGIONS = ["trails straight", "sits at an angle", "fishtails"]  # 0, 1, 2
REGION_COLORS = [COLORS[2], COLORS[0], COLORS[1]]  # aqua, blue, orange (categorical slots, fixed order)


def fig01_simos_comparison(summary):
    """Our Hopf boundary in the (L_H/L, C_F) plane against Simos Fig. 5.3 (digitised)."""
    L_H = np.linspace(0.02, 2.0, 100)
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8))
    summary["fig01"] = {}
    for ax, color, key in zip(axes[:2], COLORS, ["100", "40"]):
        hull, label = LOADS[key]
        x_A = SIMOS_FIG53[key]["x_A_L"]
        s = mr.Rigid(hull, x_A, 0.5, 0.03)
        ours = co.hopf_curve(s, L_H)
        ax.fill_between(L_H, 0, ours, color=COLORS[1], alpha=0.08, lw=0)
        ax.plot(L_H, ours, color=color, lw=2.2, label="ours (linearised full model)")
        ax.plot(SIMOS_FIG53["L_H_L"], SIMOS_FIG53[key]["C_F"], "o", color=INK, ms=6, mfc="white", mew=1.4,
                label="Simos Fig. 5.3 (digitised)")
        ours_at = co.hopf_curve(s, SIMOS_FIG53["L_H_L"])
        rel = (ours_at - np.array(SIMOS_FIG53[key]["C_F"])) / np.array(SIMOS_FIG53[key]["C_F"])
        summary["fig01"][key] = {"x_A_L": x_A, "ours_at_simos_points": list(ours_at),
                                 "relative_difference": list(rel), "max_abs_rel_diff": float(np.abs(rel).max()),
                                 "mean_rel_diff": float(rel.mean())}
        ax.text(1.0, 0.012, "unstable: fishtails", ha="center", color=INK_2, fontsize=9)
        ax.text(1.0, 0.055, "stable: trails straight", ha="center", color=INK_2, fontsize=9)
        ax.set_xlim(0, 2)
        ax.set_ylim(0, 0.06)
        ax.set_xlabel("hawser length  $L_H/L$")
        ax.set_ylabel("friction coefficient $C_F$")
        ax.set_title(f"({'a' if key == '100' else 'b'}) {label}, $x_A/L$ = {x_A}: mean difference "
                     f"{100 * rel.mean():+.1f} %", loc="left", fontsize=10)
        ax.legend(loc="lower right", fontsize=8.5)
    # (c) What Annex B gives if the lift from the sideways velocity is left out, as eq. B.13 is printed
    ax = axes[2]
    for color, key in zip(COLORS, ["100", "40"]):
        hull, label = LOADS[key]
        s = mr.Rigid(hull, SIMOS_FIG53[key]["x_A_L"], 0.5, 0.03)
        ax.semilogy(L_H, co.hopf_curve(s, L_H), color=color, lw=2.2, label=f"{label}: full linearisation")
        ax.semilogy(L_H, co.hopf_curve(s, L_H, simos_annex_b=True), color=color, lw=1.4, ls=(0, (4, 2.5)),
                    label=f"{label}: eq. B.13 as printed")
        ax.semilogy(SIMOS_FIG53["L_H_L"], SIMOS_FIG53[key]["C_F"], "o", color=INK, ms=5, mfc="white", mew=1.2)
    ax.plot([], [], "o", color=INK, ms=5, mfc="white", mew=1.2, label="Simos Fig. 5.3")
    ax.set_xlim(0, 2)
    ax.set_xlabel("hawser length  $L_H/L$")
    ax.set_ylabel("friction coefficient $C_F$ (log scale)")
    ax.set_title("(c) Leaving out the sway-velocity lift (as B.13\nis printed) does NOT give Simos's curve",
                 loc="left", fontsize=10)
    ax.legend(loc="lower right", fontsize=7.5)
    fig.suptitle("Validation: stability boundary (Hopf curve) of trailing straight, Vidal de Negreiros; "
                 "stable above the curve", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig01_simos_fig53_comparison")


def fig02_speeds_on_map(summary):
    """Like Simos Fig. 5.4: where tank tests and the real ship sit on the map."""
    L_H = np.linspace(0.02, 2.0, 100)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), sharey=True)
    summary["fig02"] = {}
    rho_ratio = FRESH_WATER.rho / SEA_WATER.rho
    for ax, key in zip(axes, ["100", "40"]):
        hull, label = LOADS[key]
        s = mr.Rigid(hull, SIMOS_FIG53[key]["x_A_L"], 0.5, 0.03)
        ax.plot(L_H, co.hopf_curve(s, L_H), color=INK, lw=2.2)
        ax.fill_between(L_H, 0, co.hopf_curve(s, L_H), color=COLORS[1], alpha=0.08, lw=0)
        model = froude_scaled(hull, rho_ratio)
        for color, U, (lo, hi) in zip(COLORS, TABLE_5_1_U_MODEL, TABLE_5_1_CT[key]):
            ax.axhspan(lo, hi, color=color, alpha=0.18, lw=0)
            ax.axhline(0.5 * (lo + hi), color=color, lw=1.4, label=f"1:90 model, {U} m/s (Table 5-1 $C_T$ range)")
        for ls, U in [((0, (4, 2.5)), 1.0), ((0, (1, 1.5)), 1.9)]:
            C = float(hydro.friction_coefficient(U, hull, SEA_WATER))
            ax.axhline(C, color=INK_2, lw=1.4, ls=ls, label=f"full-scale ship, {U} m/s ($C_F$ = {C:.4f})")
            summary["fig02"][f"{key} full {U}"] = C
        ax.set_xlim(0, 2)
        ax.set_ylim(0, 0.06)
        ax.set_xlabel("hawser length  $L_H/L$")
        ax.set_title(f"({'a' if key == '100' else 'b'}) {label}, $x_A/L$ = {SIMOS_FIG53[key]['x_A_L']}",
                     loc="left")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=7.5)
    axes[0].set_ylabel("friction coefficient $C_F$")
    fig.suptitle("Why tank tests looked stable: the model's friction (bands) sits much higher than the real "
                 "ship's (grey lines). Below the black curve the ship fishtails.", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig02_speeds_on_map")


def classify_point(hull, x_A_L, L_H_L, C_F, heading):
    """0 = trails straight, 1 = sits at an angle, 2 = fishtails (no stable rest position: the straight position
    and, if they exist, the angled positions are all unstable). `heading` = angled heading for this x_A (or None)."""
    s = mr.Rigid(hull, x_A_L, L_H_L, C_F)
    if co.is_stable(co.eigenvalues(np.zeros(4), s)):
        return 0
    if heading is not None and co.is_stable(co.eigenvalues(co.equilibrium_state(heading, s), s)):
        return 1
    return 2


def region_map(hull, xs, ys, fixed, which):
    """Classify a grid. which = 'xA-LH' (fixed C_F) or 'xA-CF' (fixed L_H)."""
    grid = np.zeros((len(ys), len(xs)), dtype=int)
    for j, xa in enumerate(xs):
        heads = co.equilibrium_headings(hull, xa)
        heading = heads[0] if heads else None
        for i, y in enumerate(ys):
            if which == "xA-LH":
                grid[i, j] = classify_point(hull, xa, y, fixed, heading)
            else:
                grid[i, j] = classify_point(hull, xa, fixed, y, heading)
    return grid


def straight_hopf_growth(hull, xs, ys, make):
    """Growth rate of the swing of trailing straight on a grid, for drawing the Hopf curve as its zero contour.

    Aft of x_Tcr the straight position is already a saddle (pitchfork), so the curve means nothing there: NaN.
    make(xa, y) builds the set-up for one grid point.
    """
    x_cr = hydro.critical_turret_position(hull)
    return np.array([[co.max_oscillatory_growth(make(xa, y)) if xa >= x_cr else np.nan for xa in xs] for y in ys])


def draw_regions(ax, xs, ys, grid):
    cmap = ListedColormap(REGION_COLORS)
    ax.pcolormesh(xs, ys, grid, cmap=cmap, vmin=-0.5, vmax=2.5, alpha=0.28, shading="nearest", lw=0)


def region_legend(ax, loc="upper right"):
    handles = [Patch(color=c, alpha=0.45, label=r) for c, r in zip(REGION_COLORS, REGIONS)]
    handles.append(plt.Line2D([], [], color=INK, lw=1.6, label="Hopf curve (trailing straight)"))
    handles.append(plt.Line2D([], [], color=INK_2, lw=1.2, ls=(0, (4, 2.5)), label="pitchfork $x_A = x_{Tcr}$"))
    ax.legend(handles=handles, loc=loc, fontsize=8)


def fig03_xA_LH_maps(summary):
    """Regime maps in the (x_A/L, L_H/L) plane at a model-like and a full-scale-like friction."""
    xs = np.linspace(0.0, 0.6, 61)
    ys = np.linspace(0.05, 1.5, 59)
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.5), sharex=True, sharey=True)
    summary["fig03"] = {}
    for row, key in enumerate(["100", "40"]):
        hull, label = LOADS[key]
        x_cr = hydro.critical_turret_position(hull)
        for col, C_F in enumerate([0.035, 0.015]):
            ax = axes[row, col]
            grid = region_map(hull, xs, ys, C_F, "xA-LH")
            draw_regions(ax, xs, ys, grid)
            growth = straight_hopf_growth(hull, xs, ys, lambda xa, LH: mr.Rigid(hull, xa, LH, C_F))
            ax.contour(xs, ys, growth, levels=[0.0], colors=[INK], linewidths=1.6)
            ax.axvline(x_cr, color=INK_2, lw=1.2, ls=(0, (4, 2.5)))
            ax.plot(x_cr, 0.05, "^", ms=10, color=INK, mec="white", mew=1.2, clip_on=False, zorder=6)
            ax.set_title(f"({'abcd'[2 * row + col]}) {label}, $C_F$ = {C_F} "
                         f"({'tank-model-like' if C_F > 0.02 else 'full-scale-like'})", loc="left", fontsize=10)
            summary["fig03"][f"{key} C_F={C_F}"] = {r: float(np.mean(grid == k)) for k, r in enumerate(REGIONS)}
            if row == 1:
                ax.set_xlabel("attachment point  $x_A/L$  (0.5 = bow)")
            if col == 0:
                ax.set_ylabel("hawser length  $L_H/L$")
    region_legend(axes[0, 1])
    fig.suptitle("What the ship does, by attachment point and hawser length (rigid hawser). "
                 "Left of the dashed line: pitchfork; black curve: Hopf", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig03_map_attachment_vs_length")


def fig04_xA_CF_maps(summary):
    """Regime maps in the (x_A/L, C_F) plane at fixed hawser length."""
    xs = np.linspace(0.0, 0.6, 61)
    ys = np.linspace(0.004, 0.08, 58)
    L_H_L = 0.5
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0), sharey=True)
    for ax, key in zip(axes, ["100", "40"]):
        hull, label = LOADS[key]
        x_cr = hydro.critical_turret_position(hull)
        grid = region_map(hull, xs, ys, L_H_L, "xA-CF")
        draw_regions(ax, xs, ys, grid)
        growth = straight_hopf_growth(hull, xs, ys, lambda xa, C: mr.Rigid(hull, xa, L_H_L, C))
        ax.contour(xs, ys, growth, levels=[0.0], colors=[INK], linewidths=1.6)
        ax.axvline(x_cr, color=INK_2, lw=1.2, ls=(0, (4, 2.5)))
        ax.set_xlabel("attachment point  $x_A/L$  (0.5 = bow)")
        ax.set_title(f"({'a' if key == '100' else 'b'}) {label}, $L_H/L$ = {L_H_L}", loc="left")
    axes[0].set_ylabel("friction coefficient $C_F$  (low = fast current, full scale)")
    region_legend(axes[1], loc="upper right")
    fig.suptitle("By attachment point and friction: attaching further forward does not always help. Ballast ship: a "
                 "narrow window just forward of $x_{Tcr}$ is stable even at low friction", x=0.01, ha="left",
                 fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig04_map_attachment_vs_friction")


def outcome_of_run(hull, x_A_L, L_H_L, C_F):
    """Time-domain check: simulate from a small disturbance and say what the ship ends up doing."""
    s = mr.Rigid(hull, x_A_L, L_H_L, C_F)
    t, y = mr.simulate(s, [np.radians(3.0), np.radians(1.0), 0, 0], t_end_star=600.0, n_out=6001, rtol=1e-8,
                       atol=1e-10)
    late = y[:, t > 500]
    swing = np.degrees(late[0].max() - late[0].min())
    if swing > 0.1:
        return 2, swing
    return (1 if abs(np.degrees(late[1, -1])) > 0.5 else 0), swing


def fig05_time_domain_check(summary):
    """Step 5 check from the brief: a handful of simulations placed in each region of the map agree with it."""
    hull, label = LOADS["100"]
    C_F = 0.035
    xs = np.linspace(0.0, 0.6, 61)
    ys = np.linspace(0.05, 1.5, 59)
    grid = region_map(hull, xs, ys, C_F, "xA-LH")
    points = [(0.15, 0.3), (0.25, 1.0), (0.3, 0.6), (0.33, 0.15), (0.38, 0.1), (0.4, 0.3), (0.45, 1.2),
              (0.5, 0.2), (0.5, 0.5), (0.5, 1.0), (0.55, 0.12), (0.05, 0.4)]
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    draw_regions(ax, xs, ys, grid)
    growth = straight_hopf_growth(hull, xs, ys, lambda xa, LH: mr.Rigid(hull, xa, LH, C_F))
    ax.contour(xs, ys, growth, levels=[0.0], colors=[INK], linewidths=1.6)
    ax.axvline(hydro.critical_turret_position(hull), color=INK_2, lw=1.2, ls=(0, (4, 2.5)))
    markers = {0: ("o", "settled straight"), 1: ("s", "settled at an angle"), 2: ("*", "kept swinging")}
    results, agree = [], 0
    for xa, LH in points:
        predicted = classify_point(hull, xa, LH, C_F, (co.equilibrium_headings(hull, xa) or [None])[0])
        got, swing = outcome_of_run(hull, xa, LH, C_F)
        agree += int(predicted == got)
        m, _ = markers[got]
        ax.plot(xa, LH, m, ms=13 if m == "*" else 10, color=REGION_COLORS[got], mec=INK, mew=1.2, zorder=6)
        results.append({"x_A_L": xa, "L_H_L": LH, "map": REGIONS[predicted], "simulation": REGIONS[got],
                        "late_swing_deg": swing})
    summary["fig05"] = {"C_F": C_F, "runs": results, "agree": f"{agree} of {len(points)}"}
    handles = [plt.Line2D([], [], marker=m, ls="none", color=REGION_COLORS[k], mec=INK, ms=10, label=f"run: {t}")
               for k, (m, t) in markers.items()]
    handles += [Patch(color=c, alpha=0.45, label=f"map: {r}") for c, r in zip(REGION_COLORS, REGIONS)]
    ax.legend(handles=handles, loc="upper left", fontsize=8, ncol=2)
    ax.set_xlabel("attachment point  $x_A/L$")
    ax.set_ylabel("hawser length  $L_H/L$")
    ax.set_title(f"Check: {agree} of {len(points)} simulations (t* = 600, released with γ = 3°, ψ = 1°) end where "
                 f"the map says\n{label}, $C_F$ = {C_F}; marker shape and colour = what the simulation did",
                 loc="left", fontsize=10)
    fig.tight_layout()
    save(fig, EXP, "fig05_time_domain_check")


def main():
    use_style()
    summary = {}
    fig01_simos_comparison(summary)
    fig02_speeds_on_map(summary)
    fig03_xA_LH_maps(summary)
    fig04_xA_CF_maps(summary)
    fig05_time_domain_check(summary)
    folder = ROOT / "data" / EXP
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / "summary.json", "w") as f:
        json.dump(summary, f, indent=1, default=float)
    print(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()
