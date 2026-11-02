"""E10 - What does a steady wind from the side do to the pitchfork? (imperfect pitchfork, Strogatz section 3.6)

Run from the project root:
    /home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E10_wind.run

Knobs: wind force (as a fraction of the steady hawser tension F_H0) and wind direction beta_W; the attachment point
x_A/L is the pitchfork knob of E03. Rigid hawser, Vidal de Negreiros 100 %, L_H/L = 0.5, C_F = 0.05 (as E03 and the
pitchfork demo: nothing fishtails, so only rest positions matter). Wind model: decision D28 (a constant force, fixed
in earth axes, acting at midship unless stated).
Writes figures to figures/E10_wind/ and numbers to data/E10_wind/summary.json.
Sign convention (Simos 2001): psi > 0 = bow turned to port, gamma > 0 = ship swung to port.
beta_W = direction the wind comes FROM: 0 = from ahead (with the current), +90 deg = from port, -90 deg = from
starboard, 180 deg = from astern. Time t* = U t / L.
"""

import json
from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np

from fishtail import continuation as co
from fishtail import hydro
from fishtail import model_rigid as mr
from fishtail.params import VIDAL_100
from fishtail.plotting import COLORS, INK, INK_2, MUTED, ROOT, current_arrow, draw_mooring, plot_branch, save, use_style

EXP = "E10_wind"
HULL = VIDAL_100
S0 = mr.Rigid(HULL, x_A_L=0.25, L_H_L=0.5, C_F=0.05)  # x_A/L = 0.25: aft of the pitchfork point, two sides possible
X_CR = hydro.critical_turret_position(HULL)  # pitchfork point without wind, 0.359


def windy(ratio, beta_deg, s=S0, **kw):
    """Set-up s with a wind of size ratio * F_H0 coming from beta_deg."""
    return replace(s, C_W=ratio * s.C_F, beta_W=np.radians(beta_deg), **kw)


def draw_diagram(ax, s, color=INK, label=None, with_folds=True, lw=2.0):
    """Heading at rest against x_A/L for the wind of s: solid stable, dashed unstable, fold points as diamonds."""
    br = co.wind_branches(s)
    first = True
    for seg in br["segments"]:
        plot_branch(ax, seg["x"], np.degrees(seg["psi"]), seg["stable"], color, label=label if first else None, lw=lw)
        first = False
    if br["straight"] is not None:
        plot_branch(ax, br["straight"]["x"], np.zeros(br["straight"]["x"].size), br["straight"]["stable"], color, lw=lw)
        ax.plot(br["pitchfork"], 0.0, "^", ms=8, color=INK, zorder=5)
    if with_folds:
        for x, p in br["folds"]:
            ax.plot(x, np.degrees(p), "D", ms=6, color=INK, mfc="white", mew=1.5, zorder=5)
    return br


def no_wind_reference(ax):
    """The perfect pitchfork (no wind) as a thin grey line."""
    draw_diagram(ax, S0, color=MUTED, lw=1.0, with_folds=False)


def style_diagram(ax, title):
    ax.axvline(S0.x_A_L, color=MUTED, lw=0.8, ls=":")
    ax.set_xlim(0, 0.6)
    ax.set_ylim(-60, 60)
    ax.set_xlabel("attachment point x_A/L (forward of midship)")
    ax.set_ylabel("heading at rest ψ [deg] (+ = bow to port)")
    ax.set_title(title)


# --- fig01: the imperfect pitchfork for a beam wind of growing size ---------------------------------------------------

def fig01(summary):
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
    out = {}
    for ax, ratio in zip(axes.flat, [0.0, 0.05, 0.10, 0.20]):
        if ratio > 0:
            no_wind_reference(ax)
        br = draw_diagram(ax, windy(ratio, 90.0))
        out[f"{ratio:.2f}"] = [[round(x, 4), round(np.degrees(p), 2)] for x, p in br["folds"]]
        style_diagram(ax, "no wind: perfect pitchfork" if ratio == 0 else
                      f"wind from port (β = 90°), force = {ratio:.2f} F_H0")
        if ratio > 0:
            ax.text(0.02, 0.04, "grey: no wind", transform=ax.transAxes, color=INK_2, fontsize=8.5)
    for ax in axes[0]:
        ax.set_xlabel("")
    for ax in axes[:, 1]:
        ax.set_ylabel("")
    axes[0, 0].plot([], [], "^", color=INK, label="pitchfork")
    axes[0, 0].plot([], [], "D", color=INK, mfc="white", label="fold (saddle-node)")
    axes[0, 0].plot([], [], color=INK, label="stable")
    axes[0, 0].plot([], [], color=INK, ls=(0, (4, 2.5)), label="unstable")
    axes[0, 0].legend(loc="upper right")
    fig.suptitle("E10 fig01 - A side wind breaks the pitchfork (rigid hawser, L_H/L = 0.5, C_F = 0.05)\n"
                 "dotted line: x_A/L = 0.25 used in figs 03 and 05", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig01_imperfect_pitchfork")
    summary["fig01_folds_beam_wind"] = out


# --- fig02: the wind direction ----------------------------------------------------------------------------------------

def fig02():
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.6), sharex=True, sharey=True)
    names = {0: "from ahead (with the current)", 30: "30° off the bow, port", 60: "60°, port side",
             90: "from port (beam)", 135: "135°, port quarter", 180: "from astern"}
    for ax, beta in zip(axes.flat, names):
        no_wind_reference(ax)
        draw_diagram(ax, windy(0.10, beta))
        style_diagram(ax, f"β = {beta}°: {names[beta]}")
    for ax in axes[0]:
        ax.set_xlabel("")
    for ax in axes[:, 1:].flat:
        ax.set_ylabel("")
    fig.suptitle("E10 fig02 - Wind direction, force = 0.10 F_H0 (grey: no wind). Wind along the current keeps the "
                 "symmetry and only moves the pitchfork;\nany sideways part unfolds it", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig02_wind_direction")


# --- fig03: hysteresis at fixed attachment point ----------------------------------------------------------------------

def signed_beam(h):
    """Beam wind of signed size h (> 0 from port, < 0 from starboard), h in units of F_H0."""
    return windy(abs(h), 90.0 if h >= 0 else -90.0)


def s_curve(s=S0, psi_max_deg=60.0, n=601):
    """Rest heading against signed beam-wind force at fixed x_A/L. The rest condition is linear in the wind force:
       C6C - x_A C2C + (x_W - x_A) Y_w* = 0  with  Y_w* = -C_W cos(psi) for wind from port,
    so every heading has one wind force: C_W = (C6C - x_A C2C) / ((x_W - x_A) cos psi)."""
    psi = np.radians(np.linspace(-psi_max_deg, psi_max_deg, n))
    alpha = np.pi - psi
    C_W = (hydro.C6C(alpha, s.hull) - s.x_A_L * hydro.C2C(alpha, s.hull)) / ((s.x_W_L - s.x_A_L) * np.cos(psi))
    h = C_W / s.C_F
    stable = []
    for p, hh in zip(psi, h):
        sk = replace(signed_beam(hh), x_A_L=s.x_A_L, x_W_L=s.x_W_L)
        stable.append(co.is_stable(co.eigenvalues(co.equilibrium_state(p, sk), sk)))
    return psi, h, np.array(stable)


def sweep(h_values, state0, t_each=120.0):
    """Change the side wind step by step, letting the ship settle each time (like slowly rising wind)."""
    psi_end, state = [], np.asarray(state0, dtype=float)
    for h in h_values:
        _, y = mr.simulate(signed_beam(h), state, t_each, n_out=3, rtol=1e-8, atol=1e-10)
        state = y[:, -1]
        psi_end.append(np.degrees(state[1]))
    return np.array(psi_end)


def fig03(summary):
    psi, h, stable = s_curve()
    keep = np.abs(h) <= 0.35
    fig, ax = plt.subplots(figsize=(8.5, 6))
    plot_branch(ax, h[keep], np.degrees(psi[keep]), stable[keep], INK, label="rest headings (stable / unstable)")
    dh = np.diff(h)
    folds = [(h[k + 1], np.degrees(psi[k + 1])) for k in range(dh.size - 1) if dh[k] * dh[k + 1] < 0]
    for fh, fp in folds:
        ax.plot(fh, fp, "D", ms=7, color=INK, mfc="white", mew=1.5, zorder=5)
    h_up = np.round(np.arange(-0.30, 0.3001, 0.01), 3)
    start = co.wind_equilibria(signed_beam(-0.30))[0]  # settled to starboard in a strong starboard wind
    up = sweep(h_up, start)
    down = sweep(h_up[::-1], co.wind_equilibria(signed_beam(0.30))[-1])
    ax.plot(h_up, up, "o", ms=4, color=COLORS[0], label="simulation, wind turning from starboard to port →")
    ax.plot(h_up[::-1], down, "s", ms=4, mfc="none", color=COLORS[1], label="simulation, back again ←")
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.set_xlabel("side wind force / F_H0   (< 0: from starboard, > 0: from port; beam wind β = ±90°)")
    ax.set_ylabel("heading at rest ψ [deg] (+ = bow to port)")
    ax.set_title("E10 fig03 - Hysteresis at x_A/L = 0.25: the ship keeps its side until the wind\n"
                 "passes a fold, then jumps (◇ = folds; L_H/L = 0.5, C_F = 0.05)")
    ax.legend(loc="upper left")
    fig.tight_layout()
    save(fig, EXP, "fig03_hysteresis")
    jump_up = h_up[np.argmax(up > 0)]
    jump_down = h_up[::-1][np.argmax(down < 0)]
    summary["fig03"] = {"folds_h_psi": [[round(a, 4), round(b, 2)] for a, b in folds],
                        "sweep_jump_up_at_h": float(jump_up), "sweep_jump_down_at_h": float(jump_down)}


# --- fig04: the cusp ----------------------------------------------------------------------------------------------------

def fold_forces(x_A_L):
    """Side-wind forces (units F_H0) at the two folds for attachment point x_A_L (empty if forward of the cusp)."""
    _, h = s_curve_shape(x_A_L)
    dh = np.diff(h)
    return [h[k + 1] for k in range(dh.size - 1) if dh[k] * dh[k + 1] < 0]


def s_curve_shape(x_A_L, n=4001):
    """As s_curve but without the stability (only the shape is needed to find the folds)."""
    psi = np.radians(np.linspace(-60.0, 60.0, n))
    alpha = np.pi - psi
    C_W = (hydro.C6C(alpha, HULL) - x_A_L * hydro.C2C(alpha, HULL)) / ((S0.x_W_L - x_A_L) * np.cos(psi))
    return psi, C_W / S0.C_F


def fig04(summary):
    xs = np.linspace(0.02, X_CR - 1e-4, 300)
    lo, hi = [], []
    for x in xs:
        f = fold_forces(x)
        lo.append(min(f) if f else np.nan)
        hi.append(max(f) if f else np.nan)
    lo, hi = np.array(lo), np.array(hi)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5.4))
    ax.fill_between(xs, lo, hi, color=COLORS[0], alpha=0.15, lw=0, label="two stable rest positions (one per side)")
    ax.plot(xs, lo, color=COLORS[0])
    ax.plot(xs, hi, color=COLORS[0], label="fold lines (saddle-node)")
    ax.plot(X_CR, 0, "^", ms=9, color=INK, label=f"cusp = the old pitchfork, x_A/L = {X_CR:.3f}")
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xlim(0, 0.45)
    ax.set_ylim(-1.0, 1.0)
    ax.set_xlabel("attachment point x_A/L")
    ax.set_ylabel("side wind force / F_H0  (+ from port, beam wind)")
    ax.set_title("Where the ship can sit on either side (C_F = 0.05)")
    ax.text(0.40, 0.55, "one rest\nposition", ha="center", color=INK_2)
    ax.legend(loc="lower left")

    # Growth law near the cusp. Strogatz's cubic normal form x' = h + r x - x^3 gives h_fold ~ r^(3/2);
    # our pitchfork has a quadratic one, x' = h + r x - x|x| (E03 fig06), which gives h_fold = r^2 / 4 ~ r^2.
    d = X_CR - xs
    near = (d > 0.003) & (d < 0.06)
    slope = np.polyfit(np.log(d[near]), np.log(hi[near]), 1)[0]
    ax2.loglog(d, hi, color=COLORS[0], label=f"model: slope {slope:.2f}")
    ref = hi[near][len(hi[near]) // 2]
    dref = d[near][len(d[near]) // 2]
    ax2.loglog(d, ref * (d / dref) ** 2, color=MUTED, ls="--", label="slope 2: x' = h + r x − x|x|")
    ax2.loglog(d, ref * (d / dref) ** 1.5, color=MUTED, ls=":", label="slope 3/2: Strogatz x' = h + r x − x³")
    ax2.set_xlabel("distance aft of the pitchfork point, x_Tcr/L − x_A/L")
    ax2.set_ylabel("fold wind force / F_H0")
    ax2.set_title("Width of the cusp near its tip")
    ax2.legend(loc="upper left")
    fig.suptitle("E10 fig04 - The cusp: the imperfect pitchfork in the plane (attachment point, side wind)",
                 fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig04_cusp")
    summary["fig04"] = {"cusp_slope": round(float(slope), 3),
                        "fold_force_at_x_0.25": [round(float(v), 4) for v in fold_forces(0.25)]}


# --- fig05: what the ship does when the wind comes up --------------------------------------------------------------------

def fig05(summary):
    start = co.equilibrium_state(co.equilibrium_headings(HULL, S0.x_A_L)[0], S0)  # resting to port, no wind
    ratios = [0.05, 0.10, 0.14, 0.16, 0.20]
    fig = plt.figure(figsize=(13, 7.5))
    ax_p = fig.add_subplot(2, 2, 1)
    ax_g = fig.add_subplot(2, 2, 3, sharex=ax_p)
    ends = {}
    for c, ratio in zip(COLORS, ratios):
        s = windy(ratio, -90.0)
        t, y = mr.simulate(s, start, 300.0, n_out=1201, rtol=1e-8, atol=1e-10)
        ax_p.plot(t, np.degrees(y[1]), color=c, label=f"{ratio:.2f} F_H0")
        ax_g.plot(t, np.degrees(y[0]), color=c)
        ends[ratio] = y[:, -1]
    ax_p.set_ylabel("heading ψ [deg]")
    ax_g.set_ylabel("hawser angle γ [deg]")
    ax_g.set_xlabel("time t* = U t / L")
    ax_p.set_title("Ship resting to port; wind from starboard switched on at t* = 0")
    ax_p.legend(title="wind force", loc="center right")
    for ax in (ax_p, ax_g):
        ax.axhline(0, color=MUTED, lw=0.8)
    for k, ratio in enumerate([0.10, 0.20]):
        ax = fig.add_subplot(2, 2, 2 + 2 * k)
        s = windy(ratio, -90.0)
        draw_mooring(ax, start, s, MUTED, alpha=0.6)
        draw_mooring(ax, ends[ratio], s, COLORS[0])
        current_arrow(ax, 0.25, 0.42)
        ax.annotate("", xy=(-0.9, 0.42), xytext=(-0.9, 0.18), arrowprops=dict(arrowstyle="-|>", color=INK_2, lw=1.4))
        ax.text(-0.87, 0.30, "wind", color=INK_2, fontsize=8.5)
        ax.set_aspect("equal")
        ax.set_xlim(-1.25, 0.45)
        ax.set_ylim(-0.6, 0.6)
        ax.set_xlabel("X / L (upstream →)")
        ax.set_ylabel("Y / L (port ↑)")
        side = "stays on the port side" if ends[ratio][1] > 0 else "swings across to starboard"
        ax.set_title(f"wind {ratio:.2f} F_H0 from starboard: {side}\n(grey: start, blue: after 300 t*)")
    fig.suptitle("E10 fig05 - A weak wind cannot push the ship across; a stronger one (past the fold, ≈ 0.15 F_H0) "
                 "can. x_A/L = 0.25, L_H/L = 0.5, C_F = 0.05", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig05_time_runs")
    summary["fig05_end_psi_deg"] = {f"{r:.2f}": round(float(np.degrees(e[1])), 2) for r, e in ends.items()}


# --- fig06: bow attachment (no pitchfork): how far the wind pushes the ship aside --------------------------------------

def fig06(summary):
    betas = np.linspace(-180, 180, 145)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), sharex=True)
    out = {}
    for c, ratio in zip(COLORS, [0.1, 0.2, 0.4]):
        g, p, F, st = [], [], [], []
        for b in betas:
            s = windy(ratio, b, x_A_L=0.5)
            eqs = co.wind_equilibria(s)
            e = min(eqs, key=lambda q: abs(q[1]))  # the one near the current line (there is only one here)
            g.append(np.degrees(e[0]))
            p.append(np.degrees(e[1]))
            F.append(mr.hawser_force_star(e, s) / s.C_F)
            st.append(co.is_stable(co.eigenvalues(e, s)))
        plot_branch(axes[0], betas, g, st, c)
        axes[0].plot([], [], color=c, label=f"{ratio:.1f} F_H0")  # legend entry (solid) only
        plot_branch(axes[1], betas, p, st, c)
        plot_branch(axes[2], betas, F, st, c)
        out[f"{ratio:.1f}"] = {"max_gamma_deg": round(max(np.abs(g)), 2), "max_psi_deg": round(max(np.abs(p)), 2),
                               "F_range": [round(min(F), 3), round(max(F), 3)], "all_stable": bool(all(st))}
    for ax, yl in zip(axes, ["hawser angle γ [deg] (+ = to port)", "heading ψ [deg] (+ = bow to port)",
                             "hawser force F_H / F_H0"]):
        ax.set_ylabel(yl)
        ax.set_xlabel("wind comes from β_W [deg]  (0 ahead, +90 port, ±180 astern)")
        ax.set_xticks([-180, -90, 0, 90, 180])
        ax.axvline(0, color=MUTED, lw=0.8)
    axes[0].plot([], [], color=INK_2, ls=(0, (4, 2.5)), label="dashed: unstable (fishtails)")
    axes[0].legend(title="wind force")
    fig.suptitle("E10 fig06 - Hawser at the bow (x_A/L = 0.5, no pitchfork): one rest position, pushed aside by the "
                 "wind. Note γ > 0 for wind from port: the ship ends UP-wind (see README).\nC_F = 0.05; a strong wind "
                 "from astern lowers the hawser tension enough to start fishtailing (dashed, see fig08)", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig06_bow_attachment")
    summary["fig06"] = out


# --- fig07: where the wind acts --------------------------------------------------------------------------------------------

def fig07():
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), sharey=True)
    for ax, x_W in zip(axes, [-0.30, 0.0, 0.5]):
        no_wind_reference(ax)
        draw_diagram(ax, windy(0.10, 90.0, x_W_L=x_W))
        ax.axvline(x_W, color=COLORS[1], lw=1.2, ls="--") if 0 <= x_W <= 0.6 else None
        where = {-0.30: "aft (superstructure)", 0.0: "at midship (our default)", 0.5: "at the bow"}[x_W]
        style_diagram(ax, f"wind acts at x_W/L = {x_W:+.2f}, {where}")
    for ax in axes[1:]:
        ax.set_ylabel("")
    fig.suptitle("E10 fig07 - Where the wind acts decides which side it favours: aft of the hawser point it turns the "
                 "bow into the wind (orange: x_A = x_W). Wind from port, 0.10 F_H0", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig07_point_of_action")


# --- fig08: does the wind move the fishtailing (Hopf) boundary? ----------------------------------------------------------

def hopf_with_wind(C_W, beta_deg, s=None):
    """Friction C_F at which the rest position (bow attachment) starts to fishtail, with a wind C_W fixed."""
    s = s or replace(S0, x_A_L=0.5)

    def growth(c):
        sc = replace(s, C_F=c, C_W=C_W, beta_W=np.radians(beta_deg))
        eqs = co.wind_equilibria(sc) if C_W > 0 else [np.zeros(4)]
        if len(eqs) != 1:
            return np.nan
        return co.oscillatory_growth(eqs[0], sc)

    grid = np.geomspace(0.004, 0.15, 50)
    roots = co.find_crossings(growth, grid)
    return roots[-1] if roots else np.nan


def fig08(summary):
    C_Ws = np.linspace(0.0, 0.006, 7)
    fig, ax = plt.subplots(figsize=(8, 5.4))
    out = {}
    for c, beta in zip(COLORS, [0, 90, 180]):
        hs = [hopf_with_wind(cw, beta) for cw in C_Ws]
        ax.plot(C_Ws, hs, "o-", color=c, label=f"wind from β = {beta}°")
        out[str(beta)] = [round(float(v), 5) for v in hs]
        for cw, hv in zip(C_Ws, hs):
            if np.isnan(hv):  # checked: the rest position stays stable for every C_F in 0.004-0.15
                ax.plot(cw, 0.022, "x", ms=9, mew=2, color=c)
                ax.annotate("no fishtailing at any C_F\nin 0.004–0.15", (cw, 0.022), xytext=(-150, -8),
                            textcoords="offset points", color=INK_2, fontsize=8.5,
                            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
        out["note"] = "NaN = no Hopf point found: stable for every C_F in 0.004-0.15"
    ax.set_xlabel("wind force C_W = F_wind / (½ ρ U² L T)")
    ax.set_ylabel("Hopf friction coefficient C_F (fishtails below it)")
    ax.set_title("E10 fig08 - Wind and the fishtailing boundary\n(bow attachment x_A/L = 0.5, L_H/L = 0.5)")
    ax.legend()
    fig.tight_layout()
    save(fig, EXP, "fig08_hopf_with_wind")
    summary["fig08_hopf_C_F"] = {"C_W": [round(float(v), 4) for v in C_Ws], **out}


def main():
    use_style()
    summary = {"x_Tcr": X_CR, "set_up": {"x_A_L": S0.x_A_L, "L_H_L": S0.L_H_L, "C_F": S0.C_F, "x_W_L": S0.x_W_L}}
    fig01(summary)
    fig02()
    fig03(summary)
    fig04(summary)
    fig05(summary)
    fig06(summary)
    fig07()
    fig08(summary)
    folder = ROOT / "data" / EXP
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / "summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
