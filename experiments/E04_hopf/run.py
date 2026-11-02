"""E04 - When does the ship start to fishtail? (Hopf bifurcation of the trailing-straight position)

Run from the project root:
    /home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E04_hopf.run

Knob: friction coefficient C_F (stands for current speed, decision D1); also the hawser length L_H/L.
Vidal de Negreiros 100 %, rigid hawser attached at the bow (x_A/L = 0.5, as in Simos's tests).
Writes figures to figures/E04_hopf/ and numbers to data/E04_hopf/summary.json.
Sign convention (Simos 2001): psi > 0 = bow turned to port, gamma > 0 = ship swung to port. Time t* = U t / L.
"""

import json
from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from fishtail import continuation as co
from fishtail import model_rigid as mr
from fishtail.params import VIDAL_100
from fishtail.plotting import COLORS, INK, INK_2, MUTED, ROOT, current_arrow, draw_mooring, save, use_style

EXP = "E04_hopf"
BASE = mr.Rigid(VIDAL_100, x_A_L=0.5, L_H_L=0.5, C_F=0.03)
BLUES = LinearSegmentedColormap.from_list("blue_ramp", ["#9ec5f4", "#2a78d6", "#0d366b"])  # one hue, light -> dark


def fig01_eigenvalue_paths(summary, C_crit):
    """Eigenvalues of trailing straight in the complex plane as the friction changes."""
    C_values = np.geomspace(0.005, 0.1, 120)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 5.0), gridspec_kw={"width_ratios": [1.5, 1]})
    for ax in (ax1, ax2):
        for C in C_values:
            eig = np.linalg.eigvals(mr.linear_jacobian_star(replace(BASE, C_F=C)))
            ax.scatter(eig.real, eig.imag, s=9, color=BLUES((np.log(C) - np.log(0.005)) / np.log(20)), lw=0)
        eig_c = np.linalg.eigvals(mr.linear_jacobian_star(replace(BASE, C_F=C_crit)))
        pair = eig_c[np.abs(eig_c.imag) > 1e-9]
        ax.plot(pair.real, pair.imag, "o", ms=9, color=INK, mfc="white", mew=1.6, zorder=5)
        ax.axvline(0, color=INK_2, lw=1.0)
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.set_xlabel("real part  [units of U/L]   (> 0: grows)")
    ax1.set_ylabel("imaginary part  [units of U/L]   (= angular frequency)")
    ax1.set_title("(a) All four eigenvalues", loc="left")
    ax2.set_xlim(-0.12, 0.16)
    ax2.set_ylim(0.2, 1.1)
    ax2.set_title("(b) Zoom on the upper member of the pair", loc="left")
    ax2.annotate("friction falls:\npair crosses into\nthe right half-plane", xy=(0.0, float(np.abs(pair.imag).max())),
                 xytext=(0.04, 0.95), fontsize=9, color=INK, arrowprops=dict(arrowstyle="->", color=INK_2, lw=0.8))
    sm = plt.cm.ScalarMappable(cmap=BLUES, norm=plt.matplotlib.colors.LogNorm(0.005, 0.1))
    cb = fig.colorbar(sm, ax=ax2, pad=0.02)
    cb.set_label("friction coefficient $C_F$")
    ax2.plot([], [], "o", ms=9, color=INK, mfc="white", mew=1.6, ls="none", label=f"Hopf point, $C_F$ = {C_crit:.4f}")
    ax2.legend(loc="lower right", fontsize=8.5)
    fig.suptitle(f"Eigenvalues of 'trailing straight' as $C_F$ goes from 0.005 (light) to 0.1 (dark):  "
                 f"100 % load, $x_A/L$ = {BASE.x_A_L}, $L_H/L$ = {BASE.L_H_L}", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig01_eigenvalue_paths")


def fig02_growth_and_period(summary):
    """Growth rate and period of the swing against friction, for three hawser lengths."""
    C_values = np.geomspace(0.004, 0.1, 150)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.6))
    summary["fig02"] = {}
    for color, L_H_L in zip(COLORS, [0.2, 0.5, 1.0]):
        s = replace(BASE, L_H_L=L_H_L)
        growth = [co.max_oscillatory_growth(replace(s, C_F=C)) for C in C_values]
        period = [2 * np.pi / co.hopf_frequency(replace(s, C_F=C)) for C in C_values]
        C_h = co.hopf_C_F(s)
        ax1.plot(C_values, growth, color=color, label=f"$L_H/L$ = {L_H_L}")
        ax1.plot(C_h, 0, "o", ms=8, color=color, mec="white", mew=1.4, zorder=5)
        ax2.loglog(C_values, period, color=color, label=f"$L_H/L$ = {L_H_L}")
        summary["fig02"][str(L_H_L)] = {"C_F_hopf": C_h, "period_at_hopf": 2 * np.pi / co.hopf_frequency(
            replace(s, C_F=C_h))}
    ax1.axhline(0, color=INK_2, lw=1.0)
    ax1.set_xscale("log")
    ax1.set_xlabel("friction coefficient $C_F$")
    ax1.set_ylabel("growth rate of the swing, Re λ  [U/L]")
    ax1.set_title("(a) Below the Hopf point (●) a small swing grows", loc="left")
    ax1.legend(loc="upper right")
    C_ref = np.array([0.006, 0.08])
    ax2.loglog(C_ref, 9.0 * np.sqrt(0.03 / C_ref), color=MUTED, ls="--", lw=1.0, label="∝ $1/\\sqrt{C_F}$ (pendulum)")
    ax2.set_yticks([3, 4, 6, 10, 20, 30])
    ax2.yaxis.set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
    ax2.yaxis.set_minor_formatter(plt.matplotlib.ticker.NullFormatter())
    ax2.set_xlabel("friction coefficient $C_F$")
    ax2.set_ylabel("period of the swing  2π/ω  [t* units, L/U]")
    ax2.set_title("(b) Period of the linear swing: roughly a pendulum\nwhose 'gravity' is the friction drag",
                  loc="left")
    ax2.legend(loc="upper right")
    fig.suptitle(f"Trailing straight, 100 % load, $x_A/L$ = {BASE.x_A_L}: eigenvalue pair of the swing",
                 x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig02_growth_and_period")


def run(s, start, t_end=300.0, n=6001):
    return mr.simulate(s, start, t_end_star=t_end, n_out=n)


def fig03_phase_portraits(summary, C_crit):
    """Phase portraits (projections) just before and just after the Hopf point."""
    cases = [(1.15 * C_crit, "(a, b) Before: $C_F$ = 1.15 × Hopf value - the swing dies out (stable focus)"),
             (0.85 * C_crit, "(c, d) After: $C_F$ = 0.85 × Hopf value - small swings grow, large ones shrink: "
                             "a limit cycle")]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.0))
    summary["fig03"] = {}
    for row, (C, title) in enumerate(cases):
        s = replace(BASE, C_F=C)
        starts = [([np.radians(2.0), 0, 0, 0], COLORS[0], "start γ = 2°"),
                  ([np.radians(30.0), 0, 0, 0], COLORS[1], "start γ = 30°")]
        for start, color, label in starts:
            t, y = run(s, start)
            axes[row, 0].plot(np.degrees(y[0]), np.degrees(y[2]), color=color, lw=0.9, label=label)
            axes[row, 1].plot(np.degrees(y[1]), np.degrees(y[3]), color=color, lw=0.9, label=label)
            summary["fig03"][f"{C:.4f} {label}"] = float(np.degrees(np.abs(y[0, -1000:]).max()))
        for col, (xl, yl) in enumerate([("hawser angle γ [deg]", "dγ/dt*  [deg per L/U]"),
                                        ("heading ψ [deg]", "dψ/dt*  [deg per L/U]")]):
            ax = axes[row, col]
            ax.plot(0, 0, "o", ms=8, color=INK, mfc=INK if row == 0 else "white", mew=1.6, zorder=5)
            ax.set_xlabel(xl)
            ax.set_ylabel(yl)
            ax.legend(loc="upper right", fontsize=8.5)
        axes[row, 0].set_title(title, loc="left", fontsize=10)
    fig.suptitle(f"Phase portraits (2-D views of the 4-D motion, so curves may cross), 100 % load, "
                 f"$L_H/L$ = {BASE.L_H_L}; ● stable, ○ unstable equilibrium", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig03_phase_portraits")


def fig04_time_series(summary, C_crit):
    """Time series before and after the Hopf point: angles and hawser force."""
    fig, axes = plt.subplots(3, 2, figsize=(12.0, 8.0), sharex=True)
    for col, C in enumerate([1.15 * C_crit, 0.85 * C_crit]):
        s = replace(BASE, C_F=C)
        t, y = run(s, [np.radians(2.0), 0, 0, 0], t_end=250.0, n=5001)
        F = np.array([mr.hawser_force_star(y[:, i], s) for i in range(t.size)]) / C
        axes[0, col].plot(t, np.degrees(y[0]), color=COLORS[0])
        axes[1, col].plot(t, np.degrees(y[1]), color=COLORS[1])
        axes[2, col].plot(t, F, color=COLORS[2])
        axes[0, col].set_title(f"{'Before' if col == 0 else 'After'} the Hopf point: $C_F$ = {C:.4f}", loc="left")
        axes[2, col].set_xlabel("time  t* = U t / L")
    axes[0, 0].set_ylabel("hawser angle γ [deg]")
    axes[1, 0].set_ylabel("heading ψ [deg]")
    axes[2, 0].set_ylabel("hawser force $F_H / F_{H0}$")
    for ax in axes.ravel():
        ax.axhline(0 if ax not in axes[2] else 1, color=MUTED, lw=0.8)
    fig.suptitle(f"Released with γ = 2°: 100 % load, $x_A/L$ = {BASE.x_A_L}, $L_H/L$ = {BASE.L_H_L}. "
                 f"$F_{{H0}}$ = tension when trailing straight", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig04_time_series")


def fig05_top_view(summary, C_crit):
    """Top view of a fully grown fishtailing cycle: where the ship goes behind the buoy."""
    s = replace(BASE, C_F=0.6 * C_crit)
    t, y = run(s, [np.radians(30.0), 0, 0, 0], t_end=400.0, n=8001)
    late = t > 300
    tl, yl = t[late], y[:, late]
    # one period: between two maxima of gamma
    idx = np.where((yl[2, :-1] > 0) & (yl[2, 1:] <= 0))[0]
    i0, i1 = idx[0], idx[1]
    period = tl[i1] - tl[i0]
    summary["fig05"] = {"C_F": s.C_F, "period": float(period),
                        "gamma_amp_deg": float(np.degrees(np.abs(yl[0]).max())),
                        "psi_amp_deg": float(np.degrees(np.abs(yl[1]).max()))}
    fig, ax = plt.subplots(figsize=(10.0, 6.4))
    bow = np.array([mr.positions(yl[:, i], s)["bow"] for i in range(i0, i1 + 1)])
    stern = np.array([mr.positions(yl[:, i], s)["stern"] for i in range(i0, i1 + 1)])
    ax.plot(bow[:, 0], bow[:, 1], color=COLORS[1], lw=1.2, label="path of the bow")
    ax.plot(stern[:, 0], stern[:, 1], color=COLORS[2], lw=1.2, label="path of the stern")
    for k, i in enumerate(np.linspace(i0, i1, 9)[:-1].astype(int)):
        draw_mooring(ax, yl[:, i], s, COLORS[0], alpha=0.35 + 0.65 * (k % 4 == 0), fill=(k % 4 == 0))
    current_arrow(ax, 0.1, 0.5)
    ax.set_aspect("equal")
    ax.set_xlim(-1.6, 0.15)
    ax.set_ylim(-0.5, 0.62)
    ax.set_xlabel("X / L  (upstream →)")
    ax.set_ylabel("Y / L  (port ↑)")
    ax.legend(loc="lower right")
    ax.set_title(f"Fishtailing seen from above: 8 snapshots over one period (T* = {period:.1f} L/U)\n"
                 f"100 % load, $x_A/L$ = {s.x_A_L}, $L_H/L$ = {s.L_H_L}, $C_F$ = {s.C_F:.4f} (0.6 × Hopf value); "
                 f"γ swings ±{np.degrees(np.abs(yl[0]).max()):.0f}°, ψ ±{np.degrees(np.abs(yl[1]).max()):.0f}°",
                 loc="left")
    fig.tight_layout()
    save(fig, EXP, "fig05_top_view")


def fig06_hawser_force(summary, C_crit):
    """Hawser force during fully grown fishtailing for three friction values."""
    fig, ax = plt.subplots(figsize=(10.0, 4.6))
    summary["fig06"] = {}
    for color, frac in zip(COLORS, [0.85, 0.6, 0.4]):
        s = replace(BASE, C_F=frac * C_crit)
        t, y = run(s, [np.radians(30.0), 0, 0, 0], t_end=360.0, n=7201)
        late = t > 300
        F = np.array([mr.hawser_force_star(y[:, i], s) for i in np.where(late)[0]]) / s.C_F
        ax.plot(t[late] - 300, F, color=color, label=f"$C_F$ = {frac} × Hopf value ({s.C_F:.4f})")
        summary["fig06"][str(frac)] = {"F_max_over_F0": float(F.max()), "F_min_over_F0": float(F.min())}
    ax.axhline(1, color=MUTED, lw=0.8)
    ax.axhline(0, color=INK_2, lw=1.0)
    ax.text(59, 0.03, "below 0: the bar would be pushed (a real rope goes slack)", ha="right", va="bottom",
            color=INK_2, fontsize=8.5)
    ax.set_xlabel("time  t* = U t / L  (after the transient)")
    ax.set_ylabel("hawser force  $F_H / F_{H0}$")
    ax.set_title(f"Hawser force during fishtailing: twice per swing it peaks, and at low friction it nearly goes "
                 f"slack\n100 % load, $x_A/L$ = {BASE.x_A_L}, $L_H/L$ = {BASE.L_H_L}; $F_{{H0}}$ = tension when "
                 f"trailing straight", loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=8.5)  # below, clear of peaks
    fig.tight_layout()
    save(fig, EXP, "fig06_hawser_force")


def main():
    use_style()
    C_crit = co.hopf_C_F(BASE)
    summary = {"C_F_hopf_base": C_crit, "omega_at_hopf": co.hopf_frequency(replace(BASE, C_F=C_crit)),
               "set_up": {"hull": BASE.hull.name, "x_A_L": BASE.x_A_L, "L_H_L": BASE.L_H_L}}
    fig01_eigenvalue_paths(summary, C_crit)
    fig02_growth_and_period(summary)
    fig03_phase_portraits(summary, C_crit)
    fig04_time_series(summary, C_crit)
    fig05_top_view(summary, C_crit)
    fig06_hawser_force(summary, C_crit)
    folder = ROOT / "data" / EXP
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
