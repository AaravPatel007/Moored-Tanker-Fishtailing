"""E07 - What changes with a real, stretchy hawser that can go slack?

Run from the project root (a few minutes):
    /home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E07_elastic_hawser.run

Rope law fitted to a real SPM nylon hawser (Huang et al. 2026, Fig. 3): tension = a * strain^n.
Ship: Vidal de Negreiros 100 %, full scale, rope at the bow (x_A/L = 0.5), unstretched length L_H/L = 0.5 (160 m).
Writes figures to figures/E07_elastic_hawser/ and numbers to data/E07_elastic_hawser/summary.json.
Sign convention (Simos 2001): psi > 0 = bow turned to port, gamma > 0 = ship swung to port. Time t* = U t / L.
"""

import json
from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from fishtail import continuation as co
from fishtail import hydro
from fishtail import model_elastic as me
from fishtail import model_rigid as mr
from fishtail.params import HUANG_HAWSER, SEA_WATER, VIDAL_100
from fishtail.plotting import COLORS, INK, INK_2, MUTED, ROOT, save, use_style

EXP = "E07_elastic_hawser"
HULL = VIDAL_100
START = (np.radians(5.0), np.radians(1.0))  # released with gamma = 5 deg, psi = 1 deg


def full_scale(U, **kw):
    """Elastic set-up at full scale, friction from the ITTC line at current speed U."""
    C = float(hydro.friction_coefficient(U, HULL, SEA_WATER))
    return me.Elastic(HULL, kw.pop("x_A_L", 0.5), kw.pop("L_H_L", 0.5), C, U=U, **kw)


def run(s, t_end=400.0, n=8001):
    t, y, T = me.simulate(s, me.from_rigid_start(*START, s), t_end, n_out=n)
    g = np.array([me.angles(y[:, i], s)[0] for i in range(t.size)])
    return t, y, T, g


def fishtail_pair(s):
    """The eigenvalue pair that continues the rigid model's swing: the complex pair with the lowest frequency."""
    eig = me.eigenvalues_star(s)
    cplx = eig[eig.imag > 1e-6]
    return cplx[np.argmin(cplx.imag)]


def fig01_rope_law(summary):
    """The measured hawser curve, our power-law fit, and where our ship sits on it."""
    e = np.array(HUANG_HAWSER["elongation_m"]) / HUANG_HAWSER["length_m"]
    T = np.array(HUANG_HAWSER["tension_kN"])
    strain = np.linspace(0, 0.22, 200)
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    ax.plot(100 * e, T, "o", color=INK, ms=7, mfc="white", mew=1.5, label="Huang et al. (2026) Fig. 3, prototype (read by eye)")
    ax.plot(100 * strain, me.ROPE_A * strain**me.ROPE_N / 1e3, color=COLORS[0],
            label=f"fit: T = {me.ROPE_A / 1e3:,.0f} kN × strain$^{{{me.ROPE_N:.2f}}}$")
    ax.plot(100 * strain, (T[-1] / e[-1]) * strain, color=MUTED, ls="--", lw=1.0, label="straight line through the last point (n = 1)")
    summary["rope"] = {"a_kN": me.ROPE_A / 1e3, "n": me.ROPE_N, "operating": {}}
    for color, U in zip(COLORS[1:], [0.5, 1.0, 2.0]):
        s = full_scale(U)
        _, e0 = me.trailing_straight(s)
        T0 = me.ROPE_A * e0**me.ROPE_N / 1e3
        ax.plot(100 * e0, T0, "s", color=color, ms=9, mec="white", mew=1.4, zorder=5,
                label=f"our ship trailing straight in {U} m/s: {T0:.0f} kN, strain {100 * e0:.1f} %")
        summary["rope"]["operating"][str(U)] = {"tension_kN": T0, "strain": e0}
    ax.set_xlabel("strain = stretch / unstretched length [%]")
    ax.set_ylabel("tension [kN]")
    ax.set_title("The rope law: a real 144 mm nylon SPM hawser (breaking load 650 t) gets stiffer as it stretches",
                 loc="left", fontsize=10)
    ax.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    save(fig, EXP, "fig01_rope_law")


def fig02_stiff_limit(summary):
    """Check from the brief: a very stiff rope must behave like the rigid hawser."""
    s = me.Elastic(HULL, 0.5, 0.5, 0.03, U=1.0)
    factors = np.geomspace(0.03, 3000, 30)
    pairs = [fishtail_pair(replace(s, rope_a=s.rope_a * f)) for f in factors]
    eig_r = co.eigenvalues(np.zeros(4), mr.Rigid(HULL, 0.5, 0.5, 0.03))
    rigid_pair = eig_r[np.abs(eig_r.imag) > 1e-9][0]
    summary["fig02"] = {"rigid_pair": [rigid_pair.real, abs(rigid_pair.imag)],
                        "elastic_x1000": [pairs[-3].real, pairs[-3].imag]}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.6))
    ax1.semilogx(factors, [p.real for p in pairs], color=COLORS[0], label="stretchy: growth rate Re λ")
    ax1.semilogx(factors, [p.imag for p in pairs], color=COLORS[1], label="stretchy: frequency Im λ")
    ax1.axhline(rigid_pair.real, color=COLORS[0], ls="--", lw=1.0, label="rigid hawser")
    ax1.axhline(abs(rigid_pair.imag), color=COLORS[1], ls="--", lw=1.0)
    ax1.axvline(1, color=MUTED, lw=0.8)
    ax1.text(1.1, 0.3, "real hawser", color=INK_2, fontsize=9)
    ax1.set_xlabel("rope stiffness ÷ real hawser's")
    ax1.set_ylabel("swing eigenvalue [U/L]")
    ax1.set_title("(a) The swing's eigenvalues approach the rigid ones", loc="left")
    ax1.legend(loc="center right", fontsize=8.5)

    t_r, y_r = mr.simulate(mr.Rigid(HULL, 0.5, 0.5, 0.03), [START[0], START[1], 0, 0], t_end_star=120, n_out=1201)
    ax2.plot(t_r, np.degrees(y_r[1]), color=INK, lw=3.0, label="rigid hawser")
    for color, f in zip(COLORS, [1000, 1]):
        s2 = replace(s, rope_a=s.rope_a * f)
        t, y, T = me.simulate(s2, me.from_rigid_start(*START, s2), 120, n_out=1201)
        ax2.plot(t, np.degrees(y[2]), color=color, lw=1.4, label=f"stretchy, {f} × real stiffness")
    ax2.set_xlabel("time t* = U t / L")
    ax2.set_ylabel("heading ψ [deg]")
    ax2.set_title("(b) Same start: stiff rope = rigid hawser", loc="left")
    ax2.legend(loc="upper left", fontsize=8.5)
    fig.suptitle("Stiff-rope check: 100 % load, $x_A/L$ = 0.5, $L_H/L$ = 0.5, $C_F$ = 0.03, U = 1 m/s",
                 x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig02_stiff_limit")


def fig03_real_ship(summary):
    """The real ship on the real rope in a 1 m/s current: fishtailing, slack rope and snatch loads."""
    s = full_scale(1.0)
    t, y, T, g = run(s, t_end=400.0)
    sr = mr.Rigid(HULL, 0.5, 0.5, s.C_F)
    t_r, y_r = mr.simulate(sr, [START[0], START[1], 0, 0], t_end_star=400, n_out=8001)
    F_r = np.array([mr.hawser_force_star(y_r[:, i], sr) for i in range(t_r.size)]) * 0.5 * 1025 * 320 * HULL.T * s.U**2
    F0 = s.C_F * 0.5 * s.fluid.rho * s.U**2 * HULL.L * HULL.T
    late = t > 300
    summary["fig03"] = {"C_F": s.C_F, "F_H0_kN": F0 / 1e3,
                        "gamma_amp_deg": float(np.degrees(0.5 * (g[late].max() - g[late].min()))),
                        "psi_amp_deg": float(np.degrees(0.5 * (y[2, late].max() - y[2, late].min()))),
                        "T_max_kN": float(T[late].max() / 1e3), "slack_fraction": float(np.mean(T[late] == 0)),
                        "rigid_F_max_kN": float(F_r[t_r > 300].max() / 1e3),
                        "rigid_F_min_kN": float(F_r[t_r > 300].min() / 1e3)}
    fig, axes = plt.subplots(3, 1, figsize=(11.0, 8.6), sharex=True)
    axes[0].plot(t_r, np.degrees(y_r[0]), color=MUTED, lw=1.0, label="rigid hawser")
    axes[0].plot(t, np.degrees(g), color=COLORS[0], label="stretchy hawser")
    axes[1].plot(t_r, np.degrees(y_r[1]), color=MUTED, lw=1.0, label="rigid hawser")
    axes[1].plot(t, np.degrees(y[2]), color=COLORS[1], label="stretchy hawser")
    axes[2].plot(t_r, F_r / 1e3, color=MUTED, lw=1.0, label="rigid hawser (negative = pushing)")
    axes[2].plot(t, T / 1e3, color=COLORS[2], label="stretchy hawser")
    slack = T == 0
    axes[2].fill_between(t, 0, 1, where=slack, color=COLORS[1], alpha=0.25, lw=0,
                         transform=axes[2].get_xaxis_transform(), label="rope slack")
    axes[2].axhline(F0 / 1e3, color=INK_2, lw=0.8, ls=":")
    axes[2].axhline(0, color=INK, lw=0.8)
    axes[0].set_ylabel("hawser angle γ [deg]")
    axes[1].set_ylabel("heading ψ [deg]")
    axes[2].set_ylabel("hawser tension [kN]")
    axes[2].set_xlabel("time t* = U t / L   (1 unit = 320 s here)")
    for ax in axes:
        ax.legend(loc="upper left", fontsize=8.5)
    fig.suptitle(f"Real ship (100 %), real rope, 1 m/s current ($C_F$ = {s.C_F:.4f}), 160 m hawser at the bow: "
                 f"it fishtails, the rope goes slack and snaps tight (dotted: steady tension {F0 / 1e3:.0f} kN)",
                 x=0.01, ha="left", fontsize=10.5)
    fig.tight_layout()
    save(fig, EXP, "fig03_real_ship_1ms")


def fig04_speed(summary):
    """Swing size and rope loads against current speed (full scale): speed now acts through friction AND rope."""
    speeds = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5]
    rows = []
    for U in speeds:
        s = full_scale(U)
        t, y, T, g = run(s, t_end=400.0, n=6001)
        late = t > 300
        F0 = s.C_F * 0.5 * s.fluid.rho * U**2 * HULL.L * HULL.T
        rows.append({"U": U, "C_F": s.C_F, "F_H0_kN": F0 / 1e3,
                     "gamma_amp": float(np.degrees(0.5 * (g[late].max() - g[late].min()))),
                     "T_max_kN": float(T[late].max() / 1e3), "slack": float(np.mean(T[late] == 0))})
        print(f"  U = {U}: {rows[-1]}")
    summary["fig04"] = rows
    U = np.array(speeds)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4))
    axes[0].plot(U, [r["gamma_amp"] for r in rows], "o-", color=COLORS[0])
    axes[0].set_ylabel("hawser-angle amplitude [deg]")
    axes[0].set_title("(a) Size of the swing", loc="left")
    axes[1].plot(U, [r["T_max_kN"] for r in rows], "o-", color=COLORS[2], label="peak tension")
    axes[1].plot(U, [r["F_H0_kN"] for r in rows], "s--", color=MUTED, label="steady tension (trailing straight)")
    axes[1].set_ylabel("tension [kN]")
    axes[1].set_title("(b) Peak rope load", loc="left")
    axes[1].legend(loc="upper left", fontsize=8.5)
    axes[2].plot(U, [100 * r["slack"] for r in rows], "o-", color=COLORS[1])
    axes[2].set_ylabel("time the rope is slack [%]")
    axes[2].set_title("(c) How often the rope goes slack", loc="left")
    for ax in axes:
        ax.set_xlabel("current speed U [m/s], full scale")
    fig.suptitle("Fully grown fishtailing on the real rope, 100 % load, $x_A/L$ = 0.5, 160 m hawser "
                 "(friction from the ITTC line at each speed)", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig04_vs_current_speed")


def hopf_elastic(s, lo=0.002, hi=0.15):
    f = lambda C: fishtail_pair(replace(s, C_F=C)).real
    grid = np.geomspace(lo, hi, 25)
    vals = [f(c) for c in grid]
    for a, b, fa, fb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
        if fa > 0 > fb:
            return brentq(f, a, b, xtol=1e-7)
    return np.nan


def fig05_boundary(summary):
    """Stability boundary (L_H/L, C_F): rigid vs real rope vs a 10x softer rope, at U = 1 m/s."""
    L = np.array([0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0])
    base = me.Elastic(HULL, 0.5, 0.5, 0.03, U=1.0)
    rigid = co.hopf_curve(mr.Rigid(HULL, 0.5, 0.5, 0.03), L)
    real = np.array([hopf_elastic(replace(base, L_H_L=float(x))) for x in L])
    soft = np.array([hopf_elastic(replace(base, L_H_L=float(x), rope_a=base.rope_a / 10)) for x in L])
    summary["fig05"] = {"L_H_L": list(L), "rigid": list(rigid), "real_rope": list(real), "soft_rope": list(soft)}
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    ax.plot(L, rigid, "o-", color=INK, lw=2.4, label="rigid hawser")
    ax.plot(L, real, "s-", color=COLORS[0], label="real nylon hawser")
    ax.plot(L, soft, "^-", color=COLORS[1], label="10 × softer rope")
    ax.set_xlabel("hawser length $L_H/L$")
    ax.set_ylabel("friction coefficient $C_F$ at the Hopf point")
    ax.set_ylim(0, None)
    ax.set_title("Stability boundary with a stretchy rope (stable above), 100 % load, $x_A/L$ = 0.5, U = 1 m/s\n"
                 "found from the eigenvalues of the 6-state model at trailing straight (rope taut there)", loc="left",
                 fontsize=10)
    ax.legend(loc="lower right")
    fig.tight_layout()
    save(fig, EXP, "fig05_stability_boundary")


def fig06_exponent(summary):
    """Does the stiffening exponent matter? Compare n = 1.86 with a linear rope of the same stiffness at rest."""
    s = full_scale(1.0)
    _, e0 = me.trailing_straight(s)
    # linear rope with the same tension AND the same slope at the trailing-straight stretch:
    # T = a n e0^(n-1) * (e - e_shift) is awkward; simplest fair comparison: n = 1 with the same tension at e0
    k_lin = s.rope_a * s.rope_n * e0 ** (s.rope_n - 1)  # slope dT/de at e0
    # linear law T = k_lin * e would give T(e0) = k_lin e0 = n T0 (too much), so match the tension instead and report
    lin_same_tension = replace(s, rope_a=s.rope_a * e0 ** (s.rope_n - 1), rope_n=1.0)
    fig, ax = plt.subplots(figsize=(11.0, 4.6))
    summary["fig06"] = {"slope_at_rest_kN": k_lin / 1e3}
    for color, (label, ss) in zip(COLORS, [(f"stiffening rope, n = {s.rope_n:.2f} (fitted)", s),
                                           ("linear rope, n = 1, same tension at rest", lin_same_tension)]):
        t, y, T, g = run(ss, t_end=400.0)
        late = t > 300
        ax.plot(t[late] - 300, T[late] / 1e3, color=color, label=label)
        summary["fig06"][label] = {"T_max_kN": float(T[late].max() / 1e3),
                                   "gamma_amp_deg": float(np.degrees(0.5 * (g[late].max() - g[late].min())))}
    ax.axhline(0, color=INK, lw=0.8)
    ax.set_xlabel("time t* (after the swing has fully grown)")
    ax.set_ylabel("hawser tension [kN]")
    ax.set_title("Rope law matters for the peak loads: real rope vs a linear rope, 1 m/s current, 100 % load, "
                 "160 m hawser", loc="left", fontsize=10)
    ax.legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    save(fig, EXP, "fig06_exponent_effect")


def main():
    use_style()
    summary = {}
    fig01_rope_law(summary)
    fig02_stiff_limit(summary)
    fig03_real_ship(summary)
    fig04_speed(summary)
    fig05_boundary(summary)
    fig06_exponent(summary)
    folder = ROOT / "data" / EXP
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / "summary.json", "w") as f:
        json.dump(summary, f, indent=1, default=float)
    print(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()
