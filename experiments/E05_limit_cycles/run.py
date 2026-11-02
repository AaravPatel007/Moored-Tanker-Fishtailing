"""E05 - How big is the fishtailing swing, and is the Hopf bifurcation super- or subcritical?

Run from the project root (takes a few minutes: every point is a long simulation):
    /home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E05_limit_cycles.run

Knob: friction coefficient C_F (swept down, then up), and hawser length L_H/L. Vidal de Negreiros 100 %, x_A/L = 0.5.
Writes figures to figures/E05_limit_cycles/ and numbers to data/E05_limit_cycles/.
Sign convention (Simos 2001): psi > 0 = bow turned to port, gamma > 0 = ship swung to port. Time t* = U t / L.
"Amplitude" = half of peak-to-peak, in degrees.
"""

import json
from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from fishtail import continuation as co
from fishtail import cycles as cy
from fishtail import model_rigid as mr
from fishtail.params import SIMOS_FIG57, SIMOS_FIG517, VIDAL_100, VIDAL_100_IPT
from fishtail.plotting import COLORS, INK, INK_2, MUTED, ROOT, save, use_style

EXP = "E05_limit_cycles"
BASE = mr.Rigid(VIDAL_100, x_A_L=0.5, L_H_L=0.5, C_F=0.03)
SMALL_KICK = np.array([0.01, 0.0, 0.0, 0.0])  # 0.6 deg hawser angle: "a tiny disturbance"
MODEL_L = 320.0 / 90.0  # length of the 1:90 tank model [m]


def sweep(s0, C_values, x_start):
    """Tank-test style sweep: for each C_F in turn, start from where the previous run ended and simulate until
    the swing stops changing. Returns a list of dicts (amplitude 0 if the ship settles)."""
    rows, x = [], np.array(x_start, float)
    for C in C_values:
        s = replace(s0, C_F=float(C))
        summ, x_end = cy.measure_cycle(s, x, t_transient=300.0)
        if summ is None:
            rows.append({"C_F": float(C), "gamma_amp_deg": 0.0, "psi_amp_deg": 0.0, "period": np.nan,
                         "F_H_max": float(C), "F_H_min": float(C)})
            x = x_end + SMALL_KICK  # a real tank always has small disturbances
        else:
            rows.append({"C_F": float(C), **{k: v for k, v in summ.items() if k != "x0"}})
            x = np.array(summ["x0"])
        print(f"  C_F = {C:.5f}: gamma amp {rows[-1]['gamma_amp_deg']:.2f} deg")
    return rows


def fig01_sweeps(summary, C_crit):
    """Amplitude against friction, sweeping down then up: no jump, no hysteresis -> supercritical Hopf."""
    fracs = np.array([1.06, 1.02, 0.98, 0.96, 0.94, 0.92, 0.9, 0.86, 0.82, 0.78, 0.74, 0.7, 0.65, 0.6, 0.55, 0.5,
                      0.45, 0.4, 0.35, 0.3])
    print("sweep down")
    down = sweep(BASE, fracs * C_crit, SMALL_KICK)
    print("sweep up")
    up = sweep(BASE, fracs[::-1] * C_crit, [np.radians(60.0), 0, 0, 0])  # start from a large swing
    summary["fig01"] = {"down": down, "up": up}

    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.6))
    for rows, marker, label, color in [(down, "v", "friction decreasing (sweep down)", COLORS[0]),
                                       (up, "^", "friction increasing (sweep up)", COLORS[1])]:
        C = np.array([r["C_F"] for r in rows])
        axes[0].plot(C, [r["gamma_amp_deg"] for r in rows], marker, color=color, ms=7, mfc="none", mew=1.4,
                     label=label)
        axes[1].plot(C, [r["psi_amp_deg"] for r in rows], marker, color=color, ms=7, mfc="none", mew=1.4,
                     label=label)
        axes[2].plot(C, np.array([r["F_H_max"] for r in rows]) / C, marker, color=color, ms=7, mfc="none", mew=1.4,
                     label=f"largest, {label.split(' (')[0]}")
        axes[2].plot(C, np.array([r["F_H_min"] for r in rows]) / C, marker, color=color, ms=5, mew=1.0,
                     label=f"smallest, {label.split(' (')[0]}")
    for ax in axes[:2]:
        ax.plot(C_crit, 0, "o", ms=9, color=INK, mec="white", mew=1.4, zorder=6)
        ax.set_ylim(bottom=-1)
    axes[2].axhline(0, color=INK_2, lw=1.0)
    axes[2].axhline(1, color=MUTED, lw=0.8)
    axes[0].set_ylabel("hawser-angle amplitude [deg]")
    axes[1].set_ylabel("heading amplitude [deg]")
    axes[2].set_ylabel("hawser force over one cycle  $F_H / F_{H0}$")
    axes[0].set_title("(a) Hawser angle γ", loc="left")
    axes[1].set_title("(b) Heading ψ", loc="left")
    axes[2].set_title("(c) Hawser force: peaks and dips", loc="left")
    for ax in axes:
        ax.set_xlabel("friction coefficient $C_F$")
    axes[0].legend(loc="upper right", fontsize=8.5)
    axes[2].legend(loc="upper right", fontsize=7.5)
    fig.suptitle(f"Fishtailing amplitude, 100 % load, $x_A/L$ = {BASE.x_A_L}, $L_H/L$ = {BASE.L_H_L}: both sweeps "
                 f"give the same curve, growing from zero at the Hopf point ● ($C_F$ = {C_crit:.4f}) - supercritical",
                 x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig01_amplitude_sweeps")


def shoot_branch(s0, fracs, C_crit):
    """Follow the cycle with the shooting method from far to near the Hopf point (precise even when the
    simulation would take very long to settle). Returns lists of (C_F, x0, T)."""
    s = replace(s0, C_F=fracs[0] * C_crit)
    summ, _ = cy.measure_cycle(s, [np.radians(30.0), 0, 0, 0])
    x0, T = np.array(summ["x0"]), summ["period"]
    out = []
    for f in fracs:
        s = replace(s0, C_F=f * C_crit)
        x0, T, ok = cy.shoot(s, x0, T)
        if not ok:
            raise RuntimeError(f"shooting failed at C_F/C_crit = {f}")
        out.append((s, x0.copy(), T))
    return out


def fig02_amplitude_law(summary, C_crit):
    """Amplitude against distance past the Hopf point: grows like the distance, not its square root - and why."""
    fracs = np.array([0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 0.98, 0.985, 0.99, 0.993, 0.995])
    branch = shoot_branch(BASE, fracs, C_crit)
    dist = 1 - fracs
    amp_g = np.array([cy.summarise(s, x0, T)["gamma_amp_deg"] for s, x0, T in branch])
    amp_p = np.array([cy.summarise(s, x0, T)["psi_amp_deg"] for s, x0, T in branch])
    slope_g = np.polyfit(np.log(dist[-6:]), np.log(amp_g[-6:]), 1)[0]
    slope_p = np.polyfit(np.log(dist[-6:]), np.log(amp_p[-6:]), 1)[0]
    summary["fig02"] = {"distance": list(dist), "gamma_amp": list(amp_g), "psi_amp": list(amp_p),
                        "loglog_slope_gamma": slope_g, "loglog_slope_psi": slope_p}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.8))
    ax1.loglog(dist, amp_g, "o", color=COLORS[0], ms=6, label=f"hawser angle γ (slope {slope_g:.2f})")
    ax1.loglog(dist, amp_p, "s", color=COLORS[1], ms=5, label=f"heading ψ (slope {slope_p:.2f})")
    ax1.loglog(dist, amp_g[-1] * dist / dist[-1], color=COLORS[0], lw=1.0)
    ax1.loglog(dist, amp_p[-1] * dist / dist[-1], color=COLORS[1], lw=1.0)
    ax1.loglog(dist, amp_g[-1] * np.sqrt(dist / dist[-1]), color=MUTED, lw=1.0, ls="--",
               label="slope ½ (Strogatz's square-root law)")
    ax1.set_xlabel("distance past the Hopf point  $(C_{F,Hopf} - C_F)/C_{F,Hopf}$")
    ax1.set_ylabel("amplitude [deg]")
    ax1.set_title("(a) Ship: amplitude ∝ distance (log-log slope 1)", loc="left")
    ax1.legend(loc="upper left", fontsize=8.5)

    # Toy model (Strogatz section 7.6, averaging): x'' - mu x' + D(x') + x = 0 with two kinds of damping D.
    mus = np.geomspace(0.003, 0.3, 9)
    for color, name, D, theory in [
        (COLORS[2], "quadratic damping  $\\dot x|\\dot x|$", lambda v: v * abs(v), lambda m: 3 * np.pi * m / 8),
        (COLORS[3], "cubic damping  $\\dot x^3$", lambda v: v**3, lambda m: np.sqrt(4 * m / 3)),
    ]:
        amps = []
        for m in mus:
            sol = solve_ivp(lambda t, z: [z[1], m * z[1] - D(z[1]) - z[0]], (0, 40 / m + 200), [0.5 * theory(m), 0],
                            rtol=1e-9, atol=1e-12, dense_output=True)
            t_late = np.linspace(sol.t[-1] - 20, sol.t[-1], 2000)
            amps.append(np.abs(sol.sol(t_late)[0]).max())
        ax2.loglog(mus, amps, "o", color=color, ms=6, label=f"{name}: simulated")
        ax2.loglog(mus, theory(mus), color=color, lw=1.2, label="  averaging theory")
    ax2.set_xlabel("distance past the Hopf point  μ")
    ax2.set_ylabel("limit-cycle amplitude")
    ax2.set_title("(b) Why: toy oscillator  $\\ddot x - μ\\dot x + D(\\dot x) + x = 0$", loc="left")
    ax2.legend(loc="upper left", fontsize=8.5)
    fig.suptitle("Strogatz's √-law needs a smooth (cubic) nonlinearity. Our hull drag is quadratic (v|v|, r|r|), "
                 "so averaging gives dr/dt = μr/2 − c r² and r ∝ μ", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig02_amplitude_law")
    return branch


def fig03_floquet(summary, C_crit):
    """Floquet multipliers of the cycle: one is 1, the others inside the unit circle -> the cycle is stable."""
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(th), np.sin(th), color=MUTED, lw=1.0)
    summary["fig03"] = {}
    x0 = np.array([np.radians(10.0), 0, 0, 0])
    for color, frac in zip(COLORS, [0.95, 0.8, 0.6, 0.4]):
        s = replace(BASE, C_F=frac * C_crit)
        summ, _ = cy.measure_cycle(s, x0)  # simulate to the cycle, then polish it by shooting
        x0, T, ok = cy.shoot(s, np.array(summ["x0"]), summ["period"])
        mu = cy.floquet_multipliers(s, x0, T)
        ax.plot(mu.real, mu.imag, "o", color=color, ms=9, mec="white", mew=1.3,
                label=f"$C_F$ = {frac} × Hopf value: |μ| = " + ", ".join(f"{abs(m):.3f}" for m in mu))
        summary["fig03"][str(frac)] = {"multipliers_abs": [float(abs(m)) for m in mu], "period": float(T),
                                       "shooting_converged": ok}
    ax.set_aspect("equal")
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)
    ax.set_xlabel("real part of multiplier μ")
    ax.set_ylabel("imaginary part")
    ax.set_title("Floquet multipliers of the fishtailing cycle\n(circle = |μ| = 1; inside = disturbance shrinks "
                 "each period)", loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=7.5)  # below, clear of the circle
    fig.tight_layout()
    save(fig, EXP, "fig03_floquet_multipliers")


def model_period_seconds(T_star, U):
    return T_star * MODEL_L / U


def fig04_period_vs_simos(summary):
    """Fishtailing period at tank-model scale against Simos's simulated values (Fig. 5.7)."""
    U = SIMOS_FIG57["U_model"]
    C_values = np.linspace(0.0220, 0.0250, 7)
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    summary["fig04"] = {}
    for color, key in zip(COLORS, ["0.36", "0.6", "0.9"]):
        s0 = mr.Rigid(VIDAL_100_IPT, 0.5, float(key), 0.023)
        x = [np.radians(30.0), 0, 0, 0]
        periods = []
        for C in C_values:
            summ, _ = cy.measure_cycle(replace(s0, C_F=float(C)), x)
            x = summ["x0"]
            periods.append(model_period_seconds(summ["period"], U))
        ax.plot(C_values, periods, color=color, label=f"ours, $L_H/L$ = {key}")
        ax.plot(SIMOS_FIG57["C_F"], SIMOS_FIG57["period_s"][key], "o", color=color, ms=8, mfc="white", mew=1.6,
                label=f"Simos Fig. 5.7, $L_H/L$ = {key}")
        summary["fig04"][key] = {"C_F": list(C_values), "period_s": periods}
    ax.set_xlabel("friction coefficient $C_F$")
    ax.set_ylabel("fishtailing period [s], 1:90 model at 0.40 m/s")
    ax.set_title("Period of the fully grown swing vs Simos's simulations (Fig. 5.7, read off the gridlines ±1 s)\n"
                 "100 % load, $x_A/L$ = 0.5, $C_Y$ = 0.86, $l_pC_Y$ = 0.044 (as Simos used)", loc="left")
    ax.legend(loc="upper right", fontsize=8, ncol=2)
    fig.tight_layout()
    save(fig, EXP, "fig04_period_vs_simos")


def fig05_amplitude_vs_simos(summary):
    """Amplitude against hawser length at tank-model scale, against Simos's theory and experiments (Fig. 5.17)."""
    L_values = np.array([0.1, 0.12, 0.14, 0.17, 0.2, 0.25, 0.3, 0.36, 0.45, 0.55, 0.65, 0.8, 0.95, 1.1, 1.25, 1.35])
    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    summary["fig05"] = {}
    for color, C, label in [(COLORS[0], SIMOS_FIG517["C_F_min"], "C_F min"),
                            (COLORS[1], SIMOS_FIG517["C_F_max"], "C_F max")]:
        g_amp, p_amp = [], []
        x = [np.radians(30.0), 0, 0, 0]
        for LH in L_values[::-1]:  # from long to short hawser, continuing from the previous cycle
            s = mr.Rigid(VIDAL_100_IPT, 0.5, float(LH), C)
            summ, x_end = cy.measure_cycle(s, x)
            if summ is None:
                g_amp.append(0.0)
                p_amp.append(0.0)
                x = x_end + SMALL_KICK
            else:
                g_amp.append(summ["gamma_amp_deg"])
                p_amp.append(summ["psi_amp_deg"])
                x = summ["x0"]
        g_amp, p_amp = g_amp[::-1], p_amp[::-1]
        ax.plot(L_values, g_amp, color=color, lw=2.0, label=f"ours γ, $C_F$ = {C} ({label})")
        ax.plot(L_values, p_amp, color=color, lw=2.0, ls=(0, (1, 1.5)), label=f"ours ψ, $C_F$ = {C}")
        summary["fig05"][str(C)] = {"L_H_L": list(L_values), "gamma_amp": g_amp, "psi_amp": p_amp}
    th = SIMOS_FIG517
    ax.plot(th["theory_L_H_L"], th["theory_gamma_CFmin_deg"], "x", color=COLORS[0], ms=7, mew=1.4,
            label="Simos theory γ, $C_F$ min (by eye)")
    ax.plot(th["theory_L_H_L"], th["theory_gamma_CFmax_deg"], "x", color=COLORS[1], ms=7, mew=1.4,
            label="Simos theory γ, $C_F$ max (by eye)")
    ax.plot(th["theory_L_H_L"], th["theory_psi_deg"], "x", color=INK_2, ms=7, mew=1.4, label="Simos theory ψ (by eye)")
    ax.plot(th["exp_L_H_L"], th["exp_gamma_deg"], "P", color=INK, ms=8, mfc="white", mew=1.4,
            label="experiment γ (Simos Fig. 5.17)")
    ax.plot(th["exp_L_H_L"], th["exp_psi_deg"], "o", color=INK, ms=7, mfc="white", mew=1.4,
            label="experiment ψ (Simos Fig. 5.17)")
    ax.set_xlabel("hawser length  $L_H/L$")
    ax.set_ylabel("amplitude [deg]")
    ax.set_xlim(0, 1.4)
    ax.set_ylim(0, 70)
    ax.set_title("Fishtailing amplitude vs hawser length: 1:90 model at 0.4 m/s, 100 % load, $x_A/L$ = 0.5\n"
                 "ours (lines) vs Simos's theory and tank tests (Fig. 5.17, read by eye ±2°)", loc="left")
    ax.legend(loc="upper right", fontsize=7.5, ncol=2)
    fig.tight_layout()
    save(fig, EXP, "fig05_amplitude_vs_simos")


def fig06_cycle_shape(summary):
    """Shape of one cycle: the hawser angle swings like a sine, the heading like a saw-tooth (Simos p. 209-210)."""
    s = mr.Rigid(VIDAL_100_IPT, 0.5, 0.6, 0.0223)
    summ, _ = cy.measure_cycle(s, [np.radians(30.0), 0, 0, 0])
    t, y = cy.one_period(s, np.array(summ["x0"]), summ["period"])
    t2 = np.concatenate([t, t[1:] + t[-1]])
    y2 = np.concatenate([y, y[:, 1:]], axis=1)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.6), gridspec_kw={"width_ratios": [1.6, 1]})
    ax1.plot(t2 / summ["period"], np.degrees(y2[0]) / summ["gamma_amp_deg"], color=COLORS[0],
             label=f"hawser angle γ ÷ {summ['gamma_amp_deg']:.0f}°")
    ax1.plot(t2 / summ["period"], np.degrees(y2[1]) / summ["psi_amp_deg"], color=COLORS[1],
             label=f"heading ψ ÷ {summ['psi_amp_deg']:.0f}°")
    ax1.plot(t2 / summ["period"], np.cos(2 * np.pi * t2 / summ["period"]), color=MUTED, lw=1.0, ls="--",
             label="pure sine, for comparison")
    ax1.set_xlabel("time ÷ period")
    ax1.set_ylabel("angle ÷ its amplitude")
    ax1.set_title("(a) Two cycles: γ is almost a pure sine; ψ is not (sharp peaks, a shoulder between)", loc="left")
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=8.5)  # below, clear of curves
    ax2.plot(np.degrees(y[0]), np.degrees(y[1]), color=INK, lw=1.6)
    ax2.axhline(0, color=MUTED, lw=0.8)
    ax2.axvline(0, color=MUTED, lw=0.8)
    ax2.set_xlabel("hawser angle γ [deg]")
    ax2.set_ylabel("heading ψ [deg]")
    ax2.set_title("(b) The cycle in the (γ, ψ) plane", loc="left")
    summary["fig06"] = {k: v for k, v in summ.items() if k != "x0"}
    fig.suptitle(f"One fishtailing cycle, 1:90 model at 0.4 m/s ($C_F$ = 0.0223), $L_H/L$ = 0.6, $x_A/L$ = 0.5: "
                 f"period {model_period_seconds(summ['period'], 0.4):.0f} s at model scale", x=0.01, ha="left",
                 fontsize=11)
    fig.tight_layout()
    save(fig, EXP, "fig06_cycle_shape")


def main():
    use_style()
    C_crit = co.hopf_C_F(BASE)
    summary = {"C_F_hopf": C_crit}
    fig01_sweeps(summary, C_crit)
    fig02_amplitude_law(summary, C_crit)
    fig03_floquet(summary, C_crit)
    fig04_period_vs_simos(summary)
    fig05_amplitude_vs_simos(summary)
    fig06_cycle_shape(summary)
    folder = ROOT / "data" / EXP
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / "summary.json", "w") as f:
        json.dump(summary, f, indent=1, default=float)
    print(json.dumps({k: v for k, v in summary.items() if k != "fig01"}, indent=1, default=float))


if __name__ == "__main__":
    main()
