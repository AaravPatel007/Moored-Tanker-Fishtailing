"""E08 - Can we stop the fishtailing? (a) a hold-back tug at the stern, (b) the tanker's own stern thruster with
PD heading control and a thrust limit. Decision D31.

Run from the project root:
    /home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E08_control.run

Base case: rigid hawser, Vidal de Negreiros 100 %, bow attachment x_A/L = 0.5, L_H/L = 0.5, C_F = 0.03 (the
Fishtailing demo default; like a 1:90 tank model at 0.17 m/s): without control it fishtails.
  Tug   : a steady pull F_tug downstream at the stern (the wind force with beta_W = 0, x_W = stern).
  Thruster: side force at the stern, N_c* = -k_p psi - k_d r*, clipped at thrust_limit * F_H0 (fishtail/control.py).
All thrust limits, tug pulls and the stern positions are PLACEHOLDERs (no cited values yet, D31).
"Optimal" (agreed with the student): the gains with the smallest thrust needed that still give a decay rate
sigma = -max Re(lambda) >= SIGMA_MIN. Thrust needed = peak |Y_T| / F_H0 when the controller is switched on during
the full uncontrolled fishtail (worst of 4 phases of the swing), measured without a limit.
Writes figures to figures/E08_control/ and numbers to data/E08_control/summary.json.
Sign convention (Simos 2001): psi > 0 = bow turned to port, gamma > 0 = ship swung to port. Time t* = U t / L.
"""

import json
from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np

from fishtail import continuation as co
from fishtail import control
from fishtail import cycles as cy
from fishtail import hydro
from fishtail import model_rigid as mr
from fishtail.params import SEA_WATER, VIDAL_100, X_TUG_L
from fishtail.plotting import COLORS, INK, INK_2, MUTED, ROOT, plot_branch, save, use_style

EXP = "E08_control"
HULL = VIDAL_100
BASE = mr.Rigid(HULL, x_A_L=0.5, L_H_L=0.5, C_F=0.03)
SIGMA_MIN = 0.05  # required decay rate [U/L]: a small swing halves within ln 2 / 0.05 = 14 t*. PROPOSED, to agree
KP = np.round(np.linspace(0.0, 0.04, 21), 4)  # proportional gains on the map
KD = np.round(np.linspace(0.0, 0.10, 26), 4)  # derivative gains on the map
NUDGE = [np.radians(3.0), np.radians(1.0), 0.0, 0.0]


def tug(s, ratio):
    """Set-up s with a hold-back tug pulling ratio * F_H0 downstream from the stern."""
    return replace(s, C_W=ratio * s.C_F, beta_W=0.0, x_W_L=X_TUG_L)


def sigma(s):
    """Decay rate of the straight position [U/L] (exact linearisation)."""
    return control.decay_rate(co.eigenvalues(np.zeros(4), s))


def base_cycle(s=BASE, n_phases=4):
    """The fully grown uncontrolled fishtail of s: summary and n_phases states spread over one period."""
    summ, _ = cy.measure_cycle(s, [np.radians(30.0), 0, 0, 0])
    t, y = cy.one_period(s, np.array(summ["x0"]), summ["period"], n=4 * n_phases + 1)
    return summ, [y[:, k] for k in range(0, 4 * n_phases, 4)]


def end_amplitude(s, start, t_end=600.0):
    """Half the peak-to-peak hawser angle over the last 50 t* of a run [deg] (0 = the swing died out)."""
    t, y = mr.simulate(s, start, t_end, n_out=int(3 * t_end) + 1, rtol=1e-8, atol=1e-10)
    return float(np.degrees(np.ptp(y[0, -150:])) / 2)


def peak_thrust(s, starts, t_end=150.0):
    """Largest |Y_T| / F_H0 over runs from each start, with no thrust limit (the thrust the controller asks for)."""
    s = replace(s, thrust_limit=np.inf)
    peak = 0.0
    for x0 in starts:
        t, y = mr.simulate(s, x0, t_end, n_out=1501, rtol=1e-8, atol=1e-10)
        peak = max(peak, float(np.abs(control.thruster_ratio(y, s)).max()))
    return peak


def sigma_map(s, kp=KP, kd=KD):
    return np.array([[sigma(replace(s, k_p=float(p), k_d=float(d))) for d in kd] for p in kp])


# --- fig01: the hold-back tug -------------------------------------------------------------------------------------

def fig01(summary):
    ratios = np.round(np.linspace(0.0, 1.0, 41), 4)
    stable = [co.is_stable(co.eigenvalues(np.zeros(4), tug(BASE, r))) for r in ratios]
    crit = co.find_crossings(lambda r: sigma(tug(BASE, r)), [0.0, 0.25, 0.5, 1.0])[0]
    amps = []
    x = [np.radians(30.0), 0, 0, 0]
    for r in ratios:  # from no tug upwards: each run starts from the previous swing
        summ, x_end = cy.measure_cycle(tug(BASE, r), x, t_transient=300.0)
        amps.append(0.0 if summ is None else summ["gamma_amp_deg"])
        x = x_end if summ is None else summ["x0"]
    amps = np.array(amps)
    C_Fs = np.round(np.linspace(0.006, 0.06, 28), 5)
    crit_C = []
    for c in C_Fs:
        f = lambda r: sigma(tug(replace(BASE, C_F=float(c)), r))
        roots = co.find_crossings(f, np.linspace(0.0, 3.0, 31))
        crit_C.append(roots[0] if roots else (0.0 if f(0.0) > 0 else np.nan))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4))
    plot_branch(a1, ratios, np.zeros(ratios.size), stable, INK, label="trailing straight")
    has = amps > 0
    xs = np.concatenate([ratios[has], [crit]])
    ys = np.concatenate([amps[has], [0.0]])
    a1.plot(xs, ys, color=COLORS[0], label="swing amplitude ±")
    a1.plot(xs, -ys, color=COLORS[0])
    a1.plot(crit, 0, "o", color=INK, ms=7, zorder=5, label=f"Hopf: F_tug = {crit:.3f} F_H0")
    a1.set_xlabel("tug pull F_tug / F_H0 (downstream, at the stern)")
    a1.set_ylabel("hawser angle γ: rest value or swing range [deg]")
    a1.set_title("Hold-back tug: the swing shrinks to zero\n(C_F = 0.03, x_A/L = 0.5, L_H/L = 0.5)")
    a1.legend()
    a2.plot(C_Fs, crit_C, "o-", color=COLORS[1], ms=3.5)
    a2.fill_between(C_Fs, 0, np.nan_to_num(crit_C), color=COLORS[1], alpha=0.12, lw=0)
    a2.text(0.012, 0.15, "fishtails", color=INK_2)
    a2.text(0.035, 1.0, "trails straight", color=INK_2)
    a2.axvline(co.hopf_C_F(BASE), color=MUTED, ls=":", lw=1)
    a2.set_xlabel("friction coefficient C_F (faster current → smaller)")
    a2.set_ylabel("tug pull needed to stop fishtailing, F_tug / F_H0")
    a2.set_title("Tug pull needed against friction\n(dotted: Hopf point without a tug)")
    fig.tight_layout()
    save(fig, EXP, "fig01_tug")
    summary["tug"] = {"hopf_ratio_at_base": round(crit, 4), "C_F": C_Fs.tolist(),
                      "ratio_needed": [None if np.isnan(v) else round(float(v), 4) for v in crit_C],
                      "amp_deg": [round(float(a), 3) for a in amps], "ratios": ratios.tolist()}


# --- fig02: root locus of the thruster gains -------------------------------------------------------------------------

def fig02(summary):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    for ax, key, vals, other in [(axes[0], "k_d", np.linspace(0, 0.3, 301), "k_p = 0"),
                                 (axes[1], "k_p", np.linspace(0, 0.1, 301), "k_d = 0")]:
        eigs = np.array([co.eigenvalues(np.zeros(4), replace(BASE, **{key: float(v)})) for v in vals])
        sc = None
        for j in range(4):
            sc = ax.scatter(eigs[:, j].real, eigs[:, j].imag, c=vals, s=4, cmap="viridis")
        ax.plot(eigs[0].real, eigs[0].imag, "x", color=INK, ms=8, mew=2, label="no control")
        ax.axvline(0, color=INK, lw=0.8)
        ax.axvspan(0, 0.2, color=COLORS[7], alpha=0.06, lw=0)
        ax.set_xlim(-1.2, 0.15)
        ax.set_xlabel("Re λ [U/L]")
        ax.set_title(f"Eigenvalues of trailing straight as {key} grows ({other})")
        fig.colorbar(sc, ax=ax, label=key)
        ax.legend(loc="lower left")
        crit = co.find_crossings(lambda v: sigma(replace(BASE, **{key: v})), list(vals[::10]))
        summary[f"hopf_{key}"] = round(crit[0], 5) if crit else None
    axes[0].set_ylabel("Im λ [U/L]")
    axes[0].annotate("this eigenvalue does not move:\nheading control cannot reach it",
                     (-0.2205, 0), xytext=(-1.1, -0.9), color=INK_2, fontsize=8.5,
                     arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    fig.tight_layout()
    save(fig, EXP, "fig02_root_locus")
    best = max(((sigma(replace(BASE, k_d=float(v))), float(v)) for v in np.linspace(0, 0.3, 301)))
    summary["fastest_decay_kp0"] = {"sigma": round(best[0], 4), "k_d": round(best[1], 4)}


# --- fig03 + fig04: the (k_p, k_d) maps and the least-thrust point ★ --------------------------------------------------

def fig03_04(summary, starts):
    S = sigma_map(BASE)
    P = np.array([[peak_thrust(replace(BASE, k_p=float(p), k_d=float(d)), starts) if S[i, j] > 0 else np.nan
                   for j, d in enumerate(KD)] for i, p in enumerate(KP)])
    ok = S >= SIGMA_MIN
    i, j = np.unravel_index(np.nanargmin(np.where(ok, P, np.nan)), P.shape)
    star = (float(KP[i]), float(KD[j]))
    fine_p, fine_d = np.linspace(0, 0.04, 81), np.linspace(0, 0.1, 101)
    SF = sigma_map(BASE, fine_p, fine_d)
    fig, ax = plt.subplots(figsize=(6.6, 5))
    m = ax.pcolormesh(fine_d, fine_p, SF, cmap="viridis", shading="auto")
    fig.colorbar(m, ax=ax, label="decay rate σ = −max Re λ [U/L]  (< 0: fishtails)")
    cs = ax.contour(fine_d, fine_p, SF, levels=[0.0, 0.02, SIGMA_MIN, 0.1, 0.15], colors=["white"] + [INK] * 4,
                    linewidths=[2.2, 0.8, 1.6, 0.8, 0.8])
    ax.clabel(cs, fmt=lambda v: "Hopf σ = 0" if v == 0 else f"σ = {v:g}", fontsize=8)
    ax.plot(star[1], star[0], "*", color="white", mec=INK, ms=16, label="★ least thrust with σ ≥ %.2f" % SIGMA_MIN)
    ax.set_xlabel("derivative gain k_d (damper on turning)")
    ax.set_ylabel("proportional gain k_p (spring to straight)")
    ax.set_title("Stern thruster: stability of trailing straight\n(C_F = 0.03, x_A/L = 0.5, L_H/L = 0.5)")
    ax.text(0.075, 0.02, "σ capped at 0.22:\none eigenvalue\ndoes not move", ha="center", fontsize=8.5, color=INK)
    ax.legend(loc="upper right", markerscale=0.6, framealpha=0.9, frameon=True)
    fig.tight_layout()
    save(fig, EXP, "fig03_gain_map")
    fig, ax = plt.subplots(figsize=(6.6, 5))
    m = ax.pcolormesh(KD, KP, P, cmap="magma_r", shading="auto")
    fig.colorbar(m, ax=ax, label="thrust needed: peak |Y_T| / F_H0 (no limit)")
    ax.contour(fine_d, fine_p, SF, levels=[0.0], colors=[INK], linewidths=2.0)
    ax.contour(fine_d, fine_p, SF, levels=[SIGMA_MIN], colors=[INK], linewidths=1.2, linestyles="--")
    ax.plot(star[1], star[0], "*", color="white", mec=INK, ms=16)
    ax.annotate(f"★ k_p = {star[0]:g}, k_d = {star[1]:g}\nthrust {P[i, j]:.3f} F_H0", (star[1], star[0]),
                xytext=(20, 25), textcoords="offset points", fontsize=8.5)
    ax.set_xlabel("derivative gain k_d")
    ax.set_ylabel("proportional gain k_p")
    ax.set_title("Thrust needed to stop the full fishtail\n(solid: Hopf σ = 0, dashed: σ = %.2f; blank: fishtails)"
                 % SIGMA_MIN)
    fig.tight_layout()
    save(fig, EXP, "fig04_thrust_needed")
    summary["star"] = {"k_p": star[0], "k_d": star[1], "sigma": round(float(S[i, j]), 4),
                       "thrust_needed": round(float(P[i, j]), 4), "sigma_min": SIGMA_MIN}
    summary["sigma_max_on_map"] = round(float(np.nanmax(SF)), 4)
    return star


# --- fig05: time runs, controller switched on during the full fishtail ----------------------------------------------

def fig05(star, x0):
    cases = [("no control", BASE), ("★ gains", replace(BASE, k_p=star[0], k_d=star[1])),
             ("4 × ★ gains", replace(BASE, k_p=4 * star[0], k_d=4 * star[1]))]
    fig, axes = plt.subplots(3, 1, figsize=(9, 7.5), sharex=True)
    for c, (name, s) in zip(COLORS, cases):
        t, y = mr.simulate(s, x0, 120.0, n_out=1201, rtol=1e-8, atol=1e-10)
        axes[0].plot(t, np.degrees(y[0]), color=c, label=name)
        axes[1].plot(t, np.degrees(y[1]), color=c)
        axes[2].plot(t, control.thruster_ratio(y, s) if control.is_on(s) else 0 * t, color=c)
    axes[0].set_ylabel("hawser angle γ [deg]")
    axes[1].set_ylabel("heading ψ [deg]")
    axes[2].set_ylabel("thruster Y_T / F_H0")
    axes[2].set_xlabel("time t* = U t / L")
    axes[0].set_title("Controller switched on at t* = 0 during the full fishtail (no thrust limit)")
    axes[0].legend(ncol=3)
    fig.tight_layout()
    save(fig, EXP, "fig05_time_runs")


# --- fig06: saturation -> bistability ---------------------------------------------------------------------------------

def fig06(summary, star, x0):
    limits = np.round(np.geomspace(0.002, 0.2, 25), 5)
    s_star = replace(BASE, k_p=star[0], k_d=star[1])
    big = [end_amplitude(replace(s_star, thrust_limit=float(l)), x0) for l in limits]
    small = [end_amplitude(replace(s_star, thrust_limit=float(l)), NUDGE) for l in limits]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    a1.semilogx(limits, big, "o-", color=COLORS[1], ms=4, label="start: full fishtail")
    a1.semilogx(limits, small, "s--", color=COLORS[0], ms=4, label="start: small nudge (3°, 1°)")
    a1.set_xlabel("thrust limit / F_H0")
    a1.set_ylabel("swing left after 600 t*, γ amplitude [deg]")
    a1.set_title("Same gains (★), same limit, different start → different end\n(where the two curves differ: "
                 "bistable)")
    a1.legend()
    for c, (lim, x) in zip(COLORS, [(0.01, x0), (0.01, NUDGE), (0.02, x0)]):
        s = replace(s_star, thrust_limit=lim)
        t, y = mr.simulate(s, x, 200.0, n_out=2001, rtol=1e-8, atol=1e-10)
        a2.plot(t, np.degrees(y[0]), color=c,
                label=f"limit {lim} F_H0, " + ("full fishtail" if x is x0 else "nudge"))
    a2.set_xlabel("time t* = U t / L")
    a2.set_ylabel("hawser angle γ [deg]")
    a2.set_title("Thruster at its limit: a big swing survives,\na small one dies out")
    a2.legend(fontsize=8)
    fig.tight_layout()
    save(fig, EXP, "fig06_saturation")
    bist = [float(l) for l, b, sm in zip(limits, big, small) if b > 0.5 and sm < 0.05]
    summary["saturation"] = {"limits": limits.tolist(), "amp_from_cycle": [round(v, 3) for v in big],
                             "amp_from_nudge": [round(v, 3) for v in small],
                             "bistable_limits": bist}


# --- fig07: the (L_H/L, C_F) stability map with and without control ---------------------------------------------------

def fig07(summary, star):
    L_H = np.linspace(0.1, 1.5, 29)
    cases = [("no control", BASE), ("tug 0.5 F_H0", tug(BASE, 0.5)),
             ("stern thruster, ★ gains", replace(BASE, k_p=star[0], k_d=star[1]))]
    fig, ax = plt.subplots(figsize=(7, 5))
    out = {}
    for c, (name, s) in zip(COLORS, cases):
        if s.C_W != 0.0:  # the tug pull is a fraction of F_H0, so it changes with C_F: keep the ratio fixed
            def h(LH):
                g = lambda cf: sigma(tug(replace(BASE, L_H_L=LH, C_F=cf), 0.5))
                roots = co.find_crossings(g, list(np.geomspace(1e-3, 0.2, 40)))
                return roots[0] if roots else np.nan
            curve = np.array([h(LH) for LH in L_H])
        else:
            curve = co.hopf_curve(s, L_H)
        ax.plot(L_H, curve, color=c, label=name)
        none = np.isnan(curve)  # no Hopf point: checked, then stable for every C_F down to 1e-4
        if none.any():
            ax.plot(L_H[none], np.zeros(none.sum()) + 0.0008, "x", color=c, ms=6, mew=1.6)
            ax.annotate("×: no fishtailing at any C_F ≥ 0.0001", (L_H[none][0], 0.0008), xytext=(10, 18),
                        textcoords="offset points", color=INK_2, fontsize=8.5)
        out[name] = [None if np.isnan(v) else round(float(v), 5) for v in curve]
    for U in (0.5, 1.0, 2.0):
        cf = hydro.friction_coefficient(U, HULL, SEA_WATER)
        ax.axhline(cf, color=MUTED, lw=0.8, ls=":")
        ax.text(1.5, cf, f" real ship {U} m/s", va="center", fontsize=8, color=INK_2)
    ax.set_xlabel("hawser length L_H / L")
    ax.set_ylabel("Hopf friction coefficient C_F (fishtails below the line)")
    ax.set_title("Fishtailing boundary with and without control\n(x_A/L = 0.5; gains and tug ratio fixed in "
                 "non-dimensional form)")
    ax.set_ylim(0, None)
    ax.legend(loc="center right")
    fig.tight_layout()
    save(fig, EXP, "fig07_stability_map")
    summary["fig07"] = {"L_H_L": L_H.round(4).tolist(), **out}


# --- fig08: stern vs bow thruster ---------------------------------------------------------------------------------------

def fig08(summary):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    out = {}
    for ax, (name, x_T) in zip(axes, [("stern, x_T/L = −0.5", -0.5), ("bow, x_T/L = +0.45", 0.45)]):
        S = sigma_map(replace(BASE, x_T_L=x_T))
        m = ax.pcolormesh(KD, KP, S, cmap="viridis", shading="auto", vmin=-0.3, vmax=0.22)
        ax.contour(KD, KP, S, levels=[0.0], colors=["white"], linewidths=2)
        ax.set_title(f"Thruster at the {name}")
        ax.set_xlabel("derivative gain k_d")
        out[name] = {"stable_fraction": round(float(np.mean(S > 0)), 3), "sigma_max": round(float(S.max()), 4)}
    axes[0].set_ylabel("proportional gain k_p")
    fig.colorbar(m, ax=axes, label="decay rate σ [U/L] (white line: Hopf σ = 0)")
    fig.suptitle("Same controller, same yaw moment: why the thruster goes at the stern (hawser at the bow)")
    save(fig, EXP, "fig08_bow_vs_stern")
    summary["bow_vs_stern"] = out


# --- fig09: the real ship, fixed thruster / tug in kN, against current speed -----------------------------------------

def fig09(summary):
    """Real ship (sea water): tune the ★-type gains once at 1 m/s, then keep the PHYSICAL gains (kN m per rad,
    kN m s per rad) and the physical tug pull (kN) fixed while the current speed changes."""
    U = np.linspace(0.25, 2.0, 36)
    U0 = 1.0
    real = lambda u, **kw: replace(BASE, U=float(u), C_F=float(hydro.friction_coefficient(u, HULL, SEA_WATER)), **kw)
    s0 = real(U0)
    S = sigma_map(s0, KP * 4, KD * 4)  # the real ship fishtails harder (lower C_F): larger gain range
    kp_grid, kd_grid = KP * 4, KD * 4
    i, j = np.unravel_index(np.argmax(np.where(S >= SIGMA_MIN, -(kp_grid[:, None] + kd_grid[None, :]), -np.inf)),
                            S.shape)  # smallest gains meeting the margin (thrust grows with the gains)
    kp0, kd0 = float(kp_grid[i]), float(kd_grid[j])
    K0 = control.force_scale(s0) * HULL.L  # (1/2) rho U0^2 L^2 T
    Kp_phys, Kd_phys = kp0 * K0, kd0 * K0 * HULL.L / U0  # N m per rad,  N m s per rad
    r_tug = co.find_crossings(lambda r: sigma(tug(s0, r)), list(np.linspace(0, 3, 31)))
    r_tug = (r_tug[0] if r_tug else 3.0) * 1.2  # 20 % above the Hopf pull at 1 m/s
    F_tug = r_tug * s0.C_F * control.force_scale(s0)  # N
    sig_c, sig_t, sig_0 = [], [], []
    for u in U:
        su = real(u)
        Ku = control.force_scale(su) * HULL.L
        sig_c.append(sigma(replace(su, k_p=Kp_phys / Ku, k_d=Kd_phys / (Ku * HULL.L / u))))
        sig_t.append(sigma(tug(su, F_tug / (su.C_F * control.force_scale(su)))))
        sig_0.append(sigma(su))
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.plot(U, sig_0, color=MUTED, label="no control")
    ax.plot(U, sig_c, color=COLORS[0], label=f"thruster, gains fixed in kN·m (tuned at {U0} m/s)")
    ax.plot(U, sig_t, color=COLORS[1], label=f"tug, pull fixed at {F_tug / 1e3:.0f} kN")
    ax.axhline(0, color=INK, lw=0.8)
    ax.axvline(U0, color=MUTED, ls=":", lw=0.8)
    ax.set_xlabel("current speed U [m/s] (real ship; enters through C_F(Re) and through the fixed kN settings)")
    ax.set_ylabel("decay rate σ [U/L]  (< 0: fishtails)")
    ax.set_title("Real ship: a controller set up for one current, used at others")
    ax.legend(fontsize=8.5)
    fig.tight_layout()
    save(fig, EXP, "fig09_real_ship_speed")
    summary["real_ship"] = {"U0": U0, "C_F_U0": round(s0.C_F, 5), "k_p": kp0, "k_d": kd0,
                            "Kp_kNm_per_rad": round(Kp_phys / 1e3), "Kd_kNms_per_rad": round(Kd_phys / 1e3),
                            "tug_kN": round(F_tug / 1e3, 1), "tug_ratio_U0": round(r_tug, 3),
                            "F_H0_kN_U0": round(s0.C_F * control.force_scale(s0) / 1e3, 1),
                            "U": U.round(3).tolist(), "sigma_thruster": [round(v, 4) for v in sig_c],
                            "sigma_tug": [round(v, 4) for v in sig_t], "sigma_none": [round(v, 4) for v in sig_0]}


def main():
    use_style()
    summ, starts = base_cycle()
    summary = {"base": {"x_A_L": BASE.x_A_L, "L_H_L": BASE.L_H_L, "C_F": BASE.C_F,
                        "cycle_gamma_amp_deg": round(summ["gamma_amp_deg"], 2),
                        "cycle_psi_amp_deg": round(summ["psi_amp_deg"], 2), "period": round(summ["period"], 3),
                        "sigma": round(sigma(BASE), 4)}}
    x0 = np.array(summ["x0"])
    fig01(summary)
    fig02(summary)
    star = fig03_04(summary, starts)
    fig05(star, x0)
    fig06(summary, star, x0)
    fig07(summary, star)
    fig08(summary)
    fig09(summary)
    folder = ROOT / "data" / EXP
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / "summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
