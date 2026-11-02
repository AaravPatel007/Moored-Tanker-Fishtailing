"""Compute the data for the slider demos (decision D8: Python computes, the page plays back).

Run from the project root (takes several minutes; every slider position is a full simulation):
    /home/aarav-patel/wec_rl/.venv/bin/python -m interactive.build_demos            # all three demos
    /home/aarav-patel/wec_rl/.venv/bin/python -m interactive.build_demos pitchfork  # just one

Writes interactive/<demo>/data.js. It is JSON wrapped in one line of JavaScript
(window.FISHTAIL_DEMOS["<demo>"] = {...}), because a web page opened straight from disk (file://) is not allowed to
read a .json file, but it may load a .js file. The page interactive/index.html loads all three.

Demos (rigid hawser, Vidal de Negreiros 100 %):
    pitchfork      slider x_A/L  (attachment point), L_H/L = 0.5, C_F = 0.05      -> E03
                   two runs per slider value: nudged to port ("runs") and to starboard ("runs_stbd")
    fishtailing    slider C_F    (friction, i.e. current speed), x_A/L = 0.5, L_H/L = 0.5   -> E04, E05
    hawser_length  slider L_H/L  (hawser length), x_A/L = 0.5, C_F = 0.03          -> E06
    wind           two sliders: wind force and wind direction; x_A/L = 0.25, L_H/L = 0.5, C_F = 0.05  -> E10
                   Settings stored force-major: index = i_force * n_directions + i_direction. The ship keeps its
                   history: "trans" maps "from_setting|side|to_setting" to a transition run (see build_wind).
Sign convention (Simos 2001): psi > 0 = bow turned to port, gamma > 0 = ship swung to port. Time t* = U t / L.
Angles are stored in degrees, rates in degrees per unit t*.
"""

import json
import multiprocessing
import sys
from dataclasses import replace

import contourpy
import numpy as np
from scipy.optimize import brentq

from fishtail import continuation as co
from fishtail import control
from fishtail import cycles as cy
from fishtail import hydro
from fishtail import model_rigid as mr
from fishtail.params import FRESH_WATER, SEA_WATER, VIDAL_100, X_T_L, X_TUG_L, froude_scaled
from fishtail.plotting import ROOT

HULL = VIDAL_100
T_END = 300.0  # length of each animated run [t*]
N_FRAMES = 901  # frames per run (one every 1/3 of t*)
START = [np.radians(3.0), np.radians(1.0), 0.0, 0.0]  # every run starts with a small nudge: gamma 3 deg, psi 1 deg
START_STBD = [-v for v in START]  # the same nudge to starboard (pitchfork demo: the page lets you pick the side)
MODEL = froude_scaled(HULL, FRESH_WATER.rho / SEA_WATER.rho)


def r(values, digits):
    """Round for a smaller file."""
    return [round(float(v), digits) for v in values]


def speed_for_friction(C_F, hull, fluid):
    """Current speed [m/s] at which the ITTC line (eq. 3.20) gives this C_F, or None if outside 0.01-20 m/s."""
    f = lambda U: hydro.friction_coefficient(U, hull, fluid) - C_F
    if f(0.01) * f(20.0) > 0:
        return None
    return float(brentq(f, 0.01, 20.0))


def eig_list(eig):
    return [[round(float(e.real), 4), round(float(e.imag), 4)] for e in eig]


def status_of(y):
    """What the ship is doing at the end of the run, from the last 30 t*."""
    late = y[:, -91:]
    swing = np.degrees(late[0].max() - late[0].min())
    if swing > 0.2:
        return "fishtails"
    if abs(np.degrees(y[1, -1])) <= 0.5:
        return "trails straight"
    return "sits at an angle to " + ("port" if y[0, -1] > 0 else "starboard")  # side = sign of the hawser angle


def one_run(s, start=START, status=None):
    """Animated run for one slider position, plus its equilibria and eigenvalues. status(y) names what the ship does
    (default status_of)."""
    status = status or status_of
    t, y = mr.simulate(s, start, t_end_star=T_END, n_out=N_FRAMES, rtol=1e-8, atol=1e-10)
    F = np.array([mr.hawser_force_star(y[:, i], s) for i in range(t.size)]) / s.C_F  # F_H / F_H0
    eq = []
    for state in co.all_equilibria(s):
        eig = co.eigenvalues(state, s)
        eq.append({"g": round(float(np.degrees(state[0])), 3), "p": round(float(np.degrees(state[1])), 3),
                   "stable": co.is_stable(eig), "eig": eig_list(eig),
                   "side": int(np.sign(round(float(np.degrees(state[0])), 3))),  # +1 port, -1 starboard, 0 straight
                   "F": round(float(mr.hawser_force_star(state, s) / s.C_F), 3)})  # hawser force at rest, F_H/F_H0
    U_model = speed_for_friction(s.C_F, MODEL, FRESH_WATER)
    U_full = speed_for_friction(s.C_F, HULL, SEA_WATER)
    return {
        "x_A_L": s.x_A_L, "L_H_L": s.L_H_L, "C_F": round(s.C_F, 5),
        "status": status(y), "equilibria": eq,
        "U_model": None if U_model is None else round(U_model, 3),
        "U_full": None if U_full is None else round(U_full, 3),
        "dt": round(float(t[1] - t[0]), 6),
        "g": r(np.degrees(y[0]), 2), "p": r(np.degrees(y[1]), 2),
        "gd": r(np.degrees(y[2]), 2), "pd": r(np.degrees(y[3]), 2), "F": r(F, 3),
    }


def crossing(kind, s_crit):
    """Eigenvalues of the straight position exactly at the bifurcation point. The page marks the ones that sit on the
    imaginary axis: one real eigenvalue at 0 (pitchfork) or a complex pair at +-i omega (Hopf)."""
    return {"type": kind, "eig": eig_list(co.eigenvalues(np.zeros(4), s_crit))}


def common(demo_id, title, knob, s0):
    """Fields every demo has. U_model / U_full: current speed [m/s] at which the tank model / the real ship has the
    fixed C_F of s0 (None if outside 0.01-20 m/s); with L_m and model_scale the page turns t* into real seconds,
    t = t* L / U (the Fishtailing tab, whose slider is C_F, uses the speeds stored with each run instead)."""
    U_model = speed_for_friction(s0.C_F, MODEL, FRESH_WATER)
    U_full = speed_for_friction(s0.C_F, HULL, SEA_WATER)
    return {"id": demo_id, "title": title, "knob": knob, "B_L": HULL.B / HULL.L, "hull": HULL.name,
            "L_m": HULL.L, "model_scale": HULL.scale, "fixed": {"x_A_L": s0.x_A_L, "L_H_L": s0.L_H_L, "C_F": s0.C_F},
            "U_model": None if U_model is None else round(U_model, 4),
            "U_full": None if U_full is None else round(U_full, 4),
            "start": {"g": 3.0, "p": 1.0}, "T_end": T_END}


def build_pitchfork():
    s0 = mr.Rigid(HULL, x_A_L=0.4, L_H_L=0.5, C_F=0.05)
    values = np.round(np.arange(0.0125, 0.501, 0.0125), 4)  # up to the bow (0.5): further forward is off the ship
    x_cr = hydro.critical_turret_position(HULL)
    # Bifurcation diagram: heading at rest against x_A/L (all equilibria), stability along each branch.
    xs = np.linspace(0.0, 0.5, 201)
    straight_stable, angled_x, angled_p, angled_stable = [], [], [], []
    for xa in xs:
        s = replace(s0, x_A_L=float(xa))
        straight_stable.append(co.is_stable(co.eigenvalues(np.zeros(4), s)))
        heads = co.equilibrium_headings(HULL, xa)
        if heads:
            state = co.equilibrium_state(heads[0], s)
            angled_x.append(float(xa))
            angled_p.append(float(np.degrees(heads[0])))
            angled_stable.append(co.is_stable(co.eigenvalues(state, s)))
    angled_x.append(x_cr)  # close the branches exactly at the pitchfork
    angled_p.append(0.0)
    angled_stable.append(angled_stable[-1])
    diagram = {
        "xlabel": "attachment point x_A/L", "ylabel": "heading at rest ψ [deg]",
        "branches": [
            {"x": r(xs, 4), "y": [0.0] * xs.size, "stable": straight_stable, "name": "trailing straight"},
            {"x": r(angled_x, 4), "y": r(angled_p, 3), "stable": angled_stable, "name": "angled (port)"},
            {"x": r(angled_x, 4), "y": r(-np.array(angled_p), 3), "stable": angled_stable, "name": "angled (starboard)"},
        ],
        "points": [{"x": round(x_cr, 4), "y": 0.0, "type": "pitchfork"}],
        "marker": "heading",
    }
    runs, runs_stbd = [], []  # same slider values, nudged to port / to starboard
    for v in values:
        print(f"pitchfork: x_A/L = {v:.4f}")
        runs.append(one_run(replace(s0, x_A_L=float(v))))
        runs_stbd.append(one_run(replace(s0, x_A_L=float(v)), start=START_STBD))
    knob = {"key": "x_A_L", "label": "Attachment point", "symbol": "x_A/L", "values": r(values, 4),
            "default": int(np.argmin(np.abs(values - 0.25))), "digits": 3,
            "help": "Where the hawser is tied, measured forward from midship (0.5 = at the bow)."}
    out = common("pitchfork", "Pitchfork", knob, s0)
    out.update({"diagram": diagram, "runs": runs, "runs_stbd": runs_stbd, "x_Tcr": round(x_cr, 4),
                "crossing": crossing("pitchfork", replace(s0, x_A_L=x_cr)),
                "story": "Slide the attachment point aft. Past x_Tcr = %.3f the straight position turns unstable and "
                         "the ship settles at an angle instead, to port or starboard." % x_cr})
    return out


def cycle_branch(s0, key, values, crit):
    """Amplitude of the fully grown swing for each slider value, sweeping from the most unstable end so that each
    run starts from the previous cycle. Returns gamma amplitudes [deg] (0 where the ship settles)."""
    amps = {}
    x = [np.radians(30.0), 0, 0, 0]
    order = np.argsort(values) if key == "C_F" else np.argsort(values)[::-1]  # most unstable first
    for i in order:
        s = replace(s0, **{key: float(values[i])})
        summ, x_end = cy.measure_cycle(s, x, t_transient=300.0)
        if summ is None:
            amps[i] = (0.0, 0.0)
            x = x_end + np.array([0.01, 0, 0, 0])
        else:
            amps[i] = (summ["gamma_amp_deg"], summ["psi_amp_deg"])
            x = summ["x0"]
        print(f"   cycle {key} = {values[i]:.4f}: gamma amp {amps[i][0]:.2f}")
    return np.array([amps[i][0] for i in range(len(values))]), np.array([amps[i][1] for i in range(len(values))])


def hopf_diagram(s0, key, values, crit, xlabel):
    """Bifurcation diagram for a Hopf knob: straight branch (stable / unstable) and +- cycle amplitude of gamma."""
    stable = [co.is_stable(co.eigenvalues(np.zeros(4), replace(s0, **{key: float(v)}))) for v in values]
    g_amp, p_amp = cycle_branch(s0, key, values, crit)
    has_cycle = g_amp > 0
    branches = [{"x": r(values, 5), "y": [0.0] * len(values), "stable": stable, "name": "trailing straight"}]
    # Cycle branch: from the Hopf point outwards, through every slider value that has a cycle
    cx = [crit] + [float(v) for v, h in zip(values, has_cycle) if h]
    cy_ = [0.0] + [float(a) for a, h in zip(g_amp, has_cycle) if h]
    idx = np.argsort(cx)
    cx, cy_ = np.array(cx)[idx], np.array(cy_)[idx]
    branches.append({"x": r(cx, 5), "y": r(cy_, 3), "stable": [True] * len(cx), "name": "swing amplitude (+)",
                     "cycle": True})
    branches.append({"x": r(cx, 5), "y": r(-cy_, 3), "stable": [True] * len(cx), "name": "swing amplitude (−)",
                     "cycle": True})
    return {"xlabel": xlabel, "ylabel": "hawser angle γ: rest value or swing range [deg]", "branches": branches,
            "points": [{"x": round(crit, 5), "y": 0.0, "type": "hopf"}], "marker": "swing",
            "amp_gamma": r(g_amp, 3), "amp_psi": r(p_amp, 3)}


def build_fishtailing():
    s0 = mr.Rigid(HULL, x_A_L=0.5, L_H_L=0.5, C_F=0.03)
    values = np.round(np.linspace(0.012, 0.06, 41), 5)
    crit = co.hopf_C_F(s0)
    diagram = hopf_diagram(s0, "C_F", values, crit, "friction coefficient C_F")
    runs = []
    for v in values:
        print(f"fishtailing: C_F = {v:.5f}")
        runs.append(one_run(replace(s0, C_F=float(v))))
    knob = {"key": "C_F", "label": "Friction coefficient", "symbol": "C_F", "values": r(values, 5),
            "default": int(np.argmin(np.abs(values - 0.03))), "digits": 4, "reverse": True,
            "help": "Skin friction of the hull. Faster current → larger Reynolds number → smaller C_F "
                    "(the only way current speed enters this model)."}
    out = common("fishtailing", "Fishtailing", knob, s0)
    out.update({"diagram": diagram, "runs": runs, "hopf": round(crit, 5), "crossing": crossing("hopf", replace(s0, C_F=crit)),
                "story": "Slide the friction down (= a faster current). Below C_F = %.4f the straight position "
                         "loses stability through a Hopf bifurcation and the ship fishtails; lower friction means a "
                         "bigger swing." % crit})
    return out


def build_hawser_length():
    s0 = mr.Rigid(HULL, x_A_L=0.5, L_H_L=0.5, C_F=0.03)
    values = np.round(np.concatenate([np.linspace(0.05, 0.4, 15), np.linspace(0.45, 1.5, 22)]), 4)
    f = lambda LH: co.max_oscillatory_growth(replace(s0, L_H_L=LH))
    crit = brentq(f, 0.05, 1.5)
    diagram = hopf_diagram(s0, "L_H_L", values, crit, "hawser length L_H/L")
    runs = []
    for v in values:
        print(f"hawser_length: L_H/L = {v:.4f}")
        runs.append(one_run(replace(s0, L_H_L=float(v))))
    knob = {"key": "L_H_L", "label": "Hawser length", "symbol": "L_H/L", "values": r(values, 4),
            "default": int(np.argmin(np.abs(values - 0.5))), "digits": 3,
            "help": "Length of the (rigid) hawser in ship lengths: 1.0 = 320 m."}
    out = common("hawser_length", "Hawser length", knob, s0)
    out.update({"diagram": diagram, "runs": runs, "hopf": round(crit, 4),
                "crossing": crossing("hopf", replace(s0, L_H_L=float(crit))),
                "story": "At this friction (C_F = %g, like a 1:90 tank model at about %.2f m/s) a short hawser keeps "
                         "the ship straight. Longer than L_H/L = %.2f, the ship fishtails."
                         % (s0.C_F, out["U_model"], crit)})
    return out


def wind_diagram(s, marker_x):
    """Bifurcation diagram (rest heading against x_A/L) for the wind of s, in the format the page draws."""
    br = co.wind_branches(s, x_hi=0.5, n_psi=401, n_straight=61)  # attachment point at most at the bow
    branches = [{"x": r(seg["x"], 4), "y": r(np.degrees(seg["psi"]), 2), "stable": [bool(v) for v in seg["stable"]]}
                for seg in br["segments"]]
    points = [{"x": round(x, 4), "y": round(float(np.degrees(p)), 2), "type": "fold"} for x, p in br["folds"]]
    if br["straight"] is not None:
        branches.append({"x": r(br["straight"]["x"], 4), "y": [0.0] * br["straight"]["x"].size,
                         "stable": [bool(v) for v in br["straight"]["stable"]]})
        points.append({"x": round(br["pitchfork"], 4), "y": 0.0, "type": "pitchfork"})
    return {"xlabel": "attachment point x_A/L", "ylabel": "heading at rest ψ [deg]", "branches": branches,
            "points": points, "marker": "heading", "marker_x": marker_x}


def wind_fold_ratio(s0, beta_deg):
    """Wind force (units F_H0) above which only ONE rest position is left at x_A/L of s0, for wind from beta_deg.
    At fixed x_A/L the rest condition is linear in the wind force (see experiments/E10_wind/run.py, s_curve), so
    every heading has one wind force; the fold is the smallest positive turning point of that curve.
    None for wind along the current (no sideways part: the symmetry survives, no fold)."""
    beta = np.radians(beta_deg)
    if abs(np.sin(beta)) < 1e-9:
        return None
    psi = np.radians(np.linspace(-60.0, 60.0, 6001))
    alpha = np.pi - psi
    moment = hydro.C6C(alpha, HULL) - s0.x_A_L * hydro.C2C(alpha, HULL)
    # rest condition: moment + (x_W - x_A) Y_w* = 0 with Y_w* = -C_W sin(beta - psi)
    with np.errstate(divide="ignore"):  # inf where psi = beta (wind along the hull: no sideways part)
        C_W = moment / ((s0.x_W_L - s0.x_A_L) * np.sin(beta - psi))
    d = np.diff(C_W)
    turns = [C_W[k + 1] for k in range(d.size - 1) if d[k] * d[k + 1] < 0 and C_W[k + 1] > 0]
    return round(float(min(turns) / s0.C_F), 4) if turns else None


def force_diagram(s):
    """Bifurcation diagram 'rest heading against wind force' for the direction of s, in the format the page draws."""
    br = co.force_branches(s, ratio_max=0.4, n_psi=2001)  # fine grid: the curve is steep near zero force
    for seg in br["segments"]:  # keep every 4th point (plus the ends) for a smaller file
        keep = np.unique(np.r_[np.arange(0, seg["ratio"].size, 4), seg["ratio"].size - 1])
        for key in ("ratio", "psi", "stable"):
            seg[key] = seg[key][keep]
    branches = [{"x": r(seg["ratio"], 4), "y": r(np.degrees(seg["psi"]), 2),
                 "stable": [bool(v) for v in seg["stable"]]} for seg in br["segments"]]
    points = [{"x": round(x, 4), "y": round(float(np.degrees(p)), 2), "type": "fold"} for x, p in br["folds"]]
    if br["straight"] is not None:
        branches.append({"x": r(br["straight"]["ratio"], 4), "y": [0.0] * br["straight"]["ratio"].size,
                         "stable": [bool(v) for v in br["straight"]["stable"]]})
    return {"xlabel": "wind force F_w/F_H0", "ylabel": "heading at rest ψ [deg]", "branches": branches,
            "points": points, "marker": "heading"}


# The wind demo keeps the ship's history (so hysteresis can be seen): moving a slider one step means "the wind changes
# to the new value and then stays", and the ship starts from where it was resting before. For every wind setting,
# every side the ship can rest on there (port / starboard), and every neighbouring setting one slider step away, we
# store one run: start = that rest position, wind = the neighbour's. The page chains these runs as the sliders move.
T_STEP = 150.0  # length of one transition run [t*]; long enough to settle, also next to a fold (slow passage)
N_STEP = 301  # frames per transition run (one every 0.5 t*)


def stable_by_side(s):
    """Stable rest positions of s, as {+1: port state, -1: starboard state} (side = sign of the hawser angle)."""
    out = {}
    for e in co.all_equilibria(s):
        if co.is_stable(co.eigenvalues(e, s)):
            side = 1 if e[0] > 0 else -1
            if side in out:
                raise RuntimeError("two stable rest positions on the same side - not expected here")
            out[side] = e
    return out


def trim_settled(run):
    """Drop the frames at the end of a run whose stored (rounded) values no longer change: the ship has settled. The
    page repeats the last frame to get back the full length (N_STEP), so nothing is lost at the stored precision.
    Keeps the wind data.js small enough to publish (most runs settle in the first third)."""
    keys = [k for k in ("g", "p", "gd", "pd", "F", "T") if k in run]
    last = len(run["g"]) - 1
    while last > 0 and all(run[k][last - 1] == run[k][-1] for k in keys):
        last -= 1
    for k in keys:
        run[k] = run[k][:last + 1]
    return run


def transition_run(s, start, targets):
    """One transition: wind of s switched on with the ship at rest in `start`. Also says on which side it ended."""
    # A swing right across to the other side next to a fold can take longer than T_STEP; then run on (same frame
    # spacing) up to 4 T_STEP, so the run really ends at a rest position and its side is named correctly.
    for k in (1, 4):
        t, y = mr.simulate(s, start, t_end_star=k * T_STEP, n_out=k * (N_STEP - 1) + 1, rtol=1e-8, atol=1e-10)
        end = y[:, -1]
        side = min(targets, key=lambda q: np.linalg.norm(targets[q][:2] - end[:2]))  # nearest stable rest position
        miss = np.degrees(np.linalg.norm(targets[side][:2] - end[:2]))
        if miss <= 0.5:
            break
        print(f"   note: after {k * T_STEP:.0f} t* the run is {miss:.2f} deg from a rest position" +
              (", running on" if k == 1 else " (still settling)"))
    F = np.array([mr.hawser_force_star(y[:, i], s) for i in range(t.size)]) / s.C_F
    return trim_settled({"dt": round(float(t[1] - t[0]), 6), "end_side": side,
                         "g": r(np.degrees(y[0]), 2), "p": r(np.degrees(y[1]), 2),
                         "gd": r(np.degrees(y[2]), 2), "pd": r(np.degrees(y[3]), 2), "F": r(F, 3)})


def nudged_runs(s, targets, size_deg=0.5):
    """Runs that start on the UNSTABLE rest position (the middle branch) of s, nudged slightly along its unstable
    direction (eigenvector of the growing eigenvalue), once each way. Returned as {end side: run}; {} if s has no
    unstable rest position (past a fold only one rest position is left)."""
    unstable = [e for e in co.all_equilibria(s) if not co.is_stable(co.eigenvalues(e, s))]
    if not unstable:
        return {}
    e = min(unstable, key=lambda q: abs(q[1]))
    w, V = np.linalg.eig(mr.jacobian_star(e, s))
    v = V[:, np.argmax(w.real)].real
    v = v / np.max(np.abs(v[:2])) * np.radians(size_deg) * np.sign(v[0])  # 0.5 deg nudge, towards port first
    out = {}
    for sign in (1, -1):
        run = transition_run(s, e + sign * v, targets)
        out[run["end_side"]] = run
    return out


def build_wind():
    """Two knobs: wind force (fraction of the steady hawser tension F_H0) and wind direction. Experiment E10."""
    s0 = mr.Rigid(HULL, x_A_L=0.25, L_H_L=0.5, C_F=0.05)  # as the pitchfork demo; aft of the pitchfork point
    ratios = np.round(np.linspace(0.0, 0.4, 31), 4)  # steps of 0.4/30 = 0.0133 F_H0 (student asked for 3x finer)
    betas = np.arange(-180, 181, 30)
    n1, n2 = ratios.size, betas.size
    setup = lambda i, j: replace(s0, C_W=float(ratios[i]) * s0.C_F, beta_W=np.radians(float(betas[j])))

    # 1. every wind setting: rest positions (with side and hawser force), diagram against x_A/L
    settings, rest = [], {}
    for i in range(n1):
        for j in range(n2):
            s = setup(i, j)
            rest[i, j] = stable_by_side(s)
            eq = []
            for e in co.all_equilibria(s):
                eig = co.eigenvalues(e, s)
                eq.append({"g": round(float(np.degrees(e[0])), 3), "p": round(float(np.degrees(e[1])), 3),
                           "stable": co.is_stable(eig), "eig": eig_list(eig), "side": 1 if e[0] > 0 else -1,
                           "F": round(float(mr.hawser_force_star(e, s) / s.C_F), 3)})
            settings.append({"equilibria": eq, "diagram": wind_diagram(s, s0.x_A_L), "x_A_L": s0.x_A_L,
                             "L_H_L": s0.L_H_L, "wind": {"ratio": float(ratios[i]), "beta": int(betas[j])}})
        print(f"wind: settings for force {ratios[i]:.2f} done")

    # 2. transitions to every neighbouring setting, from every side the ship can rest on
    runs, trans = [], {}
    for i in range(n1):
        for j in range(n2):
            for (a, b) in [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]:
                if not (0 <= a < n1 and 0 <= b < n2):
                    continue
                for side, start in rest[i, j].items():
                    trans[f"{i * n2 + j}|{side}|{a * n2 + b}"] = len(runs)
                    runs.append(transition_run(setup(a, b), start, rest[a, b]))
        print(f"wind: transitions from force {ratios[i]:.2f} done ({len(runs)} runs)")

    # 3. start on the unstable middle position, nudged to either side (lets the page show every branch)
    for i in range(n1):
        for j in range(n2):
            nudges = nudged_runs(setup(i, j), rest[i, j])
            settings[i * n2 + j]["nudge"] = {}
            for side, run in nudges.items():
                settings[i * n2 + j]["nudge"][str(side)] = len(runs)
                runs.append(run)
    print(f"wind: nudged runs done ({len(runs)} runs in total)")

    reference = wind_diagram(s0, s0.x_A_L)  # no wind: the perfect pitchfork, drawn in grey behind
    knob = {"key": "wind_ratio", "label": "Wind force", "symbol": "F_w/F_H0", "values": r(ratios, 4),
            "default": 0, "digits": 3,
            "help": "Size of the steady wind force on the ship, as a fraction of the hawser tension F_H0 when the "
                    "ship trails straight with no wind. It acts at midship and keeps its direction while the ship "
                    "turns (decision D28)."}
    knob2 = {"key": "beta_W", "label": "Wind from", "symbol": "β_W", "values": [int(b) for b in betas],
             "default": int(np.argmin(np.abs(betas - (-90)))), "digits": 0,
             "help": "Direction the wind comes from: 0° = from ahead (with the current), +90° = from port, "
                     "−90° = from starboard, ±180° = from astern."}
    out = common("wind", "Wind", knob, s0)
    out["knob"]["default"] = 0  # start with no wind ...
    out.update({"T_end": T_STEP, "n_frames": N_STEP,  # runs are trimmed once settled (trim_settled); page pads them
                "knob2": knob2, "runs": runs, "trans": trans, "settings": settings, "reference": reference,
                "wind": True, "start_side": 1,  # ... and the ship resting to port
                "force_diagrams": [force_diagram(setup(0, j)) for j in range(n2)],
                "fold_ratio": [wind_fold_ratio(s0, b) for b in betas],
                "xrange_xa": [0.0, 0.5], "yrange_xa": [-60.0, 60.0],
                "xrange_force": [0.0, 0.4], "yrange_force": [-20.0, 20.0],
                "story": "The attachment point is fixed at x_A/L = 0.25, aft of the pitchfork point, so without wind "
                         "the ship can rest on either side. Each slider step changes the wind and then holds it, and "
                         "the ship carries on from where it was, so it remembers its side: it rides its branch until "
                         "the branch ends at a fold (◇), then swings across, and it does not swing back when the wind "
                         "drops again. Use 'Put the ship on' (or click a branch in the diagram) to place the ship on "
                         "any rest position that exists for the current wind, including the unstable middle one, "
                         "from which a small nudge sends it to one side or the other. Wind along the current "
                         "(0°, ±180°) keeps the symmetry: both sides stay possible at any force."})
    return out


# --- Control tabs (step 8, D31, experiment E08) -----------------------------------------------------------------------
# Both tabs have two sliders but, unlike the wind tab, no history: every slider position is one run that starts from
# the FULL UNCONTROLLED FISHTAIL (the control is switched on at t* = 0 while the ship is swinging).
# Settings are stored knob-1-major (index = i1 * n2 + i2), like the wind tab. Runs are trimmed once settled.
T_CTRL = 200.0  # length of each control run [t*]
N_CTRL = 601  # frames per control run (3 per t*, as the other tabs)
# Slider ranges, chosen from the E08 results (data/E08_control/summary.json). PLACEHOLDER thruster/tug sizes (D31).
TUG_MAX = 2.4  # largest tug pull on the slider [F_H0]
KD_MAX, KP_MAX = 0.06, 0.03  # largest thruster gains on the sliders (the Hopf line and ★ lie well inside)
LIMITS = (np.inf, 0.02, 0.01)  # thrust-limit switch [F_H0]: none / medium / small (PLACEHOLDER)


def equilibria_of(s):
    """Rest positions of s with stability, eigenvalues, side and hawser force (the format of one_run)."""
    eq = []
    for state in co.all_equilibria(s):
        eig = co.eigenvalues(state, s)
        eq.append({"g": round(float(np.degrees(state[0])), 3), "p": round(float(np.degrees(state[1])), 3),
                   "stable": co.is_stable(eig), "eig": eig_list(eig),
                   "side": int(np.sign(round(float(np.degrees(state[0])), 3))),
                   "F": round(float(mr.hawser_force_star(state, s) / s.C_F), 3)})
    return eq


def control_run(s, start):
    """One control run from `start`: angles, rates, hawser force and (thruster tab) thruster force Y_T / F_H0."""
    t, y = mr.simulate(s, start, t_end_star=T_CTRL, n_out=N_CTRL, rtol=1e-8, atol=1e-10)
    F = np.array([mr.hawser_force_star(y[:, i], s) for i in range(t.size)]) / s.C_F
    status = status_of(y)
    if status == "fishtails" and co.is_stable(co.eigenvalues(np.zeros(4), s)):
        # Still swinging at the end of the shown run although the straight position is stable: is it a lasting
        # fishtail (thrust limit, bistable), or a swing still slowly dying out (near the Hopf line it decays slowly)?
        # Run on for 300 t* (not stored) and compare the last two 100 t* windows. Where the straight position is
        # unstable the swing cannot die out, so no check is needed (the student's suggestion: saves time).
        t2, y2 = mr.simulate(s, y[:, -1], t_end_star=300.0, n_out=901, rtol=1e-8, atol=1e-10)
        a_last, a_prev = np.ptp(y2[0, -300:]), np.ptp(y2[0, -600:-300])
        if not (np.degrees(a_last) > 0.5 and a_last > 0.95 * a_prev):
            status = "swing dying out"
    run = {"dt": round(float(t[1] - t[0]), 6), "status": status,
           "g": r(np.degrees(y[0]), 2), "p": r(np.degrees(y[1]), 2),
           "gd": r(np.degrees(y[2]), 2), "pd": r(np.degrees(y[3]), 2), "F": r(F, 3)}
    if control.is_on(s):
        run["T"] = r(control.thruster_ratio(y, s), 4)
    return trim_settled(run)


def uncontrolled_cycle(s):
    """A state on the fully grown fishtail of s (no control), or the small nudge START if s does not fishtail."""
    summ, _ = cy.measure_cycle(s, [np.radians(30.0), 0, 0, 0])
    return START if summ is None else list(summ["x0"])


def newtons_per_F_H0(C_F):
    """F_H0 [N] for the 1:90 tank model and for the real ship at the speed that gives C_F (None where none does)."""
    out = {}
    for name, hull, fluid in (("model", MODEL, FRESH_WATER), ("ship", HULL, SEA_WATER)):
        U = speed_for_friction(C_F, hull, fluid)
        out[name] = None if U is None else round(C_F * 0.5 * fluid.rho * U**2 * hull.L * hull.T, 4)
    return out


def cycle_sweep(setups):
    """Swing amplitude of gamma [deg] for each set-up, in the order given (most unstable first: each run starts from
    the previous cycle). 0 where the ship settles."""
    amps, x = [], [np.radians(30.0), 0, 0, 0]
    for s in setups:
        summ, x_end = cy.measure_cycle(s, x, t_transient=300.0)
        amps.append(0.0 if summ is None else summ["gamma_amp_deg"])
        x = x_end + np.array([0.01, 0, 0, 0]) if summ is None else summ["x0"]
    return np.array(amps)


def tug_setup(s, ratio):
    """Hold-back tug: a steady pull ratio * F_H0 downstream at the stern (the wind force with beta_W = 0, D31)."""
    return replace(s, C_W=float(ratio) * s.C_F, beta_W=0.0, x_W_L=X_TUG_L)


TUG_RATIOS = np.round(np.concatenate([np.linspace(0.0, 1.0, 41), np.linspace(1.2, TUG_MAX, 7)]), 4)
TUG_C_F = np.round(np.linspace(0.010, 0.058, 13), 4)  # 0.010 = real ship at about 1 m/s; 0.030 = the base case
TUG_BASE = mr.Rigid(HULL, x_A_L=0.5, L_H_L=0.5, C_F=0.03)
WORKERS = 10  # parallel processes (the student allowed up to 10 cores)


def tug_column(c):
    """Everything for one friction value C_F = c (one column of the tug grid): diagram, Hopf pull, runs, settings.
    Columns do not depend on each other, so they are computed in parallel (multiprocessing)."""
    s0, ratios = TUG_BASE, TUG_RATIOS
    n1 = ratios.size
    sc = replace(s0, C_F=float(c))
    start = uncontrolled_cycle(sc)
    sig = lambda q: -co.max_oscillatory_growth(tug_setup(sc, q))
    roots = co.find_crossings(sig, list(np.linspace(0.0, 4 * TUG_MAX, 41)))
    crit = roots[0] if roots else None  # tug pull at the Hopf point (None: no fishtailing, or beyond 4x range)
    stable = [co.is_stable(co.eigenvalues(np.zeros(4), tug_setup(sc, q))) for q in ratios]
    amps = cycle_sweep([tug_setup(sc, q) for q in ratios])
    has = amps > 0
    branches = [{"x": r(ratios, 4), "y": [0.0] * n1, "stable": stable, "name": "trailing straight"}]
    if has.any():
        cx = [float(v) for v in ratios[has]] + ([crit] if crit is not None and crit <= ratios[-1] else [])
        cyv = [float(a) for a in amps[has]] + ([0.0] if crit is not None and crit <= ratios[-1] else [])
        branches += [{"x": r(cx, 4), "y": r(cyv, 3), "stable": [True] * len(cx), "name": "swing (+)", "cycle": True},
                     {"x": r(cx, 4), "y": r(-np.array(cyv), 3), "stable": [True] * len(cx), "name": "swing (−)",
                      "cycle": True}]
    points = [] if crit is None or crit > ratios[-1] else [{"x": round(crit, 4), "y": 0.0, "type": "hopf"}]
    U_model, U_full = speed_for_friction(float(c), MODEL, FRESH_WATER), speed_for_friction(float(c), HULL, SEA_WATER)
    runs, settings = [], []
    for q in ratios:
        s = tug_setup(sc, q)
        run = control_run(s, start)
        run.update({"U_model": None if U_model is None else round(U_model, 3),
                    "U_full": None if U_full is None else round(U_full, 3)})
        runs.append(run)
        settings.append({"equilibria": equilibria_of(s), "x_A_L": s0.x_A_L, "L_H_L": s0.L_H_L,
                         "tug": {"ratio": float(q)}})
    print(f"tug: C_F = {c:.4f} done (Hopf pull {crit})", flush=True)
    return {"diagram": {"xlabel": "tug pull F_tug/F_H0", "ylabel": "hawser angle γ [deg]", "branches": branches,
                        "points": points, "marker": "swing"},
            "hopf": None if crit is None else round(crit, 4),
            "crossing": None if crit is None else crossing("hopf", tug_setup(sc, crit)),
            "scale": newtons_per_F_H0(float(c)), "runs": runs, "settings": settings}


def build_tug():
    """Two knobs: tug pull (fraction of F_H0) and friction C_F. Diagram: swing amplitude against the tug pull."""
    s0, ratios, C_Fs = TUG_BASE, TUG_RATIOS, TUG_C_F
    # fine steps up to 1 F_H0 (base case: Hopf at 0.30), coarser above (a real ship, C_F ~ 0.01, needs about 2.2)
    n1, n2 = ratios.size, C_Fs.size
    with multiprocessing.Pool(WORKERS) as pool:
        cols = pool.map(tug_column, [float(c) for c in C_Fs], chunksize=1)
    diagrams = [c["diagram"] for c in cols]
    hopfs = [c["hopf"] for c in cols]
    crossings = [c["crossing"] for c in cols]
    scale = [c["scale"] for c in cols]
    runs = [cols[j]["runs"][i] for i in range(n1) for j in range(n2)]  # knob-1-major: index = i * n2 + j
    settings = [cols[j]["settings"][i] for i in range(n1) for j in range(n2)]
    knob = {"key": "tug_ratio", "label": "Tug pull", "symbol": "F_tug/F_H0", "values": r(ratios, 4),
            "default": 0, "digits": 3,
            "help": "Steady pull of a hold-back tug made fast at the stern, pulling downstream, as a fraction of the "
                    "hawser pull F_H0 when the ship trails straight without a tug (PLACEHOLDER range, D31)."}
    knob2 = {"key": "C_F", "label": "Friction", "symbol": "C_F", "values": r(C_Fs, 4),
             "default": int(np.argmin(np.abs(C_Fs - 0.03))), "digits": 4,
             "help": "Skin friction of the hull: faster current = smaller C_F (the only way current speed enters)."}
    out = common("tug", "Tug", knob, s0)
    out.update({"T_end": T_CTRL, "n_frames": N_CTRL, "knob2": knob2, "grid": True, "runs": runs,
                "settings": settings, "diagrams": diagrams, "hopfs": hopfs, "crossings": crossings,
                "F_H0_N": scale, "tug": True,
                "yrange_swing": [-1.1 * max(max(abs(v) for b in dg["branches"] for v in b["y"]) for dg in diagrams) - 1,
                                 1.1 * max(max(abs(v) for b in dg["branches"] for v in b["y"]) for dg in diagrams) + 1],
                "story": "A hold-back tug at the stern pulls steadily downstream; every run starts while the ship is "
                         "fishtailing fully and the tug takes up its pull at t* = 0. When the bow swings off, the "
                         "tug's line no longer runs along the ship: its sideways part pulls the stern back, like the "
                         "tail of a weathervane, and the extra tension in the hawser pulls the ship back to the "
                         "middle. Above the Hopf pull ● the swing dies out."})
    return out


def contour_lines(xs, ys, Z, level):
    """Polylines where Z = level on the grid (xs across, ys up), for the page to draw: [[x...], [y...]] pieces."""
    gen = contourpy.contour_generator(xs, ys, Z)
    return [[r(seg[:, 0], 5), r(seg[:, 1], 5)] for seg in gen.lines(level)]


def thruster_column(job):
    """Runs for one thrust limit and one k_d, over all k_p (a column of the thruster grid)."""
    lim, d, kps, start = job
    s0 = TUG_BASE  # same base case: x_A/L = 0.5, L_H/L = 0.5, C_F = 0.03
    runs = [control_run(replace(s0, k_p=p, k_d=d, thrust_limit=lim), start) for p in kps]
    print(f"thruster: limit {lim}, k_d = {d:.4f} done", flush=True)
    return runs


def build_thruster():
    """Two knobs: derivative gain k_d (across) and proportional gain k_p; a thrust-limit switch (3 run sets).
    Diagram: map of the decay rate sigma over (k_d, k_p) with the Hopf line, the margin line and the least-thrust
    point ★ of E08."""
    s0 = mr.Rigid(HULL, x_A_L=0.5, L_H_L=0.5, C_F=0.03)
    e08 = json.loads((ROOT / "data" / "E08_control" / "summary.json").read_text())  # ★ and sigma_min from E08
    kds, kps = np.round(np.linspace(0.0, KD_MAX, 31), 5), np.round(np.linspace(0.0, KP_MAX, 16), 5)
    n1, n2 = kds.size, kps.size
    start = uncontrolled_cycle(s0)
    fd, fp = np.linspace(0, KD_MAX, 101), np.linspace(0, KP_MAX, 81)
    sig = np.array([[control.decay_rate(co.eigenvalues(np.zeros(4), replace(s0, k_p=float(p), k_d=float(d))))
                     for d in fd] for p in fp])
    settings = []
    for i, d in enumerate(kds):
        for j, p in enumerate(kps):
            s = replace(s0, k_p=float(p), k_d=float(d))
            settings.append({"equilibria": equilibria_of(s), "x_A_L": s0.x_A_L, "L_H_L": s0.L_H_L,
                             "gains": {"k_p": float(p), "k_d": float(d)}})
    # one job per (thrust limit, k_d): the runs for every k_p; jobs are independent, so run in parallel
    jobs = [(float(lim), float(d), [float(p) for p in kps], list(start)) for lim in LIMITS for d in kds]
    with multiprocessing.Pool(WORKERS) as pool:
        cols = pool.map(thruster_column, jobs, chunksize=1)
    run_sets = [[run for col in cols[k * n1:(k + 1) * n1] for run in col] for k in range(len(LIMITS))]
    knob = {"key": "k_d", "label": "Damping gain", "symbol": "k_d", "values": r(kds, 5), "digits": 3,
            "default": int(np.argmin(np.abs(kds - e08["star"]["k_d"]))),
            "help": "How hard the stern thruster pushes against the turning rate (a damper on yaw)."}
    knob2 = {"key": "k_p", "label": "Spring gain", "symbol": "k_p", "values": r(kps, 5), "digits": 3,
             "default": int(np.argmin(np.abs(kps - e08["star"]["k_p"]))),
             "help": "How hard the thruster turns the bow back towards straight per degree off (a spring)."}
    out = common("thruster", "Thruster", knob, s0)
    out.update({"T_end": T_CTRL, "n_frames": N_CTRL, "knob2": knob2, "grid": True, "runs": run_sets[0],
                "settings": settings, "thruster": True, "x_T_L": X_T_L,
                "limits": [None if not np.isfinite(l) else l for l in LIMITS], "limit_default": 0,
                "F_H0_N": newtons_per_F_H0(s0.C_F),
                "map": {"kd": r(fd, 5), "kp": r(fp, 5), "sigma": [r(row, 4) for row in sig],
                        "hopf": contour_lines(fd, fp, sig, 0.0),
                        "margin": contour_lines(fd, fp, sig, e08["star"]["sigma_min"]),
                        "sigma_min": e08["star"]["sigma_min"], "star": e08["star"]},
                "story": "The tanker's own stern thruster, run by a PD controller: N = −k_p ψ − k_d r. Every run "
                         "starts while the ship is fishtailing fully and the controller is switched on at t* = 0. "
                         "The map shows how fast a small swing dies out for each pair of gains; ★ is the pair that "
                         "needs the least thrust while still meeting the margin (E08). With a thrust limit the "
                         "thruster may be too weak for the big swing even where the map says 'stable'."})
    out["_extra"] = {f"limit{k}": runs for k, runs in enumerate(run_sets) if k > 0}  # written to their own files
    return out


BUILDERS = {"pitchfork": build_pitchfork, "fishtailing": build_fishtailing, "hawser_length": build_hawser_length,
            "wind": build_wind, "tug": build_tug, "thruster": build_thruster}


def main():
    names = sys.argv[1:] or list(BUILDERS)
    for name in names:
        data = BUILDERS[name]()
        folder = ROOT / "interactive" / name
        folder.mkdir(parents=True, exist_ok=True)
        # Extra run sets (thruster tab: one per thrust limit) go to their own files, so each stays below 16 MB:
        # interactive/<name>/<key>.js sets window.FISHTAIL_EXTRA["<name>"]["<key>"].
        for key, value in data.pop("_extra", {}).items():
            text = json.dumps(value, separators=(",", ":"))
            with open(folder / f"{key}.js", "w") as f:
                f.write("window.FISHTAIL_EXTRA = window.FISHTAIL_EXTRA || {};\n")
                f.write(f'window.FISHTAIL_EXTRA["{name}"] = window.FISHTAIL_EXTRA["{name}"] || {{}};\n')
                f.write(f'window.FISHTAIL_EXTRA["{name}"]["{key}"] = {text};\n')
            print(f"wrote interactive/{name}/{key}.js ({len(text) / 1e6:.2f} MB)")
        text = json.dumps(data, separators=(",", ":"))
        with open(folder / "data.js", "w") as f:
            f.write("window.FISHTAIL_DEMOS = window.FISHTAIL_DEMOS || {};\n")
            f.write(f'window.FISHTAIL_DEMOS["{name}"] = {text};\n')
        print(f"wrote interactive/{name}/data.js ({len(text) / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
