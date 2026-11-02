"""Equilibria of the rigid-hawser system, their eigenvalues, and the bifurcation points (steps 4 and 5).

Sign convention (Simos 2001): x forward, y to port; psi > 0 = bow turned to port; gamma > 0 = ship swung to port.
All rates and eigenvalues are in units of U/L (time t* = U t / L).

Equilibria ("the ship at rest", gamma_dot = psi_dot = 0). Midship then sits still, so the water passes it at the
current speed: u = U cos psi, v = -U sin psi, r = 0, angle of attack alpha = pi - psi. The three equations (5.1) become
    surge: X(psi) + F_H cos(gamma + psi) = 0
    sway : Y(psi) - F_H sin(gamma + psi) = 0
    yaw  : N(psi) - x_A F_H sin(gamma + psi) = 0
Yaw minus x_A times sway: N(psi) - x_A Y(psi) = 0, i.e. the water's moment about the ATTACHMENT POINT is zero.
This is one equation for psi alone: in coefficients  C6C(pi - psi) - (x_A/L) C2C(pi - psi) = 0.
So the equilibrium heading depends only on x_A/L (not on the hawser length or friction). The surge and sway equations
then give the hawser force and direction: F_H = sqrt(X^2 + Y^2), gamma + psi = atan2(Y, -X).
"""

from dataclasses import replace

import numpy as np
from scipy.optimize import brentq

from fishtail import control, hydro
from fishtail import model_rigid as mr


# --- Equilibria -----------------------------------------------------------------------------------------------------

def moment_about_attachment(psi, hull, x_A_L):
    """Static water moment about the attachment point, in coefficient form: C6C - (x_A/L) C2C at alpha = pi - psi.

    > 0 turns the bow to port. At an equilibrium it is zero.
    """
    alpha = np.pi - psi
    return hydro.C6C(alpha, hull) - x_A_L * hydro.C2C(alpha, hull)


def equilibrium_headings(hull, x_A_L, n_grid=3601):
    """All headings psi in (0, pi) where the moment about the attachment point is zero (sign changes on a grid,
    then refined with Brent's method). Mirror images -psi are also equilibria; psi = 0 and psi = pi always are."""
    grid = np.linspace(1e-6, np.pi - 1e-6, n_grid)
    values = moment_about_attachment(grid, hull, x_A_L)
    roots = []
    for a, b, fa, fb in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if fa * fb < 0:
            roots.append(brentq(moment_about_attachment, a, b, args=(hull, x_A_L), xtol=1e-13))
    return roots


def equilibrium_state(psi, s):
    """Full state (gamma, psi, 0, 0) of the equilibrium with heading psi, with the hawser in tension.

    (With compression allowed there is a second solution, the ship upstream of the buoy, pushed by the bar;
    we do not follow it - it is not what a moored ship does.)  Includes the wind if s.C_W > 0 (D28).
    """
    X, Y, _ = hydro.static_forces(s.U * np.cos(psi), -s.U * np.sin(psi), s.hull, s.fluid, s.C_F)
    if s.C_W != 0.0:
        X_w, Y_w, _ = mr.wind_forces(psi, s)
        X, Y = X + X_w, Y + Y_w
    if control.is_on(s):
        Y = Y + control.thruster_force(psi, 0.0, s)
    gamma = np.arctan2(Y, -X) - psi
    gamma = (gamma + np.pi) % (2 * np.pi) - np.pi  # wrap into (-pi, pi]
    return np.array([gamma, psi, 0.0, 0.0])


def all_equilibria(s, include_reversed=False):
    """List of equilibrium states for set-up s: trailing straight, and each angled pair (+psi, -psi).

    include_reversed adds psi = pi (ship trailing stern-first); it is never stable here and is left out of plots.
    With wind there is no mirror symmetry and no exactly straight position: see wind_equilibria.
    """
    if s.C_W != 0.0 or control.is_on(s):
        if not wind_is_symmetric(s):
            return wind_equilibria(s)
        # Tug or head-on wind and/or thruster (D31): still mirror-symmetric, so straight + angled pairs (+psi, -psi)
        states = [np.zeros(4)]
        for st in wind_equilibria(s):
            if st[1] > 1e-9:
                states += [st, equilibrium_state(-st[1], s)]
        return states
    states = [np.zeros(4)]
    for psi in equilibrium_headings(s.hull, s.x_A_L):
        states.append(equilibrium_state(psi, s))
        states.append(equilibrium_state(-psi, s))
    if include_reversed:
        states.append(equilibrium_state(np.pi, s))
    return states


# --- Equilibria with wind (the imperfect pitchfork, E10) ---------------------------------------------------------------
# The trick above still works: yaw minus x_A times sway removes the unknown hawser force, so the heading at rest
# solves ONE equation  C6C + (x_W/L) Y_w* - (x_A/L) (C2C + Y_w*) = 0,  with Y_w* = -C_W sin(beta_W - psi) the wind's
# sideways force in coefficient form. This is linear in x_A/L, so every heading psi has exactly one attachment point:
#     x_A/L = (C6C + (x_W/L) Y_w*) / (C2C + Y_w*).
# Drawing (x_A/L, psi) for all psi gives the whole bifurcation diagram, folds included, with no root-finding.

def wind_side_coefficient(psi, s):
    """Sideways wind force in body axes, divided by (1/2) rho U^2 L T."""
    return -s.C_W * np.sin(s.beta_W - psi)


def moment_about_attachment_wind(psi, s):
    """Water + wind (+ stern thruster, D31) moment about the attachment point, coefficient form (> 0 turns the bow
    to port)."""
    alpha = np.pi - psi
    Yw = wind_side_coefficient(psi, s)
    m = hydro.C6C(alpha, s.hull) + s.x_W_L * Yw - s.x_A_L * (hydro.C2C(alpha, s.hull) + Yw)
    if control.is_on(s):  # the thruster's side force Y_T* at x_T, about the attachment point
        m = m + (s.x_T_L - s.x_A_L) * np.vectorize(control.thruster_side_coefficient)(psi, s)
    return m


def attachment_for_heading(psi, s):
    """The attachment point x_A/L at which heading psi is an equilibrium (inf where the side force is zero)."""
    alpha = np.pi - psi
    Yw = wind_side_coefficient(psi, s)
    num = hydro.C6C(alpha, s.hull) + s.x_W_L * Yw
    den = hydro.C2C(alpha, s.hull) + Yw
    with np.errstate(divide="ignore", invalid="ignore"):
        return num / den


def wind_equilibria(s, psi_max=np.radians(89.0), n_grid=3601):
    """All equilibrium states with the bow within psi_max of the current (sorted by psi), wind included."""
    grid = np.linspace(-psi_max, psi_max, n_grid)
    values = moment_about_attachment_wind(grid, s)
    states = []
    for a, b, fa, fb in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if fa == 0.0:
            states.append(equilibrium_state(a, s))
        elif fa * fb < 0:
            states.append(equilibrium_state(brentq(moment_about_attachment_wind, a, b, args=(s,), xtol=1e-13), s))
    return states


def wind_is_symmetric(s):
    """True when the wind has no sideways part (none, from dead ahead or from dead astern): then the mirror symmetry
    survives, trailing straight is still an equilibrium, and the bifurcation is still a perfect pitchfork."""
    return s.C_W == 0.0 or abs(np.sin(s.beta_W)) < 1e-12


def wind_branches(s, x_lo=0.0, x_hi=0.6, psi_max_deg=85.0, n_psi=681, n_straight=121):
    """Bifurcation diagram 'heading at rest psi against attachment point x_A/L' for the wind of s (E10, demo).

    Returns a dict:
      segments  list of dicts (x, psi [rad], stable): pieces of the curve x_A/L(psi) inside [x_lo, x_hi]
      folds     list of (x_A/L, psi): turning points of the curve (saddle-node bifurcations), where dx_A/dpsi = 0
      straight  dict (x, stable) for trailing straight, only when the wind is symmetric, else None
      pitchfork x_A/L where the angled branches meet trailing straight (symmetric wind only), else None
    """
    u = np.linspace(-1.0, 1.0, n_psi)
    psi = np.radians(psi_max_deg) * u**3  # dense near psi = 0, where the branches flatten out towards large x_A/L
    xa = attachment_for_heading(psi, s)
    ok = np.isfinite(xa) & (xa >= x_lo) & (xa <= x_hi)
    segments, folds = [], []
    i = 0
    while i < psi.size:  # cut the curve into connected pieces that stay inside the plotted range
        if not ok[i]:
            i += 1
            continue
        j = i
        while j + 1 < psi.size and ok[j + 1] and abs(xa[j + 1] - xa[j]) < 0.05:
            j += 1
        idx = np.arange(i, j + 1)
        stable = []
        for k in idx:
            sk = replace(s, x_A_L=float(xa[k]))
            stable.append(is_stable(eigenvalues(equilibrium_state(psi[k], sk), sk)))
        segments.append({"x": xa[idx], "psi": psi[idx], "stable": np.array(stable)})
        dx = np.diff(xa[idx])
        for k in range(dx.size - 1):  # a fold is where x_A/L stops growing and turns back
            at_pitchfork = wind_is_symmetric(s) and abs(psi[idx[k + 1]]) < 1e-3  # symmetric: that turn is the pitchfork
            if dx[k] * dx[k + 1] < 0 and not at_pitchfork:
                folds.append((float(xa[idx[k + 1]]), float(psi[idx[k + 1]])))
        i = j + 1
    straight, pitchfork = None, None
    if wind_is_symmetric(s):
        xs = np.linspace(x_lo, x_hi, n_straight)
        stable = [is_stable(eigenvalues(np.zeros(4), replace(s, x_A_L=float(x)))) for x in xs]
        straight = {"x": xs, "stable": np.array(stable)}
        pitchfork = float(attachment_for_heading(1e-7, s))  # where the angled branch reaches psi = 0
    return {"segments": segments, "folds": folds, "straight": straight, "pitchfork": pitchfork}


def force_branches(s, ratio_max=0.4, psi_max_deg=60.0, n_psi=801):
    """Bifurcation diagram 'heading at rest psi against wind force' for fixed x_A/L and wind direction of s.

    Wind force as the ratio F_wind / F_H0 = C_W / C_F. At fixed x_A/L the rest condition
        C6C - x_A C2C + (x_W - x_A) Y_w* = 0,   Y_w* = -C_W sin(beta_W - psi)
    is linear in C_W, so every heading has exactly one wind force: C_W = (C6C - x_A C2C) / ((x_W - x_A) sin(beta_W - psi)).
    Returns a dict like wind_branches: segments (ratio, psi, stable), folds (ratio, psi), straight (symmetric wind only:
    trailing straight is then a rest position for every force).
    """
    u = np.linspace(-1.0, 1.0, n_psi)
    psi = np.radians(psi_max_deg) * u**3  # dense near psi = 0
    alpha = np.pi - psi
    moment = hydro.C6C(alpha, s.hull) - s.x_A_L * hydro.C2C(alpha, s.hull)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = moment / ((s.x_W_L - s.x_A_L) * np.sin(s.beta_W - psi)) / s.C_F
    ok = np.isfinite(ratio) & (ratio >= 0.0) & (ratio <= ratio_max)
    segments, folds = [], []
    i = 0
    while i < psi.size:
        if not ok[i]:
            i += 1
            continue
        j = i
        while j + 1 < psi.size and ok[j + 1] and abs(ratio[j + 1] - ratio[j]) < 0.05:
            j += 1
        idx = np.arange(i, j + 1)
        stable = []
        for k in idx:
            sk = replace(s, C_W=float(ratio[k]) * s.C_F)
            stable.append(is_stable(eigenvalues(equilibrium_state(psi[k], sk), sk)))
        segments.append({"ratio": ratio[idx], "psi": psi[idx], "stable": np.array(stable)})
        d = np.diff(ratio[idx])
        for k in range(d.size - 1):
            if d[k] * d[k + 1] < 0:
                folds.append((float(ratio[idx[k + 1]]), float(psi[idx[k + 1]])))
        i = j + 1
    straight = None
    if abs(np.sin(s.beta_W)) < 1e-12:  # wind along the current: trailing straight stays a rest position
        rs = np.linspace(0.0, ratio_max, 41)
        stable = []
        for rr in rs:
            sk = replace(s, C_W=float(rr) * s.C_F)
            stable.append(is_stable(eigenvalues(np.zeros(4), sk)))
        straight = {"ratio": rs, "stable": np.array(stable)}
    return {"segments": segments, "folds": folds, "straight": straight}


# --- Stability ------------------------------------------------------------------------------------------------------

def eigenvalues(state, s):
    """Eigenvalues (units U/L) of the Jacobian at an equilibrium. Exact formula at trailing straight without wind,
    otherwise finite differences."""
    if np.allclose(state, 0.0) and wind_is_symmetric(s):
        J = mr.linear_jacobian_star(s)
    else:
        J = mr.jacobian_star(state, s)
    eig = np.linalg.eigvals(J)
    return eig[np.argsort(-eig.real)]  # largest real part first


def is_stable(eig):
    return bool(np.all(eig.real < 0))


def classify(eig, tol=1e-9):
    """Name the equilibrium type the way Strogatz does (for a 4-D system: by its leading eigenvalues)."""
    lead = eig[0]
    if lead.real < -tol:
        return "stable focus" if np.any(np.abs(eig.imag) > tol) else "stable node"
    n_unstable = int(np.sum(eig.real > tol))
    if np.abs(lead.imag) > tol:
        return "unstable focus (oscillatory growth)"
    return "saddle" if n_unstable < 4 else "unstable node"


def max_oscillatory_growth(s):
    """Largest real part among the COMPLEX eigenvalues of trailing straight (the Hopf pair); -inf if none."""
    eig = np.linalg.eigvals(mr.linear_jacobian_star(s))
    complex_ones = eig[np.abs(eig.imag) > 1e-9]
    return complex_ones.real.max() if complex_ones.size else -np.inf


def oscillatory_growth(state, s):
    """Largest real part among the COMPLEX eigenvalues at any equilibrium state; -inf if all are real."""
    eig = eigenvalues(state, s)
    complex_ones = eig[np.abs(eig.imag) > 1e-7]
    return complex_ones.real.max() if complex_ones.size else -np.inf


def find_crossings(f, grid):
    """All points in [grid[0], grid[-1]] where f changes sign (scan the grid, then Brent's method)."""
    values = [f(x) for x in grid]
    roots = []
    for a, b, fa, fb in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if np.isfinite(fa) and np.isfinite(fb) and fa * fb < 0:
            roots.append(brentq(f, a, b, xtol=1e-12))
    return roots


# --- Bifurcation points ---------------------------------------------------------------------------------------------

def hopf_C_F(s, C_lo=1e-4, C_hi=0.2, n_grid=60):
    """Friction coefficient at which trailing straight loses stability through a Hopf bifurcation (complex pair
    crossing the imaginary axis), for the hull, x_A/L and L_H/L of s. Returns None if no crossing in [C_lo, C_hi].

    Scans a grid for the first sign change of the pair's real part (from unstable at low friction to stable at
    high friction), then refines with Brent's method.
    """
    f = lambda c: max_oscillatory_growth(replace(s, C_F=c))
    grid = np.geomspace(C_lo, C_hi, n_grid)
    values = [f(c) for c in grid]
    for a, b, fa, fb in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if np.isfinite(fa) and np.isfinite(fb) and fa > 0 > fb:
            return brentq(f, a, b, xtol=1e-12)
    return None


def hopf_frequency(s):
    """Imaginary part (units U/L) of the eigenvalue pair with the largest real part: the angular frequency of the
    swing near a Hopf point. The period in t* is 2 pi / omega."""
    eig = np.linalg.eigvals(mr.linear_jacobian_star(s))
    complex_ones = eig[np.abs(eig.imag) > 1e-9]
    return float(np.abs(complex_ones[np.argmax(complex_ones.real)].imag))


def hopf_curve(s, L_H_values, simos_annex_b=False):
    """Hopf friction C_F for each hawser length (the stability boundary of Simos Fig. 5.3)."""
    out = []
    for LH in L_H_values:
        s2 = replace(s, L_H_L=LH)
        if simos_annex_b:
            g = lambda c: _growth_simos(replace(s2, C_F=c))
            grid = np.geomspace(1e-4, 0.5, 80)
            vals = [g(c) for c in grid]
            root = next((brentq(g, a, b) for a, b, fa, fb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:])
                         if fa > 0 > fb), None)
            out.append(np.nan if root is None else root)
        else:
            root = hopf_C_F(s2)
            out.append(np.nan if root is None else root)
    return np.array(out)


def _growth_simos(s):
    eig = np.linalg.eigvals(mr.linear_jacobian_star(s, simos_annex_b=True))
    return eig.real.max()
