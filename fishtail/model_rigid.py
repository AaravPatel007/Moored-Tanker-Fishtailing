"""Ship held by a RIGID hawser (always straight, fixed length): 4 states (Simos 2001, Ch. 5 and Annex B).

Sign convention (Simos 2001 thesis, used everywhere in this package):
    body axes at MIDSHIP, x forward, y to PORT; u, v = hull velocity RELATIVE TO THE WATER; r = yaw rate.
Earth axes (top view):
    X points UPSTREAM (into the current), Y to port of a ship pointing upstream. The buoy is at the origin and the
    water flows in the -X direction at speed U. The ship trails downstream of the buoy (at negative X).
    psi   = heading, angle from +X to the bow; psi = 0: bow straight into the current; psi > 0: bow turned to port.
    gamma = hawser angle to the current line; the attachment point A is at L_H (-cos gamma, sin gamma),
            so gamma > 0 means the ship has swung to the +Y (port) side.
States: (gamma, psi, gamma_dot, psi_dot).

The equations come from derivations/rigid_hawser.py (sympy); results in derivations/rigid_hawser_results.md.
    kinematics (eq. 5.2):  u = U cos psi + L_H gamma_dot sin(gamma + psi)
                           v = -U sin psi + L_H gamma_dot cos(gamma + psi) - x_A psi_dot,     r = psi_dot
    hawser loads (eq. B.3; the '+' printed in eq. 5.3 is a misprint, open issue 4):
                           F_X = F_H cos(gamma + psi),  F_Y = -F_H sin(gamma + psi),  N = -x_A F_H sin(gamma + psi)
    equations of motion (eq. 5.1) with the Munk term inside C6C (decision D4).
The hawser force F_H is whatever keeps the hawser length fixed, so it is an unknown like the two accelerations:
at every instant we solve 3 linear equations for (gamma_ddot, psi_ddot, F_H).

Wind (decision D28, experiment E10): a steady force of fixed size and fixed direction in EARTH axes, acting at a
point x_W on the centreline (default: midship). Its size is given as C_W = F_wind / ((1/2) rho U^2 L T), the same
scale as the current forces, so trailing straight F_wind / F_H0 = C_W / C_F. beta_W is the direction the wind comes
FROM, measured from +X (upstream) towards +Y (port): 0 = from ahead, like the current; 90 deg = from port; 180 deg =
from astern. The wind does not change with the ship's heading (no wind coefficients): our simplification.

Two ways to measure time:
    rhs()       seconds (used to show that current speed only matters through C_F, experiment E02);
    rhs_star()  t* = U t / L, rates in units of U/L (used for all bifurcation work: then U drops out, D1).
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from fishtail import hydro
from fishtail import control
from fishtail.params import SEA_WATER, X_T_L, Fluid, Hull


@dataclass(frozen=True)
class Rigid:
    """One mooring set-up: the hull, the three knobs, and the current speed."""

    hull: Hull
    x_A_L: float  # attachment point, forward of midship, as a fraction of L [-]   (pitchfork knob)
    L_H_L: float  # hawser length as a fraction of L [-]                            (Hopf knob)
    C_F: float  # surge friction coefficient C0 of eq. (3.20) [-]                 (stands for current speed, D1)
    fluid: Fluid = SEA_WATER
    U: float = 1.0  # current speed [m/s]. Results in t* do not depend on it while C_F is held fixed.
    C_W: float = 0.0  # wind force / ((1/2) rho U^2 L T) [-]; 0 = no wind                 (imperfection size, E10)
    beta_W: float = 0.0  # direction the wind comes FROM, from +X towards +Y (port) [rad]  (imperfection direction)
    x_W_L: float = 0.0  # where the wind force acts, forward of midship, fraction of L [-]
    # Stern thruster with PD heading control (step 8, D31; fishtail/control.py). Zero gains = no control.
    k_p: float = 0.0  # proportional gain: yaw moment per radian of heading, in (1/2) rho U^2 L^2 T [-]
    k_d: float = 0.0  # derivative gain: yaw moment per unit of r* = r L/U, same scale [-]
    thrust_limit: float = np.inf  # largest thruster force, as a fraction of F_H0 [-]  (PLACEHOLDER values)
    x_T_L: float = X_T_L  # where the thruster pushes, forward of midship, fraction of L [-] (the stern)

    @property
    def x_A(self):
        return self.x_A_L * self.hull.L  # [m]

    @property
    def L_H(self):
        return self.L_H_L * self.hull.L  # [m]

    @property
    def time_unit(self):
        return self.hull.L / self.U  # L/U [s]: one unit of t*


# --- Equations of motion (dimensional, SI units, time in seconds) ---------------------------------------------------

def velocities(state, s):
    """Hull velocity relative to the water (u, v [m/s]) and yaw rate r [rad/s], eq. (5.2)."""
    gamma, psi, gamma_dot, psi_dot = state
    u = s.U * np.cos(psi) + s.L_H * gamma_dot * np.sin(gamma + psi)
    v = -s.U * np.sin(psi) + s.L_H * gamma_dot * np.cos(gamma + psi) - s.x_A * psi_dot
    return u, v, psi_dot


def wind_forces(psi, s):
    """Wind force on the ship in BODY axes: F_X [N], F_Y [N], N [N m] about midship.

    In earth axes the wind blows towards -(cos beta_W, sin beta_W). Projecting on the body axes
    e_x = (cos psi, sin psi) and e_y = (-sin psi, cos psi) gives the two components below.
    """
    F_w = s.C_W * 0.5 * s.fluid.rho * s.U**2 * s.hull.L * s.hull.T  # size of the wind force [N]
    X_w = -F_w * np.cos(s.beta_W - psi)
    Y_w = -F_w * np.sin(s.beta_W - psi)
    return X_w, Y_w, s.x_W_L * s.hull.L * Y_w


def accelerations(state, s):
    """Solve the 3 equations of motion for (gamma_ddot [rad/s^2], psi_ddot [rad/s^2], F_H [N])."""
    gamma, psi, gamma_dot, psi_dot = state
    h = s.hull
    m11, m22 = h.M + h.M11, h.M + h.M22  # mass + added mass in surge, sway
    m26 = h.M * h.x_CG + h.M26  # sway-yaw coupling (the CG is ahead of midship, plus added mass)
    i66 = h.I_Zo + h.M66  # yaw inertia about midship + added inertia
    L_H, x_A = s.L_H, s.x_A
    sn, cs = np.sin(gamma + psi), np.cos(gamma + psi)  # direction of the hawser seen from the ship

    u, v, r = velocities(state, s)
    X_R, Y_R, N_R = hydro.current_forces(u, v, r, h, s.fluid, s.C_F)  # water force on the hull
    if s.C_W != 0.0:  # add the wind (D28); with C_W = 0 the model is exactly the one of steps 3-6
        X_w, Y_w, N_w = wind_forces(psi, s)
        X_R, Y_R, N_R = X_R + X_w, Y_R + Y_w, N_R + N_w
    if control.is_on(s):  # add the stern thruster (D31): a side force Y_T at x_T, so also a yaw moment x_T Y_T
        Y_T = control.thruster_force(psi, psi_dot, s)
        Y_R, N_R = Y_R + Y_T, N_R + s.x_T_L * h.L * Y_T

    # Time derivatives of u and v, split into a part with the unknown accelerations and a known rest a_u, a_v:
    #   du/dt = L_H sn gamma_ddot + a_u,       dv/dt = L_H cs gamma_ddot - x_A psi_ddot + a_v
    a_u = -s.U * np.sin(psi) * psi_dot + L_H * cs * gamma_dot * (gamma_dot + psi_dot)
    a_v = -s.U * np.cos(psi) * psi_dot - L_H * sn * gamma_dot * (gamma_dot + psi_dot)

    # Equations (5.1) with unknowns x = [gamma_ddot, psi_ddot, F_H]; known terms moved to the right.
    #   surge: m11 du/dt - m22 v r - m26 r^2           = X_R + F_H cs
    #   sway : m22 dv/dt + m11 u r + m26 dr/dt         = Y_R - F_H sn
    #   yaw  : i66 dr/dt + m26 (dv/dt + u r)           = N_R - x_A F_H sn
    A = np.array([
        [m11 * L_H * sn, 0.0, -cs],
        [m22 * L_H * cs, m26 - m22 * x_A, sn],
        [m26 * L_H * cs, i66 - m26 * x_A, x_A * sn],
    ])
    b = np.array([
        X_R + m22 * v * r + m26 * r**2 - m11 * a_u,
        Y_R - m11 * u * r - m22 * a_v,
        N_R - m26 * (a_v + u * r),
    ])
    gamma_ddot, psi_ddot, F_H = np.linalg.solve(A, b)
    return gamma_ddot, psi_ddot, F_H


def rhs(t, state, s):
    """d(state)/dt in seconds. t is unused: the current is steady, so the system is autonomous."""
    gamma_ddot, psi_ddot, _ = accelerations(state, s)
    return np.array([state[2], state[3], gamma_ddot, psi_ddot])


def hawser_force(state, s):
    """Hawser force F_H [N] (tension > 0, compression < 0) at a dimensional state."""
    return accelerations(state, s)[2]


# --- Non-dimensional time t* = U t / L (all bifurcation work) ---------------------------------------------------------

def to_dimensional(state_star, s):
    """(gamma, psi, dgamma/dt*, dpsi/dt*) -> (gamma, psi, dgamma/dt, dpsi/dt) in rad/s."""
    g, p, gd, pd = state_star
    return np.array([g, p, gd / s.time_unit, pd / s.time_unit])


def rhs_star(t_star, state_star, s):
    """d(state*)/dt* with t* = U t / L; angular rates in units of U/L."""
    T = s.time_unit
    g_dd, p_dd, _ = accelerations(to_dimensional(state_star, s), s)
    return np.array([state_star[2], state_star[3], g_dd * T**2, p_dd * T**2])


def hawser_force_star(state_star, s):
    """Hawser force divided by the friction drag scale (1/2) rho U^2 L T. Trailing straight it equals C_F."""
    h = s.hull
    return hawser_force(to_dimensional(state_star, s), s) / (0.5 * s.fluid.rho * s.U**2 * h.L * h.T)


def jacobian_star(state_star, s, h=1e-6):
    """4x4 Jacobian of rhs_star by central finite differences. Its eigenvalues are in units of U/L."""
    J = np.zeros((4, 4))
    x0 = np.asarray(state_star, dtype=float)
    for k in range(4):
        dx = np.zeros(4)
        dx[k] = h
        J[:, k] = (rhs_star(0, x0 + dx, s) - rhs_star(0, x0 - dx, s)) / (2 * h)
    return J


# --- Exact linearisation about trailing straight (gamma = psi = 0) ---------------------------------------------------

def linear_matrices(s, simos_annex_b=False):
    """Mass, damping and stiffness matrices A, B, C of  A q'' + B q' + C q = 0,  q = (gamma, psi), about trailing
    straight. Formulas from derivations/rigid_hawser.py (section 4); dimensional SI units.

    Water forces near head-on, to first order (E01): Y = Y_v v + Y_r r, N = N_v v + N_r r, X = -F_H0, with
        Y_v, N_v from the head-on slopes of C2C and C6C (Clarke et al. 1982, eqs. 3.28, 3.31),
        Y_r, N_r from the 'wing' yaw-rate terms of eq. (4.37) (the other yaw-rate terms are second order).
    simos_annex_b=True drops the lift caused by the sideways velocity (L_H gamma' - x_A psi') from the damping,
    as Simos's eq. B.13 does (he writes the lift with the heading psi only). Used only to compare with his Fig. 5.3.
    Also covered (step 8, D31): a steady force along the current line (wind or a hold-back tug, sin beta_W = 0), and
    the stern thruster while it is below its limit (near straight it always is).
    """
    if s.C_W != 0.0 and abs(np.sin(s.beta_W)) > 1e-12:
        raise ValueError("linear_matrices: wind with a sideways part has no straight equilibrium; use jacobian_star")
    h, f, U = s.hull, s.fluid, s.U
    rho, L, T, B = f.rho, h.L, h.T, h.B
    dC2, dC6 = hydro.clarke_slopes(h)  # dC2C/dalpha, dC6C/dalpha at alpha = 180 deg (both < 0)
    Y_v = 0.5 * rho * L * T * U * dC2  # dF_Y/dv  [N s/m]
    N_v = 0.5 * rho * L**2 * T * U * dC6  # dN/dv    [N s]
    Y_r = (0.25 * rho * np.pi * T**2 * L * (1 - 4.4 * B / L + 0.160 * B / T) + h.M11) * U  # dF_Y/dr [N s]
    N_r = (-(1 / 8) * rho * np.pi * T**2 * L**2 * (1 + 0.16 * B / T - 2.2 * B / L) + h.M26) * U  # dN/dr [N m s]
    F_H0 = s.C_F * 0.5 * rho * U**2 * L * T  # hawser tension = friction drag when trailing straight [N]

    m22 = h.M + h.M22
    m26 = h.M * h.x_CG + h.M26
    i66 = h.I_Zo + h.M66
    L_H, x_A = s.L_H, s.x_A
    A = np.array([[L_H * m22, m26 - x_A * m22],
                  [L_H * m26, i66 - m26 * x_A]])
    if simos_annex_b:
        Bm = np.array([[0.0, (h.M11 - h.M22) * U - Y_r],
                       [0.0, -N_r]])
    else:
        Bm = np.array([[-L_H * Y_v, (h.M11 - h.M22) * U - Y_r + Y_v * x_A],
                       [-L_H * N_v, -N_r + N_v * x_A]])
    # A steady force along the current line (hold-back tug or head-on wind): F_w at x_W, from the direction beta_W,
    # c = cos(beta_W) = +1 (pulls downstream) or -1. It adds c F_w to the hawser tension, and once the ship turns by
    # psi it gets a sideways part Y_w = c F_w psi (moment x_W Y_w) about midship.
    c = np.cos(s.beta_W)
    F_w = s.C_W * 0.5 * rho * U**2 * L * T
    F_H = F_H0 + c * F_w  # hawser tension when trailing straight [N]
    x_W = s.x_W_L * L
    C = np.array([[F_H, F_H + U * Y_v - c * F_w],
                  [F_H * x_A, F_H * x_A + U * N_v - x_W * c * F_w]])
    if control.is_on(s):
        # Stern thruster, linear part (D31): Y_T = N_c / x_T with N_c = -K (k_p psi + k_d (L/U) psi'), K = (1/2) rho
        # U^2 L^2 T. Moved to the left-hand side it adds to the psi column of C (k_p) and of B (k_d).
        K = 0.5 * rho * U**2 * L**2 * T
        x_T = s.x_T_L * L
        lever = np.array([1.0 / x_T, 1.0])  # sway gets Y_T, yaw gets x_T Y_T = N_c
        C[:, 1] += K * s.k_p * lever
        Bm[:, 1] += K * s.k_d * (L / U) * lever
    return A, Bm, C


def linear_jacobian_star(s, simos_annex_b=False):
    """Exact 4x4 Jacobian at trailing straight, in t* units (eigenvalues in units of U/L)."""
    A, Bm, C = linear_matrices(s, simos_annex_b)
    Ainv = np.linalg.inv(A)
    T = s.time_unit
    # q'' = -A^-1 C q - A^-1 B q'.  In t*: d/dt* = T d/dt, so the C block gets T^2 and the B block gets T.
    J = np.zeros((4, 4))
    J[0:2, 2:4] = np.eye(2)
    J[2:4, 0:2] = -Ainv @ C * T**2
    J[2:4, 2:4] = -Ainv @ Bm * T
    return J


# --- Simulation and drawing helpers -----------------------------------------------------------------------------------

def simulate(s, state0_star, t_end_star, n_out=2001, rtol=1e-9, atol=1e-11):
    """Integrate rhs_star from t* = 0 to t_end_star. Returns t* (n_out,) and states (4, n_out).

    Tight tolerances as the brief asks: near a bifurcation transients are slow and errors would pile up.
    """
    t_eval = np.linspace(0.0, t_end_star, n_out)
    sol = solve_ivp(rhs_star, (0.0, t_end_star), np.asarray(state0_star, dtype=float), args=(s,),
                    t_eval=t_eval, rtol=rtol, atol=atol, method="RK45")
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.t, sol.y


def positions(state, s):
    """Top-view positions in units of L (earth axes: X upstream, Y to port; buoy at the origin).

    Returns a dict with the attachment point A, midship O, bow and stern, each an (X, Y) pair.
    """
    gamma, psi = state[0], state[1]
    e_x = np.array([np.cos(psi), np.sin(psi)])
    A = s.L_H_L * np.array([-np.cos(gamma), np.sin(gamma)])
    O = A - s.x_A_L * e_x
    return {"A": A, "O": O, "bow": O + 0.5 * e_x, "stern": O - 0.5 * e_x}
