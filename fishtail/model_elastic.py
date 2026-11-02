"""Ship on a STRETCHY hawser that can go slack: 6 states (step 7).

Sign convention (Simos 2001 thesis, used everywhere in this package):
    body axes at MIDSHIP, x forward, y to PORT; u, v = hull velocity RELATIVE TO THE WATER; r = yaw rate.
Earth axes as in model_rigid.py: X points UPSTREAM, Y to port of a ship pointing upstream, buoy at the origin,
the water flows in the -X direction at speed U.

States (SI units, time in seconds): (X, Y, psi, u, v, r)
    X, Y  = position of midship [m];  psi = heading [rad];  u, v [m/s], r [rad/s] relative to the water.

Kinematics (ground velocity = velocity relative to the water + the water's velocity (-U, 0)):
    dX/dt = u cos psi - v sin psi - U,    dY/dt = u sin psi + v cos psi,    dpsi/dt = r
Equations of motion: Simos eq. (5.1), as in the rigid model, but now the hawser force comes from a ROPE LAW instead of
being an unknown:
    strain  e = (distance buoy -> attachment point) / L_H - 1
    tension F_H = a e^n  if e > 0 (taut),   F_H = 0 if e <= 0 (slack)
a and n are fitted to a real SPM nylon hawser (Huang et al. 2026, Fig. 3; see rope_from_huang()).
The rope pulls the attachment point straight towards the buoy; it never pushes.
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from fishtail import hydro
from fishtail.params import HUANG_HAWSER, SEA_WATER, Fluid, Hull


def rope_from_huang():
    """Fit tension = a * strain^n to the Huang et al. (2026) hawser curve (log-log least squares).

    Returns (a [N], n [-]). Strain = elongation / 70 m.
    """
    e = np.array(HUANG_HAWSER["elongation_m"][1:]) / HUANG_HAWSER["length_m"]
    T = np.array(HUANG_HAWSER["tension_kN"][1:]) * 1e3
    n, log_a = np.polyfit(np.log(e), np.log(T), 1)
    return float(np.exp(log_a)), float(n)


ROPE_A, ROPE_N = rope_from_huang()


@dataclass(frozen=True)
class Elastic:
    """One set-up with a stretchy hawser. Same knobs as model_rigid.Rigid, plus the rope law."""

    hull: Hull
    x_A_L: float  # attachment point forward of midship / L
    L_H_L: float  # UNSTRETCHED hawser length / L
    C_F: float  # friction coefficient
    U: float = 1.0  # current speed [m/s] - matters now: the rope does not scale with U^2
    rope_a: float = ROPE_A  # tension at 100 % strain [N] (fitted, see rope_from_huang)
    rope_n: float = ROPE_N  # exponent of the rope law [-]
    fluid: Fluid = SEA_WATER

    @property
    def x_A(self):
        return self.x_A_L * self.hull.L

    @property
    def L_H(self):
        return self.L_H_L * self.hull.L

    @property
    def time_unit(self):
        return self.hull.L / self.U


def attachment_point(state, s):
    X, Y, psi = state[0], state[1], state[2]
    return np.array([X + s.x_A * np.cos(psi), Y + s.x_A * np.sin(psi)])


def tension(state, s):
    """Rope tension [N] and strain [-] at a state. Zero tension when the rope is slack."""
    A = attachment_point(state, s)
    strain = np.hypot(A[0], A[1]) / s.L_H - 1.0
    return (s.rope_a * strain**s.rope_n if strain > 0 else 0.0), strain


def rhs(t, state, s):
    """d(state)/dt in seconds."""
    X, Y, psi, u, v, r = state
    h = s.hull
    m11, m22 = h.M + h.M11, h.M + h.M22
    m26 = h.M * h.x_CG + h.M26
    i66 = h.I_Zo + h.M66
    ex = np.array([np.cos(psi), np.sin(psi)])
    ey = np.array([-np.sin(psi), np.cos(psi)])

    # Rope force: towards the buoy, size from the rope law
    A = attachment_point(state, s)
    dist = np.hypot(A[0], A[1])
    T, _ = tension(state, s)
    F = -T * A / dist  # earth axes
    FX, FY = F @ ex, F @ ey
    N = s.x_A * FY  # moment about midship of a force at (x_A, 0)

    X_R, Y_R, N_R = hydro.current_forces(u, v, r, h, s.fluid, s.C_F)

    # Equations (5.1): unknowns du/dt, dv/dt, dr/dt
    M = np.array([[m11, 0.0, 0.0],
                  [0.0, m22, m26],
                  [0.0, m26, i66]])
    b = np.array([X_R + FX + m22 * v * r + m26 * r**2,
                  Y_R + FY - m11 * u * r,
                  N_R + N - m26 * u * r])
    du, dv, dr = np.linalg.solve(M, b)
    return np.array([u * np.cos(psi) - v * np.sin(psi) - s.U,
                     u * np.sin(psi) + v * np.cos(psi),
                     r, du, dv, dr])


def trailing_straight(s):
    """Equilibrium: ship straight behind the buoy, rope stretched until its tension equals the friction drag."""
    F_H0 = s.C_F * 0.5 * s.fluid.rho * s.U**2 * s.hull.L * s.hull.T
    strain0 = (F_H0 / s.rope_a) ** (1.0 / s.rope_n)
    return np.array([-s.L_H * (1 + strain0) - s.x_A, 0.0, 0.0, s.U, 0.0, 0.0]), strain0


def angles(state, s):
    """Hawser angle gamma and heading psi [rad], defined as in the rigid model (for comparisons)."""
    A = attachment_point(state, s)
    return np.arctan2(A[1], -A[0]), state[2]


def from_rigid_start(gamma, psi, s):
    """State at rest with the rope just taut at the trailing-straight stretch, at hawser angle gamma and heading psi
    (same kind of start as the rigid runs)."""
    eq, strain0 = trailing_straight(s)
    ell = s.L_H * (1 + strain0)
    A = ell * np.array([-np.cos(gamma), np.sin(gamma)])
    O = A - s.x_A * np.array([np.cos(psi), np.sin(psi)])
    # at rest over the ground -> relative velocity = (U, 0) in earth axes, put into body axes
    return np.array([O[0], O[1], psi, s.U * np.cos(psi), -s.U * np.sin(psi), 0.0])


def jacobian(state, s):
    """6x6 Jacobian by central differences (valid where the rope is taut, e.g. at trailing straight)."""
    scales = np.array([s.hull.L, s.hull.L, 1.0, s.U, s.U, s.U / s.hull.L]) * 1e-6
    J = np.zeros((6, 6))
    for k in range(6):
        dx = np.zeros(6)
        dx[k] = scales[k]
        J[:, k] = (rhs(0, state + dx, s) - rhs(0, state - dx, s)) / (2 * scales[k])
    return J


def eigenvalues_star(s):
    """Eigenvalues of trailing straight in units of U/L (largest real part first)."""
    eq, _ = trailing_straight(s)
    eig = np.linalg.eigvals(jacobian(eq, s)) * s.time_unit
    return eig[np.argsort(-eig.real)]


def simulate(s, state0, t_end_star, n_out=3001, rtol=1e-8, atol=None):
    """Integrate from state0 for t_end_star units of L/U. Returns t* (n,), states (6, n), tension [N] (n,).

    max_step keeps the solver from stepping over a short slack-taut switch.
    """
    T = s.time_unit
    if atol is None:
        atol = np.array([1e-6 * s.hull.L, 1e-6 * s.hull.L, 1e-9, 1e-9 * s.U, 1e-9 * s.U, 1e-12 * s.U / s.hull.L])
    t_eval = np.linspace(0.0, t_end_star * T, n_out)
    sol = solve_ivp(rhs, (0.0, t_end_star * T), state0, args=(s,), t_eval=t_eval, rtol=rtol, atol=atol,
                    max_step=0.05 * T)
    if not sol.success:
        raise RuntimeError(sol.message)
    tens = np.array([tension(sol.y[:, i], s)[0] for i in range(sol.t.size)])
    return sol.t / T, sol.y, tens
