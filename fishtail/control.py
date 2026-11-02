"""Heading control of the moored tanker with its own stern thruster (step 8, experiment E08, decision D31).

Sign convention (Simos 2001, as in the whole package): body axes at midship, x forward, y to PORT;
psi > 0 = bow turned to port; r = psi_dot > 0 = turning to port; a yaw moment N > 0 turns the bow to port.

The control law is a PD ("proportional + derivative") law on the heading, in non-dimensional form:
    N_c* = -k_p psi - k_d r*,      N_c* = N_c / ((1/2) rho U^2 L^2 T),   r* = r L / U.
  k_p: a spring back to straight (more moment the further the bow is off),
  k_d: a damper on the turning (more moment the faster the ship turns).
The thruster is a SIDE FORCE Y_T at x_T (the stern, x_T < 0), so it gives the moment N_c = x_T Y_T and also
pushes the ship sideways:  Y_T = N_c / x_T.
Saturation: the thruster cannot push harder than |Y_T| <= thrust_limit * F_H0, with F_H0 = C_F (1/2) rho U^2 L T the
hawser pull when trailing straight. Beyond that the force is simply clipped (held at the limit).
With k_p = k_d = 0 the thruster does nothing and the model is exactly the one of steps 3-6.
"""

import numpy as np


def is_on(s):
    """True when the controller has a non-zero gain."""
    return s.k_p != 0.0 or s.k_d != 0.0


def force_scale(s):
    """(1/2) rho U^2 L T [N]: the scale of all the water forces."""
    return 0.5 * s.fluid.rho * s.U**2 * s.hull.L * s.hull.T


def thruster_force(psi, r, s):
    """Side force of the stern thruster Y_T [N] for heading psi [rad] and yaw rate r [rad/s] (PD law, then clip)."""
    r_star = r * s.hull.L / s.U  # yaw rate in units of U/L
    N_c = -(s.k_p * psi + s.k_d * r_star) * force_scale(s) * s.hull.L  # wanted yaw moment [N m]
    Y_T = N_c / (s.x_T_L * s.hull.L)  # side force that gives it [N]
    Y_max = s.thrust_limit * s.C_F * force_scale(s)  # thrust limit [N]
    return np.clip(Y_T, -Y_max, Y_max)


def thruster_side_coefficient(psi, s):
    """Thruster force at REST (r = 0) with heading psi, divided by (1/2) rho U^2 L T; used for the equilibria."""
    return thruster_force(psi, 0.0, s) / force_scale(s)


def decay_rate(eig):
    """How fast a small disturbance dies out: sigma = -(largest real part) [U/L]. sigma < 0 means it grows."""
    return float(-np.max(np.real(eig)))


def thruster_ratio(state_star, s):
    """Thruster force as a fraction of F_H0 at a non-dimensional state (gamma, psi, gamma', psi'); for plots."""
    r = state_star[3] * s.U / s.hull.L
    return thruster_force(state_star[1], r, s) / (s.C_F * force_scale(s))
