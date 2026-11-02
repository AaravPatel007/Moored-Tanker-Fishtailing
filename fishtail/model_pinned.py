"""Bare-hull yaw check (experiment E01): the hull is held by a vertical pin at MIDSHIP and can only turn.

Sign convention (Simos 2001 thesis, used everywhere in this package):
    x forward, y to PORT, z up, origin at MIDSHIP; heading psi = 0 means the bow points straight into the current,
    psi > 0 means the bow has turned to port. Velocities u, v are relative to the water.

There is no rope here: the pin takes all the sideways and fore-aft force, so only the turning motion is left.
States: psi (heading) [rad] and r (yaw rate) [rad/s]. Two states, so the motion can be drawn as a phase portrait.

Physics:
    midship does not move, so the water passes it at the current speed U:   u = U cos psi,  v = -U sin psi
    yaw equation, third row of (3.9). Its coupling term (M x_CG + M26)(dv/dt + u r) is zero here,
    because dv/dt + u r is the sideways acceleration of midship, and midship is held still:
        (I_Zo + M66) dr/dt = N(u, v, r)
"""

import numpy as np

from fishtail.hydro import current_forces


def rhs(t, state, hull, fluid, U):
    """Right-hand side d(psi, r)/dt of the pinned hull (t is unused; the system is autonomous)."""
    psi, r = state
    u, v = U * np.cos(psi), -U * np.sin(psi)
    # C0 (friction) only changes the fore-aft force, which the pin carries; it does not enter the yaw moment.
    _, _, N = current_forces(u, v, r, hull, fluid, C0=0.0)
    return np.array([r, N / (hull.I_Zo + hull.M66)])


def jacobian(state, hull, fluid, U, h=1e-7):
    """2x2 Jacobian of rhs at `state`, by central finite differences.

    Step sizes are scaled so that they are small compared with 1 rad and with the yaw-rate scale U/L.
    """
    steps = [h, h * U / hull.L]
    J = np.zeros((2, 2))
    for k in range(2):
        dx = np.zeros(2)
        dx[k] = steps[k]
        J[:, k] = (rhs(0, state + dx, hull, fluid, U) - rhs(0, state - dx, hull, fluid, U)) / (2 * steps[k])
    return J
