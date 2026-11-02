"""Water (current) forces on the hull: Leite et al. (1998) static model + Simos (2001) yaw-rate extension.

Sign convention (Simos 2001 thesis, used everywhere in this package):
    x forward, y to PORT, z up, origin at MIDSHIP.
    u, v = velocity of the hull RELATIVE TO THE WATER, in body axes; r = yaw rate (positive turns the bow to port).
    F_X, F_Y = force on the hull [N]; N = yaw moment about midship [N m], positive turns the bow to port.
    Angle of attack alpha: 180 deg = current from dead ahead, 90 deg = current hitting the starboard side.

Equation numbers are those of the Simos (2001) thesis; summary in notes/02_hydrodynamic_model.md.
The Munk moment is already inside C6C (eq. 3.29), so it is never added a second time.
"""

import numpy as np


# --- Static part: the hull is not turning (r = 0) --------------------------------------------------------------------

def angle_of_attack(u, v):
    """Angle of attack alpha [rad] in [0, 2 pi), eq. (3.13).

    (3.13) is written alpha = pi + arctan(v/u). arctan2 gives the same answer when u > 0, and also gives the right
    angle when the hull moves backwards through the water (u < 0, current from astern -> alpha = 0).
    """
    return np.mod(np.pi + np.arctan2(v, u), 2 * np.pi)


def friction_coefficient(U, hull, fluid):
    """Surge friction coefficient C0 of eq. (3.20); Simos calls it C_F(Re) in Ch. 5.

    ITTC-1957 friction line 0.075/(log10 Re - 2)^2, scaled from wetted area S to the reference area T L,
    and multiplied by (1 + k) for the extra drag caused by the 3-D hull shape.
    """
    Re = U * hull.L / fluid.nu  # Reynolds number
    return 0.075 / (np.log10(Re) - 2) ** 2 * hull.S / (hull.T * hull.L) * (1 + hull.k_form)


def C1C(alpha, hull, C0):
    """Surge force coefficient, eq. (3.14): friction + small-aspect-ratio wing terms."""
    TL = hull.T / hull.L
    return C0 * np.cos(alpha) + (np.pi * TL / 8) * (np.cos(3 * alpha) - np.cos(alpha))


def C2C(alpha, hull):
    """Sway force coefficient, eq. (3.21)."""
    TL = hull.T / hull.L
    s, c = np.sin(alpha), np.cos(alpha)
    cross_flow = (hull.C_Y - np.pi * TL / 2) * s * np.abs(s)  # sideways drag, strongest at alpha = 90 deg
    induced_drag = (np.pi * TL / 2) * s**3  # wing induced drag
    lift = np.pi * TL * (1 + 0.4 * hull.C_B * hull.B / hull.T) * s * np.abs(c)  # wing lift, Clarke et al. (1982)
    return cross_flow + induced_drag + lift


def C6C_terms(alpha, hull):
    """The three parts of the yaw moment coefficient, eq. (3.29): (cross_flow, munk, lift)."""
    TL = hull.T / hull.L
    s, c = np.sin(alpha), np.cos(alpha)
    cross_flow = -hull.lp_L * hull.C_Y * s * np.abs(s)  # sideways drag acting l_p aft of midship
    munk = -np.pi * TL * s * c  # Munk moment: turns the hull broadside to the flow
    lift = -((1 + np.abs(c)) / 2) ** 2 * np.pi * TL * (0.5 - 2.4 * TL) * s * np.abs(c)  # hull lift + rudder
    return cross_flow, munk, lift


def C6C(alpha, hull):
    """Yaw moment coefficient about midship, eq. (3.29)."""
    cross_flow, munk, lift = C6C_terms(alpha, hull)
    return cross_flow + munk + lift


def static_forces(u, v, hull, fluid, C0):
    """Forces when the hull is not turning, eq. (3.12): returns F_X [N], F_Y [N], N [N m]."""
    alpha = angle_of_attack(u, v)
    q = 0.5 * fluid.rho * (u**2 + v**2) * hull.L * hull.T  # (1/2) rho U^2 L T
    return q * C1C(alpha, hull, C0), q * C2C(alpha, hull), q * hull.L * C6C(alpha, hull)


# --- Sideways drag along the hull ----------------------------------------------------------------------------------

def sectional_drag(x_L, hull):
    """Sideways drag coefficient C_D of the hull slice at x/L (-0.5 = stern, +0.5 = bow).

    Our simplification (decision D6): a straight line, chosen so that it matches the two measured numbers exactly:
        integral of C_D d(x/L)          = C_Y               (total sideways drag)
        integral of C_D (x/L) d(x/L)    = -(l_p/L) C_Y      (it acts l_p AFT of midship)
    """
    return hull.C_Y * (1 - 12 * hull.lp_L * x_L)


def drag_integral(j, v, r, hull):
    """The integral I_j(v, r) of eq. (4.38), computed exactly (j = 0, 1, 2, 3).

    I_j = integral from -1/2 to 1/2 of C_D(xi) * sign(v + r L xi) * xi^j  d(xi),   with xi = x/L.

    v + r L xi is the sideways speed of the slice at xi; its sign says from which side the water crosses that slice.
    It changes sign at most once, at xi0 = -v/(r L). On each side of xi0 the integrand is a polynomial in xi,
    which we integrate by hand (primitive below).
    """
    slope = -12 * hull.lp_L  # C_D(xi) = C_Y (1 + slope * xi)

    def primitive(xi):  # integral of C_Y (1 + slope xi) xi^j
        return hull.C_Y * (xi ** (j + 1) / (j + 1) + slope * xi ** (j + 2) / (j + 2))

    points = [-0.5, 0.5]  # the stern and the bow
    if r != 0:
        xi0 = -v / (r * hull.L)  # slice with zero sideways speed
        if -0.5 < xi0 < 0.5:
            points = [-0.5, xi0, 0.5]

    total = 0.0
    for left, right in zip(points[:-1], points[1:]):
        middle = 0.5 * (left + right)
        side = np.sign(v + r * hull.L * middle)  # +1, -1 (or 0 if the hull is not moving sideways at all)
        total += side * (primitive(right) - primitive(left))
    return total


# --- Yaw-rate part: extra forces because the hull is turning (Simos 2001, eq. 4.37) ---------------------------------

def yaw_rate_terms(u, v, r, hull, fluid):
    """Extra forces caused by the yaw rate r, split into their physical parts, eq. (4.37).

    Returns a dict (forces in N, moments in N m). Every part is zero when r = 0, so nothing is counted twice with
    the static forces above (thesis p. 126).
    """
    rho, L, T, B = fluid.rho, hull.L, hull.T, hull.B
    lpCY = hull.lp_L * hull.C_Y  # (4.37) prints "l_p C_Y"; the units only work if l_p means l_p/L (open issue 6)

    I0 = drag_integral(0, v, r, hull)
    I1 = drag_integral(1, v, r, hull)
    I2 = drag_integral(2, v, r, hull)
    I3 = drag_integral(3, v, r, hull)
    I3_turn_only = drag_integral(3, 0.0, r, hull)  # I3(0, r): the same hull turning on the spot

    # Surge, eq. (4.19), from wing theory. alpha_a = effective angle of attack of the hull seen as a wing.
    # Simos wants cos(alpha_a) to have the sign of u (p. 119), so alpha_a is the two-argument angle atan2(rL/2 - v, u)
    # (decision D12). For the moored ship u > 0 and this is the same as his formula.
    size = np.hypot(u, r * L / 2 - v)
    cos_alpha_a = u / size if size > 0 else 0.0
    X = -0.25 * rho * np.pi * T**2 * L * v * r - (1 / 16) * rho * np.pi * T**2 * L**2 * cos_alpha_a * r**2

    # Wing (lift) effects of turning, eq. (4.14): Clarke et al. (1982) Y_r and N_r. The "+ M11 u r" and "+ M26 u r"
    # take out the added-mass (ideal-fluid) part of Clarke's measured values, because that part is already on the
    # left side of the equations of motion (3.9) (thesis p. 116).
    Y_wing = 0.25 * rho * np.pi * T**2 * L * (1 - 4.4 * B / L + 0.160 * B / T) * u * r + hull.M11 * u * r
    N_wing = -(1 / 8) * rho * np.pi * T**2 * L**2 * (1 + 0.16 * B / T - 2.2 * B / L) * abs(u) * r + hull.M26 * u * r

    # Cross-flow drag of each slice with its own sideways speed v + r x, minus what the static part already has,
    # eq. (4.29). The square brackets are zero when r = 0.
    Y_cross_flow = (-0.5 * rho * T * L * (I0 * v**2 - hull.C_Y * v * abs(v))
                    - rho * T * L**2 * I1 * v * r
                    - 0.5 * rho * T * L**3 * I2 * r**2)
    N_cross_flow = (-0.5 * rho * T * L**2 * (I1 * v**2 + lpCY * v * abs(v))
                    - rho * T * L**3 * I2 * v * r
                    - 0.5 * rho * T * L**4 * I3 * r**2)

    # Vortices shed at bow and stern when the hull turns, eq. (4.33). Built so that a hull turning on the spot feels
    # the measured moment N = -(1/32) rho T L^4 C_Y r|r| (eq. 4.31, from the IPT yaw-rotation tests, Fig. 4.6).
    N_vortex = -0.5 * rho * T * L**4 * (hull.C_Y / 16 * r * abs(r) - I3_turn_only * r**2)

    return {"X": X, "Y_wing": Y_wing, "Y_cross_flow": Y_cross_flow,
            "N_wing": N_wing, "N_cross_flow": N_cross_flow, "N_vortex": N_vortex}


def yaw_rate_forces(u, v, r, hull, fluid):
    """Sum of the yaw-rate parts: returns dF_X [N], dF_Y [N], dN [N m]."""
    t = yaw_rate_terms(u, v, r, hull, fluid)
    return t["X"], t["Y_wing"] + t["Y_cross_flow"], t["N_wing"] + t["N_cross_flow"] + t["N_vortex"]


def current_forces(u, v, r, hull, fluid, C0):
    """Total water force on the hull, eq. (4.34): static part + yaw-rate part. Returns F_X [N], F_Y [N], N [N m]."""
    FX, FY, N = static_forces(u, v, hull, fluid, C0)
    dFX, dFY, dN = yaw_rate_forces(u, v, r, hull, fluid)
    return FX + dFX, FY + dFY, N + dN


# --- Small-angle checks near head-on current (alpha = 180 deg) -------------------------------------------------------

def head_on_slopes(hull, h=1e-8):
    """dC2C/dalpha and dC6C/dalpha at alpha = 180 deg (current from dead ahead), by central differences.

    The step h is very small because the cross-flow terms sin(alpha)|sin(alpha)| make the error of a central
    difference proportional to h (not h^2) at this point.
    """
    dC2 = (C2C(np.pi + h, hull) - C2C(np.pi - h, hull)) / (2 * h)
    dC6 = (C6C(np.pi + h, hull) - C6C(np.pi - h, hull)) / (2 * h)
    return dC2, dC6


def clarke_slopes(hull):
    """The same two slopes written out by hand from Clarke et al. (1982), eqs. (3.28) and (3.31)."""
    TL = hull.T / hull.L
    dC2 = -np.pi * TL * (1 + 0.4 * hull.C_B * hull.B / hull.T)
    dC6 = -np.pi * TL * (0.5 + 2.4 * TL)
    return dC2, dC6


def critical_turret_position(hull):
    """x_Tcr / L, eq. (3.37): a hull held by a pivot at x_T points steadily into the current only if x_T > x_Tcr.

    About the pivot, the water's turning moment is N - x_T F_Y. For a small heading change this stops pushing the
    bow away from the current when x_T = (dC6C/dalpha) / (dC2C/dalpha) * L.
    """
    dC2, dC6 = head_on_slopes(hull)
    return dC6 / dC2
