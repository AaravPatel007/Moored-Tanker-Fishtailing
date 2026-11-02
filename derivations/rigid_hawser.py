"""Step 3 derivation (sympy): the ship on a RIGID hawser -> 4 states (gamma, psi and their rates).

Run from the project root:
    /home/aarav-patel/wec_rl/.venv/bin/python -m derivations.rigid_hawser
It prints each result and writes them all to derivations/rigid_hawser_results.md.

What it does, in order:
  1. Geometry -> velocity of the hull relative to the water. Must give Simos (2001) eq. 5.2 (= B.2).
  2. Direction of the hawser pull -> force and moment on the ship. Settles open issue 4 (eq. 5.3 vs eq. B.3).
  3. Equations of motion (Simos eq. 5.1, Munk term left out because it is inside C6C, decision D4)
     -> a 3x3 LINEAR system for the unknowns (gamma_ddot, psi_ddot, F_H).
     fishtail/model_rigid.py codes exactly this system; tests/test_rigid.py checks that the two agree.
  4. Linearisation about "trailing straight" (gamma = psi = 0) -> A q'' + B q' + C q = 0 with q = (gamma, psi),
     to compare with Simos Annex B (eq. B.18), and the pitchfork condition det C = 0 <=> x_A = x_Tcr.

Sign convention (Simos 2001 thesis, used everywhere in this project):
    body axes at MIDSHIP, x forward, y to PORT; u, v = hull velocity RELATIVE TO THE WATER; r = yaw rate.
Earth axes (top view), our choice for this derivation:
    X points UPSTREAM (into the current), Y to port of a ship that points upstream. The buoy is at the origin.
    The water flows in the -X direction at speed U.
    psi   = heading: angle from +X to the bow, positive = bow turned to port (counter-clockwise seen from above).
    gamma = hawser angle: angle between the hawser and the current line. The attachment point A is at
            L_H (-cos gamma, sin gamma): gamma > 0 means the ship has swung to the +Y (port) side.
"""

from pathlib import Path

import sympy as sp

t = sp.symbols("t")
U, L_H, x_A = sp.symbols("U L_H x_A", positive=True)  # current speed, hawser length, attachment point
M, M11, M22 = sp.symbols("M M_11 M_22", positive=True)  # mass, added masses in surge and sway
m26, I66 = sp.symbols("m_26 I_66", positive=True)  # m26 = M x_CG + M26,  I66 = I_Zo + M66 (thesis eq. 3.9)
F_H = sp.symbols("F_H")  # hawser force, tension positive
X_R, Y_R, N_R = sp.symbols("X_R Y_R N_R")  # water forces on the hull (any model; here Leite + Simos)

gamma = sp.Function("gamma")(t)
psi = sp.Function("psi")(t)

# Plain symbols for the states, used once the time derivatives have been taken
g, p, gd, pd, gdd, pdd = sp.symbols("gamma psi gamma_dot psi_dot gamma_ddot psi_ddot")


def to_symbols(expr):
    """Replace gamma(t), d gamma/dt, ... by the plain symbols g, gd, gdd (highest derivatives first)."""
    return expr.subs({gamma.diff(t, 2): gdd, psi.diff(t, 2): pdd,
                      gamma.diff(t): gd, psi.diff(t): pd,
                      gamma: g, psi: p})


# --- 1. Kinematics: velocity of midship relative to the water, in body axes -----------------------------------------

def kinematics():
    """Return u, v, r (as functions of gamma(t), psi(t)) and the earth-frame unit vectors of the body axes."""
    e_x = sp.Matrix([sp.cos(psi), sp.sin(psi)])  # body x axis (towards the bow) in earth axes
    e_y = sp.Matrix([-sp.sin(psi), sp.cos(psi)])  # body y axis (to port)
    A = L_H * sp.Matrix([-sp.cos(gamma), sp.sin(gamma)])  # attachment point; the buoy is at the origin
    O = A - x_A * e_x  # midship is x_A behind the attachment point, along the ship
    water = sp.Matrix([-U, 0])  # current velocity
    V_rel = O.diff(t) - water  # velocity of midship relative to the water, earth axes
    u = sp.simplify(V_rel.dot(e_x))
    v = sp.simplify(V_rel.dot(e_y))
    r = psi.diff(t)
    return u, v, r, e_x, e_y, A


# --- 2. Hawser force on the ship ------------------------------------------------------------------------------------

def hawser_loads(e_x, e_y, A):
    """Force (F_X, F_Y) in body axes and moment N about midship from a hawser with tension F_H.

    A taut hawser pulls the attachment point A straight towards the buoy, i.e. along -A / |A|.
    """
    pull = F_H * (-A / L_H)  # earth axes
    FX = sp.simplify(pull.dot(e_x))
    FY = sp.simplify(pull.dot(e_y))
    N = x_A * FY  # moment about midship of a force applied at (x_A, 0): x F_Y - y F_X with y = 0
    return FX, FY, N


# --- 3. Equations of motion -> 3x3 linear system --------------------------------------------------------------------

def linear_system():
    """Equations (5.1) with the rigid-hawser kinematics, written as  Mat @ [gamma_ddot, psi_ddot, F_H] = rhs.

    Left-hand sides of (5.1) (added-mass terms included, Munk term not added - it is inside C6C, decision D4):
        surge: (M + M11) du/dt - (M + M22) v r - m26 r^2
        sway : (M + M22) dv/dt + (M + M11) u r + m26 dr/dt
        yaw  : I66 dr/dt + m26 (dv/dt + u r)
    Right-hand sides: water force (X_R, Y_R, N_R) + hawser load (F_X, F_Y, N).
    """
    u, v, r, e_x, e_y, A = kinematics()
    FX, FY, N = hawser_loads(e_x, e_y, A)
    du, dv, dr = u.diff(t), v.diff(t), r.diff(t)
    surge = (M + M11) * du - (M + M22) * v * r - m26 * r**2 - (X_R + FX)
    sway = (M + M22) * dv + (M + M11) * u * r + m26 * dr - (Y_R + FY)
    yaw = I66 * dr + m26 * (dv + u * r) - (N_R + N)
    eqs = [sp.expand(to_symbols(e)) for e in (surge, sway, yaw)]
    Mat, rhs = sp.linear_eq_to_matrix(eqs, [gdd, pdd, F_H])
    return sp.simplify(Mat), sp.simplify(rhs)


# --- 4. Linearisation about trailing straight -----------------------------------------------------------------------

def linearisation():
    """Linear equations A q'' + B q' + C q = 0, q = (gamma, psi), about gamma = psi = 0 (u = U, v = r = 0).

    Water forces near this point, to first order (E01: the static coefficients have Clarke's slopes there, and every
    yaw-rate term except the 'wing' term is of second order):
        X_R = X0 (the friction drag, constant),   Y_R = Y_v v + Y_r r,   N_R = N_v v + N_r r.
    The hawser tension is then F_H = -X0 + (second-order terms), so in the sway and yaw equations we may put
    F_H = F_H0 = -X0: a change of F_H only enters multiplied by sin(gamma + psi), which is already small.
    """
    Y_v, Y_r, N_v, N_r, F_H0 = sp.symbols("Y_v Y_r N_v N_r F_H0")
    u, v, r, e_x, e_y, A = kinematics()
    FX, FY, N = hawser_loads(e_x, e_y, A)
    du, dv, dr = u.diff(t), v.diff(t), r.diff(t)
    sway = (M + M22) * dv + (M + M11) * u * r + m26 * dr - (Y_v * v + Y_r * r + FY)
    yaw = I66 * dr + m26 * (dv + u * r) - (N_v * v + N_r * r + N)
    eqs = [to_symbols(e).subs(F_H, F_H0) for e in (sway, yaw)]

    zero = {g: 0, p: 0, gd: 0, pd: 0, gdd: 0, pdd: 0}
    A_mat = sp.Matrix([[sp.diff(e, s).subs(zero) for s in (gdd, pdd)] for e in eqs])
    B_mat = sp.Matrix([[sp.diff(e, s).subs(zero) for s in (gd, pd)] for e in eqs])
    C_mat = sp.Matrix([[sp.diff(e, s).subs(zero) for s in (g, p)] for e in eqs])
    return sp.simplify(A_mat), sp.simplify(B_mat), sp.simplify(C_mat), (Y_v, Y_r, N_v, N_r, F_H0)


def main():
    lines = ["# Rigid-hawser model: sympy results", "",
             "Written by `derivations/rigid_hawser.py` (do not edit by hand). Symbols: `m_26 = M x_CG + M26`, "
             "`I_66 = I_Zo + M66`; `X_R, Y_R, N_R` = water force and moment on the hull.", ""]

    def show(title, expr):
        text = f"**{title}**\n\n```\n{sp.pretty(expr, use_unicode=True, wrap_line=False)}\n```\n"
        print(title)
        print(sp.pretty(expr, use_unicode=True, wrap_line=False))
        print()
        lines.append(text)

    # 1. Kinematics vs Simos (5.2)
    u, v, r, e_x, e_y, A = kinematics()
    u_simos = U * sp.cos(psi) + L_H * gamma.diff(t) * sp.sin(gamma + psi)
    v_simos = -U * sp.sin(psi) + L_H * gamma.diff(t) * sp.cos(gamma + psi) - x_A * psi.diff(t)
    show("1. u (relative to the water)", to_symbols(u))
    show("   v (relative to the water)", to_symbols(v))
    ok_u = sp.simplify(sp.expand_trig(u - u_simos)) == 0
    ok_v = sp.simplify(sp.expand_trig(v - v_simos)) == 0
    lines.append(f"Same as Simos eq. 5.2 / B.2: u {'YES' if ok_u else 'NO'}, v {'YES' if ok_v else 'NO'}.\n")
    print("matches Simos (5.2):", ok_u, ok_v, "\n")

    # 2. Hawser loads vs (5.3) and (B.3)
    FX, FY, N = hawser_loads(e_x, e_y, A)
    show("2. Hawser force along the ship F_X", to_symbols(FX))
    show("   Hawser force sideways F_Y", to_symbols(FY))
    show("   Hawser moment about midship N", to_symbols(N))
    plus = sp.simplify(sp.expand_trig(FY - F_H * sp.sin(gamma + psi))) == 0
    minus = sp.simplify(sp.expand_trig(FY + F_H * sp.sin(gamma + psi))) == 0
    verdict = (f"F_Y = +F_H sin(gamma+psi) (eq. 5.3): {'YES' if plus else 'NO'};  "
               f"F_Y = -F_H sin(gamma+psi) (eq. B.3): {'YES' if minus else 'NO'}.")
    lines.append(verdict + " **Open issue 4 settled: B.3 is right, the '+' in 5.3 is a misprint.**\n")
    print(verdict, "\n")

    # 3. The 3x3 system solved at every time step
    Mat, rhs = linear_system()
    show("3. Matrix of the linear system  Mat @ [gamma_ddot, psi_ddot, F_H] = rhs", Mat)
    show("   right-hand side rhs", rhs)

    # 4. Linearisation
    A_mat, B_mat, C_mat, (Y_v, Y_r, N_v, N_r, F_H0) = linearisation()
    show("4. Linearised about gamma = psi = 0:  A q'' + B q' + C q = 0.   Mass matrix A", A_mat)
    show("   Damping matrix B", B_mat)
    show("   Stiffness matrix C", C_mat)
    detC = sp.factor(C_mat.det())
    show("   det C", detC)
    lines.append("With F_H0 > 0 (tension) and U > 0, det C > 0 exactly when N_v - x_A Y_v > 0. Writing "
                 "Y_v = (1/2) rho L T U C2', N_v = (1/2) rho L^2 T U C6' (C2', C6' = slopes of the static "
                 "coefficients at head-on, both negative) this is x_A > L C6'/C2' = x_Tcr (eq. 3.37). "
                 "det C changes sign there, so one real eigenvalue crosses zero: **the pitchfork is at "
                 "x_A = x_Tcr for every hawser length and friction coefficient.**\n")
    lines.append("Compared with Simos eq. B.18: the mass matrix and the stiffness matrix have the same structure. "
                 "Our damping matrix has extra terms -Y_v L_H (first column) and +Y_v x_A: the sideways lift "
                 "force caused by the sideways velocity L_H gamma' - x_A psi' of the hull. Simos's B.13 writes the "
                 "lift with the heading psi only, so his B.15/B.18 has zeros in the first column. See "
                 "experiments/E04 and E06 for what this changes.\n")

    out = Path(__file__).with_name("rigid_hawser_results.md")
    out.write_text("\n".join(lines))
    print(f"written to {out}")


if __name__ == "__main__":
    main()
