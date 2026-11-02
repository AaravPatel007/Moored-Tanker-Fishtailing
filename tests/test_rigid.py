"""Step 3 checks: the rigid-hawser model (fishtail/model_rigid.py) against its sympy derivation and basic physics."""

from dataclasses import replace

import numpy as np
import pytest
import sympy as sp

from derivations import rigid_hawser as der
from fishtail import hydro
from fishtail import model_rigid as mr
from fishtail.params import VIDAL_100, VIDAL_40

SETUPS = [mr.Rigid(VIDAL_100, x_A_L=0.5, L_H_L=0.5, C_F=0.02),
          mr.Rigid(VIDAL_40, x_A_L=0.4, L_H_L=0.36, C_F=0.035),
          mr.Rigid(VIDAL_100, x_A_L=0.25, L_H_L=1.0, C_F=0.015)]
STATES = [np.array([0.3, -0.1, 0.002, -0.001]), np.array([-0.6, 0.25, -0.004, 0.003]),
          np.array([1.0, 0.5, 0.0, 0.002])]  # dimensional: rad, rad, rad/s, rad/s


@pytest.fixture(scope="module")
def sympy_system():
    """The 3x3 system from the sympy derivation, turned into numerical functions (done once: sympy is slow)."""
    Mat, rhs = der.linear_system()
    args = [der.g, der.p, der.gd, der.pd, der.U, der.L_H, der.x_A, der.M, der.M11, der.M22, der.m26, der.I66,
            der.X_R, der.Y_R, der.N_R]
    return sp.lambdify(args, Mat, "numpy"), sp.lambdify(args, rhs, "numpy")


@pytest.mark.parametrize("s", SETUPS)
@pytest.mark.parametrize("state", STATES)
def test_code_matches_sympy_derivation(sympy_system, s, state):
    """The hand-coded equations in model_rigid.accelerations() equal the sympy-derived ones."""
    Mat_f, rhs_f = sympy_system
    h = s.hull
    u, v, r = mr.velocities(state, s)
    X_R, Y_R, N_R = hydro.current_forces(u, v, r, h, s.fluid, s.C_F)
    args = [*state, s.U, s.L_H, s.x_A, h.M, h.M11, h.M22, h.M * h.x_CG + h.M26, h.I_Zo + h.M66, X_R, Y_R, N_R]
    expected = np.linalg.solve(np.array(Mat_f(*args), dtype=float), np.array(rhs_f(*args), dtype=float).ravel())
    assert np.allclose(mr.accelerations(state, s), expected, rtol=1e-10, atol=1e-14)


def test_sympy_settles_hawser_sign():
    """Open issue 4: the hawser's sideways force is F_Y = -F_H sin(gamma + psi) (Simos eq. B.3, not 5.3)."""
    u, v, r, e_x, e_y, A = der.kinematics()
    FX, FY, N = der.hawser_loads(e_x, e_y, A)
    assert sp.simplify(sp.expand_trig(FY + der.F_H * sp.sin(der.gamma + der.psi))) == 0
    assert sp.simplify(sp.expand_trig(N + der.x_A * der.F_H * sp.sin(der.gamma + der.psi))) == 0


@pytest.mark.parametrize("s", SETUPS)
def test_trailing_straight_is_an_equilibrium(s):
    """gamma = psi = 0 at rest: nothing moves, and the hawser carries exactly the friction drag."""
    assert np.allclose(mr.rhs_star(0, np.zeros(4), s), 0.0, atol=1e-14)
    assert mr.hawser_force_star(np.zeros(4), s) == pytest.approx(s.C_F)


@pytest.mark.parametrize("s", SETUPS)
def test_jacobian_finite_difference_matches_exact(s):
    """Step 3 check from the brief: the finite-difference Jacobian equals the exact (sympy) linearisation."""
    J_fd = mr.jacobian_star(np.zeros(4), s)
    J_exact = mr.linear_jacobian_star(s)
    assert np.allclose(J_fd, J_exact, rtol=1e-4, atol=1e-5 * np.abs(J_exact).max())


@pytest.mark.parametrize("s", SETUPS)
def test_current_speed_drops_out_at_fixed_friction(s):
    """Decision D1: with C_F held fixed, the equations in t* = U t / L are the same for every current speed."""
    x = np.array([0.4, -0.2, 0.05, 0.03])
    reference = mr.rhs_star(0, x, s)
    for U in [0.2, 0.7, 2.5]:
        assert np.allclose(mr.rhs_star(0, x, replace(s, U=U)), reference, rtol=1e-10, atol=1e-14)


@pytest.mark.parametrize("s", SETUPS)
def test_port_starboard_mirror(s):
    """Mirror image of a motion is also a motion: flipping the signs of all four states flips the rates."""
    for x in STATES:
        x_star = np.array([x[0], x[1], x[2] * s.time_unit, x[3] * s.time_unit])
        assert np.allclose(mr.rhs_star(0, -x_star, s), -mr.rhs_star(0, x_star, s), rtol=1e-10, atol=1e-14)


def test_kinematics_from_simulated_positions():
    """Independent check of eq. (5.2): differentiate the simulated midship position numerically and compare with
    the relative velocity used by the model (relative velocity = ground velocity minus current)."""
    s = SETUPS[0]
    t, y = mr.simulate(s, [0.3, 0.1, 0.0, 0.0], t_end_star=2.0, n_out=4001)
    dt = t[1] - t[0]
    O = np.array([mr.positions(y[:, i], s)["O"] for i in range(t.size)])  # in units of L
    V_ground = np.gradient(O, dt, axis=0)  # in units of U (since time is t* = U t / L)
    for i in [500, 2000, 3500]:
        psi = y[1, i]
        e_x = np.array([np.cos(psi), np.sin(psi)])
        e_y = np.array([-np.sin(psi), np.cos(psi)])
        V_rel = V_ground[i] - np.array([-1.0, 0.0])  # the water moves at -X with unit speed
        u, v, _ = mr.velocities(mr.to_dimensional(y[:, i], s), s)
        assert V_rel @ e_x == pytest.approx(u / s.U, abs=1e-5)
        assert V_rel @ e_y == pytest.approx(v / s.U, abs=1e-5)


@pytest.mark.parametrize("hull", [VIDAL_100, VIDAL_40])
def test_pitchfork_at_critical_attachment_point(hull):
    """det C changes sign at x_A = x_Tcr (eq. 3.37), for any hawser length and friction: one real eigenvalue of the
    trailing-straight Jacobian passes through zero there."""
    x_cr = hydro.critical_turret_position(hull)
    for L_H_L, C_F in [(0.3, 0.02), (1.0, 0.04)]:
        below = mr.Rigid(hull, x_cr - 0.02, L_H_L, C_F)
        above = mr.Rigid(hull, x_cr + 0.02, L_H_L, C_F)
        at = mr.Rigid(hull, x_cr, L_H_L, C_F)
        assert np.linalg.det(mr.linear_matrices(below)[2]) < 0 < np.linalg.det(mr.linear_matrices(above)[2])
        eig = np.linalg.eigvals(mr.linear_jacobian_star(at))
        assert np.min(np.abs(eig)) < 1e-6
