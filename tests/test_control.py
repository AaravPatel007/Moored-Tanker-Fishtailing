"""Control checks (step 8, decision D31): the stern thruster with PD heading control, and the hold-back tug
(a steady pull downstream at the stern, built from the wind force with beta_W = 0 and x_W = stern)."""

from dataclasses import replace

import numpy as np

from fishtail import continuation as co
from fishtail import control
from fishtail import model_rigid as mr
from fishtail.params import VIDAL_100, X_TUG_L

BASE = mr.Rigid(VIDAL_100, x_A_L=0.5, L_H_L=0.5, C_F=0.03)  # the Fishtailing tab default: it fishtails
STATE = np.array([0.3, 0.1, 0.02, -0.01])  # an arbitrary moving state (rad, rad per t*)


def tug(s, ratio):
    """Set-up s with a hold-back tug pulling ratio * F_H0 downstream from the stern."""
    return replace(s, C_W=ratio * s.C_F, beta_W=0.0, x_W_L=X_TUG_L)


def test_zero_gains_is_the_old_model():
    """k_p = k_d = 0 must give exactly the model of steps 3-6, whatever the thrust limit."""
    off = replace(BASE, thrust_limit=0.1)
    assert np.array_equal(mr.rhs_star(0, STATE, off), mr.rhs_star(0, STATE, BASE))


def test_thruster_law_and_sign():
    """Bow turned to port (psi > 0) -> the moment turns it back to starboard (N_c < 0): with the thruster at the stern
    (x_T < 0) the side force pushes the stern to PORT (Y_T > 0)."""
    s = replace(BASE, k_p=0.01)
    Y_T = control.thruster_force(0.1, 0.0, s)
    assert Y_T > 0
    N_c = s.x_T_L * s.hull.L * Y_T
    expected = -0.01 * 0.1 * control.force_scale(s) * s.hull.L
    assert np.isclose(N_c, expected, rtol=1e-12)


def test_thrust_never_exceeds_the_limit():
    s = replace(BASE, k_p=0.05, k_d=0.05, thrust_limit=0.2)
    for psi, r in [(0.5, 0.0), (-0.5, 0.0), (0.0, 0.01), (0.3, -0.02)]:
        assert abs(control.thruster_force(psi, r, s)) <= 0.2 * s.C_F * control.force_scale(s) * (1 + 1e-12)
    t, y = mr.simulate(s, [0.3, 0.2, 0.0, 0.0], 60.0, n_out=601, rtol=1e-8, atol=1e-10)
    ratios = np.array([control.thruster_ratio(y[:, k], s) for k in range(y.shape[1])])
    assert np.all(np.abs(ratios) <= 0.2 + 1e-12)


def test_linear_jacobian_with_control_matches_finite_differences():
    s = replace(BASE, k_p=0.004, k_d=0.003)
    J_exact = mr.linear_jacobian_star(s)
    J_fd = mr.jacobian_star(np.zeros(4), s)
    assert np.allclose(J_exact, J_fd, rtol=1e-4, atol=1e-5 * np.abs(J_fd).max())


def test_linear_jacobian_with_tug_matches_finite_differences():
    s = tug(BASE, 0.5)
    J_fd = mr.jacobian_star(np.zeros(4), s)
    assert np.allclose(mr.linear_jacobian_star(s), J_fd, rtol=1e-4, atol=1e-5 * np.abs(J_fd).max())


def test_mirror_symmetry_with_control_and_tug():
    """Mirror image (gamma, psi) -> (-gamma, -psi) of a solution is a solution: thruster law and tug are both odd."""
    s = tug(replace(BASE, k_p=0.01, k_d=0.01, thrust_limit=0.3), 0.3)
    assert np.allclose(mr.rhs_star(0, -STATE, s), -mr.rhs_star(0, STATE, s), rtol=1e-12, atol=1e-14)


def test_straight_is_still_an_equilibrium():
    for s in (replace(BASE, k_p=0.01, k_d=0.01), tug(BASE, 0.5)):
        assert np.allclose(mr.rhs_star(0, np.zeros(4), s), 0.0, atol=1e-14)
        assert np.allclose(co.all_equilibria(s)[0], 0.0)


def test_thrust_limit_does_not_move_the_eigenvalues_at_straight():
    """At straight the thruster force is zero, below any limit, so the eigenvalues (and the Hopf gain) do not depend
    on the limit."""
    a = co.eigenvalues(np.zeros(4), replace(BASE, k_p=0.004, k_d=0.003))
    b = co.eigenvalues(np.zeros(4), replace(BASE, k_p=0.004, k_d=0.003, thrust_limit=0.05))
    assert np.allclose(a, b)


def test_base_case_fishtails_without_control():
    assert co.max_oscillatory_growth(BASE) > 0
