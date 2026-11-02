"""Steps 4 and 5 checks: equilibria, their stability, and the bifurcation points (fishtail/continuation.py)."""

from dataclasses import replace

import numpy as np
import pytest

from fishtail import continuation as co
from fishtail import hydro
from fishtail import model_rigid as mr
from fishtail.params import SIMOS_FIG53, VIDAL_100, VIDAL_40

HULLS = [VIDAL_100, VIDAL_40]


@pytest.mark.parametrize("hull", HULLS)
@pytest.mark.parametrize("x_A_L, L_H_L, C_F", [(0.1, 0.5, 0.03), (0.25, 1.0, 0.015), (0.5, 0.3, 0.05)])
def test_equilibria_satisfy_the_equations(hull, x_A_L, L_H_L, C_F):
    """Every state returned as an equilibrium really is one: all rates of change are zero there."""
    s = mr.Rigid(hull, x_A_L, L_H_L, C_F)
    for state in co.all_equilibria(s, include_reversed=True):
        assert np.allclose(mr.rhs_star(0, state, s), 0.0, atol=1e-9)


@pytest.mark.parametrize("hull", HULLS)
def test_angled_heading_depends_only_on_attachment_point(hull):
    """The rest heading comes from 'zero moment about the attachment point', so the hawser length and the friction
    cannot change it (they do change the hawser angle)."""
    x_A_L = 0.6 * hydro.critical_turret_position(hull)
    a = co.equilibrium_state(co.equilibrium_headings(hull, x_A_L)[0], mr.Rigid(hull, x_A_L, 0.3, 0.015))
    b = co.equilibrium_state(co.equilibrium_headings(hull, x_A_L)[0], mr.Rigid(hull, x_A_L, 1.2, 0.05))
    assert a[1] == pytest.approx(b[1], abs=1e-12)
    assert abs(a[0] - b[0]) > np.radians(1.0)


@pytest.mark.parametrize("hull", HULLS)
def test_angled_positions_exist_only_aft_of_x_Tcr(hull):
    """Pitchfork: angled rest positions exist just aft of x_Tcr and not just forward of it."""
    x_cr = hydro.critical_turret_position(hull)
    assert len(co.equilibrium_headings(hull, x_cr - 0.01)) == 1
    assert len(co.equilibrium_headings(hull, x_cr + 0.01)) == 0


@pytest.mark.parametrize("hull", HULLS)
def test_branches_grow_linearly_not_as_square_root(hull):
    """Doubling the distance past the pitchfork doubles the angle (normal form x' = r x - x|x|), E03 fig06."""
    x_cr = hydro.critical_turret_position(hull)
    psi_1 = co.equilibrium_headings(hull, x_cr - 1e-3)[0]
    psi_2 = co.equilibrium_headings(hull, x_cr - 2e-3)[0]
    assert psi_2 / psi_1 == pytest.approx(2.0, rel=0.01)


@pytest.mark.parametrize("hull, x_A_L", [(VIDAL_100, 0.5), (VIDAL_40, 0.4)])
def test_hopf_point_is_a_complex_pair_crossing(hull, x_A_L):
    """At the Hopf friction the swing's eigenvalue pair sits on the imaginary axis with a nonzero frequency;
    a little more friction is stable, a little less is not."""
    s = mr.Rigid(hull, x_A_L, 0.5, 0.03)
    C_h = co.hopf_C_F(s)
    assert co.max_oscillatory_growth(replace(s, C_F=C_h)) == pytest.approx(0.0, abs=1e-9)
    assert co.hopf_frequency(replace(s, C_F=C_h)) > 0.1
    assert co.is_stable(co.eigenvalues(np.zeros(4), replace(s, C_F=1.02 * C_h)))
    assert not co.is_stable(co.eigenvalues(np.zeros(4), replace(s, C_F=0.98 * C_h)))


@pytest.mark.parametrize("key, hull, tol", [("100", VIDAL_100, 0.06), ("40", VIDAL_40, 0.10)])
def test_hopf_curve_matches_simos_fig_5_3(key, hull, tol):
    """Validation target: our Hopf boundary against Simos Fig. 5.3 (digitised, +-0.0005).
    100 %: mean difference about +0.5 %; 40 %: about +5.6 % (open issue 18)."""
    s = mr.Rigid(hull, SIMOS_FIG53[key]["x_A_L"], 0.5, 0.03)
    ours = co.hopf_curve(s, SIMOS_FIG53["L_H_L"])
    assert np.allclose(ours, SIMOS_FIG53[key]["C_F"], rtol=tol)


def test_eigenvalue_crossing_matches_time_domain():
    """Step 4 check from the brief: either side of the eigenvalue crossing, a small disturbance really dies out or
    grows in a simulation."""
    s = mr.Rigid(VIDAL_100, 0.5, 0.5, 0.03)
    C_h = co.hopf_C_F(s)
    start = [np.radians(1.0), 0.0, 0.0, 0.0]
    for factor, should_grow in [(1.1, False), (0.9, True)]:
        t, y = mr.simulate(replace(s, C_F=factor * C_h), start, t_end_star=120.0, n_out=1201)
        late_swing = np.abs(y[0, -200:]).max()
        assert (late_swing > np.radians(1.0)) == should_grow


def test_pitchfork_crossing_matches_time_domain():
    """Either side of x_Tcr (high friction, so no fishtailing): the ship ends straight or at the predicted angle."""
    hull = VIDAL_100
    x_cr = hydro.critical_turret_position(hull)
    start = [np.radians(3.0), np.radians(1.0), 0.0, 0.0]
    t, y = mr.simulate(mr.Rigid(hull, x_cr + 0.05, 0.5, 0.05), start, t_end_star=200.0, n_out=401)
    assert abs(y[1, -1]) < np.radians(0.05)
    x_aft = x_cr - 0.1
    t, y = mr.simulate(mr.Rigid(hull, x_aft, 0.5, 0.05), start, t_end_star=200.0, n_out=401)
    assert y[1, -1] == pytest.approx(co.equilibrium_headings(hull, x_aft)[0], abs=1e-4)
