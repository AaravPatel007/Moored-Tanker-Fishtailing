"""Wind checks (decision D28, experiment E10): a steady wind force of fixed size and direction on the rigid model."""

from dataclasses import replace

import numpy as np
import pytest

from fishtail import continuation as co
from fishtail import hydro
from fishtail import model_rigid as mr
from fishtail.params import VIDAL_100

S0 = mr.Rigid(VIDAL_100, x_A_L=0.25, L_H_L=0.5, C_F=0.05)
STATE = np.array([0.3, 0.1, 0.02, -0.01])  # an arbitrary moving state (rad, rad per t*)


def test_zero_wind_is_the_old_model():
    """C_W = 0 must give exactly the model of steps 3-6, whatever the wind direction."""
    windless = replace(S0, C_W=0.0, beta_W=1.0, x_W_L=-0.3)
    assert np.array_equal(mr.rhs_star(0, STATE, windless), mr.rhs_star(0, STATE, S0))


def test_wind_force_is_fixed_in_earth_axes():
    """Rotating the body components back to earth axes gives the same force for any heading: size C_W q and
    pointing where the wind blows, -(cos beta, sin beta)."""
    s = replace(S0, C_W=0.01, beta_W=np.radians(70.0))
    q = 0.5 * s.fluid.rho * s.U**2 * s.hull.L * s.hull.T
    for psi in np.radians([-40.0, 0.0, 25.0, 90.0]):
        X, Y, _ = mr.wind_forces(psi, s)
        earth = X * np.array([np.cos(psi), np.sin(psi)]) + Y * np.array([-np.sin(psi), np.cos(psi)])
        assert np.allclose(earth, -s.C_W * q * np.array([np.cos(s.beta_W), np.sin(s.beta_W)]), rtol=1e-12)


@pytest.mark.parametrize("ratio, beta_deg, x_W_L", [(0.05, 90, 0.0), (0.2, 45, 0.0), (0.3, -120, -0.3),
                                                     (0.1, 0, 0.0)])
def test_wind_equilibria_satisfy_the_equations(ratio, beta_deg, x_W_L):
    """Every state returned as a wind equilibrium really is one (all rates of change zero) and pulls the hawser."""
    s = replace(S0, C_W=ratio * S0.C_F, beta_W=np.radians(beta_deg), x_W_L=x_W_L)
    states = co.wind_equilibria(s)
    assert states
    for state in states:
        assert np.allclose(mr.rhs_star(0, state, s), 0.0, atol=1e-9)
        assert mr.hawser_force_star(state, s) > 0


def test_mirror_symmetry():
    """Wind from starboard is the mirror image of wind from port: (gamma, psi) -> (-gamma, -psi)."""
    port = co.wind_equilibria(replace(S0, C_W=0.004, beta_W=np.radians(60.0)))
    stbd = co.wind_equilibria(replace(S0, C_W=0.004, beta_W=np.radians(-60.0)))
    a = np.array(sorted(st[1] for st in port))
    b = np.array(sorted(-st[1] for st in stbd))
    assert np.allclose(a, b, atol=1e-9)


def test_aligned_wind_keeps_the_pitchfork():
    """Wind from dead ahead has no sideways part: trailing straight stays an equilibrium (perfect pitchfork)."""
    s = replace(S0, C_W=0.01, beta_W=0.0)
    assert co.wind_is_symmetric(s)
    assert np.allclose(mr.rhs_star(0, np.zeros(4), s), 0.0, atol=1e-12)


def test_side_wind_breaks_the_pitchfork():
    """Any sideways wind: trailing straight is no longer an equilibrium, and forward of the old pitchfork point the
    single rest position is turned slightly off the current line (the 'imperfect' pitchfork)."""
    s = replace(S0, x_A_L=0.5, C_W=0.002, beta_W=np.radians(90.0))
    assert not np.allclose(mr.rhs_star(0, np.zeros(4), s), 0.0, atol=1e-6)
    states = co.wind_equilibria(s)
    assert len(states) == 1
    assert 0 < abs(np.degrees(states[0][1])) < 5


def test_parametric_curve_matches_root_finding():
    """The closed form x_A/L(psi) and the root search give the same rest headings."""
    s = replace(S0, C_W=0.003, beta_W=np.radians(-90.0))
    for state in co.wind_equilibria(s):
        assert np.isclose(co.attachment_for_heading(state[1], s), s.x_A_L, atol=1e-9)


def test_wind_without_side_part_recovers_pitchfork_point():
    """With no wind the curve x_A/L(psi) meets psi = 0 exactly at x_Tcr of eq. (3.37)."""
    branches = co.wind_branches(replace(S0, C_W=0.0), n_psi=41, n_straight=5)
    assert np.isclose(branches["pitchfork"], hydro.critical_turret_position(VIDAL_100), atol=1e-6)


def test_folds_exist_aft_of_the_pitchfork():
    """A weak side wind turns the pitchfork into a smooth branch plus a separate branch that ends in a fold."""
    branches = co.wind_branches(replace(S0, C_W=0.0025, beta_W=np.radians(90.0)), n_psi=201)
    assert len(branches["folds"]) == 1
    x_fold, _ = branches["folds"][0]
    assert x_fold < hydro.critical_turret_position(VIDAL_100)


@pytest.mark.parametrize("beta_deg", [-150, -90, -30, 60])
def test_demo_fold_force_matches_counting_rest_positions(beta_deg):
    """The fold force marked on the demo's slider: just below it there are 3 rest positions, just above only 1."""
    from interactive.build_demos import wind_fold_ratio

    ratio = wind_fold_ratio(S0, beta_deg)
    below = co.wind_equilibria(replace(S0, C_W=0.98 * ratio * S0.C_F, beta_W=np.radians(beta_deg)))
    above = co.wind_equilibria(replace(S0, C_W=1.02 * ratio * S0.C_F, beta_W=np.radians(beta_deg)))
    assert (len(below), len(above)) == (3, 1)


@pytest.mark.parametrize("beta_deg", [-120, 90])
def test_force_curve_points_are_rest_positions(beta_deg):
    """Every point of the 'heading against wind force' curve is a rest position for that wind force."""
    s = replace(S0, beta_W=np.radians(beta_deg))
    br = co.force_branches(s, n_psi=81)
    for seg in br["segments"]:
        for ratio, psi in list(zip(seg["ratio"], seg["psi"]))[::7]:
            sk = replace(s, C_W=float(ratio) * s.C_F)
            assert np.allclose(mr.rhs_star(0, co.equilibrium_state(psi, sk), sk), 0.0, atol=1e-9)


def test_linear_matrices_refuse_sideways_wind():
    """Wind with a sideways part has no straight equilibrium, so there is no linearisation about it. (Head-on wind
    and the hold-back tug, sin beta_W = 0, are covered since step 8: tests/test_control.py.)"""
    with pytest.raises(ValueError):
        mr.linear_matrices(replace(S0, C_W=0.01, beta_W=np.radians(90.0)))
