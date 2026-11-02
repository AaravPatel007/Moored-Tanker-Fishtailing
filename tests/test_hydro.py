"""Step 2 checks: the water-force model (fishtail/hydro.py) and the bare-hull yaw check (fishtail/model_pinned.py)."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import quad

from fishtail import hydro, model_pinned
from fishtail.params import (FRESH_WATER, FULL_SCALE_CF_CHECK, SEA_WATER, TABLE_5_1_CF, TABLE_5_1_U_MODEL,
                             VIDAL_100, VIDAL_100_IPT, VIDAL_40, X_TCR_40_CHECK, froude_scaled)

HULLS = [VIDAL_100, VIDAL_40]
C0 = 0.02  # any friction value will do for these checks


# --- Basic physics ----------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("hull", HULLS)
def test_no_relative_speed_no_force(hull):
    """A hull at rest in still water (u = v = r = 0) feels no water force."""
    assert hydro.current_forces(0.0, 0.0, 0.0, hull, SEA_WATER, C0) == (0.0, 0.0, 0.0)


@pytest.mark.parametrize("hull", HULLS)
def test_angle_of_attack_conventions(hull):
    """alpha = 180 deg: hull moving forward through the water (current from ahead); 0 deg: from astern."""
    assert hydro.angle_of_attack(1.0, 0.0) == pytest.approx(np.pi)
    assert hydro.angle_of_attack(-1.0, 0.0) == pytest.approx(0.0, abs=1e-12)
    assert hydro.angle_of_attack(0.0, -1.0) == pytest.approx(np.pi / 2)  # hull sliding to starboard


@pytest.mark.parametrize("hull", HULLS)
def test_special_angles(hull):
    """At the special angles the formulas reduce to the measured numbers they were built from."""
    assert hydro.C2C(np.pi / 2, hull) == pytest.approx(hull.C_Y)  # side-on: sideways drag = C_Y
    assert hydro.C6C(np.pi / 2, hull) == pytest.approx(-hull.lp_L * hull.C_Y)  # acting l_p aft
    assert hydro.C1C(np.pi, hull, C0) == pytest.approx(-C0)  # head-on: pure friction, pointing aft
    assert hydro.C1C(0.0, hull, C0) == pytest.approx(C0)


@pytest.mark.parametrize("hull", HULLS)
def test_drag_opposes_motion(hull):
    """Sliding sideways to port (v > 0), the hull is pushed to starboard (F_Y < 0); moving forward, pushed aft."""
    FX, FY, _ = hydro.current_forces(0.0, 1.0, 0.0, hull, SEA_WATER, C0)
    assert FY < 0
    FX, FY, _ = hydro.current_forces(1.0, 0.0, 0.0, hull, SEA_WATER, C0)
    assert FX < 0 and FY == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize("hull", HULLS)
def test_port_starboard_mirror_symmetry(hull):
    """Mirror image (v -> -v, r -> -r): the hull is symmetric, so F_X stays and F_Y, N change sign."""
    for u, v, r in [(1.0, 0.3, 0.002), (0.8, -0.5, -0.004), (-0.6, 0.2, 0.003), (0.1, 0.9, 0.01)]:
        FX, FY, N = hydro.current_forces(u, v, r, hull, SEA_WATER, C0)
        FXm, FYm, Nm = hydro.current_forces(u, -v, -r, hull, SEA_WATER, C0)
        assert FXm == pytest.approx(FX, rel=1e-12)
        assert FYm == pytest.approx(-FY, rel=1e-12)
        assert Nm == pytest.approx(-N, rel=1e-12)


# --- Sideways drag along the hull and the yaw-rate terms ------------------------------------------------------------

@pytest.mark.parametrize("hull", HULLS)
def test_sectional_drag_matches_measured_numbers(hull):
    """Our straight-line C_D(x) reproduces C_Y and puts the centre of the sideways drag l_p aft of midship."""
    total, _ = quad(lambda xi: hydro.sectional_drag(xi, hull), -0.5, 0.5)
    first_moment, _ = quad(lambda xi: hydro.sectional_drag(xi, hull) * xi, -0.5, 0.5)
    assert total == pytest.approx(hull.C_Y)
    assert first_moment == pytest.approx(-hull.lp_L * hull.C_Y)
    assert hydro.sectional_drag(0.5, hull) > 0  # drag is positive everywhere along the hull


@pytest.mark.parametrize("j", [0, 1, 2, 3])
@pytest.mark.parametrize("v, r", [(0.5, 0.0), (-0.5, 0.0), (0.2, 0.001), (0.05, 0.004), (-0.1, -0.003), (0.0, 0.002)])
def test_drag_integral_exact_vs_numerical(j, v, r):
    """The hand-integrated I_j agrees with brute-force numerical integration (including when the sign flips)."""
    hull = VIDAL_100
    xi0 = -v / (r * hull.L) if r != 0 else 1.0
    f = lambda xi: hydro.sectional_drag(xi, hull) * np.sign(v + r * hull.L * xi) * xi**j
    numerical, _ = quad(f, -0.5, 0.5, points=[xi0] if abs(xi0) < 0.5 else None)
    assert hydro.drag_integral(j, v, r, hull) == pytest.approx(numerical, abs=1e-10)


@pytest.mark.parametrize("hull", HULLS)
def test_yaw_rate_terms_vanish_without_turning(hull):
    """No double counting: with r = 0 every yaw-rate term is zero (thesis p. 126), up to round-off.
    Compared with the size of the static force (1/2) rho U^2 L T."""
    for u, v in [(1.0, 0.0), (0.7, 0.4), (0.0, -1.0), (-0.5, 0.3)]:
        scale = 0.5 * SEA_WATER.rho * (u**2 + v**2) * hull.L**2 * hull.T  # force x length, so it covers N too
        for value in hydro.yaw_rate_terms(u, v, 0.0, hull, SEA_WATER).values():
            assert abs(value) < 1e-12 * scale


@pytest.mark.parametrize("hull", HULLS)
def test_cross_flow_terms_equal_strip_theory(hull):
    """Independent check of the algebra in (4.29): it must equal the slice-by-slice drag integral
    -1/2 rho T  integral C_D(x) (v + r x)|v + r x| dx   minus the cross-flow part already in the static model."""
    rho, L, T = SEA_WATER.rho, hull.L, hull.T
    for v, r in [(0.3, 0.002), (-0.2, 0.004), (0.05, -0.003), (0.0, 0.001)]:
        speed = lambda x: v + r * x
        xs = [-v / r] if abs(v / r) < L / 2 else None
        Y_strip, _ = quad(lambda x: -0.5 * rho * T * hydro.sectional_drag(x / L, hull) * speed(x) * abs(speed(x)),
                          -L / 2, L / 2, points=xs)
        N_strip, _ = quad(lambda x: -0.5 * rho * T * hydro.sectional_drag(x / L, hull) * speed(x) * abs(speed(x)) * x,
                          -L / 2, L / 2, points=xs)
        Y_static = -0.5 * rho * T * L * hull.C_Y * v * abs(v)  # cross-flow part of (3.21) at this v
        N_static = 0.5 * rho * T * L**2 * hull.lp_L * hull.C_Y * v * abs(v)  # cross-flow part of (3.29)
        t = hydro.yaw_rate_terms(0.0, v, r, hull, SEA_WATER)
        assert t["Y_cross_flow"] == pytest.approx(Y_strip - Y_static, rel=1e-8, abs=1e-6)
        assert t["N_cross_flow"] == pytest.approx(N_strip - N_static, rel=1e-8, abs=1e-3)


@pytest.mark.parametrize("hull", HULLS)
def test_turning_on_the_spot_matches_ipt_tests(hull):
    """Hull turning in still water (u = v = 0): N = -(1/32) rho T L^4 C_Y r|r|, eq. (4.31)."""
    for r in [0.001, -0.003]:
        _, _, N = hydro.current_forces(0.0, 0.0, r, hull, SEA_WATER, C0)
        assert N == pytest.approx(-(1 / 32) * SEA_WATER.rho * hull.T * hull.L**4 * hull.C_Y * r * abs(r))


@pytest.mark.parametrize("hull", HULLS)
def test_turning_on_the_spot_sideways_force(hull):
    """Turning on the spot, bow and stern slide opposite ways. Their sideways drags cancel only if the drag is the
    same at both ends (l_p = 0). Our C_D is larger aft, so the stern wins and the hull is pushed sideways."""
    r = 0.002
    symmetric = replace(hull, lp_L=0.0)
    _, FY_sym, _ = hydro.current_forces(0.0, 0.0, r, symmetric, SEA_WATER, C0)
    assert FY_sym == pytest.approx(0.0, abs=1e-9 * SEA_WATER.rho * hull.T * hull.L**3 * r**2)
    _, FY, _ = hydro.current_forces(0.0, 0.0, r, hull, SEA_WATER, C0)
    # r > 0 turns the bow to port, so the stern slides to starboard; the water pushes the stern back to port: F_Y > 0
    assert FY > 0


@pytest.mark.parametrize("hull", HULLS)
def test_yaw_damping_head_on(hull):
    """Ship facing the current and starting to turn: the water resists the turn (N has the opposite sign to r)."""
    U = 1.0
    for r in [0.0005, -0.0005]:
        _, _, N = hydro.current_forces(U, 0.0, r, hull, SEA_WATER, C0)
        assert N * r < 0


# --- Bare-hull yaw instability (the sanity check required by the project brief) ------------------------------------

@pytest.mark.parametrize("hull", HULLS)
def test_head_on_slopes_match_clarke(hull):
    """Near head-on current the static coefficients reduce to Clarke's linear derivatives, eqs. (3.28), (3.31)."""
    assert hydro.head_on_slopes(hull) == pytest.approx(hydro.clarke_slopes(hull), rel=1e-6)


@pytest.mark.parametrize("hull", HULLS)
def test_bare_hull_turns_away_from_current(hull):
    """Hull held at midship, bow turned a little to port (psi > 0): the water turns it further to port (N > 0).
    This is the Munk moment beating the hull/rudder lift. If this failed, the hydrodynamics would be wrong."""
    U = 1.0
    for psi in [np.radians(1.0), np.radians(5.0)]:
        _, _, N = hydro.static_forces(U * np.cos(psi), -U * np.sin(psi), hull, SEA_WATER, C0)
        assert N > 0
        _, _, N_mirror = hydro.static_forces(U * np.cos(-psi), -U * np.sin(-psi), hull, SEA_WATER, C0)
        assert N_mirror < 0


@pytest.mark.parametrize("hull", HULLS)
def test_pinned_hull_head_on_is_a_saddle(hull):
    """Linearised pinned hull at psi = 0: one eigenvalue > 0 and one < 0 (a saddle point), so it is unstable."""
    eig = np.linalg.eigvals(model_pinned.jacobian(np.zeros(2), hull, SEA_WATER, U=1.0))
    assert np.all(np.isreal(eig))
    assert eig.real.max() > 0 > eig.real.min()


def test_critical_turret_position():
    """x_Tcr/L from the slopes equals eq. (3.37); for the 40 % ship it is ~0.19 as stated in the thesis (p. 198)."""
    for hull in HULLS:
        TL, CBBT = hull.T / hull.L, hull.C_B * hull.B / hull.T
        assert hydro.critical_turret_position(hull) == pytest.approx((0.5 + 2.4 * TL) / (1 + 0.4 * CBBT), rel=1e-6)
    assert hydro.critical_turret_position(VIDAL_40) == pytest.approx(X_TCR_40_CHECK, abs=0.005)


# --- Friction coefficient: the only way current speed enters the rigid model ----------------------------------------

def test_friction_full_scale_matches_thesis():
    """40 % ship, full scale, U = 1.90 m/s: Simos gets C_F = 0.015 (p. 195)."""
    cf = hydro.friction_coefficient(FULL_SCALE_CF_CHECK["U"], VIDAL_40, SEA_WATER)
    assert cf == pytest.approx(FULL_SCALE_CF_CHECK["C_F"], rel=0.02)


def test_friction_model_scale_matches_table_5_1_slowest_speed():
    """1:90 model at 0.2 m/s: Table 5-1 values within 5 % (the water temperature of the tests is not known).
    The faster rows of Table 5-1 do NOT match k = 0.39; see notes/open_issues.md."""
    for key, hull in [("100", VIDAL_100), ("40", VIDAL_40)]:
        model = froude_scaled(hull, FRESH_WATER.rho / SEA_WATER.rho)
        cf = hydro.friction_coefficient(TABLE_5_1_U_MODEL[0], model, FRESH_WATER)
        assert cf == pytest.approx(TABLE_5_1_CF[key][0], rel=0.05)


def test_friction_falls_with_speed():
    """Faster current -> larger Reynolds number -> smaller friction coefficient."""
    speeds = np.linspace(0.2, 2.5, 20)
    cf = hydro.friction_coefficient(speeds, VIDAL_100, SEA_WATER)
    assert np.all(np.diff(cf) < 0)


# --- Parameters ------------------------------------------------------------------------------------------------------

def test_froude_scaling_keeps_pure_numbers():
    model = froude_scaled(VIDAL_100_IPT, 1.0)
    assert model.L == pytest.approx(320 / 90)
    assert model.C_Y == VIDAL_100_IPT.C_Y and model.lp_L == VIDAL_100_IPT.lp_L
    assert np.sqrt(model.I_z / model.M) / model.L == pytest.approx(np.sqrt(VIDAL_100.I_z / VIDAL_100.M) / 320)


@pytest.mark.parametrize("hull", HULLS)
def test_radius_of_gyration_is_quarter_length(hull):
    """Why we read Annex A's I_z as 'about the CG' (decision D13): sqrt(I_z/M) = 0.25 L."""
    assert np.sqrt(hull.I_z / hull.M) / hull.L == pytest.approx(0.25, abs=0.002)
