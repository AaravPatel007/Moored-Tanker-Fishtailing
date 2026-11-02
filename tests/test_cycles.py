"""Step 6 checks: limit cycles, shooting and Floquet multipliers (fishtail/cycles.py). Each takes a few seconds."""

from dataclasses import replace

import numpy as np
import pytest

from fishtail import continuation as co
from fishtail import cycles as cy
from fishtail import model_rigid as mr
from fishtail.params import VIDAL_100

BASE = mr.Rigid(VIDAL_100, x_A_L=0.5, L_H_L=0.5, C_F=0.03)


@pytest.fixture(scope="module")
def cycle_at_80_percent():
    C_h = co.hopf_C_F(BASE)
    s = replace(BASE, C_F=0.8 * C_h)
    summ, _ = cy.measure_cycle(s, [np.radians(10.0), 0, 0, 0])
    x0, T, ok = cy.shoot(s, np.array(summ["x0"]), summ["period"])
    return s, x0, T, ok, summ


def test_shooting_closes_the_orbit(cycle_at_80_percent):
    s, x0, T, ok, summ = cycle_at_80_percent
    assert ok
    assert np.allclose(cy.flow(s, x0, T), x0, atol=1e-8)
    assert T == pytest.approx(summ["period"], rel=1e-3)  # simulation and shooting agree


def test_floquet_multipliers(cycle_at_80_percent):
    """One multiplier is exactly 1 (a kick along the cycle only shifts the timing); the cycle is stable."""
    s, x0, T, ok, _ = cycle_at_80_percent
    mu = cy.floquet_multipliers(s, x0, T)
    assert np.abs(mu[0] - 1.0) < 1e-4
    assert np.all(np.abs(mu[1:]) < 1.0)


def test_cycle_is_symmetric(cycle_at_80_percent):
    """Fishtailing about trailing straight swings equally to port and starboard."""
    s, x0, T, ok, _ = cycle_at_80_percent
    summ = cy.summarise(s, x0, T)
    assert abs(summ["gamma_mean_deg"]) < 0.01 and abs(summ["psi_mean_deg"]) < 0.01


def test_period_near_hopf_matches_linear_frequency():
    """Close to the Hopf point the cycle's period tends to 2 pi / omega of the crossing eigenvalue pair."""
    C_h = co.hopf_C_F(BASE)
    s = replace(BASE, C_F=0.98 * C_h)
    summ, _ = cy.measure_cycle(s, [np.radians(3.0), 0, 0, 0])
    assert summ["period"] == pytest.approx(2 * np.pi / co.hopf_frequency(replace(BASE, C_F=C_h)), rel=0.02)


def test_no_cycle_above_hopf_point():
    """Supercritical Hopf: above the Hopf friction even a large swing dies out (no hidden large cycle)."""
    C_h = co.hopf_C_F(BASE)
    summ, _ = cy.measure_cycle(replace(BASE, C_F=1.05 * C_h), [np.radians(60.0), 0, 0, 0])
    assert summ is None
