"""Step 7 checks: the stretchy, slack-able hawser model (fishtail/model_elastic.py)."""

from dataclasses import replace

import numpy as np
import pytest

from fishtail import continuation as co
from fishtail import model_elastic as me
from fishtail import model_rigid as mr
from fishtail.params import VIDAL_100

S = me.Elastic(VIDAL_100, x_A_L=0.5, L_H_L=0.5, C_F=0.03, U=1.0)


def test_rope_fit_is_stiffening():
    """The fitted exponent of the real hawser (Huang et al. 2026) is about 1.9: stiffer the more it stretches."""
    assert 1.7 < me.ROPE_N < 2.0


def test_trailing_straight_is_an_equilibrium():
    eq, strain0 = me.trailing_straight(S)
    assert strain0 > 0
    rates = me.rhs(0, eq, S)
    assert np.abs(rates[:3]).max() < 1e-9 and np.abs(rates[3:]).max() < 1e-12


def test_slack_rope_carries_no_force():
    """Ship closer to the buoy than the rope length: no tension at all (a rope cannot push)."""
    state = np.array([-0.5 * S.L_H - S.x_A, 0.0, 0.0, S.U, 0.0, 0.0])
    T, strain = me.tension(state, S)
    assert strain < 0 and T == 0.0


def test_stiff_rope_gives_rigid_eigenvalues():
    """Step 7 check from the brief: a very stiff rope behaves like the rigid hawser."""
    stiff = replace(S, rope_a=S.rope_a * 1000)
    eig = me.eigenvalues_star(stiff)
    pair = min((e for e in eig if e.imag > 1e-6), key=lambda e: e.imag)
    eig_r = co.eigenvalues(np.zeros(4), mr.Rigid(VIDAL_100, 0.5, 0.5, 0.03))
    rigid = eig_r[np.abs(eig_r.imag) > 1e-9][0]
    assert pair.real == pytest.approx(rigid.real, abs=2e-3)
    assert pair.imag == pytest.approx(abs(rigid.imag), rel=2e-3)
