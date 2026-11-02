"""Limit cycles (fishtailing) of the rigid-hawser system (step 6).

Sign convention (Simos 2001): psi > 0 = bow turned to port; gamma > 0 = ship swung to port. Time t* = U t / L.

Three tools, from simple to less simple:
  1. measure_cycle(): just simulate until the swing stops changing, then measure amplitude and period.
     This is what a towing-tank test does. It only finds STABLE cycles.
  2. shoot(): the shooting method (Strogatz glossary in notes/glossary.md). Guess a point on the cycle and the period,
     simulate one period, and correct the guess with Newton's method until the orbit closes on itself.
     Also finds UNSTABLE cycles if the guess is close.
  3. floquet_multipliers(): how a small kick to the cycle grows or shrinks after one period. One multiplier is always
     1 (a kick along the cycle just shifts the timing); the cycle is stable if all the others are inside the unit circle.

Where a cycle "starts": we use the moment the hawser angle gamma is at its maximum (d gamma/dt* = 0, decreasing).
"""

import numpy as np
from scipy.integrate import solve_ivp

from fishtail import model_rigid as mr

TIGHT = dict(rtol=1e-10, atol=1e-12)  # for shooting and Floquet multipliers (finite differences need accuracy)
NORMAL = dict(rtol=1e-9, atol=1e-11)  # for long simulations


def _gamma_max_event(t, x, s):
    """Zero when d gamma/dt* = 0; with direction = -1 only maxima of gamma count."""
    return x[2]


_gamma_max_event.direction = -1


def flow(s, x0, T, tol=TIGHT):
    """State after simulating for a time T (t* units) from x0."""
    sol = solve_ivp(mr.rhs_star, (0.0, T), x0, args=(s,), **tol)
    return sol.y[:, -1]


def gamma_maxima(s, x0, t_end, tol=NORMAL):
    """Simulate from x0 for t_end and return the times and states at every maximum of gamma."""
    sol = solve_ivp(mr.rhs_star, (0.0, t_end), x0, args=(s,), events=_gamma_max_event, **tol)
    return sol.t_events[0], sol.y_events[0], sol.y[:, -1]


def one_period(s, x0, T, n=801, tol=TIGHT):
    """Densely sampled states over one period starting at x0: returns t (n,), states (4, n)."""
    t = np.linspace(0.0, T, n)
    sol = solve_ivp(mr.rhs_star, (0.0, T), x0, args=(s,), t_eval=t, **tol)
    return sol.t, sol.y


def summarise(s, x0, T):
    """Amplitudes (half of peak-to-peak, degrees), period and hawser force range over one period."""
    t, y = one_period(s, x0, T)
    FH = np.array([mr.hawser_force_star(y[:, i], s) for i in range(t.size)])
    return {
        "period": float(T),  # t* units
        "gamma_amp_deg": float(np.degrees(0.5 * (y[0].max() - y[0].min()))),
        "psi_amp_deg": float(np.degrees(0.5 * (y[1].max() - y[1].min()))),
        "gamma_mean_deg": float(np.degrees(0.5 * (y[0].max() + y[0].min()))),
        "psi_mean_deg": float(np.degrees(0.5 * (y[1].max() + y[1].min()))),
        "F_H_max": float(FH.max()),  # divided by (1/2) rho U^2 L T, like C_F
        "F_H_min": float(FH.min()),
        "x0": [float(v) for v in x0],
    }


def measure_cycle(s, x0, t_transient=400.0, t_chunk=150.0, max_chunks=12, tol_amp=1e-4, amp_zero_deg=0.05):
    """Simulate like a tank test: run t_transient, then keep going in chunks until two successive gamma maxima
    agree to tol_amp (relative). Returns (summary dict, final state), or (None, final state) if the motion died out
    (amplitude below amp_zero_deg degrees): then the ship trails straight or sits still at an angle.
    """
    x = flow(s, np.asarray(x0, float), t_transient, tol=NORMAL)
    T_guess = None
    for _ in range(max_chunks):
        times, states, x_end = gamma_maxima(s, x, t_chunk)
        if len(times) < 3:  # no swing, or a very slow one
            x = x_end
            if np.degrees(np.abs(x[2]) + np.abs(x[3])) < 1e-6:
                break
            continue
        g_max = states[:, 0]
        T_guess = float(times[-1] - times[-2])
        if np.abs(g_max[-1] - g_max[-2]) <= tol_amp * max(np.abs(g_max[-1]), 1e-12):
            summary = summarise(s, states[-1], T_guess)
            if summary["gamma_amp_deg"] < amp_zero_deg and summary["psi_amp_deg"] < amp_zero_deg:
                return None, states[-1]
            return summary, states[-1]
        x = states[-1]
    # Not converged: either it decayed to an equilibrium, or (close to a Hopf point) the swing is still creeping
    # towards the cycle very slowly - "critical slowing down". In the second case the shooting method finds the
    # cycle directly, because Newton's method does not care how slow the approach is.
    t, y = one_period(s, x, 20.0, tol=NORMAL)
    swing = np.degrees(y[0].max() - y[0].min())
    if swing < 2 * amp_zero_deg:
        return None, y[:, -1]
    if T_guess is not None:
        x_c, T_c, ok = shoot(s, x, T_guess)
        if ok and np.all(np.abs(floquet_multipliers(s, x_c, T_c)[1:]) < 1.0):  # a STABLE cycle
            return summarise(s, x_c, T_c), x_c
    raise RuntimeError(f"cycle did not settle (swing {swing:.3f} deg still changing) for {s}")


def monodromy(s, x0, T, h=1e-6):
    """4x4 matrix: how a small change of the start state changes the state one period later (finite differences)."""
    Phi = np.zeros((4, 4))
    for k in range(4):
        dx = np.zeros(4)
        dx[k] = h
        Phi[:, k] = (flow(s, x0 + dx, T) - flow(s, x0 - dx, T)) / (2 * h)
    return Phi


def shoot(s, x0, T, max_iter=8, tol=1e-9):
    """Newton's method for a closed orbit: find x0 (with d gamma/dt* = 0) and T such that flow(x0, T) = x0.

    Unknowns: gamma, psi, psi_dot at the start (gamma_dot is fixed at 0 - the 'phase condition') and the period T.
    Returns x0, T, converged (bool).
    """
    x0 = np.array(x0, float)
    x0[2] = 0.0
    for _ in range(max_iter):
        xT = flow(s, x0, T)
        residual = xT - x0
        if np.max(np.abs(residual)) < tol:
            return x0, T, True
        Phi = monodromy(s, x0, T)
        # d(residual)/d(unknowns): columns for gamma, psi, psi_dot (Phi - I), and for T (the rate at the end point)
        J = np.zeros((4, 4))
        J[:, 0:2] = Phi[:, 0:2] - np.eye(4)[:, 0:2]
        J[:, 2] = Phi[:, 3] - np.eye(4)[:, 3]
        J[:, 3] = mr.rhs_star(0, xT, s)
        step = np.linalg.solve(J, -residual)
        x0[0] += step[0]
        x0[1] += step[1]
        x0[3] += step[2]
        T += step[3]
    return x0, T, bool(np.max(np.abs(flow(s, x0, T) - x0)) < tol)


def floquet_multipliers(s, x0, T):
    """Eigenvalues of the monodromy matrix, largest magnitude first."""
    mu = np.linalg.eigvals(monodromy(s, x0, T))
    return mu[np.argsort(-np.abs(mu))]
