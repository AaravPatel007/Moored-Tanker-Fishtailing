"""Swing amplitude against hawser length for the Hawser length tab (x_A/L = 0.5, C_F = 0.03, loaded ship).

The hawser length is swept up (each run starting from a small nudge, gamma = 3 deg, psi = 1 deg, as in the lab) and
down (starting from a large swing), each run continuing from where the previous one ended. Matching curves mean the
Hopf bifurcation has no hysteresis (supercritical). The grid is fine near the Hopf point.

Run from the project root (a few minutes):
    python3 -m report.hawser_sweep
Writes report/data/hawser_sweep.json.
"""

import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from fishtail import continuation as co
from fishtail import cycles as cy
from fishtail import model_rigid as mr
from fishtail.params import VIDAL_100

OUT = Path(__file__).resolve().parent / "data" / "hawser_sweep.json"
C_F = 0.03
NUDGE = np.radians([3.0, 1.0, 0.0, 0.0])


def critical_length():
    return brentq(lambda L: co.max_oscillatory_growth(mr.Rigid(VIDAL_100, 0.5, L, C_F)), 0.1, 0.4)


def grid(L_cr):
    fine = L_cr + np.array([-0.02, -0.01, -0.005, 0.0025, 0.005, 0.01, 0.02, 0.03, 0.05, 0.08])
    coarse = np.array([0.05, 0.1, 0.15, 0.35, 0.45, 0.6, 0.75, 0.9, 1.05, 1.2, 1.35, 1.5])
    return np.unique(np.round(np.concatenate([fine, coarse]), 4))


def one_sweep(direction):
    L_cr = critical_length()
    Ls = grid(L_cr) if direction == "up" else grid(L_cr)[::-1]
    x = NUDGE.copy() if direction == "up" else np.radians([50.0, 0.0, 0.0, 0.0])
    rows = []
    for L in Ls:
        s = mr.Rigid(VIDAL_100, 0.5, float(L), C_F)
        try:
            summ, x_end = cy.measure_cycle(s, x, t_transient=600.0, max_chunks=30)
        except RuntimeError:
            summ, x_end = None, x
        if summ is None:
            rows.append({"L_H_L": float(L), "gamma_amp_deg": 0.0, "psi_amp_deg": 0.0})
            x = np.asarray(x_end) + NUDGE
        else:
            mu = cy.floquet_multipliers(s, np.array(summ["x0"]), summ["period"])
            rows.append({"L_H_L": float(L), "gamma_amp_deg": summ["gamma_amp_deg"], "psi_amp_deg": summ["psi_amp_deg"],
                         "period": summ["period"], "second_multiplier_abs": float(np.abs(mu[1]))})
            x = np.array(summ["x0"])
        print(f"{direction} L={L:.4f}: gamma {rows[-1]['gamma_amp_deg']:.2f} deg", flush=True)
    return {"direction": direction, "L_cr": L_cr, "C_F": C_F, "rows": rows}


def main():
    with Pool(2) as pool:
        results = pool.map(one_sweep, ["up", "down"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=1))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
