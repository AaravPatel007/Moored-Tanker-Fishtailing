"""Fine up/down sweep of the hawser length around the Hopf point: is the onset smooth, or is there a jump?

E05 fig05 uses a coarse grid in L_H/L, and the amplitude goes from 0 to about 30 deg between two neighbouring grid
points. This script resolves that step, with the same hull and friction values as E05 fig05, sweeping L_H/L up (from
a small disturbance) and down (from a large swing). Each run continues from where the previous one ended.

Run from the project root (a few minutes, four processes):
    python3 -m report.hawser_onset_sweep
Writes report/data/hawser_onset.json.
"""

import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from fishtail import continuation as co
from fishtail import cycles as cy
from fishtail import model_rigid as mr
from fishtail.params import VIDAL_100_IPT

OUT = Path(__file__).resolve().parent / "data" / "hawser_onset.json"
SMALL_KICK = np.array([0.01, 0.0, 0.0, 0.0])
OFFSETS = np.array([-0.02, -0.01, -0.005, 0.0025, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08, 0.12, 0.18])


def critical_length(C_F):
    return brentq(lambda L: co.max_oscillatory_growth(mr.Rigid(VIDAL_100_IPT, 0.5, L, C_F)), 0.05, 0.3)


def one_sweep(job):
    C_F, direction = job
    L_cr = critical_length(C_F)
    grid = L_cr + (OFFSETS if direction == "up" else OFFSETS[::-1])
    x = SMALL_KICK.copy() if direction == "up" else np.array([np.radians(50.0), 0.0, 0.0, 0.0])
    rows = []
    for L in grid:
        s = mr.Rigid(VIDAL_100_IPT, 0.5, float(L), C_F)
        try:
            summ, x_end = cy.measure_cycle(s, x, t_transient=600.0, max_chunks=30)
            status = "ok"
        except RuntimeError:
            summ, x_end, status = None, x, "not settled"
        if summ is None:
            rows.append({"L_H_L": float(L), "gamma_amp_deg": 0.0, "psi_amp_deg": 0.0, "status": status})
            x = np.asarray(x_end) + SMALL_KICK
        else:
            mu = cy.floquet_multipliers(s, np.array(summ["x0"]), summ["period"])
            rows.append({"L_H_L": float(L), "gamma_amp_deg": summ["gamma_amp_deg"], "psi_amp_deg": summ["psi_amp_deg"],
                         "period": summ["period"], "status": status,
                         "second_multiplier_abs": float(np.abs(mu[1]))})
            x = np.array(summ["x0"])
        print(f"C_F={C_F} {direction} L={L:.4f}: gamma {rows[-1]['gamma_amp_deg']:.2f} deg ({status})", flush=True)
    return {"C_F": C_F, "direction": direction, "L_cr": L_cr, "rows": rows}


def main():
    jobs = [(C_F, d) for C_F in (0.0223, 0.0247) for d in ("up", "down")]
    with Pool(4) as pool:
        results = pool.map(one_sweep, jobs)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=1))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
