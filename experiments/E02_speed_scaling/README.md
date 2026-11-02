# E02 — Does the current speed change what the moored ship does?

**Run:** `/home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E02_speed_scaling.run` (about 20 s)
**Code:** `fishtail/model_rigid.py` (`rhs` in seconds, `rhs_star` in t*), `fishtail/hydro.py` (`friction_coefficient`)
**Checks:** `tests/test_rigid.py::test_current_speed_drops_out_at_fixed_friction`
**Output:** `figures/E02_speed_scaling/fig01…fig03` (PNG), `data/E02_speed_scaling/summary.json`

## Question
Does a faster current make the ship behave differently, or just faster? (Decision D1,
`notes/01_why_current_speed_matters.md`.)

## What is varied
Current speed U, first with the friction coefficient C_F held fixed, then with C_F taken from the ITTC friction line
(eq. 3.20), which depends on U through the Reynolds number. Vidal de Negreiros 100 %, rigid hawser, x_A/L = 0.5.

## What we expected
- C_F fixed: the same motion, only faster; in t* = U t / L the runs are identical.
- C_F from the friction line: faster current → smaller C_F → less stable.

## What we got

| Check | Result |
|---|---|
| fig01: three speeds (0.5, 1, 2 m/s), C_F = 0.035 fixed, integrated **in seconds** | Different in real time; in t* the largest difference is **4 × 10⁻⁷ degrees** (round-off) |
| fig02: same set-up (L_H/L = 0.15), C_F from the friction line | Hopf boundary there: C_F = **0.0254**. 1:90 model at 0.2 m/s: C_F = 0.0288 → **settles**. Model at 0.4 / 0.6 m/s: 0.0247 / 0.0227 → **fishtails**, more at 0.6. Real ship at 1.9 m/s: C_F = 0.0092 → fishtails strongly (±8° heading) |
| fig03: growth rate vs full-scale speed, 0.25–3 m/s | x_A/L = 0.5: unstable for every hawser length ≥ 0.1 L at every speed; L_H/L = 0.05 is stable only below about 0.35 m/s |

## What it means physically
- Every force (water, inertia, hawser) grows like U², so a faster current just plays the same "film" faster. This is
  why all our bifurcation results are in t* and do not mention U.
- The one exception is skin friction, which gets relatively smaller at large Reynolds number. That friction drag is
  what keeps the hawser tight and pulls the bow back into line (the "gravity" of the pendulum), so less friction means
  less stability. Current speed enters the rigid model **only** through C_F.
- A 1:90 tank model has about 3× the friction coefficient of the real ship at the same Froude speed (E01 fig08). That is
  why Simos's model at 0.2 m/s was stable while the real ship would fishtail (fig02, Simos p. 195).

## Things that look odd (flagged)
- fig02 uses the 100 % ship; Simos's p. 195 statement is for the 40 % ship. The conclusion is the same, E06 fig02 shows
  both loads.
