# E01 — Water forces on the hull, and the bare-hull yaw check

**Run:** `/home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E01_hydrodynamics.run` (about 10 s)
**Code:** `fishtail/hydro.py` (forces), `fishtail/model_pinned.py` (pinned hull), `fishtail/params.py` (numbers)
**Checks:** `tests/test_hydro.py` (21 checks, 59 cases with the parameter variations; all pass)
**Output:** `figures/E01_hydrodynamics/fig01…fig08` (PNG), `data/E01_hydrodynamics/summary.json`

## Question
1. Does our code reproduce the published current-force model (Leite et al. 1998 + Simos 2001)?
2. Is the bare hull unstable in yaw when it faces the current? (If not, the forces are wrong.)

## What is varied
Angle of attack α (0–180°), heading ψ, yaw rate r, position along the hull x/L, pivot position x_T/L, current speed U.
Hull: Vidal de Negreiros VLCC, 100 % (loaded) and 40 % (ballast).

## What we expected
- Coefficient curves with the same shape and size as the published theory curves (Simos Figs. 3.7, 3.10, 3.11).
- Near head-on current, the **Munk moment** (turns the hull broadside) beats the hull/rudder lift, so a small heading
  error grows.
- A pivot further forward than x_Tcr (eq. 3.37) makes the hull point into the current; Simos gives x_Tcr ≈ 0.19 for 40 %.

## What we got

| Check | Result | Reference |
|---|---|---|
| C1C, C2C, C6C shapes | Same shape and size as the published theory curves, including the small dip of C2C at exactly 90° (it is in Leite's formula and in his Fig. 3.10) | Figs. 3.7, 3.10, 3.11 (visual comparison) |
| Head-on slopes | Exactly Clarke's linear derivatives (eqs. 3.28, 3.31) | test |
| Moment near head-on, 100 % | dC6C/dψ = +0.211 (Munk) − 0.071 (lift + rudder) + 0 (cross-flow) = **+0.139 > 0 → unstable** | fig02 |
| Pinned hull, head-on | eigenvalues **+1.47 and −2.98** (units U/L) → **saddle** | fig03 |
| Pinned hull, end state | settles at **ψ = ±74.8°** (nearly broadside), eigenvalues −2.23, −2.41 (stable node) | fig03 |
| x_Tcr/L | **0.359** (100 %), **0.193** (40 %; thesis ≈ 0.19, p. 198 ✔) | fig04 |
| Yaw-rotating test (yaw-rate terms) | sway force: RMS difference **0.85 N** from Simos's theory points; yaw moment: **1.6 N m** (curves ±20 N, ±20 N m) | fig07, Simos Fig. 4.9 |
| Friction C_F, full scale 40 %, 1.90 m/s | **0.0150** (thesis 0.015, p. 195 ✔) | fig08 |
| Friction C_F, 1:90 model, 0.20 m/s | 0.0288 / 0.0469 (Table 5-1: 0.0279 / 0.0459, within 3 %) | fig08 |

Figures:
- **fig01** static coefficients C1C, C2C, C6C vs α, both loads; C1C at model scale and full scale.
- **fig02** the three parts of the yaw moment; zoom near head-on showing the total pushes the bow further away.
- **fig03** bare hull pinned at midship: time history and phase portrait (saddle at head-on, stable nodes at ±75°).
- **fig04** static yaw stiffness vs pivot position x_T/L; zero crossing = critical turret position (▲).
- **fig05** our straight-line sideways drag distribution C_D(x).
- **fig06** extra force and moment caused by turning, ship facing the current.
- **fig07** reproduction of Simos Fig. 4.9 (yaw-rotating test, model scale).
- **fig08** friction coefficient vs current speed, model and full scale, with Simos's Table 5-1.

## What it means physically
- **Why the bare hull turns away.** When the bow is turned a little off the current, the water flowing round the
  hull creates the Munk moment, which turns the hull further broadside. The hull + rudder act like a short wing whose
  lift pushes back, but only about a third as strongly (0.071 vs 0.211). Net result: head-on is a **saddle point**
  (Strogatz §5.2, §6.4): one direction grows (e-folding time ≈ 0.7 L/U), the other dies out.
- **Where the pinned hull ends up.** Held at midship, it swings to ±75°, nearly broadside, with the bow slightly
  upstream. There the Munk moment is balanced by the sideways drag, which acts behind midship (l_p aft), so the stern
  is pushed downstream. Strong yaw damping makes it creep in without overshooting (both eigenvalues real).
- **Why the hawser attachment point matters (preview of E03).** Holding the hull further forward adds a restoring
  moment: the sideways force F_Y acts behind the pivot like a weathervane's tail. Past x_Tcr the weathervane effect
  beats the Munk moment. x_Tcr is 0.36 L for the loaded ship but only 0.19 L in ballast. In ballast the hull is
  shallow and wide (B/T = 6.3 instead of 2.5), so its sideways lift force (the weathervane's "tail") is relatively
  much larger: in eq. 3.37 the bottom term 1 + 0.4 C_B B/T is 2.92 instead of 1.84.
- **Turning resists turning.** Ship facing the current, a yaw rate always produces a moment against it (fig06b):
  yaw damping, dN'/dr* ≈ −0.05 for small r*. It is mostly the wing (lift) term at fishtailing yaw rates.
- **Speed enters through friction only.** fig08: the model's C_F is 3× the full-scale value at the same Froude speed.
  This is why model tests are more stable than real ships (see `notes/01_why_current_speed_matters.md`).

## Things that look odd (flagged)
1. **Yaw moment in the yaw-rotating check is 1.5–2 N m too negative near beam-on** (ψ ≈ 60–110°, 250–300°), and
   fine at head-on and stern-on. That pattern fits the `v·r` cross-flow term, which depends on the *shape* of C_D(x):
   our straight line puts more drag near the ends than Simos's curve (his goes to zero at bow and stern, Fig. 4.4).
   This term is small at the small angles of fishtailing. It will be tested by the planned D6 comparison run.
2. **Table 5-1 is not self-consistent.** Its C_F values at 0.4–0.6 m/s lie on the k = 0.25 friction line, not k = 0.39
   as its caption says; only the 0.2 m/s row matches k = 0.39 (fig08). The full-scale 0.015 (p. 195) matches
   k = 0.39. Simos's simulations used measured C_T values anyway, so this matters only when we compare with his runs.
   Recorded in `notes/open_issues.md`.
3. **x_Tcr for 100 %:** eq. 3.37 gives 0.359; Leite's theory curve in Simos Fig. 3.21 seems to reach zero nearer
   0.38–0.39 (read by eye). The cause is unknown (perhaps slightly different hull data in Leite 1998, which we cannot
   read). To revisit in E03.
4. The kinks in our curves at exactly 90° and 270° (fig07) come from the |cos α| and |u| in the published formulas,
   not from the code.

## Not done
- Experimental data points (the ■ marks in Simos Figs. 3.7, 3.10, 3.11) were **not digitised**; the comparison with
  those figures is visual. Fig. 4.9's theory points were read by eye (±1.5) and are marked as such on fig07.
