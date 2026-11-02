# E04 — When does the ship start to fishtail? (Hopf bifurcation)

**Run:** `/home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E04_hopf.run` (about 1 minute)
**Code:** `fishtail/model_rigid.py` (`linear_jacobian_star`), `fishtail/continuation.py` (`hopf_C_F`, `hopf_frequency`)
**Checks:** `tests/test_continuation.py::test_hopf_point_is_a_complex_pair_crossing`, `…matches_time_domain`
**Output:** `figures/E04_hopf/fig01…fig06` (PNG), `data/E04_hopf/summary.json`

## Question
As the friction falls (= faster current, D1), when does the straight position stop damping out a swing, and what
does the ship do then?

## What is varied
Friction coefficient C_F (0.004–0.1), for hawser lengths L_H/L = 0.2, 0.5, 1.0. Vidal 100 %, x_A/L = 0.5 (bow).

## What we expected
A complex pair of eigenvalues crosses the imaginary axis (Hopf, Strogatz §8.2); the damped swing becomes a
self-sustained swing (limit cycle). Simos found Hopf boundaries (Fig. 5.3) with limit cycles in the unstable region.

## What we got

| Check | Result | Figure |
|---|---|---|
| Hopf friction, L_H/L = 0.2 / 0.5 / 1.0 | C_F = **0.0298 / 0.0425 / 0.0478** | fig02 |
| Frequency at the Hopf point (L_H/L = 0.5) | ω = 0.738 U/L, period 2π/ω = **8.51 L/U** | fig01, fig02 |
| Period vs friction | close to ∝ 1/√C_F, like a pendulum | fig02b |
| Before (C_F = 1.15 × Hopf) | stable focus: swings from 2° and 30° both spiral in | fig03, fig04 |
| After (C_F = 0.85 × Hopf) | from 2° the swing grows, from 30° it shrinks: both reach the **same** cycle (γ ±18°) | fig03, fig04 |
| Fully grown cycle, C_F = 0.6 × Hopf | γ ±43°, ψ ±13°, period 12.0 L/U; bow moves on a circle, stern traces a figure-8 | fig05 |
| Hawser force over a cycle | 0.85×: 0.87–1.22 F_H0; 0.6×: 0.42–2.19 F_H0; 0.4×: **−0.015**–3.42 F_H0 | fig06 |
| Time-domain check either side of the crossing | decays at 1.1×, grows at 0.9× | test |

Figures:
- **fig01** the four eigenvalues in the complex plane as C_F changes (colour = C_F); the pair crosses at ○.
- **fig02** growth rate and period of the swing against C_F, three hawser lengths.
- **fig03** phase portraits before and after the Hopf point (two 2-D views of the 4-D motion).
- **fig04** time series of γ, ψ and hawser force before and after.
- **fig05** top view of one fishtailing cycle (8 snapshots, path of bow and stern).
- **fig06** hawser force during fishtailing for three frictions.

## What it means physically
- **The pendulum.** Trailing straight, the hawser tension equals the friction drag. When the ship swings sideways,
  that tension pulls it back, like gravity on a pendulum (Simos §5.1). Less friction → weaker pull → slower swing
  (period ∝ 1/√C_F, fig02b) and less damping.
- **What drives the swing.** Each sideways swing turns the bow a little off the current; the Munk moment then pushes
  the bow further, and the lift on the hull pushes the ship sideways. Below the Hopf friction this push puts more energy
  into the swing per cycle than the water's damping takes out, so a small swing grows. It stops growing when the
  quadratic cross-flow drag (which grows faster than the push) balances it: a limit cycle.
- **Supercritical.** Small and large starting swings reach the same cycle and the cycle grows from zero at the Hopf
  point (E05 confirms: no hysteresis). So the model gives a supercritical Hopf, unlike the subcritical one seen in
  Aghamohammadi & Thompson's flume tests (R6). We report what the model gives (brief §3).
- **Snatch loads.** During fishtailing the hawser force swings between a fraction and several times its
  straight-line value, twice per cycle. At low friction (0.4 × Hopf) it briefly goes negative: the rigid bar would
  push. A real rope would go slack there and then snatch tight — the stretchy-hawser model (E07) is needed for that.

## Things that look odd (flagged)
1. Simos's Annex B (eq. B.13) writes the hull lift with the heading ψ only, leaving out the lift caused by the sideways
   velocity of the hull. Done that way, the Hopf boundary is 3–6 times too high (E06 fig01c). With the full
   linearisation (sympy, `derivations/rigid_hawser_results.md`) we reproduce his Fig. 5.3, so the omission is only in
   the printed text, not in his computation.
2. F_H < 0 at C_F = 0.4 × Hopf (L_H/L = 0.5): flagged as decision D10 asks. Simos's experiments used rigid bars that can
   push (p. 222), so the comparison with him is still fair.
