# E03 — When does the ship stop trailing straight and sit at an angle? (pitchfork)

**Run:** `/home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E03_pitchfork.run` (about 1 minute)
**Code:** `fishtail/continuation.py` (equilibria, eigenvalues), `fishtail/model_rigid.py`
**Checks:** `tests/test_continuation.py` (equilibria, pitchfork position, linear growth, time-domain check)
**Output:** `figures/E03_pitchfork/fig01…fig07` (PNG), `data/E03_pitchfork/summary.json`

## Question
Where must the hawser be attached for the ship to trail straight behind the buoy? What happens if it is attached
too far aft?

## What is varied
Attachment point x_A/L from 0 (midship) to 0.6 (0.5 = bow). Rigid hawser L_H/L = 0.5. Both loads.
Friction: C_F = 0.05 (high, like a small tank model, so nothing fishtails) and C_F = 0.015 (like the real ship).

## What we expected
The straight position loses stability at x_Tcr/L (eq. 3.37) = 0.36 (100 %) / 0.19 (40 %) for any hawser length and
friction (derivation: `derivations/rigid_hawser_results.md`, det C changes sign there). Two mirror-image angled
positions appear.

## What we got

| Check | Result | Figure |
|---|---|---|
| Pitchfork position | x_Tcr/L = **0.359** (100 %), **0.193** (40 %), exactly as eq. 3.37 | fig01, fig02 |
| Type | **Supercritical**: the angled positions exist aft of x_Tcr and are stable (C_F = 0.05) | fig01 |
| Eigenvalue at the pitchfork | one real eigenvalue passes through 0 on both branches | fig02 |
| How the angle grows | **linearly**: log-log slope 1.00 for both loads, not ½ | fig06 |
| Time runs (100 %, C_F = 0.05) | x_A/L = 0.25: settles at ψ = 10.3°, γ = 53.4° (= predicted); x_A/L = 0.45: back to straight | fig05 |
| Low friction (C_F = 0.015) | 100 %: straight position fishtails everywhere forward of x_Tcr; angled branch loses stability by Hopf at x_A/L = 0.042 and 0.347. 40 %: straight is stable only for 0.193 < x_A/L < 0.215 | fig03 |
| Simos's static tests, 40 %, x_A/L = 0.20 (Table 5-2) | our model: ψ = 0 (x_A just forward of x_Tcr); measured ≈ 1° at 0.2 m/s ✔ but 13–17° at 0.4–0.6 m/s ✘ | fig07 |
| Same with Simos's wave correction k_w (Table 5-3) | 10.9°, 11.5°, 13.8° vs measured 13°, 14°, 17° | fig07 |

Figures:
- **fig01** pitchfork diagram: rest heading vs x_A/L, both loads (stable solid, unstable dashed, ▲ pitchfork).
- **fig02** real parts of the eigenvalues along both branches (100 %).
- **fig03** the same diagram at full-scale friction, with Hopf points ● on the branches.
- **fig04** hawser angle and hawser force at the angled position for three frictions (heading does not change).
- **fig05** time runs and top views either side of x_Tcr.
- **fig06** why the branches leave in a straight line (normal form ẋ = rx − x|x|).
- **fig07** comparison with Simos's static-bifurcation tests and his wave correction k_w.

## What it means physically
- **Weathervane.** The water's sideways force acts near midship. If the hawser is tied well forward of it, that force
  is like a weathervane's tail and swings the ship back into line. Tied aft of x_Tcr, the Munk moment (E01) wins and the
  ship sheers off to one side, like a kite.
- **Where it sits.** At rest the water's turning moment *about the attachment point* must be zero. That one equation
  fixes the heading, so the heading depends only on x_A/L, not on hawser length or friction (fig04a). The hawser angle
  does depend on friction: at low friction the sideways lift pulls the ship far out to the side (γ up to 70°, fig04b),
  and the hawser force can be many times the straight-line value (fig04c).
- **Strogatz's normal form, with a twist.** Strogatz's pitchfork ẋ = rx − x³ grows like √r. Here the restoring
  term is the cross-flow drag, which grows like ψ|ψ| (quadratic, but odd), so the normal form is ẋ = rx − x|x| and the
  angle grows like r: ψ ≈ c₂Δ/q with c₂ = (πT/L)(1 + 0.4 C_B B/T) and q = l_p C_Y + x_A(C_Y − πT/2L) (fig06; good to
  about 20°). Same symmetry, same stability exchange, different exponent.
- **Low friction (real ship).** The pitchfork is still there, but the branches lose stability to oscillations (Hopf):
  for the loaded ship there is no attachment point at L_H/L = 0.5 where the real ship trails straight at C_F = 0.015.

## Things that look odd (flagged)
1. **Simos's static tests at higher speeds are not reproduced** (fig07). Simos himself explains this with wave-making
   (Froude) effects that the hydrodynamic model does not have, and fits a correction k_w (Annex C). With his k_w our
   static model gives angles 2–3° below the measured ones. Since k_w was fitted to these same tests, the remaining
   difference may come from how he fitted it (not described in a form we could reproduce). Open issue 17.
2. Open issue 11 (Leite's curve in Simos Fig. 3.21 seeming to cross zero at 0.38–0.39 for 100 %) is not resolved by
   this experiment: eq. 3.37 and our model agree on 0.359.
3. The angled positions only exist with the hawser in tension in our code; a rigid bar could also hold the ship
   upstream of the buoy (in compression). Those are not followed (not what a moored ship does).

## Not done
- Imperfect pitchfork (wind from one side). Optional; would need a wind-force model.
