# E05 — How big is the fishtailing swing, and is the Hopf bifurcation super- or subcritical?

**Run:** `/home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E05_limit_cycles.run` (about 8 minutes)
**Code:** `fishtail/cycles.py` (simulate-until-settled, shooting, Floquet multipliers)
**Checks:** `tests/test_cycles.py`
**Output:** `figures/E05_limit_cycles/fig01…fig06` (PNG), `data/E05_limit_cycles/summary.json`

## Question
Past the Hopf point (E04), how large does the swing become, how long does one swing take, and does the swing appear
smoothly (supercritical) or suddenly with hysteresis (subcritical)?

## What is varied
Friction C_F swept **down** (starting from a tiny disturbance each time) and then **up** (starting from a large swing),
each run continuing from where the previous one ended, like a tank test. Hawser length L_H/L for the comparison with
Simos. Vidal 100 %, x_A/L = 0.5.

## What we expected
Strogatz §8.2: supercritical → amplitude grows from zero like √(distance past the Hopf point), same curve both ways;
subcritical → a jump to a large swing and hysteresis.

## What we got

| Check | Result | Figure |
|---|---|---|
| Sweep down vs sweep up (L_H/L = 0.5) | **the same curve**: no jump, no hysteresis → **supercritical** | fig01 |
| Size of the swing | grows from 0 at C_F = 0.0425 to γ ±72°, ψ ±15° at C_F = 0.013 | fig01 |
| Hawser force over a cycle | peaks up to 4.2 F_H0; the dip goes **below zero** (bar pushing) for C_F ≲ 0.017 | fig01c |
| Growth law near the Hopf point (shooting, 0.5–20 % past it) | log-log slope **0.98** (γ) and **0.97** (ψ): amplitude ∝ distance, **not** √distance | fig02a |
| Toy oscillator ẍ − μẋ + D(ẋ) + x = 0 | D = ẋ\|ẋ\|: amplitude = 3πμ/8 (∝ μ); D = ẋ³: amplitude = √(4μ/3). Simulations match both averaging formulas | fig02b |
| Floquet multipliers at 0.95, 0.8, 0.6, 0.4 × Hopf value | one multiplier = 1.000 each time; the others inside the unit circle → stable cycles | fig03 |
| Period vs Simos Fig. 5.7 (model at 0.4 m/s, L_H/L = 0.36 / 0.6 / 0.9) | within **±1 s** for 0.36 and 0.6, **+2 to +3 s** for 0.9 (periods 96–141 s) | fig04 |
| Amplitude vs hawser length vs Simos Fig. 5.17 (C_F = 0.0223 and 0.0247) | ours within ~2–4° of Simos's theory curves everywhere; within a few degrees of his tank tests for L_H/L ≥ 0.5. Fishtailing starts at L_H/L ≈ 0.12–0.14 (tests: "about 0.10", p. 214). Largest γ swing ±53° near L_H/L = 0.3 | fig05 |
| Shape of one cycle (L_H/L = 0.6, C_F = 0.0223) | γ ±48°, almost a pure sine; ψ ±15°, clearly not a sine (sharp peaks, a shoulder in between). Period 124 s at model scale | fig06 |

Figures:
- **fig01** amplitude of γ and ψ and the hawser-force range against C_F, sweeping down and up.
- **fig02** the growth law near the Hopf point (log-log), and the toy oscillator that explains it.
- **fig03** Floquet multipliers in the complex plane.
- **fig04** fishtailing period at tank-model scale against Simos's simulated periods (Fig. 5.7).
- **fig05** amplitude against hawser length against Simos's theory and his tank tests (Fig. 5.17).
- **fig06** one cycle: γ(t) and ψ(t) normalised, and the cycle in the (γ, ψ) plane.

## What it means physically
- **Supercritical.** The ship does not jump suddenly into a large swing: just past the Hopf point the swing is small
  and it grows steadily as the friction falls (faster current). Coming back, it shrinks along the same curve.
- **Why "∝ distance" instead of "∝ √distance".** Near the Hopf point the motion is almost a sine wave of slowly
  changing amplitude r. Averaging the forces over one swing (Strogatz §7.6) gives dr/dt = (growth rate) r − (what
  limits the swing). The limiting term comes from the cross-flow drag, which is quadratic in the sideways speed
  (v|v|), so it averages to c r², not c r³. Setting dr/dt = 0 gives r ∝ growth rate ∝ distance past the Hopf point.
  Strogatz's √ law assumes a smooth cubic term. Both versions are "supercritical Hopf"; only the exponent differs.
- **Floquet multipliers.** After one full swing a small kick to the cycle is multiplied by these numbers. One is always
  1 (a kick *along* the cycle only shifts the timing). The others are below 1, so the swing is self-correcting: that is
  what "stable limit cycle" means. Near the Hopf point the second multiplier approaches 1 (slow recovery = critical
  slowing down).
- **Pendulum period.** The period grows as the friction falls (∝ roughly 1/√C_F at small amplitude) and with hawser
  length, like a pendulum whose "gravity" is the friction drag. Matching Simos's periods checks the inertia terms
  (masses, added masses) of our model.

## Things that look odd (flagged)
1. **Critical slowing down made the plain simulation fail** within 1 % of the Hopf point (the swing was still
   creeping after 2000 t*). `measure_cycle` now falls back to the shooting method there, which finds the cycle
   directly (and checks that it is stable).
2. **The rigid bar has to push** (F_H < 0) at low friction (fig01c, below C_F ≈ 0.017 at L_H/L = 0.5). A real rope
   would go slack; open issue 19 → step 7.
3. The amplitude law relies on the cross-flow model, the least certain part of the hydrodynamics at small angles
   (open issue 21).
4. **Short hawsers (L_H/L = 0.19, 0.36): the tank tests swing less than both models** (fig05). Simos saw the same
   and explains it by the tests being close to the stability boundary, where the swing had not finished growing in
   the 3–4 cycles the tank length allowed (p. 209). Our model agrees with his theory there, not with the tests.
5. Simos (p. 209–210) describes ψ(t) as "saw-tooth-like". Ours is clearly non-sinusoidal (fig06) but with equal rise
   and fall times; we cannot tell from his printed time series how close the shapes are.
