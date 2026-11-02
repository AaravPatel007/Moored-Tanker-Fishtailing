# E07 — What changes with a real, stretchy hawser that can go slack?

**Run:** `/home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E07_elastic_hawser.run` (about 3 minutes)
**Code:** `fishtail/model_elastic.py` (6 states: midship X, Y, heading ψ, and u, v, r relative to the water)
**Checks:** `tests/test_elastic.py`
**Output:** `figures/E07_elastic_hawser/fig01…fig06` (PNG), `data/E07_elastic_hawser/summary.json`

## Question
A real nylon hawser stretches and cannot push. Does that change the stability boundary, the swing, and the rope loads?

## The rope law (decision D25)
Fitted to the measured tension–elongation curve of a real SPM nylon hawser (Huang et al. 2026, R11, Fig. 3;
70 m, 144 mm, breaking load 650 t, points read by eye):

    tension = a · strain^n,   a = 10 580 kN,  n = 1.86    (taut);    tension = 0  (slack, strain ≤ 0)

The same strain law is used for our 160 m hawser (L_H/L = 0.5), single line, no rope damping.

## What we got

| Check | Result | Figure |
|---|---|---|
| Rope fit | n = 1.86: the rope gets stiffer as it stretches. Our ship trailing straight: 10 / 35 / 129 kN at 0.5 / 1 / 2 m/s, strain 2.3 / 4.7 / 9.4 % | fig01 |
| Stiff-rope limit (1000 × stiffer) | swing eigenvalue 0.0495 ± 0.6217i vs rigid 0.0495 ± 0.6221i ✔ | fig02 |
| Stability boundary, real rope (U = 1 m/s) | almost the rigid curve: Hopf C_F up to 5 % higher at short hawsers (slightly less stable), within 1 % for L_H/L ≥ 1 | fig05 |
| 10 × softer rope | up to 16 % higher at short hawsers, 3 % lower at L_H/L = 2 | fig05 |
| Real ship, 1 m/s (C_F = 0.0100) | fishtails: γ ±44° (rigid: ±80°), ψ ±15.5° (saw-tooth), period about half the rigid one | fig03 |
| Rope loads at 1 m/s | rope **slack 74 % of the time**; snaps tight to **501 kN** = 14 × the steady 35 kN (rigid model: 172 kN, and −15 kN = pushing) | fig03 |
| Against current speed (0.25–2.5 m/s) | peak tension grows from 67 kN to 2 890 kN (≈ 45 % of the breaking load at 2.5 m/s); slack 54–79 % of the time | fig04 |
| Rope exponent | linear rope with the same tension at rest: peak 387 kN instead of 501 kN, swing larger (61° vs 44°) | fig06 |

## What it means physically
- **Boundaries hardly move.** At trailing straight the rope is taut and its stretch is small, so the linear
  stability (eigenvalues) is nearly the rigid one. The rigid model is a good tool for *where* fishtailing starts.
- **Loads change completely.** Once the ship fishtails, it swings towards the buoy, the rope goes slack, the ship
  drifts back with the current and the rope snaps tight. These snatch loads are several times what the rigid model
  predicts and grow roughly like U² (both water forces and the stretch grow). A rigid-rope model cannot give them.
- **Current speed now acts twice:** through friction (C_F) and through the rope (the same rope is "softer" relative to
  the water forces in a slow current).
- **The heading becomes a saw-tooth** (fig03): it turns slowly while the rope is slack and quickly when it snaps tight.
  Simos describes the measured ψ(t) of his tests as saw-tooth-like (p. 209–210), although his tests used rigid bars.

## Things that look odd (flagged)
1. **0.75 m/s does not follow the trend** (fig04a: γ ±37° between ±63° and ±44°). The swing there may settle on a
   different snap pattern; there could be more than one stable cycle. Not investigated (open issue 22).
2. **No damping in the rope.** Real nylon has hysteresis (energy lost every stretch cycle, see Michel & Davies 2023),
   which would lower the snatch peaks. Our peaks are therefore an upper-side estimate.
3. Single line with a strain law taken from a hawser for a smaller (220 m) tanker. Real VLCC terminals may use different
   hawsers; fig05 and fig06 show how sensitive results are to stiffness and exponent.
4. Rope mass and drag in the water are ignored (a straight, massless spring).
