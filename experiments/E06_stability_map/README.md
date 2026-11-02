# E06 — Where in the design space is the moored ship safe? (stability maps, main result)

**Run:** `/home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E06_stability_map.run` (about 3 minutes)
**Code:** `fishtail/continuation.py`, `fishtail/model_rigid.py`
**Checks:** `tests/test_continuation.py::test_hopf_curve_matches_simos_fig_5_3`; fig05 is the step-5 check
**Output:** `figures/E06_stability_map/fig01…fig05` (PNG), `data/E06_stability_map/summary.json`

## Question
For which combinations of attachment point, hawser length and friction (current speed) does the ship trail straight,
sit at an angle, or fishtail? Does our model reproduce Simos's stability boundary (Fig. 5.3)?

## What is varied
Two knobs at a time: (L_H/L, C_F) at fixed x_A/L; (x_A/L, L_H/L) at fixed C_F; (x_A/L, C_F) at fixed L_H/L.
Both loads.

## What we expected
Simos Fig. 5.3: stable above a Hopf curve in the (L_H/L, C_F) plane. A pitchfork line at x_A = x_Tcr (E03).

## What we got

| Check | Result | Figure |
|---|---|---|
| Hopf curve vs Simos Fig. 5.3, 100 %, x_A/L = 0.5 | mean difference **+0.5 %**, largest 3 % (17 digitised points) | fig01a |
| Hopf curve vs Simos Fig. 5.3, 40 %, x_A/L = 0.4 | mean difference **+5.6 %** (ours a little more stable everywhere) | fig01b |
| Eq. B.13 as printed (no lift from sideways velocity) | boundary 3–6 times too high: not what Simos computed | fig01c |
| Where tests and ship sit (cf. Simos Fig. 5.4) | 40 %: model at 0.2 m/s stable for every hawser length; real ship (C_F ≈ 0.015) unstable for L_H/L > 0.05 | fig02 |
| Regime maps | see below | fig03, fig04 |
| Time-domain check (step 5) | **12 of 12** simulations end in the regime the map predicts | fig05 |

Fraction of the (x_A/L 0–0.6, L_H/L 0.05–1.5) plane in each regime (fig03):

| | trails straight | sits at an angle | fishtails |
|---|---|---|---|
| 100 %, C_F = 0.035 (tank-model-like) | 11 % | 56 % | 32 % |
| 100 %, C_F = 0.015 (full-scale-like) | 1 % | 53 % | 46 % |
| 40 %, C_F = 0.035 | 39 % | 33 % | 28 % |
| 40 %, C_F = 0.015 | 5 % | 32 % | 62 % |

"Fishtails" here means *no stable rest position*: the straight position and any angled ones are all unstable.

Figures:
- **fig01** our Hopf boundary against Simos Fig. 5.3 (digitised), both loads; and the "B.13 as printed" variant.
- **fig02** the boundary with the tank-test friction bands (Table 5-1) and full-scale friction lines.
- **fig03** regime maps in (x_A/L, L_H/L), two loads × two frictions.
- **fig04** regime maps in (x_A/L, C_F) at L_H/L = 0.5.
- **fig05** the map with 12 time-domain runs on it (marker = what the simulation did).

## What it means physically
- **Short hawsers are safer** (low curves at small L_H/L): with a short hawser the pendulum is stiff and the swing is
  well damped. For the loaded ship the boundary rises steeply up to L_H/L ≈ 0.5 and then flattens near C_F ≈ 0.048.
- **The real ship is on the wrong side.** Full-scale friction (0.009–0.015) is far below the boundary for practical
  hawser lengths, so the real ship on a rigid hawser at the bow fishtails, while 1:90 tank models (C_F 0.02–0.05) can
  look stable. This is the scale effect Simos describes (p. 194–195).
- **Forward is not always better** (fig04). For the ballast ship there is a narrow window just forward of x_Tcr
  (x_A/L ≈ 0.19–0.21) where the straight position is stable down to very low friction; further forward it fishtails.
  For the loaded ship the Hopf curve is highest near x_A/L ≈ 0.49, so the bow is close to the *least* stable
  attachment point at a given friction. We have not yet found a simple physical explanation for the shape of these
  curves; it comes out of the eigenvalues and is confirmed by the simulations in fig05.
- For the ballast ship at tank-model friction the fishtailing region is a closed "island" (fig03c): attaching right at
  the bow is stable again. This fits Simos's footnote (p. 192) that x_A/L = 0.5 gave only a weak instability for 40 %,
  which is why he tested at 0.4.
- **Cross-flow coefficients do not move the boundaries.** In the linearisation C_Y and l_p appear only at second order,
  so the Hopf curve is exactly independent of the shape of C_D(x) (our straight-line simplification, D6). They only
  affect the size of the swing (E05). This answers part of open issue 16.

## Things that look odd (flagged)
1. **40 % boundary 5.6 % above Simos's.** Not from the I_z reference point (D13; changing it moves the curve by
   < 0.3 %). A 3 % change in the draft T would remove the difference, so it is within the uncertainty of the hull
   data and of how the lift slopes are normalised. Cause unknown. New open issue 18.
2. The narrow vertical strips of "fishtails" at x_A/L < 0.08 (fig03a, b) are where the angled position itself loses
   stability by a Hopf bifurcation (E03 fig03): the ship would swing about an angled position. Not studied further.
3. The regions are computed from eigenvalues on a 61 × 59 grid, so their edges are blocky; the black Hopf curve is
   the accurate boundary.
