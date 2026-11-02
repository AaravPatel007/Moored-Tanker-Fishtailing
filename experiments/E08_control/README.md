# E08 — Can we stop the fishtailing? Hold-back tug and stern-thruster PD control

**Run:** `/home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E08_control.run` (about 20 minutes)
**Code:** `fishtail/control.py` (PD law, clip), thruster and tug terms in `fishtail/model_rigid.py` (equations of
motion and `linear_matrices`), rest positions in `fishtail/continuation.py`
**Checks:** `tests/test_control.py` (9 checks)
**Output:** `figures/E08_control/fig01…fig09` (PNG), `data/E08_control/summary.json`. Demo: the **Tug** and
**Thruster** tabs (D5).

## Question
Without control the base case (rigid hawser, bow attachment x_A/L = 0.5, L_H/L = 0.5, C_F = 0.03, like a 1:90 tank
model at 0.17 m/s) fishtails: γ swings ±33°, ψ ±11°, period 10.7 t*. Can (a) a hold-back tug or (b) the tanker's own
stern thruster stop it, and what does a thrust limit change? Decision D31 (choices agreed with the student).

## The two options
- **Hold-back tug:** a steady pull F_tug downstream at the stern (real SPM/FPSO practice to keep the hawser taut;
  secondary sources, see D31). Code: the wind force with β_W = 0 at x_W = stern. Not feedback, but its sideways part
  grows with ψ by geometry (a weathervane tail) and it adds to the hawser tension.
- **Stern thruster:** side force at x_T/L = −0.5, PD law N_c* = −k_p ψ − k_d r* (units ½ρU²L²T), clipped at
  thrust_limit × F_H0. Least-thrust gains ★: smallest peak |Y_T| (switched on during the full fishtail, worst of 4
  phases, no limit) with decay rate σ = −max Re λ ≥ σ_min = 0.05 U/L (PROPOSED; a small swing halves in 14 t*).
- All forces given as fractions of F_H0 (PLACEHOLDER sizes; open issue 25). Real-ship numbers: fig09 only, at the
  real ship's own friction.

## What we got

| Check | Result | Figure |
|---|---|---|
| Tug at the base case | **Hopf at F_tug = 0.301 F_H0**; the swing shrinks linearly to zero (same linear law as E05) | fig01 |
| Tug against friction | needed pull grows fast as C_F falls: ≈0.3 at C_F 0.03, ≈2.2 F_H0 at C_F 0.01 (real ship ≈1 m/s) | fig01 |
| Thruster root locus | Hopf at **k_d = 0.0064** (k_p = 0) or **k_p = 0.0131** (k_d = 0). k_d pulls the pair left until it meets the real axis; more k_d **overdamps** (one real λ moves back right) | fig02 |
| A limit of heading control | one eigenvalue, **−0.2205 U/L**, does not move with the gains (a slow sideways drift with almost no heading change), so σ ≤ 0.22 whatever the gains | fig02, fig03 |
| Gain map and ★ | ★ = **k_p 0.004, k_d 0.012**: σ = 0.055, peak thrust **0.082 F_H0**. Thrust needed grows mainly with k_p | fig03, fig04 |
| Time runs | ★ stops the full fishtail in ~70 t*; 4 × ★ in ~15 t* but needs 0.22 F_H0 at switch-on | fig05 |
| Saturation (★ gains) | limit ≥ 0.0136 F_H0: always stops. **0.0063–0.0113 F_H0: bistable**: the full fishtail survives (20–28°), a 3° nudge dies out. ≤ 0.0052: even the nudge grows into a swing | fig06 |
| Stability map | tug 0.5 F_H0 lowers the Hopf C_F from 0.0425 to 0.025 (L_H/L = 0.5); the ★ thruster much more, and for L_H/L ≥ 0.6 removes fishtailing at every C_F ≥ 0.0001 | fig07 |
| Bow vs stern thruster | stern: 98 % of the gain map stable; bow (x_T/L = +0.45): 83 %, and the damping gain hardly acts (the bow is next to the hawser point) | fig08 |
| Real ship, fixed kN settings | tuned at 1 m/s (thruster) or set at 94 kN = 2.7 F_H0 (tug): the tug fails above ≈1.1 m/s, the thruster above ≈1.4 m/s (forces scale with U², the settings do not) | fig09 |

## What it means (Strogatz)
- Both controls move the same Hopf bifurcation: the knob is now the tug pull or the gain.
- The eigenvalues describe small swings only. With a thrust limit the straight position can be stable while the
  big fishtail is also stable (bistability, Ch. 8): which one the ship ends in depends on where it starts.
- The −0.22 eigenvalue shows that what you measure limits what you can control: the heading does not see the drift.

## Caveats
- PLACEHOLDER sizes for the thrust limit, tug pull and positions (D31, open issue 25); σ_min to be agreed.
- The tug line is assumed exactly downstream; a real tug pulls at an angle and moves.
- Real thrusters lose effect when the ship moves through the water (not modelled).
