# Progress summary — read this when you come back

Last updated: 2026-11-02 (step 8 control code; see the 2026-11-02 section below). Original deadline target
~2026-10-31 has passed — check the real hand-in date with the student.
New to the project? Start with [../README.md](../README.md).
Plain-words version of all the work so far: [project_overview.md](project_overview.md).
All the equations, explained simply, plus the parameter list: [03_mathematics_explained.md](03_mathematics_explained.md).

## Where we are
Steps 1–7 are **done** (core model, bifurcations, limit cycles, stability maps, stretchy hawser), plus the wind
extension (E10) and slider demos D1–D4. Step 8 (control) is **half done**: the code and its checks exist
(2026-11-02), the experiment E08 and demo D5 do not. Then the report. All 141 automatic checks pass.
Nothing is committed to git yet.

## Decisions you made
- **Hull:** Vidal de Negreiros VLCC (320 m), 100 % and 40 % loading (D3).
- **Sideways drag along the hull:** a straight line, more drag at the stern (D6).
- **Slider demos:** Python computes, a web page plays back (D8).
- **Figures:** PNG only (D21).

Decisions taken today that you have not looked at yet: D17–D23 in [00_decision_log.md](00_decision_log.md).
D23 (typical parameter values) is marked PROPOSED.

## What steps 3–6 produced

| Step | Code | Experiment | Key result |
|---|---|---|---|
| 3 rigid-hawser model | `derivations/rigid_hawser.py`, `fishtail/model_rigid.py` | E02 | Hawser sign settled (B.3, D17). Current speed matters only through C_F (difference 4e-7° at fixed C_F) |
| 4 equilibria | `fishtail/continuation.py` | E03 | Pitchfork at x_A/L = 0.359 (100 %) / 0.193 (40 %), supercritical; the angle grows **linearly** past it |
| 5 bifurcations, maps | `continuation.py` | E04, E06 | Hopf at C_F = 0.0425 (L_H/L = 0.5, bow); **Simos Fig. 5.3 reproduced** (+0.5 % / +5.6 %); 12/12 simulations agree with the maps |
| 6 limit cycles | `fishtail/cycles.py` | E05 | Supercritical Hopf, no hysteresis; amplitude ∝ distance past the Hopf point; periods match Simos Fig. 5.7 within 1–3 s |
| demos | `interactive/build_demos.py`, `interactive/index.html` | D1–D3 | Slider + ship animation + phase plane + bifurcation diagram + eigenvalues + time trace. Online: https://aaravpatel007.github.io/Moored-Tanker-Fishtailing/interactive/ |

Checks: `tests/test_rigid.py` (sympy vs code, Jacobian vs finite differences, kinematics, pitchfork),
`tests/test_continuation.py` (equilibria, Hopf, Fig. 5.3, time-domain checks), `tests/test_cycles.py` (shooting,
Floquet multipliers, period near the Hopf point).

## Things found today that the report should explain
1. **Simos's Annex B (eq. B.13) as printed leaves out a lift force**; his Fig. 5.3 was clearly computed with it (D19).
2. **Linear, not square-root, growth** of both bifurcations. Reason: the restoring nonlinearity is quadratic drag
   (v|v|), so the normal forms are ẋ = rx − x|x| and dr/dt = μr/2 − c r². Shown with a toy oscillator (E05 fig02).
3. **The hawser can need to push** at low friction during fishtailing (E04 fig06): a real rope goes slack → step 7.
4. **Real ship vs tank model:** full-scale friction is far below the stability boundary; tank models look stable
   because of their higher friction (E02, E06 fig02).
5. **Attaching further forward is not always more stable** (E06 fig04).

## Things still open (details in [open_issues.md](open_issues.md))
- 40 % Hopf boundary 5.6 % above Simos's (issue 18). Water ρ, ν still PLACEHOLDER (issue 10).
- Our model gives a supercritical Hopf; a flume study (R6) saw a subcritical one (issue 20).
- Wave (Froude) effects in Simos's 40 % tests (issue 17); effect of the C_D(x) shape on swing size (issue 16).
- The 8 old PDF figures in `figures/E01_hydrodynamics/` can be deleted if you agree.

## Step 7 (done after the above)
Stretchy hawser, `fishtail/model_elastic.py`, E07. Rope law from a real SPM nylon hawser (Huang et al. 2026):
T = 10 580 kN · strain^1.86 (D25). Stiff rope reproduces the rigid eigenvalues. Boundaries move < 5 %. Real ship at
1 m/s: rope slack 74 % of the time, snatch peak 501 kN (14 × steady; rigid model 172 kN). 113 checks pass.

## 2026-10-18: water properties, wind (E10, demo D4)
- **Water properties** now cited from ITTC (2011), R15 (D27). No numbers changed.
- **D23 agreed. No stretchy-rope demo tab** (D26).
- **Wind** (student's request, D28 PROPOSED): constant force, fixed direction, at midship; two knobs, force
  F_w/F_H0 and direction β_W. `fishtail/model_rigid.py` (C_W, beta_W, x_W_L; zero wind = old model exactly),
  `fishtail/continuation.py` (wind equilibria, whole diagram from a formula), `tests/test_wind.py` (17 checks),
  E10 (8 figures), **Wind tab** in the demo page (D4). 130 checks pass.
- E10 results: imperfect pitchfork with folds; at x_A/L = 0.25 the far side vanishes at ±0.153 F_H0 (beam wind),
  hysteresis confirmed by sweeps; cusp width ∝ distance^2.12 (x|x| form predicts 2, Strogatz's cubic 3/2).
  **Surprise, flagged:** the ship ends on the UP-wind side (bow weathervanes into the wind, then hull lift wins).
  Wind from astern makes fishtailing start sooner, from ahead or the side later.

## 2026-10-18 (later): demo page feedback (D30)
- Pitchfork tab: "Start the ship" bar (nudge to port or starboard, or rest on an angled branch; click a branch too).
- Eigenvalue panel redone: axes fitted to the data, fast eigenvalue in a table instead of off-scale labels, ▲/● mark
  the crossing points, a table + one sentence say why the position is (un)stable and whether it is pitchfork- or
  Hopf-type. Bifurcation-diagram axes and colours explained per tab. Pitchfork, fishtailing, hawser-length data
  rebuilt (starboard runs are exact mirrors of the port runs). 132 checks pass. Same private link, version 6.
- Second round (D30, end): panels less wordy (one-line verdict, short legends, text in "How to read it"), three
  equal detail graphs, wind force slider 3× finer (31 values, runs trimmed once settled → 10 MB), slow wind runs
  extended to 600 t*, playback 1× = 1 t* per second with 1/2/4/8× plus real-time line, floating mini top view.

- Last change (artifact version 10): the floating ship view is **off by default**; the "Show ship view" button in the
  Detail heading turns it on; it hides once more than 60 % of the Detail graphs are scrolled past. Not yet tested by
  scrolling in a real browser (headless Chrome cannot scroll). Playback is back to 10 t*/s with ½×/1×/2×/4×.

## Where to pick up (end of 2026-10-18 session; superseded by the 2026-11-02 sections below)
- Demo page (D30) is finished unless the student reports problems after trying it.
- Student still to check: scrolling/drag/resize of the ship view; D28 (wind at midship or aft, open issue 24).
- Next real work: step 8 (yaw thruster control) or the report. Do not start without the student's go-ahead.

## 2026-11-02: step 8 control — code and checks (D31)
- **`fishtail/control.py` (new):** PD heading control by a stern thruster, N_c* = −k_p ψ − k_d r*, made by a side
  force Y_T = N_c / x_T at the stern (so it also pushes sideways), clipped to |Y_T| ≤ thrust_limit · F_H0.
- **`fishtail/model_rigid.py`:** new fields `k_p`, `k_d`, `thrust_limit` (default ∞), `x_T_L` (default −0.5).
  Thruster force/moment added to the equations of motion and, as linear terms, to `linear_matrices`.
  `linear_matrices` now also accepts a steady force along the current line (sin β_W = 0).
- **Hold-back tug:** no new code; it is the wind force with β_W = 0 at x_W = stern (`params.X_TUG_L = −0.5`).
- **`fishtail/continuation.py`:** equilibria include the thruster's force at rest; symmetric cases (tug, thruster,
  head-on wind) still give straight + angled mirror pairs.
- **`fishtail/params.py`:** `X_T_L`, `X_TUG_L` — both PLACEHOLDER (open issue 25).
- **`tests/test_control.py` (new, 9 checks); `tests/test_wind.py`** docstring points to it. **141 checks pass.**
- E08 and demo D5 (Tug and Thruster tabs) were completed later the same day; results in
  [../experiments/README.md](../experiments/README.md) and the E08 README.

## 2026-11-02 (later): E08 and demo D5 done
- **Choices (student, D31):** both a hold-back tug and a stern thruster (no bow thruster); non-dimensional gains plus
  N / kN conversion; optimal = least thrust with σ ≥ σ_min (0.05 U/L PROPOSED); thruster tab = map of (k_d, k_p) +
  3-button thrust limit (none / 0.02 / 0.01 F_H0, PLACEHOLDER).
- **E08 (9 figures):** tug Hopf 0.301 F_H0 at C_F 0.03, ≈2.2 F_H0 at real-ship friction; thruster Hopf k_d 0.0064 /
  k_p 0.0131; ★ k_p 0.004, k_d 0.012 (0.082 F_H0); eigenvalue −0.22 cannot be moved (σ ≤ 0.22); overdamping;
  bistable for limits 0.0063–0.0113 F_H0; bow thruster much weaker; real ship with fixed kN settings fails at
  faster currents (tug above ≈1.1 m/s, thruster above ≈1.4 m/s).
- **Demo:** two new tabs. Data built in parallel (10 cores): tug 6.3 MB; thruster 3.5 + 6.8 + 10.3 MB (one file per
  limit). Same private link, version 11.
- **Student to check:** σ_min = 0.05; the tabs in a real browser. Still PLACEHOLDER: thrust/tug sizes (issue 25).

## 2026-11-03: demo page text and range fixes (student's feedback)
- Pitchfork slider and diagram now end at x_A/L = 0.5 (the bow); further forward is off the ship. Wind tab's
  "ψ vs attachment point" axis also ends at 0.5. Pitchfork data rebuilt; wind data only had its axis range changed.
- Number fixes in the page text: Hawser-length tab C_F = 0.03 ↔ 1:90 model at 0.17 m/s (was "0.3 m/s"); control
  check runs on for 300 t* (was "800"); fast eigenvalue below −1 U/L (was "−1.7"; Tug tab reaches −1.02); real-ship
  time scale shown only for C_F < 0.021 (ITTC line at 1 cm/s).
- "How to read it" shortened in every tab.

## Next steps
1. **Student to check:** σ_min = 0.05 (D31), the Tug and Thruster tabs, D28 (wind at midship or aft? issue 24).
2. **Re-run E08** (`$PY -m experiments.E08_control.run`, ~20 min): only to redraw the labels of fig03 (★ legend) and
   fig07 (× = no fishtailing at any C_F); the fixes are already in `run.py`, the numbers do not change. Stopped
   on 2026-11-02 at the student's request after fig02; summary.json and all figures from the first run are intact.
3. **Optional:** cited thrust / tug values (issue 25) to replace the PLACEHOLDERs.
4. **Report and slides** (step 10). Step 9 (forcing) stays cut unless time is left over.
5. Optional: open issue 22.

## How to run things (from the project root)
```bash
PY=/home/aarav-patel/wec_rl/.venv/bin/python
$PY -m experiments.E03_pitchfork.run           # any experiment: E01 ... E07, E10 (E05, E06 take a few minutes)
$PY -m derivations.rigid_hawser                # sympy derivation -> derivations/rigid_hawser_results.md
$PY -m interactive.build_demos                 # demo data (several minutes); then open interactive/index.html
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 $PY -m pytest # all checks
```
