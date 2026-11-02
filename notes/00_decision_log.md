# Decision log

Newest at the bottom. Each entry: what was decided, why, and where the evidence is.
**Where this log disagrees with `docs/PROJECT_BRIEF.md`, this log wins.**

Status words: **AGREED** (student said yes), **PROPOSED** (waiting for the student), **DONE** (carried out).

---

### D1 — Current speed is not a free knob in the rigid model (2026-10-17, AGREED)
- **Brief said:** use current speed `U_c` as the main bifurcation parameter for the rigid-hawser model.
- **Problem:** every force in the rigid model scales with `U²`, so after measuring time in units of `L/U` the speed
  cancels out — *except* through the friction coefficient `C_F(Re)`, which depends weakly on speed via the
  Reynolds number. Full explanation: [01_why_current_speed_matters.md](01_why_current_speed_matters.md).
- **Decision:** rigid-model knobs are `x_A/L` (attachment point), `L_H/L` (hawser length) and `C_F` (which stands for
  current speed through `Re`). Plots "against current speed" are allowed only if the x-axis is mapped through `C_F(Re)`
  and this is said on the plot. Current speed also enters through the stretchy hawser (stiffness ratio) and the
  controller (fixed gains).
- **Evidence:** Simos (2001) thesis, §5.1, Annex B (non-dimensional equations B.4–B.7), Fig. 5.3.

### D2 — Hydrodynamic model: Leite (1998) + Simos (2001) "quasi-explicit heuristic model" (2026-10-17, AGREED)
- Option B of the brief. Static forces from Leite et al. (1998); yaw-rate terms from Simos (2001).
- Needs only three hull numbers: `C_F` (friction), `C_Y` (sideways drag), `l_p/L` (where the sideways drag acts).
- All equations transcribed in [02_hydrodynamic_model.md](02_hydrodynamic_model.md) with equation numbers.
- Option A (MMG manoeuvring derivatives, KVLCC2) is **not** used for the main results: the KVLCC2 captive tests
  only cover drift angles up to ±16° (Dai & Li 2019, Table 2), and the full derivative table is in a paywalled paper.

### D3 — Hull choice (2026-10-17, AGREED: Vidal de Negreiros)
- Student approved KVLCC2 "if dimensions check out". They do (L=320 m, B=58 m, T=20.8 m, ∇=312 622 m³,
  `x_G`=11.2 m, `C_b`=0.810; added masses m'x=0.022, m'y=0.223, J'z=0.011 — Dai & Li 2019, Tables 1 and 5,
  quoting Yasukawa & Yoshimura 2015).
- **But:** the chosen hydrodynamic model (D2) needs `C_Y`, `l_p/L`, `M26`, `M66`, yaw inertia — none of which are
  published for KVLCC2 in a source we can read.
- **Proposal:** use the **Vidal de Negreiros** VLCC (same size class: L=320 m) from Simos (2001) Annex A:
  every number is published, `C_Y` and `l_p/L` were measured on this hull, and Simos's Fig. 5.3 gives a stability
  boundary to validate against. Keep KVLCC2 numbers in `sources.md` for reference only.
- **Student agreed (2026-10-17)** after asking what the ship is. Cross-flow coefficients used: 100 % → C_Y = 0.84,
  l_p/L = 0.059 (Leite, Tables 3-1/3-2); 40 % → C_Y = 0.52, l_p C_Y = 0.052 (IPT 1998, Table 4-5; Leite did not test
  this condition). For reproducing Simos's yaw-rotating figures, 100 % uses the IPT values C_Y = 0.86,
  l_p C_Y = 0.044 (Table 4-5), as he did.

### D4 — Equations of motion in velocities relative to the water (2026-10-17, AGREED)
- Use Simos (2001) eq. 3.9 / 5.1: added-mass terms on the left side, written with relative velocities `u, v`.
- The Munk moment `(M11 − M22) u v` is **already inside** the static yaw coefficient `C6C(α)` (footnote to eq. 3.8),
  so it is **not** added separately. This resolves review point 2b.

### D5 — No double-counting of cross-flow drag (2026-10-17, AGREED)
- In Simos's model the yaw-rate cross-flow terms are built to be **zero when r = 0** (eq. 4.29), so the static
  coefficients are recovered exactly. This resolves review point 2c.

### D6 — Sectional drag distribution `C_D(x)` (2026-10-17, AGREED, sign corrected)
- Simos builds `C_D(x)` from a digitised numerical curve (Aarsnes 1984, via Faltinsen 1990) times two step
  "weighting functions" chosen so that ∫C_D = `C_Y` and ∫C_D·x = `l_p C_Y` (eqs. 4.27–4.28, Fig. 4.3–4.4).
- **Decision:** use the simplest shape that satisfies the same two conditions exactly — a straight line,
  labelled as *our simplification*. Simos states (p. 190, §5.4) that the cross-flow coefficients have very little
  effect on stability and limit-cycle amplitude; we will check that with one comparison run (step 5).
- **Sign correction (step 2):** the line is `C_D(x) = C_Y [1 − 12 (l_p/L)(x/L)]` (more drag **aft**), not `1 + …` as
  first proposed. Evidence: (i) Simos's own C_D in Fig. 4.4 is larger aft (weighting K2 ≈ 1.45 aft, 1/K2 forward);
  (ii) the bracket in (4.29) only vanishes at r = 0 if ∫C_D x dx = −l_p C_Y L²; (iii) eq. 3.33 uses the arm
  (x_T + l_p/L). So the "+" printed in (4.28) is a sign slip. Resolves open issues 5 and 6.
- **Seen in E01:** the straight line gives a yaw moment 1.5–2 N m larger than Simos's near beam-on in the yaw-rotating
  test (fig07), consistent with it putting more drag near bow and stern than his curve.

### D7 — Python environment (2026-10-17, DONE)
- Use the existing `/home/aarav-patel/wec_rl/.venv` (student's request). Do not install into it.

### D8 — Interactive demos are "Python computes, HTML plays back" (2026-10-17, AGREED)
- Python runs the model over a grid of slider values → `interactive/<demo>/data.json`.
- A standalone HTML page shows slider + bifurcation diagram with a moving marker + time trace + top-view ship animation.
- Reason: one source of physics (no second model in JavaScript to keep in sync); works offline for the presentation.

### D9 — Scope order (2026-10-17, AGREED from brief)
- Steps 1–6 are the core. Cut step 9 (forcing) first — it is also close to the forced Van der Pol topic a classmate
  already has. Then cut 8, then 7.

### D10 — Rigid hawser can push; real rope cannot (2026-10-17, AGREED)
- The rigid model can give a negative (compressive) hawser force. Simos's experiments used rigid bars that *can* push
  (p. 222), so for comparing with him compression is allowed. For statements about real ropes, we compute the hawser
  force along every branch and cycle and flag any `F_H < 0`.

### D11 — Angle of attack computed with atan2 (2026-10-17, DONE — small, reversible)
- (3.13) writes α = π + arctan(v/u). We use α = π + atan2(v, u): identical for u > 0 (the moored ship always has
  u > 0), and also correct when the hull moves backwards through the water (current from astern → α = 0).

### D12 — cos α_a in the surge yaw-rate term (2026-10-17, DONE — small, reversible)
- (4.19) has cos[atan(−v/u + rL/2u)] and Simos says (p. 119) this must have the sign of u. A one-argument arctan
  cannot do that, so we use the two-argument angle: cos α_a = u / √(u² + (rL/2 − v)²). Same value for u > 0.

### D13 — Yaw inertia about midship (2026-10-17, DONE — small, reversible)
- Annex A gives I_z without saying about which point. √(I_z/M) = 80 m = 0.25 L exactly for both loads (the usual
  rule-of-thumb radius of gyration about the CG), so we read it as about the CG and use
  I_Zo = I_z + M x_CG² (+1.5 % for 100 %, +0.3 % for 40 %). Resolves open issue 7.

### D14 — "+M11 u r" and "+M26 u r" kept exactly as printed (2026-10-17, DONE)
- Thesis p. 116: Clarke's measured Y_r, N_r include the added-mass (ideal-fluid) force, which (3.9) already has on the
  left-hand side; these terms remove it. Confirmed numerically: our reproduction of the yaw-rotating test (E01 fig07)
  matches Simos's theory points (RMS 0.85 N on ±20 N). Resolves open issue 8.

### D15 — Bare-hull check = hull pinned at midship (2026-10-17, DONE)
- A free hull with no rope just drifts with the current, so "unstable in yaw" needs the hull held somewhere.
  We hold it by a vertical pin at midship (fishtail/model_pinned.py, 2 states). Result (E01): head-on is a saddle
  (eigenvalues +1.47, −2.98 in units U/L). Resolves open issue 9.

### D16 — Pitchfork position depends on the loading condition (2026-10-17, correction)
- The "x_A/L ≈ 0.19" quoted earlier is for the **40 % (ballast)** condition only (thesis p. 197–198, eq. 3.37).
  For the **100 %** condition eq. 3.37 gives **0.359**. E03's expected value is updated accordingly.

### D17 — Sign of the hawser side force: Annex B is right (2026-10-17, DONE — settles open issue 4)
- Our sympy derivation from the geometry (`derivations/rigid_hawser.py`) gives F_Y = **−**F_H sin(γ+ψ) and
  N = −x_A F_H sin(γ+ψ), as in eq. B.3. The "+" in eq. 5.3 is a misprint. The same derivation reproduces the
  kinematics (5.2) exactly. Results: `derivations/rigid_hawser_results.md`; check: `tests/test_rigid.py`.

### D18 — How the rigid-hawser equations are solved (2026-10-17, DONE — small, reversible)
- The hawser force F_H is an unknown, like the two angular accelerations. At every time step we solve the three
  equations of motion (5.1) as a 3×3 linear system for (γ̈, ψ̈, F_H) (`fishtail/model_rigid.py`). Simple, and it
  gives the hawser force for free. Equivalent to Simos's approach of rewriting (5.1) in γ and ψ.
- Bifurcation work is done in t* = U t / L. The code evaluates the (dimensional) forces at U = 1 m/s with C_F given
  explicitly; `tests/test_rigid.py` checks that the result in t* is the same for any U (D1).

### D19 — Stability boundaries from the full linearisation, not eq. B.13 as printed (2026-10-17, DONE)
- Simos's eq. B.13 writes the hull's lift force with the heading ψ only. The full model's lift depends on the
  angle of the *relative* velocity, which also contains the sideways velocity L_H γ̇ − x_A ψ̇ of the hull; that adds
  damping terms (−Y_v L_H, +Y_v x_A, −N_v L_H, +N_v x_A) missing from his B.15/B.18.
- Evidence: with the full linearisation our Hopf curve matches his Fig. 5.3 to +0.5 % (100 %) and +5.6 % (40 %);
  with B.13 as printed the curve is 3–6× too high (E06 fig01). So his figure was computed with the full terms.
- Our linearisation is derived symbolically and matches finite differences of the nonlinear model (tests).

### D20 — Only equilibria with the hawser in tension are followed (2026-10-17, DONE — small)
- A rigid bar could also hold the ship upstream of the buoy, in compression. That is not what a moored ship does, so
  these equilibria are not computed or plotted (E03). Compression *during* motion is still flagged (D10).

### D21 — Figures are saved as PNG only (2026-10-17, AGREED — student's request)
- `fishtail/plotting.save` no longer writes PDF copies; figure conventions updated. The 8 PDFs already made for E01 are
  still in `figures/E01_hydrodynamics/` (not deleted without asking).

### D22 — Demo data as `data.js`, one page for all demos (2026-10-17, DONE — small, reversible)
- A web page opened from disk (`file://`) is not allowed to read a `.json` file, but it may load a `.js` file. So
  `interactive/build_demos.py` writes `interactive/<demo>/data.js` = one line of JavaScript wrapping the JSON.
- One page `interactive/index.html` with three tabs (pitchfork, fishtailing, hawser length) instead of one page per
  demo: same layout, less code. D4 (control) will be added as a fourth tab after E08.

### D24 — Demo page layout (2026-10-17, AGREED — student's request)
- Light mode only. The first screen (no scrolling, for presenting) shows only the slider, the ship from above and
  the bifurcation diagram. Phase plane (limit cycle), eigenvalues and time trace are below; all explanation text is
  at the bottom. Plain technical look: white, black rules, square corners, one orange accent meaning "the ship now".

### D25 — Stretchy-rope law from a real SPM hawser (2026-10-17, AGREED: student asked for n > 1 from a real source)
- Rope law tension = a · strain^n, fitted (log-log least squares) to the measured "prototype" curve of a 144 mm nylon
  SPM hawser, 70 m, breaking load 650 t (Huang et al. 2026, R11, Fig. 3, read by eye): **a = 10 580 kN, n = 1.86**.
- Searched first for the power-law values used by Bernitsas and co-workers (1980s–90s): the papers are paywalled and
  no value could be confirmed, so none is quoted. Fibre-rope research models (Michel & Davies 2023) use hysteresis
  models that are too detailed for this project.
- Same strain law applied to our 160 m hawser; one line; no rope damping (so snatch peaks are on the high side).

### D23 — Parameter choices for experiments and demos (2026-10-17, AGREED 2026-10-18)
- Hull for the main runs: Vidal 100 %, attachment at the bow x_A/L = 0.5 (as in Simos's 100 % tests), L_H/L = 0.5.
- Friction values used as "typical": C_F = 0.015 ≈ real ship (E01 fig08: 0.009–0.016 at 1–2 m/s), C_F = 0.035 and
  0.05 ≈ 1:90 tank model (0.02–0.05). For comparisons with Simos (E05) we use his C_Y = 0.86, l_p C_Y = 0.044 and his
  C_F range 0.0223–0.0247, as he did.
- Demos: pitchfork at C_F = 0.05 (clean pitchfork, nothing fishtails); fishtailing and hawser-length at x_A/L = 0.5;
  hawser-length at C_F = 0.03. Every run starts from γ = 3°, ψ = 1°.

### D26 — No stretchy-rope tab in the demo page (2026-10-18, AGREED — student's request)
- The stretchy rope barely moves the stability boundaries (E07), so a slider tab would look almost the same as the
  rigid tabs. Its real result (slack/snatch loads) stays in E07 and the report.

### D27 — Water properties from ITTC (2011) (2026-10-18, AGREED — student's request)
- Sea water ρ = 1025 kg/m³ (student's choice; ITTC Table 3 gives 1026.02 at 15 °C), ν = 1.19e-6 m²/s
  (ITTC Table 3, 15 °C: 1.1892e-6). Fresh water (tank) 1000 kg/m³, 1.14e-6 m²/s (ITTC Table 1, 15 °C: 999.10,
  1.1386e-6). Source R15. No number in the code changed, so no result changes. Settles most of open issue 10.

### D28 — Wind model: constant force, fixed direction, acting at midship (2026-10-18, PROPOSED — easy to change)
- Student asked for wind with two knobs, force and direction, and a 4th demo tab, to see the imperfect pitchfork.
- Model: F_wind of fixed size and fixed earth direction β_W (where it comes FROM; 0 = from ahead, +90° = from port),
  acting at x_W = midship. Size in coefficient form C_W = F_wind / (½ ρ U² L T); shown as F_wind / F_H0 = C_W / C_F.
- Why so simple: wind coefficients for this hull are not available in a source we can read (OCIMF data are
  paywalled); a constant force is the cleanest "imperfection" h of Strogatz §3.6. Heading-dependence ignored.
- x_W is a parameter (default 0). E10 fig07 shows what changes if it moves: aft → same picture; at/ahead of the
  hawser point → favoured side flips. **To confirm with the student:** keep midship, or use aft (superstructure)?
- With C_W = 0 the model is exactly the old one (checked in tests/test_wind.py).
- Demo D4 (Wind tab) is added as the 4th tab; control becomes D5.

### D29 — Wind demo keeps the ship's history (2026-10-18, DONE — student's feedback)
- First version restarted every run from the same start, so hysteresis could not be seen, and its diagram had x_A/L
  on the axis, which no slider moves. Now: each slider step = "wind changes, then holds"; the run starts from the rest
  position the ship was in (side carried over). Python precomputes one run per (setting, side, neighbouring setting);
  big slider jumps are walked one step at a time, only the last step is animated. Main diagram: rest heading vs wind
  force (the slider moves along it). Still "Python computes, page plays back" (D8).
- Second round of feedback: "Put the ship on" buttons (or clicking a branch in the diagram) place the ship on any
  stable rest position of the current wind, or on the unstable middle one with a 0.5° nudge along its growing
  eigenvector (precomputed, 190 runs), so every branch can be reached and then followed with the sliders.
- Playback now runs at 10 t* per second at 1× for every tab (runs have different frame spacing). The water streaks
  move right to left at a fixed on-screen speed: tied to t* they aliased and seemed to run backwards.

### D30 — Demo page: pitchfork start choice, clearer eigenvalue panel (2026-10-18, DONE — student's feedback)
- **Pitchfork tab** gets the "Start the ship" bar: straight with the 3°/1° nudge to port (as before) or to starboard
  (new Python runs `runs_stbd`, the mirrored start), or sitting at rest on the port / starboard angled branch (only aft
  of x_Tcr; rest values and hawser force from Python). Clicking a branch in the diagram does the same.
- **Eigenvalue panel:** axes fitted per tab to every slow eigenvalue (the old fixed range −0.6…0.3 cut off the
  pitchfork's growing eigenvalue, up to +0.67); the fast real eigenvalue (< −1.7) is listed, not plotted, instead of
  the overlapping "◀ −2.41" labels; filled = stable / open = unstable rest position; ▲ / ● mark where the eigenvalues
  cross at the pitchfork / Hopf point (eigenvalues at the critical value, from Python: `crossing` in data.js); a table
  and a sentence under the plot say which eigenvalue makes the position unstable and what that means. Wind tab: grey
  trail = the other forces of the current direction only (was all 143 settings: blurry clouds).
- Bifurcation-diagram axes and colours are now explained per tab in "How to read it" (blue = edges of the limit cycle).
- Second round (student): panels less wordy. Legends cut to a few words, eigenvalue panel ends with one line
  ("Straight unstable: a real λ crossed into the RHP → pitchfork"), all longer text moved to "How to read it". Mirror
  port/starboard positions share one table row. The three detail graphs are now the same size.
- **Wind force slider 3× finer:** 31 values, step 0.0133 F_H0 (was 11, step 0.04). To stay under the 16 MB publish
  limit, each transition run stops being stored once its rounded values no longer change (`trim_settled`); the page
  repeats the last frame, so nothing is lost at the stored precision.
- **Playback speed:** student asked for 1× = real time with 2×/4×/8×. Real time per unit t* is L/U: 148 s (Pitchfork,
  Wind: 1:90 model at 0.024 m/s), 21 s (Hawser length), up to 7 h for the real ship, so real-time playback would be
  unwatchable. Student chose "keep t*, label honestly". 1× = 1 t* per second was tried and was too choppy (runs store
  only 2-3 frames per t*), so the speeds went back to the old ones: 1× = 10 t* per second, ½×-4×. Kept: the line under
  the speed menu giving 1 t* in real time for the model and the ship (speeds from Python, `common()`).
- **Slow swings across in the wind demo:** with the finer force grid, 8 transition runs (next to folds, ship swinging
  right across) had not settled after 150 t* (80-100° away from any rest position), so their end side could be named
  wrongly. Such runs now run on to 600 t* (checked: all settle by then). The page sizes the time axis per run.
- **Mini top view** (student chose it from three options): a small floating top view, bottom right, appears when the
  big top view is scrolled out of sight, so the detail graphs can be related to the ship; click it to scroll back.
  It is **off by default**: the "Show ship view" button in the Detail heading turns it on (the same button or ×
  turns it off). It can be dragged by its bar and resized from its bottom-right corner; position and size are
  remembered in the viewer's browser only. It shows only while the Detail graphs are on screen and hides once more
  than 60 % of them have been scrolled past.

### D31 — Control: PD stern thruster with saturation, plus a hold-back tug (2026-11-02, AGREED; numbers PLACEHOLDER)
- **Agreed with the student (2026-11-02):** option C = both a hold-back tug (tab "Tug") and the tanker's own stern
  thruster with PD heading control (tab "Thruster"); no bow thruster (with the hawser at the bow it has almost no
  lever arm; E08 fig08 compares). Gains non-dimensional, with forces also given in N (1:90 model) and kN (real ship)
  where a current speed gives the friction. "Optimal" = **least thrust that still meets a margin** σ ≥ σ_min
  (σ = −max Re λ; σ_min = 0.05 U/L PROPOSED). Thruster tab: sliders k_d and k_p (a map of the two gains instead of a
  1-D diagram) plus a 3-button thrust-limit switch (none / 0.02 / 0.01 F_H0, PLACEHOLDER), so the nonlinear effect
  of saturation (bistability) is visible.
- **Why the tug pulls backwards** (student's question): once the bow swings off, the downstream pull at the stern gets
  a sideways part ≈ F_tug·ψ, a restoring moment with a full ship length of lever arm (a weathervane tail, a "free"
  P controller), and it adds to the hawser tension. Sources (secondary, OCIMF not read): the tug keeps the hawser
  taut and helps the tanker weathervane (Marine Insight; Offshore Magazine; Rigzone). No source found that names it a
  fishtailing cure: in the report that is a model result, not an assumption.
- **Stern thruster.** The ship's own side thruster at the stern, x_T = −0.5 L (PLACEHOLDER position), driven by a PD
  law on the heading, written non-dimensionally so the gains do not depend on current speed in the rigid model:
  N_c* = −k_p ψ − k_d r*, with N_c* = N_c / (½ρU²L²T) and r* = rL/U. This is the brief's N_c = −K_p ψ − K_d r (§5).
- **It is a side force, not a pure moment.** The thruster pushes sideways with Y_T = N_c / x_T, so it also adds a sway
  force. Both enter the equations of motion (`model_rigid.py`), the straight-trailing linearisation
  (`linear_matrices`: k_p adds to the ψ column of the stiffness matrix C, k_d to that of the damping matrix B) and the
  equilibrium finder (`continuation.py`: thruster force at rest, r = 0).
- **Saturation:** |Y_T| ≤ thrust_limit × F_H0 (fraction of the straight-trailing hawser pull); beyond that the force
  is clipped. Default limit = ∞. At trailing straight the thrust is zero, so the limit cannot move the eigenvalues or
  the Hopf gain there; it matters only for large swings (time runs, E08).
- **Hold-back tug** (second control option): a steady pull downstream from the stern. No new code: it is the wind
  force with β_W = 0 (pulls along the current) and x_W = stern (`X_TUG_L = −0.5`, PLACEHOLDER). It keeps the mirror
  symmetry, so straight trailing is still an equilibrium; `linear_matrices` now accepts any force along the current
  line (sin β_W = 0), not only zero wind.
- **Zero gains = the model of steps 3–6 exactly** (brief's check for step 8). Checked in `tests/test_control.py`
  (9 checks: zero gains, sign of the law, the limit is never exceeded in a simulation, linear Jacobian vs finite
  differences with thruster and with tug, mirror symmetry, straight is still an equilibrium, limit-independence of the
  eigenvalues, the test case fishtails without control).
- **Still open:** cited values for the thrust limit, thruster position and tug pull (open issue 25); σ_min to be
  confirmed by the student.
- **First results (E08):** tug Hopf at 0.301 F_H0 (base case C_F = 0.03); at real-ship friction (C_F ≈ 0.01) about
  2.2 F_H0. Thruster: Hopf at k_d = 0.0064 (k_p = 0) or k_p = 0.0131 (k_d = 0). One eigenvalue (−0.2205 U/L) does not
  move with the gains (a sideways drift a heading controller cannot reach), so σ ≤ 0.22. Too much k_d overdamps (the
  pair splits into two real λ and one moves back right). ★ = k_p 0.004, k_d 0.012 (thrust needed 0.082 F_H0).
  Saturation gives bistability: with the ★ gains, for limits 0.0063–0.0113 F_H0 the full fishtail survives while a 3°
  nudge dies out. Full table: experiments/E08_control/README.md.
- **Demo D5 (two tabs, 2026-11-02).** Every run starts from the full uncontrolled fishtail with the control switched
  on at t* = 0 (200 t* shown). Tug tab: 48 pulls (0–1 in steps of 0.025, then to 2.4) × 13 C_F (0.010–0.058).
  Thruster tab: 31 k_d (0–0.06) × 16 k_p (0–0.03) × 3 limits (none / 0.02 / 0.01 F_H0, one data file each, < 16 MB).
  Run label: a run still swinging at 200 t* while the straight position is stable is run on for 300 t* more (not
  stored) and called "fishtails" only if its swing no longer shrinks (last 100 t* ≥ 95 % of the 100 t* before),
  otherwise "swing dying out"; where the straight position is unstable no check is needed (student's suggestion,
  shortened from 800 t*). Builds run in parallel (multiprocessing, 10 cores, student's OK).
