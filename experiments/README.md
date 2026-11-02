# Experiments

Each experiment answers **one question** by changing **one knob** (sometimes two).
Each lives in its own folder `EXX_name/` with:
- `run.py` — runs on its own from the project root: `$PY -m experiments.EXX_name.run`
- `README.md` — the question, what is varied, what we expected, what we got, and what it means physically.

Figures go to `figures/EXX_name/`, saved data to `data/EXX_name/`.
Knobs and conventions: see the top-level `README.md`. Status: ☐ not started, ◐ in progress, ☑ done.

| ID | Question | Knob(s) | Expected (from theory / Simos 2001) | Figures (planned) | Demo | Status |
|---|---|---|---|---|---|---|
| E01 | Do our water forces match the published model and data? Is the bare hull unstable in yaw? | angle of attack α (0–180°), yaw rate r | Curves match Leite/Simos; Munk moment wins near α = 180° | 8 figures, see `E01_hydrodynamics/README.md` (coefficients vs α; moment parts; pinned-hull phase portrait; x_Tcr; C_D(x); yaw-rate terms; Simos Fig. 4.9 check; C_F vs U) | – | ☑ |
| E02 | Does current speed change the behaviour? | U (two values), with `C_F` fixed vs `C_F(Re)` | Fixed `C_F`: same motion, faster. With `C_F(Re)`: faster current → less stable | 3 figures: ψ(t) in seconds and in t*; model vs full-scale friction; growth rate vs speed. **Got: as expected** | – | ☑ |
| E03 | When does the ship stop trailing straight and sit at an angle? (pitchfork) | `x_A/L` | Straight position loses stability below x_Tcr/L (eq. 3.37) = **0.36 for 100 %, 0.19 for 40 %**; two mirror-image angled positions appear | 7 figures. **Got: 0.359 / 0.193, supercritical, angle grows linearly (ψ\|ψ\| normal form)** | D1 | ☑ |
| E04 | When does the ship start to fishtail? (Hopf) | `C_F` (↔ current speed), `L_H/L` | A complex eigenvalue pair crosses zero; damped swing → sustained swing | 6 figures. **Got: Hopf at C_F = 0.0425 (L_H/L = 0.5); period ≈ 8.5 L/U; hawser force can go negative at low friction** | D2, D3 | ☑ |
| E05 | How big is the fishtailing swing, and is the Hopf super- or subcritical? | same as E04, swept up then down | Supercritical: amplitude ∝ √(distance past threshold), no hysteresis. Subcritical: jump + hysteresis | 6 figures. **Got: supercritical, no hysteresis, amplitude ∝ distance (not √); periods match Simos Fig. 5.7** | D2 | ☑ |
| E06 | Where in the design space is the ship safe? (main result) | (`L_H/L`, `C_F`) and (`x_A/L`, `L_H/L`) | Reproduce Simos Fig. 5.3; regions "trails straight", "sits at an angle", "fishtails" | 5 figures. **Got: Fig. 5.3 reproduced (+0.5 % / +5.6 %); 12/12 time runs agree with the map** | D3 | ☑ |
| E07 | What changes with a real, stretchy rope that can go slack? | rope stiffness ratio (↔ current speed), `L_H/L` | Stiff rope → same boundaries as rigid; snap loads during fishtailing | 6 figures. **Got: boundaries within 5 % of rigid; real ship at 1 m/s: rope slack 74 % of the time, snatch peak 501 kN = 14 × steady (rigid: 172 kN)** | – | ☑ |
| E08 | Can a hold-back tug or a stern thruster stop the fishtailing? | tug pull F_tug/F_H0; gains `k_p, k_d`, thrust limit | Zero control → original map; enough control → stable region grows; saturation limits it | 9 figures. **Got: tug Hopf at 0.30 F_H0 (C_F 0.03), ≈2.2 F_H0 at real-ship friction; thruster Hopf at k_d = 0.0064; ★ least-thrust gains k_p 0.004, k_d 0.012 (0.082 F_H0); one eigenvalue (−0.22) cannot be moved; bistable for limits 0.0063–0.0113 F_H0; bow thruster much weaker (D31)** | D5 | ☑ |
| E09 | (Optional, cut first) Slow forcing near the fishtailing frequency | forcing frequency/amplitude | Frequency locking (close to forced Van der Pol — a classmate's topic) | response curve | – | ☐ |
| E10 | What does a steady wind do to the pitchfork? (imperfect pitchfork) | wind force F_w/F_H0, wind direction β_W, with `x_A/L` | Strogatz §3.6: symmetry broken → one smooth branch + one ending at a fold; hysteresis; cusp | 8 figures. **Got: folds at ±0.153 F_H0 at x_A/L = 0.25, sweeps jump there; cusp width ∝ distance^2.12 (x\|x\| form predicts 2); ship ends on the UP-wind side** | D4 | ☑ |

## Interactive demos
All in one page, `interactive/index.html` (open it in a browser). Data: `python -m interactive.build_demos`.
D1–D3 are built (2026-10-17), D4 (2026-10-18), D5 (2026-11-02, two tabs: Tug and Thruster). `index.html#wind` opens a tab directly.
- **D1 Pitchfork:** slider `x_A/L` → marker moves on the pitchfork diagram; ship settles straight or to one side.
- **D2 Fishtailing:** slider for current speed (through `C_F`) → eigenvalues cross; ship's swing grows to a limit cycle.
- **D3 Hawser length:** slider `L_H/L` → same as D2 along the other axis of the stability map.
- **D4 Wind:** two sliders, wind force and direction; x_A/L = 0.25 fixed. The ship keeps its history: it starts
  resting to port with no wind, and each slider step changes the wind and holds it, starting from where the ship was
  resting (precomputed transition runs), so hysteresis is visible. "Put the ship on" (or a click on the diagram)
  places it on any branch that exists for the current wind, including the unstable middle one with a nudge. Main diagram: rest heading vs wind force for the
  chosen direction (the force slider moves along it; ◇ = fold, also marked on the slider). Second view: rest heading vs
  x_A/L for the current wind, over the grey no-wind pitchfork.
- **D5 Control (two tabs):** *Tug*: sliders tug pull and C_F; swing amplitude vs pull with the Hopf point.
  *Thruster*: sliders k_d and k_p, a map of the decay rate over the two gains (Hopf line, margin line, ★),
  and a 3-button thrust limit; red × = gains where the full fishtail survives at that limit (bistable).
  Every run starts from the full uncontrolled fishtail with the control switched on at t* = 0.
