# Project brief: fishtailing of a single-point-moored tanker

The original project plan. Later decisions are recorded in `notes/00_decision_log.md`, which takes precedence.

---

## 1. Who and why

- **Student project** for the course *Non Linear Problems in Ocean Engineering*. Text: Strogatz, *Nonlinear Dynamics and Chaos*. Syllabus: bifurcations, fixed points, limit cycles, nonlinear RAOs.
- **Worth 10% of the course grade. Time available: about 2 weeks.** Code is not the bottleneck; understanding and clear results are.
- Wanted: a topic on **manoeuvring, control and hydrodynamics** that is not too complicated but matches the standard of the topics other students chose.
- Topics already chosen by classmates: forced Van der Pol, forced Duffing, forced Duffing–Van der Pol, lid-driven cavity, bifurcations in flow past a cylinder, Lagrangian coherent structures in surface flow.
- **Chosen topic:** bifurcations and limit cycles ("fishtailing") of a tanker moored to a buoy by one hawser in a steady current, plus an optional heading-control stage.

### Why this topic
- Fills a gap: none of the six existing topics covers manoeuvring, mooring or control.
- Covers much of the syllabus in one physical system: pitchfork, imperfect pitchfork, supercritical and subcritical Hopf, limit cycle, cycle fold and hysteresis.
- Autonomous system (steady current), so the limit cycle is genuinely self-excited, not a forced response.
- Easy to run: ODEs only, no mesh or CFD. The work is in the modelling and bifurcation analysis.
- Published results exist to validate against.

### Honesty note (important for the report)
This topic is **not new research**. It has been studied since the 1980s (see references). The report must **not claim novelty**. Its value is a clean from-scratch reproduction of the known bifurcation scenario, presented through Strogatz's framework, with a control extension. Classmates' topics are equally textbook.

---

## 2. The physical problem, in simple terms

A tanker is tied by one rope (hawser) from its bow to a buoy. A steady current flows past. Three things can happen as conditions change:

1. **Trails straight:** the ship lines up with the current like a weathervane (stable).
2. **Sits at an angle:** the ship settles steadily to one side of the current.
3. **Fishtails:** the ship swings side to side in a sustained cycle, with large snatch loads in the rope.

Why it goes unstable: if the ship yaws slightly, the current hits it at an angle and creates a turning moment (the **Munk moment**) that turns it further broadside. The rope tension at the bow pulls it back. Which wins depends on current speed, rope length and where the rope is attached.

A ship moored in a current is the same problem as a ship towed through still water, seen from the tug.

---

## 3. What we expect to find (bifurcation story)

Slowly vary a parameter and watch the "trailing straight" equilibrium (heading aligned with the current).

| Event | Eigenvalue picture | What the ship does | Strogatz |
|---|---|---|---|
| Pitchfork | One real eigenvalue crosses 0 | Straight position loses stability; two mirror-image angled positions appear, **both stable** (supercritical case), left and right | §3.4, §8.1 |
| Imperfect pitchfork | Symmetry broken (e.g. wind at an angle) | One side preferred | §3.6 |
| Hopf | Complex pair α ± iω crosses the imaginary axis | Damped swing stops decaying, a self-sustained swing (limit cycle) is born, period about 2π/ω | §8.2 |
| Subcritical Hopf + cycle fold | | Sudden jump to large swings; hysteresis when the parameter is reversed | §8.2, §8.4 |

Corrections agreed during planning:
- In the pitchfork, **both branches are stable**, not one. The ship picks one depending on tiny disturbances.
- The limit cycle is **not a separate event** from the Hopf: the limit cycle is what the Hopf creates.
- The pitchfork and the Hopf need not appear in the same sweep. Which comes first depends on rope length, fairlead position and hull coefficients. This is why the main result is a **two-parameter stability map**.
- Published work found a **subcritical** Hopf with hysteresis in a flume experiment. Our model may instead give a supercritical Hopf. Report whichever the model gives; do not force a result.

---

## 4. The model

### 4.1 Full model (6 states)

States: surge speed `u`, sway speed `v`, yaw rate `r` (body-fixed axes), position `x, y`, heading `ψ`.

```
(m - X_udot) du/dt - (m - Y_vdot) v r  = X_H + T_x
(m - Y_vdot) dv/dt + (m - X_udot) u r  = Y_H + T_y
(I_z - N_rdot) dr/dt                   = N_H + x_F * T_y
dx/dt = u cos ψ - v sin ψ
dy/dt = u sin ψ + v cos ψ
dψ/dt = r
```
(Sign conventions for added mass must be checked and made consistent; verify against the cited papers.)

- `X_H, Y_H, N_H`: hydrodynamic forces and yaw moment, functions of the velocity **relative to the current** (`u_r, v_r`) and of `r`.
- `T_x, T_y`: hawser tension resolved along ship axes. `x_F`: distance from centre of gravity forward to the fairlead (bow attachment point), so the rope also creates the moment `x_F * T_y`.
- **Hawser law** (nonlinear spring that goes slack):
  `T = k (ℓ - ℓ0)^n` if `ℓ > ℓ0`, else `T = 0`. Here `ℓ` is the current fairlead-to-buoy distance, `ℓ0` the unstretched length, `n` the stiffening exponent (nylon).

### 4.2 Reduced model (rigid hawser, 4 states) — do the bifurcation analysis on this first

Treat the hawser as a rigid, always-taut link of fixed length `ℓ`. The ship then has two degrees of freedom, e.g. the hawser angle `γ` and the heading `ψ`, giving 4 states `(γ, γ̇, ψ, ψ̇)`.

- Derive it from the 6-state model with the length constraint (Lagrangian approach or constrained Newton–Euler). **Use sympy** to derive and to produce the Jacobian, and keep the derivation in the repo.
- **Validation:** the equilibria of the 4-state model must match the equilibria of the 6-state model in the limit of a very stiff rope (large `k`).
- Simos et al. (2001) used the same rigid-hawser simplification to avoid slack-line complications.

### 4.3 Hydrodynamic model (the main modelling decision)

Recommended: **Option B**
- Current-load coefficients `C_X(β), C_Y(β), C_N(β)` valid at any flow angle `β` (OCIMF-type), e.g. `Y = ½ ρ L T V_r² C_Y(β)`.
- Plus a **cross-flow drag integral** along the hull for yaw damping.
- The destabilising Munk moment must be present (a `sin 2β`-type term).

Option A (manoeuvring derivatives, Abkowitz or MMG style) is accurate near the aligned position but unreliable at large swing angles. May be used for a cross-check near equilibrium.

**Hull:** a benchmark tanker, preferably **KVLCC2** (published MMG coefficients, Yasukawa & Yoshimura 2015). KCS (a container ship) was also suggested by the student; it works as a fallback but is not a tanker. Hull dimensions and coefficients **must be taken from the published sources and checked**; do not invent or recall numbers from memory. KVLCC2 full scale is roughly L = 320 m, B = 58 m, T = 20.8 m (verify). If coefficients cannot be obtained, use clearly labelled placeholder non-dimensional parameters and say so in the report.

**Sanity check that must pass:** the bare hull (no rope) aligned with the current must be **unstable in yaw** (the Munk moment). If it is not, the hydrodynamics are wrong.

---

## 5. Parameters

### Knobs to vary (bifurcation parameters)
| Parameter | Meaning |
|---|---|
| Current speed `U_c` | Main knob. Faster current, stronger flow forces. |
| Hawser length `ℓ` | Changes how much the rope can straighten the ship. |
| Fairlead position `x_F` | How far forward of the centre of gravity the rope is tied. |
| Controller gains `K_p`, `K_d`; thruster limit | Control stage only: yaw thruster `N_c = -K_p ψ - K_d r`, with saturation. |
| Wind angle (optional) | Breaks symmetry, giving the imperfect pitchfork. |
| Rope stiffness `k`, exponent `n` (optional) | Stiff vs stretchy. |
| Water depth (optional) | Changes hydrodynamic forces. |

### Fixed ship properties (from the benchmark hull)
`m, I_z`, added masses `X_udot, Y_vdot, N_rdot`, `L, B, T`, current coefficients, `ρ`.

**Non-dimensionalise** (lengths by `L`, speeds by a reference speed, etc.). The axes of the results are then ratios such as `ℓ/L`, `x_F/L` and a current-speed parameter.

---

## 6. Work to do (in order)

Each step has a clear check. Do not move on if the check fails.

1. **Setup.** Repo structure, parameter file (one place for all numbers, with source and units for each), nondimensionalisation.
   *Check:* units consistent; every parameter has a cited source or a "placeholder" label.
2. **Hydrodynamics.** Implement current forces and moment, cross-flow drag.
   *Check:* bare-hull yaw instability (section 4.3); forces vanish at zero relative speed.
3. **Rigid-hawser model (4 states).** Derive with sympy, implement the right-hand side and an exact or finite-difference Jacobian.
   *Check:* Jacobian matches finite differences; equilibria satisfy the equations.
4. **Equilibria and stability.** Find the trailing equilibrium and any angled equilibria; compute eigenvalues along each parameter sweep. Natural-parameter continuation is enough; use pseudo-arclength near folds if needed.
   *Check:* the eigenvalue crossings found match time-domain behaviour either side.
5. **Locate bifurcations.** Detect pitchfork (real eigenvalue crossing 0) and Hopf (complex pair crossing the imaginary axis) points; trace them in two parameters to build the **stability map**.
   *Check:* map regions agree with a handful of time-domain runs placed in each region.
6. **Limit cycles.** Simulate just past the Hopf point; measure swing amplitude versus parameter, sweeping up and then down. Compute the cycle's stability (Floquet multipliers via shooting). Classify the Hopf as super- or subcritical; check for hysteresis.
   *Check:* near a supercritical Hopf, amplitude ∝ √(distance from threshold).
7. **Elastic hawser model (6 states).** Time-domain only (slack makes the equations non-smooth). Compare boundaries with the rigid model; record hawser tension, snap loads.
   *Check:* boundaries agree with the rigid model for a stiff rope.
8. **Control.** Add the PD yaw thruster with saturation; recompute the stability map with control.
   *Check:* with zero gains the original map is recovered.
9. **Optional extension.** Slow drift forcing near the fishtailing frequency (nonlinear response curve / RAO, frequency locking). **Cut this first if time is short.**
10. **Report and presentation.**

**If time is short:** drop steps 7 and 8 and keep the rest. Steps 1–6 already give the pitchfork, Hopf, limit cycle and stability map, which are the core of the course material.

---

## 7. Results to produce (save as separate scripts, one per figure)

1. Schematic of the mooring (plan view).
2. **Pitchfork diagram:** steady heading versus current speed, stable branches solid, unstable dashed.
3. **Eigenvalue plot:** eigenvalues moving with current speed, showing the real crossing (pitchfork) and the complex-pair crossing (Hopf).
4. **Swing amplitude versus current speed**, up and down sweeps (hysteresis loop if subcritical).
5. **Phase portraits** before and after the Hopf, plus the ship's path in plan view behind the buoy.
6. **Two-parameter stability map** (for example hawser length versus current speed): regions "trails straight", "sits at an angle", "fishtails". Main result and a design chart.
7. **Hawser tension time series** during fishtailing.
8. **Stability map with and without control.**
9. (Optional) nonlinear response curve with slow forcing.

Figures should be clear, labelled with units or non-dimensional definitions, and readable in a report.

---

## 8. Technical guidance

- **Python**: numpy, scipy (`solve_ivp`, `fsolve`, `eig`), matplotlib, sympy. Cross-check against MatCont or AUTO-07p only if available.
- Use tight tolerances in `solve_ivp` (for example `rtol=1e-9`, `atol=1e-12`) near bifurcations; transients near a Hopf point are long (critical slowing down), so simulate long enough.
- Keep the code readable and commented; the student must be able to explain every line. Prefer simple methods over clever ones.
- Fixed random seeds not needed (system is deterministic); every figure script must run on its own and write to `figures/`.
- Suggested layout:
  ```
  params.py        all parameters with sources
  hydro.py         current forces, cross-flow drag
  model_rigid.py   4-state model + Jacobian
  model_elastic.py 6-state model with slack hawser
  continuation.py  equilibria, bifurcation detection
  cycles.py        shooting, Floquet multipliers
  control.py       PD thruster
  fig_*.py         one script per figure
  figures/  notes/
  ```

---

## 9. Known risks

- **Hydrodynamic coefficients at large angles are uncertain.** Present results as **qualitative and parametric**, not predictions for a specific real ship. Even careful validation work needed an ad hoc correction to match static bifurcation experiments.
- **Slack hawser is non-smooth.** Do formal bifurcation analysis on the rigid (always-taut) model; use the slack model only in simulation.
- **Limit-cycle continuation is fiddly.** Use up/down sweeps and shooting; do not over-engineer.
- **Do not fabricate** parameter values, results, or citations. If something cannot be sourced, say so plainly.

---

## 10. How to work with the student

- The student is new to this topic and wants **simple explanations**; explain what each step means physically before writing code.
- Ask before major modelling decisions (hull choice, hydrodynamic model, scope cuts).
- After each step, show the result and what it means, and flag anything that looks physically wrong.
- Keep a running `notes/` log of decisions, parameter sources and open issues, so the report can be written from it.

---

## 11. References

- Strogatz, S. H. *Nonlinear Dynamics and Chaos*, 2nd ed. (§3.4, §3.6, §8.1, §8.2, §8.4, §8.7).
- Bernitsas, M. M. & Papoulias, F. A. (1986). Stability of single point mooring systems. *Applied Ocean Research* 8.
- Experimental study of the large-amplitude fish-tailing instabilities of a tanker at a single point mooring (1990). *Applied Ocean Research* 12(1). (Observed a subcritical Hopf with hysteresis and a cycle fold.)
- Papoulias, F. A. & Bernitsas, M. M. (1986). Stability of motions of moored ocean vehicles. *Dynamical Systems* 1(4).
- Garza-Rios, L. O. & Bernitsas, M. M. (1996). Analytical expressions of the stability and bifurcation boundaries for general spread mooring systems. *J. Ship Research* 40(4).
- Bernitsas, M. M. & Garza-Rios, L. O. (1997). Mooring system design based on analytical expressions of catastrophes of slow-motion dynamics. *J. Offshore Mech. Arct. Eng.* 119(2).
- Simos, A. N. et al. (2001). Experimental validation of a quasi-explicit hydrodynamic model: fishtailing instability of a single-point moored tanker in rigid-hawser configuration.
- Paton, C. G., Carra, C. & Sincock, P. (2006). Investigation of sway-yaw motions of deepwater FPSOs. OTC-18039.
- Papoulias, F. A. (1991). Analysis of straight-line motion of towed ships. *J. Ship Research* 35(4). (Related towing problem.)
- Yasukawa, H. & Yoshimura, Y. (2015). Introduction of MMG standard method for ship maneuvering predictions. *J. Mar. Sci. Technol.* (KVLCC2 coefficients.)

Backup topics considered if this one fails: towed-barge fishtailing (same maths), porpoising of planing boats, submarine dive-plane stability, line-of-sight path-following autopilot, nonlinear course stability with rudder saturation.
