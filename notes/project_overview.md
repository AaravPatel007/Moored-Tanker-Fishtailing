# What has been done — the short, simple version

Last updated 2026-11-02. A plain-words map of the whole project from the start. The detailed versions are in
`progress_summary.md`, `00_decision_log.md` (D1–D31) and each experiment's `README.md`.

## The big picture in one paragraph

A big oil tanker is tied by one rope (the hawser) to a buoy, with a current flowing past. Depending on where the rope
is tied, how long it is and how much friction the water gives, the ship either trails straight behind the buoy, sits
at a steady angle to one side, or swings from side to side forever ("fishtailing"). The project builds a computer
model of this (following Simos's 2001 thesis), finds exactly where each behaviour starts, and explains the switches
with the ideas of Strogatz's book: a **pitchfork** (straight → angled) and a **Hopf bifurcation** (straight →
swinging). Extras: a stretchy rope that can go slack, a steady wind, and (in progress) a thruster or tug to stop the
swinging. A web page with sliders plays back the results.

## Timeline

| Date | Step | What was done |
|---|---|---|
| 2026-10-17 | 1 setup | Folders, notes, sources. Ship chosen: Vidal de Negreiros, a 320 m tanker with every number published (D3). |
| 2026-10-17 | 2 water forces | `hydro.py`: the forces of the current on the hull. E01: matches Simos's thesis; a hull held at midship turns away from the current (unstable). |
| 2026-10-17 | 3 rigid rope | `model_rigid.py` + a sympy derivation that settled a sign misprint in the thesis (D17). E02: current speed matters only through friction. |
| 2026-10-17 | 4–5 rest positions, bifurcations | `continuation.py`. E03 pitchfork at x_A/L = 0.359 (loaded) / 0.193 (ballast). E04 Hopf. E06: Simos's stability chart (Fig. 5.3) reproduced within 0.5 % / 5.6 %. |
| 2026-10-17 | 6 the swing | `cycles.py`. E05: the swing grows smoothly from zero (supercritical), periods match Simos within a few seconds. |
| 2026-10-17 | demos | `interactive/index.html` with tabs Pitchfork, Fishtailing, Hawser length. |
| 2026-10-17 | 7 stretchy rope | `model_elastic.py`, rope law from a real SPM hawser (D25). E07: stability lines barely move, but the rope goes slack and snaps tight — peak force 14 × the steady pull. |
| 2026-10-18 | wind (your request) | Constant wind force (D28). E10: the pitchfork becomes "imperfect" — folds, hysteresis, a cusp. New Wind tab. Several rounds of demo-page polish (D29, D30). |
| 2026-11-02 | 8 control (part 1) | `control.py`: a stern thruster that pushes back against heading and turning (PD control), with a maximum force. Also a hold-back tug option. 9 new checks. Experiment E08 not yet run. |

## Files, and what each one does

**The model (`fishtail/` folder)**

| File | What it does, in simple words |
|---|---|
| `params.py` | Every number (ship size, masses, water properties, rope law…) with its source. Unsourced ones say PLACEHOLDER. |
| `hydro.py` | The water's push and turning force on the hull for any angle and turning speed. |
| `model_pinned.py` | The hull held by a pin at midship: the simplest check that the bare hull is unstable. |
| `model_rigid.py` | Ship + rigid rope + water (+ optional wind, tug, thruster). Gives how the ship moves next and how hard the rope pulls. |
| `model_elastic.py` | Same, with a stretchy rope that can go slack. |
| `continuation.py` | Finds where the ship can sit still and whether that position is stable; finds the exact knob value where stability is lost. |
| `cycles.py` | Measures the fishtailing swing (size, period) and checks the swing itself is stable. |
| `control.py` | The stern thruster: how hard it pushes for a given heading and turning speed, clipped at its maximum. |
| `plotting.py` | Shared drawing style (stable solid, unstable dashed), little ship drawings, saving PNGs. |

**Other folders**

| Folder | What is in it |
|---|---|
| `derivations/` | sympy maths that checks the equations from the geometry; results printed to `rigid_hawser_results.md`. |
| `experiments/` | E01–E07, E10: each has `run.py` (makes the figures) and `README.md` (question, answer, meaning). E08 is next. |
| `figures/`, `data/` | Output of the experiments (49 PNG figures so far). |
| `interactive/` | The slider web page (`index.html`), the data it plays back (`*/data.js`) and the script that makes the data. |
| `tests/` | 141 automatic checks; all pass. |
| `notes/` | Decision log, glossary, sources, open issues, explanations — the report will be written from these. |

## The most important decisions (and why)

1. **Ship: Vidal de Negreiros** rather than KVLCC2 — every number our water-force model needs is published for it,
   and Simos's chart gives something to check against (D3).
2. **Current speed is not a free knob** — it cancels out except through friction (D1).
3. **A misprinted sign and a missing lift force in the thesis** were found and handled; with the lift force kept we
   match Simos's chart almost exactly (D17, D19).
4. **The web page only plays back Python results**, so it can never disagree with the report (D8, D22).
5. **Rope law from a real hawser** (Huang et al. 2026), not a made-up spring (D25).
6. **Wind is a constant force in a fixed direction** — the simplest way to break the symmetry (D28).
7. **Control: thruster written in non-dimensional gains**, so a gain means the same at every current speed; with zero
   gain the model is exactly the old one (D31).

## What came out (headlines)

- Tie the rope aft of 0.36 L (loaded) or 0.19 L (ballast) from midship and the ship sits at an angle (pitchfork).
- Low friction (= fast current, or a full-size ship) makes the straight position start swinging (Hopf). The swing
  grows smoothly, with no jumps (supercritical).
- Surprise for the report: the angle and the swing grow **in proportion** to the distance past the bifurcation, not
  as a square root as in Strogatz, because water drag grows like speed × |speed|.
- A full-size ship has far less friction than a tank model, so it sits much deeper in the swinging region.
- A real rope goes slack during the swing and then snaps tight: rope peaks far above the rigid-rope estimate.
- A steady wind turns the pitchfork into an imperfect one with hysteresis; the ship ends up on the up-wind side.

## Not done yet / needs your decision

- Step 8: run E08 and add the control tab D5. Check the D31 choices and find real thruster / tug numbers (issues 25–26).
- Step 10: report and slides. Step 9 (forcing) is cut unless time is left over.
- Open questions: `open_issues.md`. The 8 old PDF figures of E01 can be deleted if you say so.
- Nothing has been committed to git.
