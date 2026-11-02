# When Does a Moored Tanker Start to Fishtail? Pitchfork and Hopf Bifurcations of a Single-Point-Moored Tanker in Steady Current

## ▶ [Open the interactive Fishtailing Lab](https://aaravpatel007.github.io/Moored-Tanker-Fishtailing/interactive/)

**https://aaravpatel007.github.io/Moored-Tanker-Fishtailing/interactive/** — runs in any browser, no installation.
Move a slider and watch the moored ship trail straight, settle at an angle, or fishtail, next to its bifurcation
diagram, phase plane and eigenvalues.

---

Course project, *Non Linear Problems in Ocean Engineering* (OE5320). Author: Aarav Patel.

A tanker moored to a buoy by a single hawser in a steady current either trails straight, rests at an angle
(pitchfork bifurcation) or swings from side to side (Hopf bifurcation, limit cycle). This project reproduces the
published model of Simos (2001) for the Vidal de Negreiros VLCC and explains the results with the tools of
Strogatz, *Nonlinear Dynamics and Chaos*. No new results are claimed.

- Experiments and results: [experiments/README.md](experiments/README.md)
- Equations and method: [notes/03_mathematics_explained.md](notes/03_mathematics_explained.md)
- Modelling decisions: [notes/00_decision_log.md](notes/00_decision_log.md) · References: [notes/sources.md](notes/sources.md)

## Conventions

Sign conventions follow Simos (2001): `x` forward, `y` to port, `ψ > 0` = bow to port, `α = 180°` = current from
dead ahead. Angles in radians in code, degrees on figures; time on figures is `t* = U t / L`. Figures: stable
branches solid, unstable dashed; ▲ pitchfork, ● Hopf, ◇ fold. Figures are PNG (200 dpi) in `figures/<experiment>/`.

## Running locally

Python 3 with NumPy, SciPy, Matplotlib, SymPy and pytest; run from the repository root.

```bash
python -m experiments.E03_pitchfork.run     # one experiment; writes to figures/ and data/
python -m interactive.build_demos           # rebuild the Lab's data (several minutes)
python -m pytest                            # test suite
```
