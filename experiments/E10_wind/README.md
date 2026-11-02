# E10 — What does a steady wind do to the pitchfork? (imperfect pitchfork, Strogatz §3.6)

**Run:** `/home/aarav-patel/wec_rl/.venv/bin/python -m experiments.E10_wind.run` (about 1 minute)
**Code:** wind force in `fishtail/model_rigid.py` (`wind_forces`); rest positions with wind in
`fishtail/continuation.py` (`wind_equilibria`, `attachment_for_heading`, `wind_branches`)
**Checks:** `tests/test_wind.py` (17 checks)
**Output:** `figures/E10_wind/fig01…fig08` (PNG), `data/E10_wind/summary.json`. Demo: the **Wind** tab (D4).

## Question
Without wind the ship is mirror-symmetric: when the attachment point is aft of x_Tcr/L = 0.359 it can rest equally
well to port or to starboard (pitchfork, E03). Wind from one side breaks that symmetry. What happens to the pitchfork,
and does the ship's history (which side it was on) matter?

## The wind model (decision D28, proposed)
A steady force of fixed size and fixed direction (it does not turn with the ship), acting at midship.
- **Size:** given as a fraction of the hawser tension when the ship trails straight without wind, F_w / F_H0.
  In the code: C_W = F_w / (½ ρ U² L T), so F_w / F_H0 = C_W / C_F.
- **Direction:** β_W = where the wind comes FROM. 0° = from ahead (same way as the current), +90° = from port,
  −90° = from starboard, ±180° = from astern.
- No real wind coefficients (they depend on the superstructure and are not published for this hull in a source we
  can read). Results are therefore given for a range of wind force, not for a real wind speed.

## Why it can still be solved exactly
At rest, "yaw moment − x_A × sway force" removes the unknown hawser force (as in E03). With wind this becomes
C6C + (x_W/L) Y_w − (x_A/L)(C2C + Y_w) = 0, which is **linear in x_A/L** (and also linear in the wind force).
So every heading ψ has exactly one attachment point (and one wind force) where it is a rest position: the whole
bifurcation diagram, folds included, is drawn from a formula, without searching for roots.

## What we got (C_F = 0.05, L_H/L = 0.5, as the pitchfork demo)

| Check | Result | Figure |
|---|---|---|
| Side wind on the pitchfork | the pitchfork splits into one **smooth branch** (the favoured side) and one **separate branch that ends at a fold** (◇). Fold moves aft as the wind grows: x_A/L = 0.288 / 0.266 / 0.239 for 0.05 / 0.10 / 0.20 F_H0 from port | fig01 |
| Wind direction | wind along the current (0°, 180°) keeps the symmetry: still a perfect pitchfork, hardly moved. Any sideways part unfolds it; the beam wind (±90°) the most | fig02 |
| Hysteresis at x_A/L = 0.25 | S-shaped curve of rest heading against side wind; folds at **±0.153 F_H0**. Simulated slow sweeps jump at **±0.16** (the next step past the fold) and follow different paths up and down | fig03 |
| The cusp | in the plane (x_A/L, side wind) the two fold lines meet at the old pitchfork point, forming a cusp. Near the tip the fold force grows like (distance)^**2.12** — the x|x| normal form predicts 2, Strogatz's cubic one 3/2 | fig04 |
| Time runs (ship resting to port, wind from starboard switched on) | 0.05, 0.10, 0.14 F_H0: ship stays on the port side (heading 9.4°, 8.2°, 6.7°). 0.16 and 0.20 F_H0: swings across to starboard (−12.5°, −13.0°) | fig05 |
| Bow attachment (x_A/L = 0.5, no pitchfork) | one rest position, pushed aside: up to γ = 9° / 16° / 28° and ψ = 2° / 4° / 6° for 0.1 / 0.2 / 0.4 F_H0. Hawser tension changes by ± the wind force (head/tail wind) | fig06 |
| Where the wind acts | aft of the hawser point (midship, superstructure): the bow turns into the wind. At the bow (x_W = x_A): no turning, the favoured side flips | fig07 |
| Wind and fishtailing (bow attachment) | wind from ahead and from the side make fishtailing start later (Hopf C_F 0.0425 → 0.038 and → 0.029 at C_W = 0.005); from astern earlier (→ 0.047). A beam wind C_W = 0.006 gave no fishtailing at any C_F in 0.004–0.15 | fig08 |

## What it means physically
- **Imperfect pitchfork.** Strogatz's ẋ = h + r x − x³ (§3.6): h is the wind's side push, r the distance aft of
  the pitchfork point. With h ≠ 0 there is no longer a "straight" position and no pitchfork point: one side is
  always available, the other appears only past a fold (saddle-node bifurcation).
- **Hysteresis = the ship remembers its side.** Between the fold lines (inside the cusp) both sides are stable. Which
  one the ship sits on depends on its history. A wind from the "wrong" side must exceed the fold force before the
  ship swings over, and once over, it stays even if the wind drops back.
- **The cusp is "sharp" in the same way the pitchfork was "linear".** The leading nonlinearity is the cross-flow drag
  ψ|ψ| (E03 fig06), so the normal form is ẋ = h + r x − x|x|. Its folds sit at h = ±r²/4, slope 2 on a log-log plot,
  instead of Strogatz's slope 3/2. Measured: 2.12 (fig04).
- **The ship ends up on the upwind side (surprising, flagged).** Wind from port pushes midship, which is aft of the
  hawser point, to starboard. That turns the bow to port, i.e. into the wind (a weathervane). With the bow to port
  the hull works like a wing and the current's lift pushes the whole ship to port, much more strongly than the wind
  pushes it the other way. If the wind acted at or ahead of the hawser point, the favoured side would flip (fig07).
  For a real tanker moored at the bow, with its superstructure aft, the wind always acts aft of the hawser point.
  Numbers (x_A/L = 0.25, wind 0.2 F_H0 from −120°, ship resting to starboard at ψ = −12.9°, γ = −57.3°): the current's
  sideways force on the hull is **2.2 F_H0 to starboard**, the wind's is 0.17 F_H0 to port. Even without wind the
  angled rest position has a sideways current force of 1.7 F_H0. So the ship is held out to the side by the current
  (that is the pitchfork); the wind only decides **which** side, through the turning moment it makes about the hawser
  point, not through its push.
- **Wind and fishtailing.** A wind from astern lowers the hawser tension that holds the ship in line, so it fishtails
  sooner. Wind from ahead adds tension. A side wind holds the ship off to one side, which damps the swing.

## Things that look odd or are uncertain (flagged)
1. **The wind force does not depend on heading** (no wind coefficients). Fine for small heading changes; for large
   angles a real wind force changes with the angle between wind and hull. Open issue 23.
2. **Point of action = midship** is our choice. A tanker's wind force acts somewhat aft (superstructure). fig07 shows
   x_W = −0.3 L: same picture, slightly stronger effect. Only an x_W ahead of the hawser point would change the story.
3. **Wind force is given relative to F_H0, which contains C_F.** The rest headings depend only on C_W. At full-scale
   friction (C_F ≈ 0.015) the same C_W is about 3 × more F_H0, so the real ship needs a smaller wind to cross over.
4. The wind is steady (no gusts). Gusts would be forcing, i.e. the cut step 9.
5. fig08: we only followed the Hopf point of the single rest position at the bow attachment, for three directions.
