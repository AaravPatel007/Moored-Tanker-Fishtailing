# Why current speed (almost) drops out of the rigid-hawser model

Written for the report. Plain-language first, then the formal version, then the correction that makes speed matter.

## 1. The idea in one sentence

If every force in the problem grows with current speed in exactly the same way, a faster current just replays the same
motion faster — it cannot change *which* motion happens.

## 2. What each force depends on

| Force | Depends on | If the current doubles |
|---|---|---|
| Water force on the hull, `½ ρ L T U² C(α)` | speed², and the flow angle α | ×4 |
| Extra water force from turning (cross-flow drag) | (local sideways water speed)², and that speed is `v + x r` | ×4 |
| Rope tension (rigid rope) | whatever is needed to keep the length fixed, so it balances the other forces | ×4 |
| Mass × acceleration | mass is fixed, acceleration changes (see below) | ×4 in the sped-up film |

## 3. The film argument

Film the ship in a current `U`. Play the film at double speed.
- Same path, same angles; every event takes half the time.
- Every velocity doubles. Every acceleration ×4 (velocity doubled, and the change happens in half the time).
- So mass × acceleration ×4. Water forces, computed from doubled velocities, ×4. Rope tension follows, ×4.

Every term in Newton's law is multiplied by 4, so the equation still balances: the sped-up film is a valid motion at
current `2U`. Therefore *if the coefficients C(α) are fixed numbers*, the boundaries between "trails straight",
"sits at an angle" and "fishtails" cannot depend on `U`.

## 4. Formal version: non-dimensionalisation

*Non-dimensionalisation* = rewriting the equations in the problem's own natural units so that only pure numbers
(no metres, no seconds) remain.

- lengths in units of ship length `L`;
- time in units of `L/U` (the time the current takes to pass one ship length): `t* = U t / L`;
- velocities in units of `U`, yaw rate in units of `U/L`;
- forces in units of `M U² / L` (Simos 2001, Annex B, eqs. B.4–B.6).

After this change `U` disappears from the equations of motion (B.7) and from the hawser kinematics (B.8).
What is left are pure numbers: `L_H/L`, `x_A/L`, mass ratios like `(M + M22)/M`, and the hull coefficients.

## 5. The correction: speed sneaks back in through friction

The hull coefficients are *not* all fixed numbers. The surge (along-ship) friction coefficient `C_F` depends on the
**Reynolds number** `Re = U L / ν` (ν = water's kinematic viscosity). The ITTC friction line (Simos eq. 3.20) gives

    C_0 ≈ 0.075 / (log10 Re − 2)² · S/(T L) · (1 + k)

so a faster current → larger `Re` → **smaller** `C_F`.

Why this matters so much: when the ship trails behind the buoy, the hawser tension is equal to the along-ship drag,
and that tension is what pulls the bow back into line — like gravity on a pendulum (Simos §5.1). Less friction →
less tension → weaker "straightening" → less stable.

Consequences:
1. **In the rigid model, current speed acts only through `C_F(Re)`.** Simos (2001) Fig. 5.3 plots the stability
   boundary with `C_F` on the vertical axis for exactly this reason.
2. **Model tests are more stable than real ships.** Model tests use small `Re` → large `C_F` → more stable.
   Simos (p. 194–195), **40 % (ballast) condition**: at full scale, 1.90 m/s gives `C_F` ≈ 0.015, unstable for any
   practical hawser length, while the 1:90 model at the matching speed (0.2 m/s) was stable for every hawser length.
   Our code reproduces both numbers: C_F = 0.0150 full scale, 0.047 model (E01, fig08) — the model has 3× the friction.
3. **Other routes for speed** (outside the rigid model): a stretchy hawser (rope stiffness vs flow force), the
   thruster controller (fixed gains), wind (its own speed), and wave-making (Froude number) — the last is small at
   mooring speeds (Simos p. 193: `0.034 ≤ U/√(gL) ≤ 0.102` in the tests).

## 6. Figure to make (experiment E02)

Run the rigid model at two current speeds with `C_F` held fixed; plot `ψ` against time in seconds (curves differ) and
against `t* = U t / L` (curves lie exactly on top of each other). Then repeat with `C_F(Re)` switched on and show the
curves separate — that is the whole argument in one figure.
