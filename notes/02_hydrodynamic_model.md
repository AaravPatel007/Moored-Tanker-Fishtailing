# The hydrodynamic model (what the water does to the hull)

Source: A. N. Simos (2001), PhD thesis, Escola Politécnica USP — *Modelo hidrodinâmico heurístico para análise de
navios petroleiros amarrados sujeitos à ação de correnteza* (in Portuguese). Local copy:
`references/Simos_2001_PhD_thesis_heuristic_current_model.pdf`. Equation numbers below are the thesis's.
"p." = printed page number in the thesis (PDF page = printed page + 23).

**Status: coded in `fishtail/hydro.py` and checked (step 2, experiment E01, `tests/test_hydro.py`).**
The ⚠ marks below were open issues; how each was settled is noted next to it (decisions D6, D11–D14).

---

## 0. Words used here

- **Heuristic model:** a model built from physical reasoning plus a few measured numbers, not solved from first
  principles. "Quasi-explicit" = needs only a handful of measured numbers.
- **Angle of attack α:** angle between the current and the ship's centreline. α = 180° → current from dead ahead.
- **Cross-flow principle:** treat each slice of the hull as a 2-D cylinder in a sideways flow; its drag ∝ (sideways
  speed)². Good when the flow hits the hull broadside.
- **Wing (lift) effects:** at small α the hull acts like a very short, stubby wing (aspect ratio 2T/L) and produces lift.
- **Munk moment:** a turning moment that exists even in a perfect (frictionless) fluid; it turns a slender body
  *broadside* to the flow. This is the destabilising moment in fishtailing.

## 1. Coordinates and equations of motion (§3.1, p. 31–37)

- Body axes at **midship** `o`: `x` forward, `y` to **port** (not ITTC!), right-handed. CG at `(x_CG, 0)`.
- `u, v` = velocity of the hull **relative to the water**; `r` = yaw rate.
- Relative-velocity kinematics (3.3): `u = u_f + U cos ψ`, `v = v_f − U sin ψ`, `r = ψ̇`,
  where `(u_f, v_f)` is the hull velocity over the ground in body axes and `ψ = π − α` (3.1).

Equations of motion, added-mass (ideal-fluid) terms moved to the left (3.9):

```
(M + M11) u̇ − (M + M22) v r − (M x_CG + M26) r²  = F_X,R(u, v, r)
(M + M22) v̇ + (M + M11) u r + (M x_CG + M26) ṙ  = F_Y,R(u, v, r)
(I_Z° + M66) ṙ + (M x_CG + M26)(v̇ + u r)         = N_Z,R(u, v, r) + (M11 − M22) u v
```
- `I_Z°` is the yaw inertia about midship. ⚠ Annex A gives "I_z" without saying about which point.
  **Settled (D13):** √(I_z/M) = 0.25 L, so read as about the CG; I_Z° = I_z + M x_CG².
- Footnote to (3.8): the Munk term `(M11 − M22) u v` is kept on the right **only to show it** — Leite's static yaw
  coefficient `C6C` already contains it. **In the code: drop the explicit `(M11 − M22) u v` and use `C6C`.**

## 2. Split into static part + yaw-rate part (3.10)

```
F_X,R(u, v, r) = F_X,R(u, v, 0) + ΔF_X,R(u, v, r),   ΔF_X,R = 0 when r = 0   (same for F_Y and N_Z)
```

## 3. Static part — Leite et al. (1998) (§3.2, p. 48–65)

```
F_X,R(u,v,0) = ½ ρ U² L T · C1C(α)
F_Y,R(u,v,0) = ½ ρ U² L T · C2C(α)                          (3.12)
N_Z,R(u,v,0) = ½ ρ U² L² T · C6C(α)

U = √(u² + v²),   α = π + arctan(v/u)                        (3.13)
```
Code uses α = π + atan2(v, u) (D11): the same for u > 0, and also right for u < 0.

**Surge** (3.14), with friction coefficient (3.20):
```
C1C(α) = C0 cos α + (1/8)(π T / L)(cos 3α − cos α)
C0 ≈ 0.075 / (log10 Re − 2)² · S / (T L) · (1 + k)
```
`k` = form factor: Leite et al. suggest 0.25; Simos uses **0.39** in Ch. 5 (Table 5-1 caption, p. 192).
`S` = wetted surface area. Simos notes the surge force is the least reliable part (p. 55).

**Sway** (3.21):
```
C2C(α) = (C_Y − πT/(2L)) sin α |sin α| + (πT/(2L)) sin³ α + (πT/L)(1 + 0.4 C_B B / T) sin α |cos α|
```
`C_Y` = sideways drag coefficient measured at α = 90° (Table 3-1).

**Yaw** (3.29):
```
C6C(α) = −(l_p/L) C_Y sin α |sin α| − (πT/L) sin α cos α − ((1 + |cos α|)/2)² (πT/L)(1/2 − 2.4 T/L) sin α |cos α|
```
- 1st term: cross-flow moment; sideways drag acts a distance `l_p` **aft** of midship (p. 62).
- 2nd term: **Munk moment** (uses `M22 ≈ ½ ρ π T² L`, p. 63). Simos notes it can be replaced by the exact
  `(M11 − M22) u v` when added masses are known (p. 64).
- 3rd term: lift of the hull + rudder, from Clarke et al. (1982) — *reduces* the Munk moment near α = 180° (p. 64–65).

## 4. Yaw-rate part — Simos (2001) extension (§4.4, eqs. 4.36–4.38, p. 133–134)

```
ΔF_X,R = −¼ ρπT² L · v r − (1/16) ρπT² L² cos α_a · r²

ΔF_Y,R = ¼ ρπT² L (1 − 4.4 B/L + 0.160 B/T) u r + M11 u r
         − ½ ρ T L [ I0(v,r) v² − C_Y v|v| ] − ρ T L² I1(v,r) v r − ½ ρ T L³ I2(v,r) r²

ΔN_Z,R = −⅛ ρπT² L² (1 + 0.16 B/T − 2.2 B/L) |u| r + M26 u r
         − ½ ρ T L² [ I1(v,r) v² + l_p C_Y v|v| ] − ρ T L³ I2(v,r) v r − ½ ρ T L⁴ I3(v,r) r²
         − ½ ρ T L⁴ [ (C_Y/16) r|r| − I3(0,r) r² ]

α_a = arctan(−v/u + r L / (2u))
I_j(v, r) = (1/L^(j+1)) ∫_{−L/2}^{L/2} C_D(x) · sign(v + r x) · x^j dx,   j = 0, 1, 2, 3
```
- The terms in square brackets vanish when `r = 0`, so there is no double counting with the static part (p. 126).
- ⚠ In `l_p C_Y v|v|` and in (4.28) `l_p` must mean `l_p/L` for the units to work. **Settled (D6):** yes.
- ⚠ The `+ M11 u r` and `+ M26 u r` terms partly cancel left-hand-side added-mass terms. **Settled (D14):** thesis
  p. 116 — Clarke's measured Y_r, N_r include the ideal-fluid inertia, which (3.9) already has on the left. Our
  reproduction of the yaw-rotating test (Fig. 4.9) agrees with Simos's theory.
- cos α_a must have the sign of u (p. 119); code uses cos α_a = u/√(u² + (rL/2 − v)²) (D12).
- Where the pieces come from: "wing" terms (Clarke Y_r, N_r, eq. 4.14); cross-flow terms (strip integral, eq. 4.29);
  vortex term ΔN^GV (eq. 4.33), built so that a hull turning on the spot feels N = −(1/32)ρTL⁴C_Y r|r| (eq. 4.31,
  from 32 C_6,rr/C_Y ≈ 2 in the IPT yaw-rotation tests, Fig. 4.6, Table 4-3).

**Yaw-rotating test (§4.5) — what the figures show.** Hull towed at U and spun about midship at constant r:
u = U cos ψ, v = −U sin ψ (4.39). From (4.40), the water's total sideways force on the hull is
F_Y,R + (M22 − M11) U r cos ψ, and the moment is N_R. Plotted at model scale (1:90).

**Sectional drag `C_D(x)`** (§4.2.1, eqs. 4.27–4.28): must satisfy
```
(1/L)  ∫ C_D(x) dx     = C_Y
(1/L²) ∫ C_D(x) x dx   = l_p C_Y          (l_p read as l_p/L)
```
Simos: digitised Aarsnes curve × step weighting functions (ΔL = 0.15 L). Ours (decision D6): a straight line
`C_D(x) = C_Y [1 − 12 (l_p/L)(x/L)]`, which satisfies both conditions exactly with l_p measured **aft**.
**Sign settled:** the second condition is really (1/L²)∫C_D x dx = **−** l_p C_Y. The printed "+" contradicts
Fig. 4.4 (more drag aft), the r = 0 cancellation in (4.29), and eq. 3.33.

## 5. Rigid hawser (Ch. 5, p. 189; Annex B, p. 267–269)

Kinematics (5.2 = B.2), with hawser length `L_H`, hawser angle `γ` (to the current), attachment point `x_A` forward of midship:
```
u = U cos ψ + L_H γ̇ sin(γ + ψ)
v = −U sin ψ + L_H γ̇ cos(γ + ψ) − x_A ψ̇
r = ψ̇
```
Hawser loads with hawser force `F_H` (tension positive):
```
F_X = F_H cos(γ + ψ)
F_Y = −F_H sin(γ + ψ)        ← Annex B (B.3). Eq. 5.3 prints "+": a misprint (settled by sympy, D17).
N_Z = −x_A F_H sin(γ + ψ)
```
Non-dimensional form (B.4–B.7): `t* = U t / L`, velocities ÷ U, `r* = r L / U`, forces ÷ `M U²/L`, and
`A11 = (M+M11)/M`, `A22 = (M+M22)/M`, `A26 = (M x_CG + M26)/(M L)`, `A66 = (I_Z° + M66)/(M L²)`.

**Coded in `fishtail/model_rigid.py` (step 3).** Linearisation: `derivations/rigid_hawser_results.md`. Note (D19):
eq. B.13 writes the hull lift with ψ only and so drops the lift from the hull's sideways velocity; the full
linearisation keeps it and is what reproduces Simos's Fig. 5.3.

## 6. What the thesis found (to compare against)

- Hopf boundary in the (`L_H/L`, `C_F`) plane, Vidal de Negreiros: 100 % load at `x_A/L = 0.50`, 40 % at
  `x_A/L = 0.40` (Fig. 5.3, p. 191). Stable above the curve. Limit cycles in the unstable region.
- Full scale, U = 1.90 m/s → `C_F ≈ 0.015` → unstable for any practical hawser length (p. 195).
- Static bifurcation (pitchfork) when the attachment point is aft of x_Tcr, eq. 3.37:
  `x_Tcr/L = (½ + 2.4 T/L) / (1 + 0.40 C_B B/T)` → **0.19 for 40 %** (thesis p. 197–198) and **0.36 for 100 %** (ours).
  Simos (Annex B): x_A > x_Tcr is necessary (not sufficient) for stability of the rigid-hawser system.
- Cross-flow coefficients `C_Y`, `l_p` hardly affect stability, amplitude or period of the cycles; `C_F` does (p. 190).
