# Open issues

Move an item to the decision log when it is settled. Last updated 2026-11-02 (step 8 control code, issue 25 added).

## Settled (kept here briefly so the numbering stays traceable)
1. Hull choice → Vidal de Negreiros (D3, agreed).
2. Sectional drag `C_D(x)` → straight line, more drag aft (D6, agreed, sign corrected).
3. Interactive demo format → "Python computes, HTML plays back" (D8, agreed); file format `data.js` (D22).
4. Sign of the hawser side force → `F_Y = −F_H sin(γ+ψ)`, eq. B.3; eq. 5.3 is a misprint (D17, sympy).
5. Direction of `l_p` → aft; the "+" in (4.28) is a sign slip (D6).
6. `l_p` vs `l_p/L` → it means `l_p/L` (D6; units of I_1 in 4.26).
7. `I_z` reference point → about the CG, shifted to midship (D13). E06 shows it moves the Hopf curve by < 0.3 %.
8. `+M11 u r`, `+M26 u r` → kept; explained on thesis p. 116 and confirmed by the yaw-rotating check (D14).
9. Bare-hull check → done with a hull pinned at midship (D15, E01).
12. Simos §5.2 read: periods (Fig. 5.7) and amplitudes (Fig. 5.17) compared in E05. §5.4 (sensitivity) still unread.

## Still open
10. **Fluid properties** — mostly settled (D27): ρ, ν now cited from ITTC (2011), R15. Left: the water temperature of
    the IPT tests is unknown; ν = 1.0e-6 (≈ 20 °C) fits the 0.2 m/s row of Table 5-1 slightly better than our 1.14e-6
    (3 % difference in C_F). One sentence in the report is enough.
11. **Pitchfork position, 100 %.** Eq. 3.37 and our model both give x_Tcr/L = 0.359. Leite's theory curve in Simos
    Fig. 3.21 seems to reach zero nearer 0.38–0.39 (read by eye). Cause unknown; Leite et al. (1998) is paywalled.
13. R8 (post-pitchfork paths depend on the hydrodynamic model) — read the abstract if it can be found; it affects how
    confidently we can talk about the angled equilibria.
14. Brief §11 references not yet checked (R13 in sources).
15. **Table 5-1 is not self-consistent** (found in E01, fig08): its C_F at 0.4–0.6 m/s match the k = 0.25 friction
    line, not k = 0.39 as the caption says. Matters only when comparing with Simos's own runs (E05 uses his C_F values
    directly).
16. **D6 comparison run (shape of C_D(x)).** Partly answered: C_Y and l_p enter the linearisation only at second order,
    so the **stability boundaries do not depend on the C_D(x) shape at all** (E06). What is left: their effect on the
    swing amplitude (E05). Simos says it is small (p. 190); could be checked by re-running E05 fig05 with his Fig. 4.4
    curve (needs digitising) or simply with C_Y changed by ±20 %.
17. **Froude (wave) effects.** Simos needed an extra moment correction k_w for the 40 % tests at 0.4–0.6 m/s
    (eq. 5.4, Table 5-3). Our model has none. E03 fig07 shows its effect on the pitchfork: with his k_w our angles are
    2–3° below his measured ones. Decide whether the report mentions it (suggestion: yes, one paragraph, as a limit of
    the model at tank speeds; it does not matter at full-scale mooring speeds).
18. **40 % Hopf boundary is 5.6 % above Simos's Fig. 5.3** (100 %: +0.5 %). Not caused by D13. A 3 % change in the
    draft would remove it, so it is within the uncertainty of the hull data/normalisation. Cause unknown.
19. **Compression in the hawser** → answered by E07: with a real rope the hawser goes slack for most of each swing
    and snaps tight to several times the rigid-model force. Remaining: rope damping/hysteresis is not modelled, so the
    peaks are high-side estimates.
22. **E07 at 0.75 m/s** breaks the trend of swing size vs speed: possibly a second stable snap pattern (coexisting
    cycles). Could be checked by starting from several initial conditions.
23. **Wind model is a constant force (D28).** No wind coefficients, so the force does not change with heading, and
    the point of action is midship. For the report: give results for a range of F_w/F_H0, not for a wind speed. If a
    readable source for tanker wind coefficients turns up, the force could depend on the wind-to-hull angle.
24. **Point of action of the wind (D28):** midship or aft? Student to decide; E10 fig07 shows it only matters if it
    is moved ahead of the hawser point.
20. **Supercritical vs subcritical Hopf.** Our model gives a supercritical Hopf (no hysteresis, E05). Aghamohammadi &
    Thompson (1990, R6) observed a subcritical one in a flume with a different set-up (not read; abstract only). Do not
    force agreement (brief §3); say in the report that the type can depend on the hydrodynamic model and set-up.
21. **Both bifurcations have "linear" normal forms** (E03 fig06, E05 fig02): branches grow ∝ distance, not √distance,
    because the leading nonlinearity is the quadratic cross-flow drag v|v|, r|r|. This is a property of the model, not
    a bug; checked with an averaging argument and a toy oscillator. For the report: present it as an explained
    difference from Strogatz's smooth normal forms. Worth one sentence of caution: it depends on the cross-flow model,
    which is the least certain part of the hydrodynamics near ψ = 0.
25. **Control numbers are PLACEHOLDERs (D31).** Thruster position x_T/L = −0.5, tug position −0.5, and the thrust
    limit (as a fraction of F_H0) have no source yet (`fishtail/params.py`). Needed: a cited stern-thruster thrust for
    a VLCC-class tanker, a cited hold-back tug pull, and ideally an OCIMF (or similar) source for tugs at SPMs. Until
    then, E08 should show results against the gains and the limit ratio, not as "a real thruster can do this".
26. **Stern thruster or hold-back tug?** Not yet checked against a source: whether VLCC-class tankers usually carry a
    stern thruster, and whether hold-back tugs are the usual practice at SPMs. Find a source, then decide with the
    student which option E08 leads with (the code supports both, D31).
