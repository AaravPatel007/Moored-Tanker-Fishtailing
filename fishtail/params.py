"""Every number the model uses, in one place, with its source and units.

Sign convention (Simos 2001 thesis, used everywhere in this package):
    x forward, y to PORT, z up, origin at MIDSHIP; angle of attack alpha = 180 deg means current from dead ahead;
    heading psi = pi - alpha. Velocities u, v are relative to the water.

Sources are numbered as in notes/sources.md. R2 = Simos (2001) PhD thesis; "p." = printed page number.
"PLACEHOLDER" marks a number that has no citable source yet.
All values are SI units: m, kg, s, N.
"""

from dataclasses import dataclass, replace


@dataclass(frozen=True)
class Fluid:
    """Water properties."""

    rho: float  # density [kg/m^3]
    nu: float  # kinematic viscosity [m^2/s]


# Source R15: ITTC (2011) procedure 7.5-02-01-03, values at 15 degC (decision D27).
# Sea water: Table 3 gives rho = 1026.02, nu = 1.1892e-6. rho = 1025 is the usual design value (student's choice, -0.1 %).
# Fresh water: Table 1 gives rho = 999.10, nu = 1.1386e-6 (rounded here).
SEA_WATER = Fluid(rho=1025.0, nu=1.19e-6)  # about 15 degC, full scale
FRESH_WATER = Fluid(rho=1000.0, nu=1.14e-6)  # about 15 degC, towing tank (model tests)


@dataclass(frozen=True)
class Hull:
    """One ship in one loading condition."""

    name: str
    L: float  # length [m]
    B: float  # beam [m]
    T: float  # draft [m]
    C_B: float  # block coefficient [-]
    S: float  # wetted surface [m^2]
    M: float  # mass [kg]
    I_z: float  # yaw moment of inertia [kg m^2]; taken as about the CG (see I_Zo below)
    x_CG: float  # centre of gravity, forward of midship [m]
    M11: float  # added mass in surge [kg]
    M22: float  # added mass in sway [kg]
    M66: float  # added moment of inertia in yaw [kg m^2]
    M26: float  # sway-yaw coupling added mass [kg m]
    C_Y: float  # sideways drag coefficient, current exactly side-on (alpha = 90 deg) [-]
    lp_L: float  # l_p / L: the sideways drag acts this far AFT of midship, as a fraction of L [-]
    k_form: float  # form factor in the friction line, eq. (3.20) [-]
    scale: float  # model scale used in the tank tests (model is 1:scale) [-]

    @property
    def I_Zo(self):
        """Yaw inertia about MIDSHIP [kg m^2], which is what the equations of motion (3.9) need.

        Annex A does not say about which point I_z is measured. sqrt(I_z / M) = 80 m = 0.25 L, the usual
        rule-of-thumb radius of gyration about the CG, so we take I_z about the CG and move it to midship with the
        parallel-axis rule: I_Zo = I_z + M x_CG^2 (decision D13; changes I by about 1.5 %).
        """
        return self.I_z + self.M * self.x_CG**2


# --- Vidal de Negreiros VLCC, full scale -------------------------------------------------------------------------
# Masses, sizes, S, C_B, x_CG, scale: R2 Annex A, p. 263-264.
# k_form = 0.39: R2 Table 5-1 caption, p. 192 (used by Simos in Ch. 5).

# 100 % (fully loaded). C_Y = 0.84: R2 Table 3-1, p. 60. l_p/L = 0.059: R2 Table 3-2, p. 63 (both from Leite et al. 1998).
VIDAL_100 = Hull(
    name="Vidal de Negreiros 100 %",
    L=320.0, B=54.5, T=21.47, C_B=0.83, S=2.734e4,
    M=3.219e8, I_z=2.06e12, x_CG=9.81,
    M11=1.848e7, M22=2.725e8, M66=1.579e12, M26=2.209e9,
    C_Y=0.84, lp_L=0.059, k_form=0.39, scale=90.0,
)

# Same hull with the cross-flow coefficients from the later IPT (1998) captive tests: C_Y = 0.86, l_p C_Y = 0.044
# (R2 Table 4-5, p. 141). Simos used these for his yaw-rotating comparisons (Figs. 4.9-4.17), so we do too.
VIDAL_100_IPT = replace(VIDAL_100, name="Vidal de Negreiros 100 % (IPT 1998 C_Y)", C_Y=0.86, lp_L=0.044 / 0.86)

# 40 % (ballast). Leite et al. did not test this condition; C_Y = 0.52, l_p C_Y = 0.052 come from R2 Table 4-5, p. 141.
VIDAL_40 = Hull(
    name="Vidal de Negreiros 40 %",
    L=320.0, B=54.5, T=8.63, C_B=0.76, S=1.791e4,
    M=1.117e8, I_z=0.72e12, x_CG=4.44,
    M11=0.254e7, M22=0.353e8, M66=0.222e12, M26=0.353e9,
    C_Y=0.52, lp_L=0.052 / 0.52, k_form=0.39, scale=90.0,
)


def froude_scaled(hull, rho_ratio):
    """The same hull at model scale 1:hull.scale (Froude scaling).

    Lengths are divided by scale, areas by scale^2, masses by scale^3, (mass x length) by scale^4 and
    moments of inertia by scale^5. Masses are also multiplied by rho_ratio = rho_model / rho_ship, because the
    model floats in fresh water at the same draft. Pure numbers (C_B, C_Y, l_p/L, k) do not change.
    """
    lam = hull.scale
    return replace(
        hull,
        name=hull.name + f" (model 1:{lam:g})",
        L=hull.L / lam, B=hull.B / lam, T=hull.T / lam, S=hull.S / lam**2,
        M=hull.M * rho_ratio / lam**3, I_z=hull.I_z * rho_ratio / lam**5, x_CG=hull.x_CG / lam,
        M11=hull.M11 * rho_ratio / lam**3, M22=hull.M22 * rho_ratio / lam**3,
        M66=hull.M66 * rho_ratio / lam**5, M26=hull.M26 * rho_ratio / lam**4,
        scale=1.0,
    )


# --- Published results used only for checking (not inputs to the model) ------------------------------------------

# R2 Table 5-1, p. 192: model-scale (1:90) towing speeds, measured total resistance C_T (min, max) and Simos's
# friction estimate C_F from eq. (3.20) "with k = 0.39".
TABLE_5_1_U_MODEL = [0.2, 0.4, 0.5, 0.6]  # [m/s]
TABLE_5_1_CF = {"100": [0.0279, 0.0216, 0.0206, 0.0199], "40": [0.0459, 0.0352, 0.0336, 0.0324]}
TABLE_5_1_CT = {
    "100": [(0.0254, 0.0336), (0.0223, 0.0247), (0.0211, 0.0227), (0.0199, 0.0207)],
    "40": [(0.0396, 0.0544), (0.0354, 0.0392), (0.0335, 0.0352), (0.0320, 0.0329)],
}

# R2 p. 195: 40 % condition at full scale, U = 1.90 m/s -> C_F = 0.015.
FULL_SCALE_CF_CHECK = {"U": 1.90, "C_F": 0.015}

# R2 p. 197-198 (with eq. 3.37, p. 87): critical turret position for the 40 % condition, x_Tcr/L ~ 0.19.
X_TCR_40_CHECK = 0.19

# R2 Table 5-2, p. 198: measured steady heading of the 1:90 model in the static-bifurcation tests,
# 40 % condition, L_H/L = 0.287, x_A/L = 0.20 (model speeds in m/s -> heading in degrees, "~1" printed for 0.20).
TABLE_5_2 = {"x_A_L": 0.20, "L_H_L": 0.287, "U_model": [0.20, 0.40, 0.50, 0.60], "psi_deg": [1.0, 13.0, 14.0, 17.0]}
# R2 Fig. 5.3, p. 191 (= Fig. B.1, p. 273): Simos's stability boundary (Hopf curve) in the (L_H/L, C_F) plane;
# stable above the curve. 100 % at x_A/L = 0.50, 40 % at x_A/L = 0.40. DIGITISED BY US from a 300-dpi render of the
# printed figure: frame corners located automatically, the curve traced column by column, checked by overlay
# (scratch work, not kept). Read accuracy about +-0.0005 in C_F.
SIMOS_FIG53 = {
    "L_H_L": [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 1.95],
    "100": {"x_A_L": 0.50, "C_F": [0.0110, 0.0186, 0.0244, 0.0289, 0.0354, 0.0395, 0.0424, 0.0443, 0.0459, 0.0470,
                                   0.0477, 0.0483, 0.0489, 0.0494, 0.0493, 0.0493, 0.0489]},
    "40": {"x_A_L": 0.40, "C_F": [0.0121, 0.0190, 0.0245, 0.0286, 0.0342, 0.0368, 0.0384, 0.0390, 0.0392, 0.0388,
                                  0.0385, 0.0379, 0.0365, 0.0351, 0.0336, 0.0321, 0.0310]},
}

# R2 Fig. 5.7, p. 201: Simos's SIMULATED fishtailing periods, Vidal 100 %, 1:90 model at U = 0.40 m/s (x_A/L = 0.5,
# C_Y = 0.86, l_p C_Y = 0.044), for four friction values inside the Table 5-1 range. Read by us off the printed
# gridlines, +-1 s. Periods in seconds at model scale.
SIMOS_FIG57 = {"U_model": 0.40, "C_F": [0.0223, 0.0231, 0.0239, 0.0247],
               "period_s": {"0.36": [103, 100, 98, 96], "0.6": [123, 120, 118, 116], "0.9": [138, 134, 132, 131]}}

# R2 Fig. 5.17, p. 211: fishtailing amplitudes (degrees, half of peak-to-peak) vs L_H/L, Vidal 100 %, model at
# U = 0.4 m/s, x_A/L = 0.5. READ BY EYE by us from the printed figure, +-2 deg.
# Experiments (+ gamma, o psi), and Simos's theory curves for C_F min = 0.0223 and C_F max = 0.0247 (Table 5-1 range).
SIMOS_FIG517 = {
    "exp_L_H_L": [0.10, 0.19, 0.36, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00],
    "exp_gamma_deg": [0.0, 24.5, 43.0, 49.0, 49.5, 44.5, 46.0, 43.0, 38.0],
    "exp_psi_deg": [0.0, 4.5, 10.5, 15.0, 16.0, 16.0, 19.0, 19.5, 15.5],
    "theory_L_H_L": [0.2, 0.35, 0.6, 0.8, 1.0, 1.35],
    "theory_gamma_CFmin_deg": [50.0, 55.0, 50.0, 46.0, 42.0, 36.0],
    "theory_gamma_CFmax_deg": [38.0, 48.5, 45.5, 41.0, 37.0, 32.5],
    "theory_psi_deg": [8.0, 12.0, 15.5, 17.0, 18.5, 19.5],
    "C_F_min": 0.0223, "C_F_max": 0.0247,
}

# --- Stretchy hawser (step 7) ---------------------------------------------------------------------------------------
# R11 = Huang et al. (2026), J. Mar. Sci. Eng. 14, 770, Fig. 3 ("prototype" curve) and Table 4: a real SPM nylon
# hawser, length 70 m, diameter 144 mm, breaking load 650 t, used for a 220 m tanker. Points READ BY EYE by us at the
# marked points of the printed figure, +-0.2 m and +-10 kN. Elongation [m] -> tension [kN].
HUANG_HAWSER = {
    "length_m": 70.0, "diameter_mm": 144.0, "breaking_load_t": 650.0,
    "elongation_m": [0.0, 4.3, 6.4, 9.35, 13.55],
    "tension_kN": [0.0, 60.0, 120.0, 248.0, 505.0],
}

# R2 Table 5-3, p. 217: Simos's ad hoc wave-making (Froude) correction of the static yaw moment, C6C (1 + k_w),
# eq. (5.4), fitted by him (Annex C) to the Table 5-2 tests. Not used in our model; only in E03 to show its effect.
TABLE_5_3_KW = {"U_model": [0.20, 0.40, 0.50, 0.60], "k_w": [0.0, 0.77, 0.82, 1.03]}

# --- Control (step 8, decision D31) ------------------------------------------------------------------------------------
# PLACEHOLDER: no cited values yet (open issue: find a cited tug pull / stern-thruster thrust for a VLCC-class tanker,
# and an OCIMF source for hold-back tugs at SPMs). All are given as fractions so they can be changed in one place.
X_T_L = -0.5  # where the stern thruster pushes, forward of midship as a fraction of L [-]: the stern.  PLACEHOLDER
X_TUG_L = -0.5  # where the hold-back tug's line is made fast: the stern [-].                           PLACEHOLDER
