# Rigid-hawser model: sympy results

Written by `derivations/rigid_hawser.py` (do not edit by hand). Symbols: `m_26 = M x_CG + M26`, `I_66 = I_Zo + M66`; `X_R, Y_R, N_R` = water force and moment on the hull.

**1. u (relative to the water)**

```
L_H⋅γ_dot⋅sin(γ + ψ) + U⋅cos(ψ)
```

**   v (relative to the water)**

```
L_H⋅γ_dot⋅cos(γ + ψ) - U⋅sin(ψ) - ψ_dot⋅x_A
```

Same as Simos eq. 5.2 / B.2: u YES, v YES.

**2. Hawser force along the ship F_X**

```
F_H⋅cos(γ + ψ)
```

**   Hawser force sideways F_Y**

```
-F_H⋅sin(γ + ψ)
```

**   Hawser moment about midship N**

```
-F_H⋅x_A⋅sin(γ + ψ)
```

F_Y = +F_H sin(gamma+psi) (eq. 5.3): NO;  F_Y = -F_H sin(gamma+psi) (eq. B.3): YES. **Open issue 4 settled: B.3 is right, the '+' in 5.3 is a misprint.**

**3. Matrix of the linear system  Mat @ [gamma_ddot, psi_ddot, F_H] = rhs**

```
⎡L_H⋅(M + M₁₁)⋅sin(γ + ψ)            0              -cos(γ + ψ)  ⎤
⎢                                                                ⎥
⎢L_H⋅(M + M₂₂)⋅cos(γ + ψ)  -M⋅x_A - M₂₂⋅x_A + m₂₆    sin(γ + ψ)  ⎥
⎢                                                                ⎥
⎣   L_H⋅m₂₆⋅cos(γ + ψ)         I₆₆ - m₂₆⋅x_A       x_A⋅sin(γ + ψ)⎦
```

**   right-hand side rhs**

```
⎡             2                           2                                                                                       2                                                          2                      2⎤
⎢- L_H⋅M⋅γ_dot ⋅cos(γ + ψ) - L_H⋅M₁₁⋅γ_dot ⋅cos(γ + ψ) - L_H⋅M₁₁⋅γ_dot⋅ψ_dot⋅cos(γ + ψ) + L_H⋅M₂₂⋅γ_dot⋅ψ_dot⋅cos(γ + ψ) - M⋅ψ_dot ⋅x_A + M₁₁⋅U⋅ψ_dot⋅sin(ψ) - M₂₂⋅U⋅ψ_dot⋅sin(ψ) - M₂₂⋅ψ_dot ⋅x_A + X_R + m₂₆⋅ψ_dot ⎥
⎢                                                                                                                                                                                                                    ⎥
⎢                                  2                                                            2                                                                                                                    ⎥
⎢                       L_H⋅M⋅γ_dot ⋅sin(γ + ψ) - L_H⋅M₁₁⋅γ_dot⋅ψ_dot⋅sin(γ + ψ) + L_H⋅M₂₂⋅γ_dot ⋅sin(γ + ψ) + L_H⋅M₂₂⋅γ_dot⋅ψ_dot⋅sin(γ + ψ) - M₁₁⋅U⋅ψ_dot⋅cos(ψ) + M₂₂⋅U⋅ψ_dot⋅cos(ψ) + Y_R                        ⎥
⎢                                                                                                                                                                                                                    ⎥
⎢                                                                                                   2                                                                                                                ⎥
⎣                                                                                          L_H⋅γ_dot ⋅m₂₆⋅sin(γ + ψ) + N_R                                                                                           ⎦
```

**4. Linearised about gamma = psi = 0:  A q'' + B q' + C q = 0.   Mass matrix A**

```
⎡L_H⋅(M + M₂₂)  m₂₆ - x_A⋅(M + M₂₂)⎤
⎢                                  ⎥
⎣   L_H⋅m₂₆        I₆₆ - m₂₆⋅x_A   ⎦
```

**   Damping matrix B**

```
⎡-L_H⋅Yᵥ  M₁₁⋅U - M₂₂⋅U - Yᵣ + Yᵥ⋅x_A⎤
⎢                                    ⎥
⎣-L_H⋅Nᵥ         -Nᵣ + Nᵥ⋅x_A        ⎦
```

**   Stiffness matrix C**

```
⎡  F_H0      F_H0 + U⋅Yᵥ  ⎤
⎢                         ⎥
⎣F_H0⋅x_A  F_H0⋅x_A + Nᵥ⋅U⎦
```

**   det C**

```
F_H0⋅U⋅(Nᵥ - Yᵥ⋅x_A)
```

With F_H0 > 0 (tension) and U > 0, det C > 0 exactly when N_v - x_A Y_v > 0. Writing Y_v = (1/2) rho L T U C2', N_v = (1/2) rho L^2 T U C6' (C2', C6' = slopes of the static coefficients at head-on, both negative) this is x_A > L C6'/C2' = x_Tcr (eq. 3.37). det C changes sign there, so one real eigenvalue crosses zero: **the pitchfork is at x_A = x_Tcr for every hawser length and friction coefficient.**

Compared with Simos eq. B.18: the mass matrix and the stiffness matrix have the same structure. Our damping matrix has extra terms -Y_v L_H (first column) and +Y_v x_A: the sideways lift force caused by the sideways velocity L_H gamma' - x_A psi' of the hull. Simos's B.13 writes the lift with the heading psi only, so his B.15/B.18 has zeros in the first column. See experiments/E04 and E06 for what this changes.
