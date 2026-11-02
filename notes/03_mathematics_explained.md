# The mathematics, explained simply

Written 2026-10-17 for the student. Every equation here is in the code; the file and function are given next to it.
Thesis equation numbers (Simos 2001) are given as "eq. 5.2" etc. Open this file in VS Code's Markdown preview (or on
GitHub) to see the equations typeset.

Contents
1. The picture and the unknowns
2. Water forces (step 2, short reminder)
3. Ship on a rigid rope: the equations of motion (step 3)
4. Why current speed almost drops out
5. Where can the ship sit still? (equilibria, step 4)
6. Is that position stable? (linearisation and eigenvalues, step 4)
7. The pitchfork (E03)
8. The Hopf bifurcation (E04)
9. The fishtailing swing: limit cycles (step 6, E05)
10. Are there only three graphs? All the parameters
11. The stretchy rope (step 7, not started yet)

---

## 1. The picture and the unknowns

Seen from above: a buoy, a rope (hawser) of length $L_H$, and a ship of length $L$. The current flows past at speed
$U$. The rope is tied to the ship a distance $x_A$ forward of midship (the middle of the ship); $x_A = L/2$ is the bow.

Two angles describe everything:
- **hawser angle $\gamma$**: how far the rope has swung away from the current line ($\gamma > 0$: swung to port);
- **heading $\psi$**: which way the bow points ($\psi = 0$: straight into the current; $\psi > 0$: bow turned to port).

Because the rope keeps its length, these two angles fix where the ship is. The **state** of the system is the
two angles and how fast they change:

$$
\mathbf{x} = (\gamma,\ \psi,\ \dot\gamma,\ \dot\psi) \qquad \text{(4 numbers = 4 states)}
$$

Everything below answers one question: *given the state now, what is it a moment later?*

## 2. Water forces (step 2, short reminder)

Full details: [02_hydrodynamic_model.md](02_hydrodynamic_model.md). Code: `fishtail/hydro.py`.

The water pushes on the hull with a force along the ship $X_R$, a sideways force $Y_R$ and a turning moment $N_R$.
They depend on how the hull moves **relative to the water**: forward speed $u$, sideways speed $v$, turning rate $r$.
Each has a "not turning" part and an extra part for turning:

$$
Y_R = \tfrac12 \rho (u^2+v^2) L T \; C_{2C}(\alpha) \;+\; \Delta Y(u,v,r), \qquad
\alpha = \pi + \operatorname{atan2}(v, u)
$$

and the same for $X_R$ (with $C_{1C}$) and $N_R$ (with $C_{6C}$ and an extra factor $L$). $\alpha$ is the angle at
which the water meets the hull. The yaw coefficient $C_{6C}$ contains the **Munk moment**, the moment that turns a long
body broadside to the flow. It is what makes the ship want to sheer away.

## 3. Ship on a rigid rope: the equations of motion (step 3)

Code: `fishtail/model_rigid.py` (`velocities`, `accelerations`). Derivation by sympy: `derivations/rigid_hawser.py`,
printed results in `derivations/rigid_hawser_results.md`.

### 3a. Where is the ship, and how fast does it move through the water?

Put the buoy at the origin, the $X$ axis pointing upstream. The rope end (attachment point) is at

$$
\mathbf{A} = L_H\,(-\cos\gamma,\ \sin\gamma)
$$

and midship is $x_A$ behind it along the ship: $\mathbf{O} = \mathbf{A} - x_A(\cos\psi, \sin\psi)$.
Differentiate $\mathbf{O}$ in time, add the water's speed (the water moves at $-U$ along $X$), and split along the ship
and across it. Sympy does this and gets exactly Simos's eq. 5.2:

$$
\boxed{\;u = U\cos\psi + L_H\dot\gamma\,\sin(\gamma+\psi), \qquad
v = -U\sin\psi + L_H\dot\gamma\,\cos(\gamma+\psi) - x_A\dot\psi, \qquad r = \dot\psi\;}
$$

In words: the water passes the ship at the current speed (first terms), plus the ship's own swinging motion
($L_H\dot\gamma$ terms), plus the sideways speed of midship because the ship turns about the rope end ($x_A\dot\psi$).

### 3b. How does the rope pull?

A tight rope pulls the attachment point straight towards the buoy with force $F_H$. Split along and across the ship:

$$
F_X = F_H\cos(\gamma+\psi), \qquad F_Y = -F_H\sin(\gamma+\psi), \qquad N = -x_A F_H \sin(\gamma+\psi)
$$

The thesis prints $F_Y$ with a "+" in eq. 5.3 and a "−" in eq. B.3. The sympy derivation shows **B.3 is right**
(decision D17).

### 3c. Newton's law for the ship (Simos eq. 5.1)

Mass times acceleration = water force + rope force, in the three directions (along, across, turning). The
"mass" includes the **added mass** of water that moves with the hull ($M_{11}$, $M_{22}$, $M_{66}$, $M_{26}$):

$$
\begin{aligned}
(M+M_{11})\,\dot u - (M+M_{22})\,v r - m_{26}\,r^2 &= X_R + F_H\cos(\gamma+\psi)\\
(M+M_{22})\,\dot v + (M+M_{11})\,u r + m_{26}\,\dot r &= Y_R - F_H\sin(\gamma+\psi)\\
I_{66}\,\dot r + m_{26}\,(\dot v + u r) &= N_R - x_A F_H\sin(\gamma+\psi)
\end{aligned}
$$

with $m_{26} = M x_{CG} + M_{26}$ and $I_{66} = I_{Z} + M_{66}$ (yaw inertia about midship).

### 3d. Turning it into something a computer can step forward

Differentiate the speeds of 3a once more. Each time derivative splits into a part containing the unknown angular
accelerations and a known rest:

$$
\dot u = L_H \sin(\gamma+\psi)\,\ddot\gamma + a_u, \qquad
\dot v = L_H \cos(\gamma+\psi)\,\ddot\gamma - x_A\ddot\psi + a_v
$$

The rope force $F_H$ is also unknown: a rigid rope pulls exactly as hard as needed to keep its length. So there are
**three unknowns** $(\ddot\gamma, \ddot\psi, F_H)$ and **three equations**. They are *linear* in the unknowns, so at
every instant we solve a 3 × 3 system (with $s = \sin(\gamma+\psi)$, $c = \cos(\gamma+\psi)$):

$$
\begin{bmatrix}
(M+M_{11})L_H s & 0 & -c\\
(M+M_{22})L_H c & m_{26} - (M+M_{22})x_A & s\\
m_{26} L_H c & I_{66} - m_{26}x_A & x_A s
\end{bmatrix}
\begin{bmatrix}\ddot\gamma\\ \ddot\psi\\ F_H\end{bmatrix}
=
\begin{bmatrix}
X_R + (M+M_{22})vr + m_{26}r^2 - (M+M_{11})a_u\\
Y_R - (M+M_{11})ur - (M+M_{22})a_v\\
N_R - m_{26}(a_v + ur)
\end{bmatrix}
$$

That is the whole model. A standard solver (`scipy.integrate.solve_ivp`) then steps it forward in time. A test
checks that this hand-written matrix equals the one sympy derived (`tests/test_rigid.py`).

## 4. Why current speed almost drops out (E02)

Measure time in "how long the current takes to pass one ship length": $t^* = Ut/L$. Then every term in 3c scales
like $U^2$, and $U$ cancels. The one exception is skin friction: its coefficient depends on the Reynolds number
$Re = UL/\nu$ (ITTC line, eq. 3.20):

$$
C_F \approx \frac{0.075}{(\log_{10} Re - 2)^2}\cdot\frac{S}{TL}\,(1+k)
$$

So in the rigid-rope model, **current speed matters only through $C_F$**: faster current → smaller $C_F$.
E02 fig01 shows the same run at three speeds collapsing onto one curve in $t^*$ (difference $4\times10^{-7}$ degrees).

## 5. Where can the ship sit still? (equilibria, step 4)

Code: `fishtail/continuation.py` (`equilibrium_headings`, `equilibrium_state`).

"Sitting still" means $\dot\gamma = \dot\psi = 0$ and no accelerations. Then $u = U\cos\psi$, $v = -U\sin\psi$, $r=0$
and $\alpha = \pi - \psi$. The three equations of 3c become

$$
X(\psi) + F_H\cos(\gamma+\psi) = 0,\qquad Y(\psi) - F_H\sin(\gamma+\psi) = 0,\qquad N(\psi) - x_A F_H\sin(\gamma+\psi) = 0
$$

Take the third minus $x_A$ times the second: the rope force drops out and only the heading is left,

$$
\boxed{\;N(\psi) - x_A\,Y(\psi) = 0 \quad\Longleftrightarrow\quad C_{6C}(\pi-\psi) - \frac{x_A}{L}\,C_{2C}(\pi-\psi) = 0\;}
$$

In words: **the ship sits still when the water's turning moment about the rope end is zero.** This depends only on
$x_A/L$: not on rope length, not on friction. $\psi = 0$ (trailing straight) always solves it. Other solutions are
found numerically (scan, then Brent's method). The first two equations then give the rope force and direction:
$F_H = \sqrt{X^2+Y^2}$, $\ \gamma + \psi = \operatorname{atan2}(Y, -X)$.

## 6. Is that position stable? (linearisation and eigenvalues, step 4)

Code: `fishtail/model_rigid.py` (`linear_matrices`, `linear_jacobian_star`, `jacobian_star`).

Give the ship a tiny nudge and keep only terms proportional to the nudge (Strogatz §6.3). Near trailing straight the
water forces become straight lines in $v$ and $r$:

$$
Y_R \approx Y_v v + Y_r r,\qquad N_R \approx N_v v + N_r r,\qquad X_R \approx -F_{H0}
$$

$$
Y_v = \tfrac12\rho L T U\,C_2',\quad N_v = \tfrac12\rho L^2 T U\,C_6',\quad
C_2' = -\tfrac{\pi T}{L}\Big(1 + 0.4\,C_B\tfrac{B}{T}\Big),\quad C_6' = -\tfrac{\pi T}{L}\Big(\tfrac12 + 2.4\tfrac{T}{L}\Big)
$$

$Y_r, N_r$ come from the "wing" turning terms, and $F_{H0} = C_F \tfrac12\rho U^2 LT$ is the rope tension when trailing
straight. Sympy then gives a **mass–damping–stiffness** system for $q = (\gamma, \psi)$:

$$
A\ddot q + B\dot q + C q = 0
$$

$$
A=\begin{bmatrix} L_H m_{22} & m_{26}-x_A m_{22}\\ L_H m_{26} & I_{66}-m_{26}x_A\end{bmatrix},\quad
B=\begin{bmatrix} -L_H Y_v & (M_{11}-M_{22})U - Y_r + Y_v x_A\\ -L_H N_v & -N_r + N_v x_A\end{bmatrix},\quad
C=\begin{bmatrix} F_{H0} & F_{H0}+U Y_v\\ F_{H0}x_A & F_{H0}x_A + U N_v\end{bmatrix}
$$

($m_{22} = M + M_{22}$.) Written as four first-order equations, the **Jacobian** is

$$
J = \begin{bmatrix} 0 & I\\ -A^{-1}C & -A^{-1}B\end{bmatrix}
$$

and its four **eigenvalues** $\lambda$ say what a small nudge does: it grows like $e^{\operatorname{Re}\lambda\, t}$ and
oscillates with angular frequency $\operatorname{Im}\lambda$. All real parts negative → stable.

*Note:* Simos's Annex B (eq. B.13) leaves out the $Y_v$, $N_v$ terms in $B$. Without them his Fig. 5.3 is *not*
reproduced (3–6× off); with them it is (D19). For other rest positions (sitting at an angle) the Jacobian is computed
numerically by finite differences; a test checks both methods agree at trailing straight.

## 7. The pitchfork (E03)

The product of the four eigenvalues is $\det J = \det C/\det A$. Sympy gives

$$
\det C = F_{H0}\,U\,(N_v - x_A Y_v)
$$

which changes sign when $x_A = L\,C_6'/C_2'$, i.e.

$$
\boxed{\;\frac{x_{Tcr}}{L} = \frac{\tfrac12 + 2.4\,T/L}{1 + 0.4\,C_B B/T}\;} \qquad (\text{Simos eq. 3.37})
$$

There one real eigenvalue passes through zero and two new rest positions appear (at an angle, to port and starboard):
a **pitchfork** (Strogatz §3.4). Numbers: $x_{Tcr}/L$ = **0.359** (loaded), **0.193** (ballast). It does not depend on
rope length or friction.

**How the angle grows.** Expanding the equation of section 5 for small $\psi$ gives

$$
c_2\,(x_{Tcr}-x_A)\,\psi \;-\; q\,\psi|\psi| = 0, \qquad
c_2 = \tfrac{\pi T}{L}\Big(1+0.4\,C_B\tfrac{B}{T}\Big),\quad q = \tfrac{l_p}{L}C_Y + \tfrac{x_A}{L}\Big(C_Y - \tfrac{\pi T}{2L}\Big)
$$

so $\psi_{rest} = c_2\,\Delta / q$ with $\Delta = (x_{Tcr}-x_A)/L$: the angle grows **in proportion to** $\Delta$.
Strogatz's pitchfork $\dot x = rx - x^3$ grows like $\sqrt r$; ours is $\dot x = rx - x|x|$ because the sideways drag
grows like $v|v|$ (quadratic). Same picture, different exponent (E03 fig06, log-log slope 1.00).

## 8. The Hopf bifurcation (E04)

Code: `fishtail/continuation.py` (`hopf_C_F`, `hopf_curve`).

When the friction $C_F$ is lowered, a **complex pair** of eigenvalues $\lambda = \sigma \pm i\omega$ moves into the right
half-plane: $\sigma$ changes sign while $\omega \neq 0$. That is a **Hopf bifurcation** (Strogatz §8.2): a damped swing
becomes a growing one. We find the crossing by solving

$$
\sigma(C_F) = \max \operatorname{Re}\lambda_{\text{complex}}(C_F) = 0
$$

with Brent's method. Repeating for many rope lengths gives the **stability boundary** (E06 fig01), which matches Simos
Fig. 5.3 to within 0.5 % (loaded) and 5.6 % (ballast). (Simos used the Hurwitz test on the characteristic
polynomial instead; same answer, different route.)

**Pendulum picture.** Trailing straight, the rope tension $F_{H0} \propto C_F$ pulls the swinging ship back, like
gravity on a pendulum. So the swing frequency grows roughly like $\omega \propto \sqrt{C_F}$ and the period like
$1/\sqrt{C_F}$ (E04 fig02b).

## 9. The fishtailing swing: limit cycles (step 6, E05)

Code: `fishtail/cycles.py`.

**Measuring the swing.** Simulate until the swing stops changing, then read the amplitude and period (like a tank test).

**Shooting** (for slow or unstable cycles). Find a start state $\mathbf{x}_0$ and period $T$ such that one period later
the ship is back where it started:

$$
\boldsymbol\phi_T(\mathbf{x}_0) - \mathbf{x}_0 = 0, \qquad \dot\gamma_0 = 0 \ (\text{start at the top of a swing})
$$

and solve with Newton's method.

**Floquet multipliers.** The **monodromy matrix** $\Phi = \partial\boldsymbol\phi_T/\partial\mathbf{x}_0$ says how a small
kick to the cycle changes after one period. Its eigenvalues $\mu$ are the Floquet multipliers. One is always 1 (a kick
*along* the cycle only shifts the timing); if all others have $|\mu| < 1$ the swing is stable. Ours: e.g. 1.000, 0.684,
0.059 at 80 % of the Hopf friction.

**Why the swing grows in proportion to the distance past the Hopf point.** Near the Hopf point the swing is almost a
sine wave $x = r\cos\theta$ with slowly changing size $r$. **Averaging** (Strogatz §7.6) over one swing gives an
equation for $r$. For a toy oscillator $\ddot x - \mu\dot x + D(\dot x) + x = 0$:

$$
\text{quadratic drag } D = \dot x|\dot x|:\quad \dot r = \frac{\mu}{2}r - \frac{4}{3\pi}r^2 \;\Rightarrow\; r^* = \frac{3\pi}{8}\mu
$$

$$
\text{cubic drag } D = \dot x^3:\quad \dot r = \frac{\mu}{2}r - \frac{3}{8}r^3 \;\Rightarrow\; r^* = \sqrt{\frac{4\mu}{3}}
$$

The ship's limiting force is the quadratic cross-flow drag, so its swing grows like $\mu$, not $\sqrt\mu$
(E05 fig02: slope 0.98). Up and down sweeps give the same curve: **supercritical**, no hysteresis.

---

## 10. Are there only three graphs? All the parameters

**The demo page has three tabs because the rigid-rope model has exactly three knobs.** Everything else is either a
fixed property of the ship or belongs to a later step:

| Parameter | What it is | In the rigid model? | Where it is studied |
|---|---|---|---|
| $x_A/L$ | where the rope is tied | **knob** | E03, E06, demo tab "Pitchfork" |
| $L_H/L$ | rope length | **knob** | E04, E06, demo tab "Hawser length" |
| $C_F$ | skin friction | **knob** | E04, E05, E06, demo tab "Fishtailing" |
| current speed $U$ | how fast the water flows | only through $C_F$ (section 4) | E02 |
| loading 100 % / 40 % | full or ballast ship | fixed hull data (two versions) | E03, E06 (demo uses 100 %) |
| $C_Y$, $l_p$ | sideways drag and where it acts | fixed; affect swing size only, not the boundaries | E05 (open issue 16) |
| masses, added masses | from Simos Annex A | fixed | – |
| water density, viscosity | $\rho$, $\nu$ | only through $C_F$ | (PLACEHOLDER, issue 10) |
| rope stiffness, stretch law | how much the rope stretches | **not in the rigid model** | step 7 (below) |
| thruster gains, thrust limit | yaw control | yes, optional (default off) — section 12 | step 8, E08 (to do) |
| hold-back tug pull | steady pull at the stern | yes, optional (= wind with β_W = 0 at the stern) | step 8 |
| wind force and direction | breaks the symmetry | yes, optional (D28) | E10, demo tab "Wind" |

*(Update 2026-11-02: a fourth tab, "Wind", was added for E10; a control tab D5 is planned after E08.)*

**There are many more graphs than the demos.** The experiments have 49 figures (E01: 8, E02: 3, E03: 7, E04: 6,
E05: 6, E06: 5, E07: 6, E10: 8). The demo shows one knob at a time; the **two-knob maps** (E06 fig03, fig04) show two at once, and
those are the main result.

Possible additions to the demo page (tell me which, if any):
- a **100 % / 40 % loading switch** on every tab (cheap: the data just needs computing twice);
- a **map tab**: click a point on the (rope length, friction) stability map and watch that ship;
- later a **stretchy-rope tab** (rope stiffness slider) and a **control tab** (thruster gain slider).

## 11. The stretchy rope (step 7 — now DONE, see E07)

**Update (done):** rope law fitted to a real SPM nylon hawser (Huang et al. 2026, Fig. 3):
$F_H = 10\,580\ \text{kN}\times\text{strain}^{1.86}$ (taut), 0 (slack). Code: `fishtail/model_elastic.py`;
results: `experiments/E07_elastic_hawser/README.md`. The text below is the plan as first written.

**What it is.** So far the rope is a rigid bar: fixed length, and it can push. A real hawser (nylon or polyester)
stretches when pulled and goes **slack** (no force at all) when the ship comes closer to the buoy than its length.

**Why we need it.** E04 and E05 found that at low friction the rigid bar has to *push* for a moment in every swing
($F_H < 0$). A real rope cannot; it goes slack and then snaps tight again (snatch loads). Only a stretchy model can say
what really happens there and how big the peak rope force is.

**The equations (planned).** The ship is now free (no length constraint), so it needs its position too: 6 states,
for example midship position $(X, Y)$, heading $\psi$, and velocities relative to the water $u, v, r$:

$$
\dot X = u\cos\psi - v\sin\psi - U,\qquad \dot Y = u\sin\psi + v\cos\psi,\qquad \dot\psi = r
$$

plus the same Newton equations as in 3c, where the rope force now comes from a **rope law** instead of being an
unknown. With $\ell$ = distance from buoy to attachment point and $\ell_0$ = unstretched length:

$$
F_H = \begin{cases} k\,(\ell - \ell_0)^n & \ell > \ell_0 \quad(\text{taut})\\ 0 & \ell \le \ell_0 \quad(\text{slack})\end{cases}
$$

$n = 1$ is a simple spring; $n > 1$ stiffens as it stretches (typical for synthetic ropes).

**Current speed comes back in.** In $t^*$ units the rope stiffness appears as a ratio "rope stiffness ÷ water force",
for example $\kappa = k\,\ell_0 / (\tfrac12\rho U^2 L T)$ for $n=1$. A faster current makes the same rope behave like a
softer one. So in E07 current speed matters through **both** $C_F$ and $\kappa$.

**Checks planned.** (1) A very stiff rope must give the same rest positions and stability boundaries as the rigid
model. (2) Rope force time series during fishtailing: slack periods and snatch peaks. (3) How the boundaries move as the
rope gets softer. Because "slack/taut" switches the equations on and off, this part is done with simulations only (the
eigenvalue method of section 6 does not apply across the switch).

**Decisions I need from you before starting** (these are modelling choices):
1. **Rope law:** simple spring ($n = 1$) or stiffening ($n > 1$)? I suggest starting with $n = 1$ (easy to explain),
   then one run with $n > 1$.
2. **Rope stiffness value:** needs a source (a real hawser's data sheet or a mooring guideline). Should I search for a
   citable value, or use a clearly labelled PLACEHOLDER range and present results against the stiffness ratio $\kappa$?
3. **Scope:** a full step 7 (about 1–2 days) or a short version (stiff-rope check + one snatch-load figure)?

## 12. Control: a stern thruster and a hold-back tug (step 8 — code done, E08 to do)

**What it is.** Fishtailing is a Hopf instability: the swing grows because the water does not damp it enough. A
controller adds the missing damping. Two options are in the code (D31):

**(a) Stern thruster with PD control.** The thruster wants a yaw moment
$$
N_c^* = -k_p\,\psi - k_d\,r^*,\qquad N_c^* = \frac{N_c}{\tfrac12\rho U^2 L^2 T},\qquad r^* = \frac{rL}{U}.
$$
$k_p$ is a spring (pull the bow back to straight), $k_d$ a damper (resist turning). Written this way the gains are
non-dimensional, so in the rigid model they mean the same thing at every current speed (like $C_F$, section 4).

The thruster can only push sideways, at $x_T = -0.5L$ (the stern, PLACEHOLDER). To make the moment $N_c$ it must push
with $Y_T = N_c / x_T$, so it also pushes the whole ship sideways. Then it is clipped to the thrust limit:
$$
|Y_T| \le \text{thrust\_limit}\times F_{H0},\qquad F_{H0} = C_F\,\tfrac12\rho U^2 L T .
$$
$Y_T$ is added to the sway force and $x_T Y_T$ to the yaw moment in Newton's law (section 3c).

**In the linearisation (section 6).** Near trailing straight the thrust is small, so the limit is never reached and
the law is linear. Moving it to the left of $A\ddot q + B\dot q + Cq = 0$:
$k_p$ adds $K k_p\,(1/x_T,\ 1)$ to the $\psi$ column of $C$ and $k_d$ adds $K k_d (L/U)\,(1/x_T,\ 1)$ to the
$\psi$ column of $B$, with $K = \tfrac12\rho U^2L^2T$. So $k_d$ adds damping (it should push the Hopf boundary) and
$k_p$ adds stiffness (it should push the pitchfork boundary). The limit does not appear at all: **the thrust limit
cannot change whether straight trailing is stable**, only what happens to large swings. This is a classic nonlinear
effect to show in E08: a controller that stabilises small disturbances may not catch a large swing.

**(b) Hold-back tug.** A steady pull $F_w$ downstream, made fast at the stern. Mathematically it is the wind force of
E10 with $\beta_W = 0$ and $x_W = -0.5L$. It keeps the mirror symmetry (straight is still a rest position) and
raises the hawser tension to $F_{H0} + F_w$. When the ship turns by $\psi$, the pull gets a sideways part
$F_w\psi$ at the stern: a weathervane tail that turns the bow back.

**Checks already passing** (`tests/test_control.py`): zero gains give the old model exactly; the sign of the law is
right; the limit is never exceeded, even in a simulation; the linear Jacobian matches finite differences with thruster
and with tug; mirror symmetry holds; straight is still a rest position; the limit does not move the eigenvalues.

**Planned for E08.** Hopf boundary vs $k_d$ (how much damping is needed), stability map with and without control,
thrust needed vs limit, time runs where saturation lets a large swing escape. Values to be shown as ratios
(gains, thrust_limit, $F_w/F_{H0}$) until cited thrust and tug numbers are found (open issue 25).
