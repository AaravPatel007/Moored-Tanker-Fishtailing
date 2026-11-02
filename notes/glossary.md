# Glossary

Plain-language meanings of every technical term used in the project. Add a term the first time it is used.

## The ship and the mooring
- **Single-point mooring (SPM):** the ship is tied at one point only (a buoy or a turret), so it can swing freely to
  face the weather.
- **Hawser:** the rope (or bar) between the buoy and the ship's bow.
- **Rigid hawser:** an idealised hawser that never stretches or goes slack — it keeps a fixed length. Easier to analyse.
- **Attachment point / fairlead (`x_A`):** where the hawser is fixed on the ship, measured forward from midship.
- **Turret:** a mooring point built *inside* the hull; moving it aft is like moving `x_A` aft.
- **Surge / sway / yaw:** moving forward-back / moving sideways / turning about a vertical axis.
- **Heading `ψ`:** which way the bow points, measured from the current direction.
- **Hawser angle `γ`:** angle between the hawser and the current direction.
- **Midship:** the middle of the ship's length; our body axes start here.
- **L, B, T:** ship length, beam (width), draft (depth below the waterline).
- **Block coefficient `C_B`:** how "boxy" the hull is: underwater volume ÷ (L·B·T). Tankers ≈ 0.8.
- **Wetted surface `S`:** area of hull skin in contact with water.
- **VLCC (Very Large Crude Carrier):** a very large crude-oil tanker. Both the Vidal de Negreiros and KVLCC2 are
  VLCCs, and both are 320 m long (sources.md).
- **Loading condition (100 % / 40 %):** fully loaded (deep draft) vs ballast (shallow draft). Same hull, very different
  water forces.
- **Pivot / turret position `x_T`:** where a hull is held so it can only turn. `x_Tcr` = the critical position: held
  further forward than this, the hull points into the current; further aft, it turns away (eq. 3.37).
- **Model scale 1:λ, Froude scaling:** the tank model is λ times smaller (λ = 90 here). Speeds and times scale with √λ,
  masses with λ³. Friction does *not* scale this way — that is the scale effect.

## The water forces
- **Current:** steady flow of water past the ship.
- **Relative velocity:** the ship's velocity *as seen by the water*. Forces depend on this, not on speed over ground.
- **Angle of attack `α`:** angle between the current and the ship's centreline; 180° = current from dead ahead.
- **Added mass:** when the hull accelerates it must also push water out of the way, so it feels heavier than it is.
  Different in each direction (`M11` surge, `M22` sway, `M66` yaw, `M26` sway–yaw coupling).
- **Munk moment:** a turning moment that pushes a long body to turn *broadside* to the flow. Exists even in a perfect
  fluid. It is what makes the ship want to sheer away — the destabilising effect in fishtailing.
- **Lift (wing effect):** at small angles the hull acts like a short stubby wing and makes a sideways force.
- **Cross-flow drag:** the drag of water flowing *across* the hull, slice by slice, like flow past a cylinder.
- **`C_Y`:** sideways drag coefficient (current hitting the ship exactly side-on).
- **`l_p`:** distance from midship to where the sideways drag acts.
- **Friction coefficient `C_F` (or `C_0`):** how much skin friction the hull feels when the water flows along it.
- **Form factor `k`:** extra resistance due to the hull's 3-D shape, on top of flat-plate friction.
- **Reynolds number `Re = U L / ν`:** ratio of "inertia" to "stickiness" of the flow. Bigger ship or faster flow →
  bigger `Re` → relatively less friction.
- **Froude number `U/√(gL)`:** measures how much the ship makes waves. Small at mooring speeds.
- **Scale effect:** a model test does not behave exactly like the real ship because `Re` cannot be matched.
- **Heuristic model:** a model built from physical reasoning plus a few measured numbers.
- **Captive model test:** a towing-tank test where the model is held at a fixed angle and the forces are measured.
- **MMG model / manoeuvring derivatives:** a standard Japanese ship-manoeuvring model; forces written as polynomials
  in sway speed and yaw rate. Accurate for small angles only.
- **Linear derivatives `Y_v, Y_r, N_v, N_r` (Clarke et al. 1982):** how much sideways force / yaw moment a small sway
  speed or yaw rate produces. Clarke's formulas estimate them from L, B, T, C_B alone. `N_r < 0` = yaw damping.
- **Yaw-rotating test:** the model is towed down the tank while being spun about a vertical axis at a steady rate,
  and the forces are measured. Tests the yaw-rate part of the force model.
- **Ideal-fluid (potential) force:** the part of the water force that would exist even without viscosity; it comes
  from the added masses (e.g. the Munk moment, or `M22 u r` when turning).
- **Sectional drag coefficient `C_D(x)`:** sideways drag coefficient of one thin slice of the hull at position x.
  Adding up all slices gives `C_Y`.
- **Strip (slice) integral:** adding up the forces on all the slices along the hull, each with its own local speed
  `v + r x`.
- **Vortex shedding (bow/stern vortices):** swirls of water that peel off the ends of a turning hull; they add yaw damping.
- **ITTC-1957 friction line:** standard formula 0.075/(log10 Re − 2)² for the skin friction of a smooth flat plate.
- **Static yaw stiffness:** how strongly the yaw moment changes with heading. Restoring (pushes back) = stable;
  negative = the hull turns further away.

## The mathematics (Strogatz)
- **Non-dimensionalisation:** rewriting equations in the problem's natural units so only pure numbers remain.
- **State / state space:** the list of numbers that fully describe the system at one instant (here `γ, γ̇, ψ, ψ̇`).
- **Equilibrium (fixed point):** a state where nothing changes (e.g. ship trailing straight behind the buoy).
- **Stable / unstable:** after a small nudge, does the system come back (stable) or move away (unstable)?
- **Linearisation:** approximating the equations near an equilibrium by straight-line (linear) ones.
- **Jacobian:** the matrix of the linearised equations — how each rate of change responds to each state.
- **Eigenvalue:** a number from the Jacobian that tells how a small disturbance grows or decays. Negative real part →
  decays (stable); positive real part → grows (unstable); imaginary part → it oscillates.
- **Bifurcation:** a sudden change in the *kind* of behaviour as a parameter is slowly changed.
- **Bifurcation parameter (knob):** the parameter we slowly change.
- **Pitchfork bifurcation:** one equilibrium becomes unstable and two new mirror-image ones appear (ship sits at an
  angle, left or right). *Supercritical*: new ones are stable. *Subcritical*: new ones are unstable (sudden jump).
- **Imperfect pitchfork:** a small asymmetry (e.g. wind from one side) makes one side preferred.
- **Hopf bifurcation:** a pair of complex eigenvalues crosses into the right half-plane; a steady state turns into a
  steady oscillation. *Supercritical*: the oscillation grows smoothly from zero. *Subcritical*: it jumps to a large size.
- **Limit cycle:** an isolated, self-sustained periodic motion (fishtailing). Nearby motions are drawn to it (stable)
  or pushed away (unstable).
- **Hysteresis:** the system's state depends on the direction the parameter was changed (up vs down).
- **Cycle fold (saddle-node of cycles):** a stable and an unstable limit cycle meet and disappear together.
- **Phase portrait:** a plot of one state against another (e.g. `ψ` vs `ψ̇`) showing how motions move.
- **Continuation:** following an equilibrium (or cycle) as the parameter changes, step by step.
- **Shooting method:** finding a periodic orbit by guessing a starting point and period, simulating one period, and
  correcting the guess until the orbit closes on itself.
- **Floquet multipliers:** numbers that tell whether a limit cycle is stable (all inside the unit circle, apart from
  the one that is always exactly 1).
- **Critical slowing down:** near a bifurcation, disturbances decay very slowly, so simulations need to run long.
- **Saddle point:** a fixed point with one growing and one decaying direction (one positive, one negative real
  eigenvalue). Almost every start leaves it. The bare hull facing the current is a saddle (Strogatz §5.2, §6.4).
- **Node:** a fixed point whose eigenvalues are real and of the same sign; stable node = motions creep in without
  oscillating.
- **Hurwitz criterion:** a test on the coefficients of a polynomial that tells whether all its roots have negative
  real part — a shortcut to stability without computing eigenvalues.
- **Focus (spiral):** a fixed point whose leading eigenvalues are a complex pair. Stable focus = motions spiral in
  (a damped swing); unstable focus = they spiral out.
- **Growth rate:** the real part of an eigenvalue. A small disturbance grows (or shrinks) like e^(growth rate × t*).
- **Mass, damping and stiffness matrices (A, B, C):** the linearised equations A q̈ + B q̇ + C q = 0 for q = (γ, ψ).
  det C changing sign means a real eigenvalue passes through zero (pitchfork).
- **Normal form:** the simplest equation that shows a bifurcation's behaviour near the bifurcation point, e.g.
  ẋ = rx − x³ for a pitchfork. Its last term sets how fast the new branches grow (here √r).
- **Non-smooth normal form:** in our model the restoring term is the quadratic drag x|x| instead of x³, so the normal
  form is ẋ = rx − x|x| and the branches grow like r, not √r. Same for the Hopf amplitude.
- **Averaging:** a method (Strogatz §7.6) for weakly nonlinear oscillators: average the slow change of the amplitude
  over one cycle. It gives the size of the limit cycle near a Hopf point.
- **Projection (of a phase portrait):** our state has 4 numbers, so a 2-D phase portrait shows only two of them. In a
  projection, paths may cross each other; in the full 4-D space they never do.
- **Monodromy matrix:** how a small change of the starting point changes the state after exactly one period. Its
  eigenvalues are the Floquet multipliers.
- **Phase condition:** the extra equation the shooting method needs to pick *which* point of the cycle to start from
  (we start where the hawser angle is at its maximum, dγ/dt* = 0).
- **Constrained motion (rigid hawser):** the hawser keeps a fixed length, so the hawser force is whatever is needed
  to keep it so. We solve for it together with the accelerations at every instant.
- **Regime map (stability map):** a chart over two knobs showing where the ship trails straight, sits at an angle,
  or fishtails. The lines between regions are bifurcation curves.
- **Weathervane effect:** the sideways water force acts behind the attachment point and swings the ship back into
  line, like the tail of a weathervane. It loses to the Munk moment when the attachment is too far aft (x_A < x_Tcr).
- **Fold (saddle-node bifurcation):** a point where a stable and an unstable rest position meet and disappear
  together. On a bifurcation diagram the branch turns back on itself there (◇ in our figures). Strogatz §3.1.
- **Cusp:** the corner where two fold lines meet in a plane of two knobs. Inside the cusp two stable rest positions
  exist (one per side); outside only one. In E10 the cusp tip is the old pitchfork point. Strogatz §3.6.
- **Bistable:** two different stable rest positions exist for the same knob values; which one the ship is in
  depends on its history (hysteresis).
- **Wind force ratio F_w/F_H0:** wind force divided by the hawser tension when the ship trails straight without wind.
- **Snatch load:** a sudden peak of hawser tension, e.g. when a slack rope comes tight again during fishtailing.
- **Digitising a figure:** reading numbers off a published plot. Done here by locating the axes on a high-resolution
  image and tracing the curve; accuracy is stated next to the numbers.
- **Bifurcation diagram:** a plot with the knob across and a measure of the ship's state up (rest heading ψ, or hawser
  angle γ). Each rest position is a line (solid = stable, dashed = unstable). Past a Hopf point the diagram also shows
  the limit cycle, drawn as its highest and lowest value (the edges of the swing; blue in the demo page).
- **Real / imaginary part of an eigenvalue, complex pair:** an eigenvalue λ = σ ± iω. σ (real part) = growth rate;
  ω (imaginary part) = how fast the disturbance swings, period 2π/ω. A real eigenvalue (ω = 0) means no swinging.
  Complex eigenvalues always come as a mirror-image pair (+ω and −ω). How the eigenvalue crosses σ = 0 decides
  the bifurcation: a real one through 0 → pitchfork (or fold); a complex pair → Hopf.
- **Right / left half-plane (RHP / LHP):** the two halves of the eigenvalue plot, Re λ > 0 / Re λ < 0. Any eigenvalue
  in the RHP makes a rest position unstable; "crossing into the RHP" is how a bifurcation shows up in the eigenvalues.
## Control (step 8)
- **Controller:** a rule that measures the ship's motion and decides how hard a thruster should push to correct it.
- **Stern thruster:** a propeller in a sideways tunnel near the stern that pushes the stern to port or starboard. In
  the code its position is `x_T_L` (−0.5 = the stern). *Careful:* the earlier symbol `x_T` / `x_Tcr` (turret / pivot
  position, above) is a different thing; only `x_T_L` in `control.py` and `model_rigid.py` means the thruster.
- **PD control (proportional + derivative):** the thruster moment is N_c* = −k_p ψ − k_d r*. The **P** part acts like
  a spring pulling the bow back to straight; the **D** part acts like a damper resisting turning. Damping is what a
  Hopf (fishtailing) instability lacks, so k_d is the gain expected to matter most.
- **Gain (k_p, k_d):** how strongly the controller reacts. Here in units of ½ρU²L²T (per radian, per unit r*).
- **Saturation / thrust limit:** a real thruster has a maximum force. If the rule asks for more, it gets only the
  maximum ("clipped"). Here the limit is a fraction of F_H0. Saturation is a nonlinearity: it does not change
  stability at trailing straight (thrust there is zero), only how large swings behave.
- **F_H0:** the hawser pull when the ship trails straight, = friction drag C_F · ½ρU²LT. Used as the scale for wind,
  tug and thruster forces.
- **Hold-back tug:** a tug boat behind the moored ship pulling it steadily downstream from the stern. It keeps the
  hawser taut and acts like a weathervane tail. In the code it is the wind force pointed along the current
  (β_W = 0) and applied at the stern.
- **Feedback:** the control force depends on what the ship is doing right now (ψ, r). The tug's pull is not feedback
  (it is constant), although its sideways part grows with ψ by geometry, which acts like a free P term.
- **Decay rate σ:** σ = −max Re λ of a rest position [U/L]. σ > 0: a small swing dies out like e^(−σ t*); σ < 0:
  it grows (fishtailing). "Margin" σ_min: the decay rate we require before calling the controller good enough.
- **Root locus:** the paths the eigenvalues trace in the complex plane as one gain is turned up (E08 fig02).
- **Overdamping:** too much damping: the swinging eigenvalue pair meets on the real axis and splits, and one of the
  two moves back towards 0, so a disturbance dies out more slowly again (like a door closer set too stiff).
  Strogatz §5 (damped oscillator).
- **Least-thrust gains (★):** the pair (k_p, k_d) that meets the margin σ ≥ σ_min with the smallest peak thruster
  force while it stops the full fishtail (E08; the student's definition of "optimal").
- **Bistability (control):** with a thrust limit, the same gains can end in two different states: a small swing dies
  out (the straight position is stable), but the full fishtail survives because the thruster is at its limit. The
  eigenvalues only see small swings, so they cannot show this; simulations do. Strogatz Ch. 8 (a stable rest point
  and a stable cycle side by side; in his planar pictures an unstable cycle separates them, not yet checked here).
- **Basin of attraction:** the set of starting states that end at a given stable state. Saturation shrinks the basin
  of the straight position: the smaller the limit, the smaller the swing it can still stop.

- **Fast mode:** in our rigid-hawser model each rest position has one real eigenvalue far to the left (about −1.7 to
  −4.3 U/L). That disturbance dies out within half a t* and never decides stability, so the demo lists it but does not
  plot it.
