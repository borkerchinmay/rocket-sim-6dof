# 🚀 6-DOF Model Rocket Flight Simulation with Airbrake Control
 
A high-fidelity, nonlinear six-degree-of-freedom (6-DOF) rocket flight simulation implemented entirely in **native MATLAB**, no toolboxes required. The model uses quaternion-based attitude kinematics, a data-driven aerodynamic coefficient table, a two-phase ODE integration strategy, and a modular airbrake control architecture designed to hit a precise target apogee altitude.
 
Built as a simulation platform for model rocketry GNC (Guidance, Navigation & Control) research and student-learning use.
 
---
 
## Example Output
 
```
Motor CSV loaded: 7106 rows | burn time = 7.1050 s
  Peak thrust       : 5082.3 N
  Propellant mass   : 9.2408 kg
  Dry mass          : 25.0000 kg
 
=== PHASE 1: Boost  (0 → 7.1050 s) ===
  Done: 444 steps, 0.10 s wall
  Burnout: alt=2192.8 m, speed=315.1 m/s
 
=== PHASE 2: Coast  (7.1050 s → landing) ===
  Done: 367 steps, 0.06 s wall
 
========== SIMULATION SUMMARY ==========
  Apogee           :   5229.2 m  at t = 30.11 s
  Max airspeed     :    450.0 m/s  (Mach 1.337)
  Total flight time:    66.97 s
  Landing downrange:   1667.8 m
=========================================
```
 
---
 
## Key Features
 
| Feature | Detail |
|---|---|
| **6-DOF dynamics** | Full translational + rotational EOM, no simplifying assumptions |
| **Quaternion attitude** | No gimbal lock, singularity-free at any orientation |
| **Data-driven aero** | Mach-interpolated coefficient table loaded from CSV — swap your own data in |
| **Time-varying CG** | CG shifts forward as propellant burns, updated every integration step |
| **Two-phase integration** | Boost and coast solved separately — eliminates thrust-discontinuity stiffness |
| **Wind model** | Constant / logarithmic boundary-layer / Dryden turbulence |
| **Airbrake control** | Modular controller interface: bang-bang, PID, open-loop, LQR template |
| **Actuator model** | Rate limiting, deadzone, and saturation on the airbrake servo |
| **Monte Carlo** | 200-run dispersion with tornado sensitivity chart and landing CEP |
| **ISA atmosphere** | Troposphere + stratosphere density, temperature, speed of sound |
 
---
 
##  File Structure
 
```
rocket_sim/
│
├── main_sim.m              ← Entry point — run this
│
│── Core simulation
│   ├── rocketEOM.m         ← 6-DOF equations of motion (plant)
│   ├── quat2DCM.m          ← Quaternion → Direction Cosine Matrix
│   ├── atmosphere.m        ← ISA atmosphere (density, sound speed)
│   ├── windModel.m         ← Wind profiles (constant / log-law / Dryden)
│   ├── motorInterp.m       ← Fast thrust & mass lookup (binary search)
│   └── aeroCoeffs.m        ← Mach-interpolated aero coefficient table
│
├── Data loading
│   ├── loadRocketData.m    ← Unified CSV loader (motor + aero + CG)
│   └── loadMotorCSV.m      ← Standalone motor-only CSV loader
│
├── Airbrake system
│   ├── airbrake_Aero.m      ← Panel drag physics (flat-plate model)
│   ├── airbrakeActuator.m  ← Servo model (rate limit, saturation)
│   └── airbrakeControl.m ← Swappable controller (bang-bang/PID/LQR)
│
├── Post-processing
│   ├── postProcess.m       ← Extract Euler angles, Mach, alpha, beta...
│   ├── plotResults.m       ← Flight kinematics, attitude, aero plots
│   ├── plotAirbrake.m      ← Airbrake deployment and drag breakdown
│   └── stabilityAnalysis.m ← Static margin vs Mach
│
└── Monte Carlo
    └── monteCarlo.m        ← Dispersion analysis, CEP, tornado chart
```
 
---
 
## Getting Started
 
### Requirements
- MATLAB R2018a or later
- No toolboxes required
### Run the simulation
 
1. Clone the repository
```bash
git clone https://github.com/yourusername/rocket-sim-6dof.git
cd rocket-sim-6dof
```
 
2. Place your motor/aero data CSV in the project folder
3. Open `main_sim.m` and set your rocket's parameters (section 2–4)
4. Run:
```matlab
>> main_sim
```
 
### Disable airbrakes (baseline trajectory)
```matlab
params.ab_enable = false;
```
 
### Change controller
```matlab
params.ab_ctrl.type = 'pid';    % 'bang_bang' | 'pid' | 'open_loop' | 'lqr'
```
 
---
 
## 📐 Physics & Coordinate Frames
 
### Reference frames
 
| Frame | x | y | z |
|---|---|---|---|
| **NED (inertial)** | North | East | Down |
| **Body** | Nose (thrust axis) | Right | Down relative to body |
 
Gravity acts in the `+z` NED direction. All aerodynamic forces are computed in the body frame and rotated to NED via the quaternion-derived Direction Cosine Matrix before integrating Newton's second law.
 
### State vector
 
```
x(1:3)   — NED position [m]
x(4:6)   — NED velocity [m/s]
x(7:10)  — Unit quaternion q = [q₀; q₁; q₂; q₃]  (scalar-first)
x(11:13) — Body angular rates [p; q; r]  [rad/s]
x(14)    — Airbrake deployment δ ∈ [0,1]  (when enabled)
```
 
### Equations of motion
 
**Translational (NED):**
 
$$\dot{\mathbf{v}} = \frac{1}{m} \mathbf{C}_{NB} \mathbf{F}_{body} + \mathbf{g}_{NED}$$
 
**Rotational — Euler's equations (body frame):**
 
$$\dot{{\omega}} = \mathbf{I}^{-1}\left(\mathbf{M}_{aero} - {\omega} \times \mathbf{I}{\omega}\right)$$
 
**Quaternion kinematics:**
 
$$\dot{\mathbf{q}} = \frac{1}{2} \Xi (\mathbf{q})\,{\omega}$$
 
where $\Xi(\mathbf{q})$ is the 4×3 kinematic matrix. The quaternion is renormalised at every integration step to prevent drift.
 
### Aerodynamic model
 
Forces and moments are computed in the body frame using small-angle lift theory and velocity-aligned drag:
 
| Term | Expression |
|---|---|
| Axial drag | $\mathbf{F}_{drag} = -C_D\, q_\infty A_{ref}\, \hat{\mathbf{V}}$ |
| Normal force | $F_N = C_{L\alpha}\, \alpha\, q_\infty A_{ref}$ |
| Side force | $F_S = C_{N\beta}\, \beta\, q_\infty A_{ref}$ |
| Pitch moment | $M_y = -F_N (x_{CP} - x_{CG})$ |
| Pitch damping | $M_y^{damp} = C_{Mq}\, \frac{\dot{\theta} L_{ref}}{2V}\, q_\infty A_{ref} L_{ref}$ |
| Roll damping | $M_x = C_{RM_p}\, \frac{p\, L_{ref}}{2V}\, q_\infty A_{ref} L_{ref}$ |
 
Drag is applied along the **velocity vector** (not body axis), so it is physically correct at any angle of attack including fully tumbled flight.
 
All coefficients are interpolated from a user-supplied Mach table using fast binary search — no `interp1` overhead in the ODE hot path.
 
---
 
##  Airbrake Control System
 
The airbrake system is architected as four independent, loosely-coupled modules:
 
```
State ──> Controller ──> Actuator ──> Aero ──? Force sum
           (swap here)    (rate lim)  (drag)
```
 
### Why this architecture?
 
The controller is a **single swappable function** (`airbrakeController.m`) with a fixed interface. Changing from bang-bang to PID to LQR requires changing **one line** in `main_sim.m` and nothing else. The plant, actuator, and aero models are completely untouched.
 
### Available controllers
 
| Controller | `params.ab_ctrl.type` | Description |
|---|---|---|
| Bang-bang | `'bang_bang'` | On/off based on predicted apogee vs target ± deadband |
| PID | `'pid'` | Proportional-integral-derivative on apogee error |
| Open-loop | `'open_loop'` | Fixed deployment at fixed time — for aero model validation |
| LQR template | `'lqr'` | State-feedback; provide your own gain matrix `K` |
 
### Apogee prediction
 
All closed-loop controllers use an energy-method apogee prediction:
 
$$h_{apo}^{pred} = h + \eta \frac{V^2}{2g}$$
 
where $\eta \approx 0.90$ is an empirical drag-loss factor. This is computed analytically at every step — no forward ODE integration required — keeping the controller computationally free.
 
### Airbrake drag model
 
Panels are modelled as flat plates perpendicular to flow:
 
$$\Delta C_D = C_{D,plate}(M) \cdot \frac{N_{panels} \cdot A_{panel}}{A_{ref}} \cdot \delta^{n}$$
 
where $\delta \in [0,1]$ is deployment fraction, $n$ is a geometry exponent (1 = linear area), and $C_{D,plate}(M)$ is interpolated from a Mach-dependent flat-plate table.
 
---
 
## Monte Carlo Dispersion Analysis
 
`monteCarlo.m` runs N independent trajectories with the following perturbed inputs:
 
| Input | 1-σ | Distribution |
|---|---|---|
| Launch elevation angle | 0.5° | Normal |
| Launch azimuth angle | 2.0° | Normal |
| Wind speed | 2.0 m/s | Normal |
| Wind direction | 15° | Normal |
| Dry mass | 0.3 kg | Normal |
| Drag coefficient offset | 0.04 | Normal |
| Thrust scale factor | 3% | Normal |
| CG location | 15 mm | Normal |
| Pitch/yaw inertia | 0.05 kg·m² | Normal |
 
**Outputs:** landing scatter with 2σ ellipse, CEP (50% circular error probable), apogee histogram, flight-time histogram, and a tornado sensitivity chart showing which inputs drive the most output variance.
 
---
 
##  Performance
 
Two design decisions make the simulation fast enough for Monte Carlo:
 
**1. Two-phase integration splits at burnout**
The thrust profile ends with a step discontinuity at burnout. Integrating through it forces `ode45`/`ode23` to shrink the step to near-machine-epsilon trying to resolve it. By stopping Phase 1 exactly at `t_burnout` and restarting Phase 2 from the burnout state, the discontinuity is never straddled. Both phases integrate cleanly.
 
**2. No `interp1` in the ODE hot path**
MATLAB's `interp1` carries significant per-call overhead from argument parsing and validation. All table lookups (motor, aero, atmosphere) use manual binary search with direct array indexing — roughly 20× faster per lookup.
 
Typical wall times on a standard laptop:
 
| Phase | Steps | Wall time |
|---|---|---|
| Boost (7 s burn) | ~450 | 0.10 s |
| Coast (60 s descent) | ~370 | 0.06 s |
| Full Monte Carlo (200 runs) | — | ~40 s |
 
---
 
## 📋 CSV Data Format
 
### `loadRocketData.m` — unified input file
 
```
col  1  time [s]
col  2  altitude [m]
col  3  total mass [grams]
col  4  motor mass [grams]
col  5  longitudinal moment
col  6  rotational moment
col  7  CP location from nose [cm]
col  8  CG location from nose [cm]
col  9  static margin [calibers]
col 10  Mach number [-]
col 11  thrust [N]
col 12  CD [-]
col 13  CN normal [-]
col 14  CM pitching moment [-]
col 15  CN beta [-]
col 16  CF side [-]
col 17  (unused)
col 18  CRM damping [-]
```
 
No header row. NaN values are filled automatically.
 
---
 
## 🔬 Numerical Methods
 
| Choice | Rationale |
|---|---|
| **`ode23` for boost** | Low-order Runge-Kutta handles the fast thrust dynamics and is less aggressive at shrinking step size near sharp features than `ode45` |
| **`ode113` for coast** | Variable-order Adams predictor-corrector — far more efficient than Runge-Kutta for smooth, slow-varying coast dynamics; uses ~3× fewer function evaluations |
| **Dynamic pressure threshold** | Below $q_\infty < 5$ Pa (≈10 m/s), aerodynamic forces are negligible. Setting them to zero near apogee eliminates the $1/V$ singularity in rate-damping terms, which would otherwise drive the ODE stiff |
| **Quaternion renormalisation** | Performed at every EOM call. Prevents slow numerical drift of $\|\mathbf{q}\| \neq 1$ that would corrupt the rotation matrix |
| **Per-state `AbsTol`** | Coast phase uses looser tolerances on position (0.01 m) and tighter on quaternion (1e-6), avoiding wasted effort on states that vary slowly |
 
---
 
## References
 
- Zipfel, P.H. (2007). *Modeling and Simulation of Aerospace Vehicle Dynamics*. AIAA Education Series.
- Barrowman, J.S. (1967). *The Practical Calculation of the Aerodynamic Characteristics of Slender Finned Vehicles*. Master's thesis, Catholic University of America.
- MIL-HDBK-1797 — Flying Qualities of Piloted Aircraft (Dryden turbulence model)
- International Standard Atmosphere — ISO 2533:1975
- Shepperd, S.W. (1978). Quaternion from rotation matrix. *Journal of Guidance and Control*, 1(3).

---
 
## License
 
See LICENSE for usage restrictions.