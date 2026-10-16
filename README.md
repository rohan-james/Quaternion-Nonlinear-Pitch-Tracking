# Quaternion Attitude Dynamics and Nonlinear Control

A Python implementation of rigid-body quaternion kinematics, nonlinear PD/PID attitude regulation and tracking, and a high-fidelity missile pitch channel simulation under aerodynamic disturbances and actuator constraints.

## Overview

This repository contains simulation code and mathematical derivations for advanced attitude determination and control systems (ADCS). The project bridges rigorous Lyapunov stability theory with practical engineering constraints such as actuator saturation, first-order lag dynamics, and aerodynamic moment uncertainties.

---

## Key Features

- **Rigorous Attitude Representation:** Uses scalar-first unit quaternions ($q = [q_0, q_v^T]^T$) with Baumgarte-type norm-constraint stabilization to prevent numerical drift. Also implements Direction Cosine Matrices (DCM), Euler angles, and Modified Rodrigues Parameters (MRP).
- **Nonlinear PD/PID Attitude Control:** 
  - Shortest-path selection via sign factor ($\text{sgn}(q_{e0})$) to prevent unwinding.
  - Lyapunov-stable tracking and regulation controllers.
  - Discrete-time integral action with conditional integration anti-windup clamping.
- **Smooth Trajectory Generation:** $C^2$-smooth quintic polynomial ramps for attitude commands, eliminating torque spikes and impulsive angular accelerations.
- **Aerospace Plant Models:**
  - Rigid body with asymmetric inertia tensor.
  - Full missile pitch channel model including US Standard Atmosphere (density, Mach number, dynamic pressure), aerodynamic moments ($C_{m\alpha}, C_{m\delta}$), and fin authority limits.
  - First-order actuator lag and torque magnitude saturation.

---

## Simulation Scenarios

1. **Scenario 1: Rest-to-Rest Regulation**
   - Regulates a rigid body from a $70^\circ$ initial rotation to identity attitude.
   - Confirms monotonic decay of the Lyapunov energy function $V(t)$ on a logarithmic scale.
2. **Scenario 2: Smooth Trajectory Tracking**
   - Tracks a $40^\circ$ quintic ramp about the body z-axis from rest using feedforward-plus-PD control.
3. **Scenario 3: Missile Pitch Tracking with Aero Disturbances**
   - Simulates a missile flying at Mach 2.5 at an altitude of $10,000\text{ m}$.
   - Evaluates performance under a 20% model mismatch ($\kappa = 0.8$) on aerodynamic coefficients, demonstrating the efficacy of integral action in rejecting steady-state errors.

---

## Getting Started

### Prerequisites
- Python 3.8+ 
- Required Python libraries: `numpy`, `scipy`, `matplotlib`

### Running Simulations
```bash
# Clone the repository
git clone https://github.com/your-username/quaternion-attitude-control.git
cd quaternion-attitude-control

# Run Demo
python3 demo.py
```

---
