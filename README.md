# Attitude Representation with Unit Quaternions

## Axis--Angle to Quaternion

The attitude is stored as a scalar-first unit quaternion
$\boldsymbol{q} = [q_0,\ q_1,\ q_2,\ q_3]^T = [q_0,\ \boldsymbol{q}_v^T]^T$.
For a rotation of angle $\phi$ about an axis $\boldsymbol{e}$
(normalised to unit length,
$\hat{\boldsymbol{e}} = \boldsymbol{e}/\lVert\boldsymbol{e}\rVert$):
$$\boldsymbol{q} = \begin{bmatrix} \cos(\phi/2) \\ \hat{\boldsymbol{e}}\,\sin(\phi/2) \end{bmatrix},
\qquad
\boldsymbol{q}^{*} = \begin{bmatrix} \cos(\phi/2) \\ -\hat{\boldsymbol{e}}\,\sin(\phi/2) \end{bmatrix}$$
where $\boldsymbol{q}^{*}$ is the conjugate (the inverse rotation for a
unit quaternion). A zero-length axis
($\lVert\boldsymbol{e}\rVert < 10^{-12}$) is rejected.

## Quaternion (Hamilton) Product

For $\boldsymbol{p} = [p_0,\ \boldsymbol{p}_v^T]^T$ and
$\boldsymbol{q} = [q_0,\ \boldsymbol{q}_v^T]^T$:
$$\boldsymbol{p} \otimes \boldsymbol{q} =
\begin{bmatrix}
p_0 q_0 - \boldsymbol{p}_v \cdot \boldsymbol{q}_v\\
p_0\,\boldsymbol{q}_v+ q_0\,\boldsymbol{p}_v + \boldsymbol{p}_v \times \boldsymbol{q}_v
\end{bmatrix}$$

## Rotation of a Point

A point $\boldsymbol{r} \in \mathbb{R}^3$ is embedded as the pure
quaternion $\boldsymbol{p} = [0,\ \boldsymbol{r}^T]^T$ and rotated by
the sandwich product; the rotated point is the vector part of the
result:
$$\boldsymbol{p}' = \boldsymbol{q} \otimes \boldsymbol{p} \otimes \boldsymbol{q}^{*}$$

## Error Quaternion

With current attitude $\boldsymbol{q}$ and desired attitude
$\boldsymbol{q}_d$ (both normalised first), the attitude error is the
rotation from the desired frame to the body frame,
$\boldsymbol{q}_e = \boldsymbol{q}_d^{*} \otimes \boldsymbol{q}$:
$$q_{e0} = q_0\, q_{d0} + \boldsymbol{q}_v\cdot \boldsymbol{q}_{dv},
\qquad
\boldsymbol{q}_{ev} = q_{d0}\,\boldsymbol{q}_v- q_0\,\boldsymbol{q}_{dv} + \boldsymbol{q}_v\times \boldsymbol{q}_{dv}$$
Perfect alignment corresponds to
$\boldsymbol{q}_e = [\pm 1, 0, 0, 0]^T$. The scalar attitude-error angle
used in all plots is
$$\theta_{\text{err}} = 2 \arccos\!\left( \lvert q_{e0} \rvert \right)$$
The absolute value accounts for the double cover ($\boldsymbol{q}$ and
$-\boldsymbol{q}$ are the same attitude), so
$\theta_{\text{err}} \in [0, \pi]$.

# Secondary Attitude Representations

## Direction Cosine Matrix

The direction cosine matrix (inertial $\to$ body) built from the unit
quaternion is $$C(\boldsymbol{q}) =
\begin{bmatrix}
q_0^2 + q_1^2 - q_2^2 - q_3^2 & 2(q_1 q_2 + q_0 q_3) & 2(q_1 q_3 - q_0 q_2) \\
2(q_1 q_2 - q_0 q_3) & q_0^2 - q_1^2 + q_2^2 - q_3^2 & 2(q_2 q_3 + q_0 q_1) \\
2(q_1 q_3 + q_0 q_2) & 2(q_2 q_3 - q_0 q_1) & q_0^2 - q_1^2 - q_2^2 + q_3^2
\end{bmatrix}$$

## Euler Angles (Yaw--Pitch--Roll, 3-2-1 Sequence)

$$\begin{aligned}
\theta &= -\arcsin\!\big( 2 (q_1 q_3 - q_0 q_2) \big) && \text{(pitch, argument clipped to } [-1, 1]\text{)} \\
\phi &= \operatorname{atan2}\!\big( 2 (q_2 q_3 + q_0 q_1),\ q_0^2 - q_1^2 - q_2^2 + q_3^2 \big) && \text{(roll)} \\
\psi &= \operatorname{atan2}\!\big( 2 (q_1 q_2 + q_0 q_3),\ q_0^2 + q_1^2 - q_2^2 - q_3^2 \big) && \text{(yaw)}
\end{aligned}$$

## Modified Rodrigues Parameters

The MRP vector is
$\boldsymbol{\sigma} = \hat{\boldsymbol{e}} \tan(\phi/4)$. To stay away
from the singularity at $\phi = 360^\circ$, the shadow set is selected
whenever $q_0 < 0$: $$\boldsymbol{\sigma} =
\begin{cases}
\dfrac{\boldsymbol{q}_v}{1 + q_0}, & q_0 \ge 0 \\[2ex]
\dfrac{-\boldsymbol{q}_v}{1 - q_0}, & q_0 < 0
\end{cases}$$ which guarantees
$\lVert \boldsymbol{\sigma} \rVert \le 1$.

# Rotational Equations of Motion

## Quaternion Kinematic Differential Equation

With body angular velocity
$\boldsymbol{\omega}= [\omega_x, \omega_y, \omega_z]^T$, the kinematics
$\dot{\boldsymbol{q}} = \tfrac{1}{2}\, \boldsymbol{q} \otimes [0,\ \boldsymbol{\omega}^T]^T$
are written in matrix form and augmented with a norm-constraint
stabilisation term:
$$\dot{\boldsymbol{q}} = \frac{1}{2}\, \Omega(\boldsymbol{\omega})\, \boldsymbol{q} \;+\; c \left( 1 - \boldsymbol{q}^T \boldsymbol{q} \right) \boldsymbol{q},
\qquad
\Omega(\boldsymbol{\omega}) =
\begin{bmatrix}
0 & -\omega_x & -\omega_y & -\omega_z \\
\omega_x & 0 & \omega_z & -\omega_y \\
\omega_y & -\omega_z & 0 & \omega_x \\
\omega_z & \omega_y & -\omega_x & 0
\end{bmatrix}$$

-   $\tfrac{1}{2}\Omega(\boldsymbol{\omega})\boldsymbol{q}$ ---
    **Kinematic Term**: Since $\Omega$ is skew-symmetric, this term
    alone preserves $\lVert \boldsymbol{q} \rVert$ exactly in continuous
    time.

-   $c\,(1 - \boldsymbol{q}^T\boldsymbol{q})\,\boldsymbol{q}$ ---
    **Constraint (Baumgarte-type) Term**: Vanishes on the unit sphere
    and pulls the numerical solution back towards
    $\lVert \boldsymbol{q} \rVert = 1$ when integration error causes
    drift. The gain is $c = 1$ in all simulations.

## Euler's Rigid-Body Equation

For inertia matrix $I$ and total body-frame torque $\boldsymbol{\tau}$:
$$I \dot{\boldsymbol{\omega}} = \boldsymbol{\tau} - \boldsymbol{\omega}\times (I \boldsymbol{\omega})
\quad\Longrightarrow\quad
\dot{\boldsymbol{\omega}} = I^{-1} \big( \boldsymbol{\tau} - \boldsymbol{\omega}\times (I \boldsymbol{\omega}) \big)$$
The term $\boldsymbol{\omega}\times (I\boldsymbol{\omega})$ is the
gyroscopic (cross-coupling) torque. The linear solve is performed either
with a pre-computed inverse $I^{-1}$ or, optionally, with a pre-computed
LU factorisation $P I = L U$ followed by forward/back substitution.

## Coupled Nonlinear State-Space Form

Stacking attitude and rate gives the 7-dimensional state
$\boldsymbol{x} = [\boldsymbol{q}^T,\ \boldsymbol{\omega}^T]^T \in \mathbb{R}^7$:
$$\dot{\boldsymbol{x}} = f(\boldsymbol{x}, \boldsymbol{\tau}) =
\begin{bmatrix}
\tfrac{1}{2}\, \Omega(\boldsymbol{\omega})\, \boldsymbol{q} + c \left( 1 - \boldsymbol{q}^T \boldsymbol{q} \right) \boldsymbol{q} \\[1ex]
I^{-1} \big( \boldsymbol{\tau} - \boldsymbol{\omega}\times (I \boldsymbol{\omega}) \big)
\end{bmatrix}$$

## Conservation-Law Checks

For verifying the integrator in torque-free motion:
$$T_{\text{rot}} = \frac{1}{2}\, \boldsymbol{\omega}^T I\, \boldsymbol{\omega},
\qquad
\boldsymbol{H}_B = I \boldsymbol{\omega},
\qquad
\boldsymbol{H}_I = C(\boldsymbol{q})^T\, \boldsymbol{H}_B$$ With
$\boldsymbol{\tau} = \boldsymbol{0}$, the rotational kinetic energy
$T_{\text{rot}}$ and the inertial angular momentum $\boldsymbol{H}_I$
must remain constant (the body-frame components $\boldsymbol{H}_B$
generally do not).

# Numerical Integration: Classical Fourth-Order Runge--Kutta

## Kinematics-Only Integrator

When $\boldsymbol{\omega}(t)$ is a prescribed function of time, only
$\boldsymbol{q}$ is integrated, with step size $h$ and
$g(\boldsymbol{q}, \boldsymbol{\omega})$ denoting the right-hand side of
the kinematic equation: $$\begin{aligned}
\boldsymbol{k}_1 &= h\, g\big(\boldsymbol{q}_n,\ \boldsymbol{\omega}(t_n)\big) \\
\boldsymbol{k}_2 &= h\, g\big(\boldsymbol{q}_n + \tfrac{1}{2}\boldsymbol{k}_1,\ \boldsymbol{\omega}(t_n + \tfrac{h}{2})\big) \\
\boldsymbol{k}_3 &= h\, g\big(\boldsymbol{q}_n + \tfrac{1}{2}\boldsymbol{k}_2,\ \boldsymbol{\omega}(t_n + \tfrac{h}{2})\big) \\
\boldsymbol{k}_4 &= h\, g\big(\boldsymbol{q}_n + \boldsymbol{k}_3,\ \boldsymbol{\omega}(t_n + h)\big) \\
\tilde{\boldsymbol{q}}_{n+1} &= \boldsymbol{q}_n + \tfrac{1}{6} \left( \boldsymbol{k}_1 + 2\boldsymbol{k}_2 + 2\boldsymbol{k}_3 + \boldsymbol{k}_4 \right),
\qquad
\boldsymbol{q}_{n+1} = \frac{\tilde{\boldsymbol{q}}_{n+1}}{\lVert \tilde{\boldsymbol{q}}_{n+1} \rVert}
\end{aligned}$$

## Coupled-State Integrator

For the full state $\boldsymbol{x}$ with torque law
$\boldsymbol{\tau}(t, \boldsymbol{x})$: $$\begin{aligned}
\boldsymbol{k}_1 &= f\big(\boldsymbol{x}_n,\ \boldsymbol{\tau}(t_n, \boldsymbol{x}_n)\big) \\
\boldsymbol{k}_2 &= f\big(\boldsymbol{x}_n + \tfrac{h}{2}\boldsymbol{k}_1,\ \boldsymbol{\tau}(t_n + \tfrac{h}{2},\ \boldsymbol{x}_n + \tfrac{h}{2}\boldsymbol{k}_1)\big) \\
\boldsymbol{k}_3 &= f\big(\boldsymbol{x}_n + \tfrac{h}{2}\boldsymbol{k}_2,\ \boldsymbol{\tau}(t_n + \tfrac{h}{2},\ \boldsymbol{x}_n + \tfrac{h}{2}\boldsymbol{k}_2)\big) \\
\boldsymbol{k}_4 &= f\big(\boldsymbol{x}_n + h\,\boldsymbol{k}_3,\ \boldsymbol{\tau}(t_n + h,\ \boldsymbol{x}_n + h\,\boldsymbol{k}_3)\big) \\
\boldsymbol{x}_{n+1} &= \boldsymbol{x}_n + \frac{h}{6} \left( \boldsymbol{k}_1 + 2\boldsymbol{k}_2 + 2\boldsymbol{k}_3 + \boldsymbol{k}_4 \right)
\end{aligned}$$ followed by re-normalisation of the quaternion part
only,
$\boldsymbol{q}_{n+1} \leftarrow \boldsymbol{q}_{n+1} / \lVert \boldsymbol{q}_{n+1} \rVert$.

**Torque handling inside a step.** Two modes are available:

-   **Continuous feedback** (default): the torque law is re-evaluated at
    every RK4 stage, as written above.

-   **Zero-order hold**: the torque is evaluated once,
    $\boldsymbol{\tau}_n = \boldsymbol{\tau}(t_n, \boldsymbol{x}_n)$,
    and the same value is used in all four stages. This models a
    sampled-data (digital) controller.

## Propagation Over a Time Interval

To propagate from $t_0$ to $t_f$ with $T = t_f - t_0 \ge 0$ and
tolerance $\varepsilon = 10^{-9}$:

1.  Number of full steps:
    $n_{\text{full}} = \left\lfloor T/h + \varepsilon \right\rfloor$.

2.  Remainder: $r = T - n_{\text{full}}\, h$. A final partial step of
    size $r$ is taken only if $r > \varepsilon h$.

3.  Take $n_{\text{full}}$ steps of size $h$. After each step the time
    is reset to $t_i = t_0 + i\,h$ (rather than accumulated as
    $t_{i-1} + h$) to avoid floating-point round-off drift.

4.  If required, take the partial step and set the final time exactly to
    $t_f$.

5.  Return the time history and the $(n+1) \times 7$ state history,
    including the initial state.

# Attitude Regulation: Quaternion-Feedback PD Law

## Control Law

For regulation to a fixed target attitude $\boldsymbol{q}_d$:
$$\boldsymbol{\tau} = -K_p\, \operatorname{sgn}(q_{e0})\, \boldsymbol{q}_{ev} \;-\; K_d\, \boldsymbol{\omega},
\qquad
\operatorname{sgn}(q_{e0}) = \begin{cases} +1, & q_{e0} \ge 0 \\ -1, & q_{e0} < 0 \end{cases}$$

-   $-K_p\, \operatorname{sgn}(q_{e0})\, \boldsymbol{q}_{ev}$ ---
    **Proportional (Attitude) Term**: Restoring torque along the error
    rotation axis. Since
    $\boldsymbol{q}_{ev} = \hat{\boldsymbol{e}} \sin(\theta_{\text{err}}/2)$,
    it behaves as a nonlinear spring.

-   $\operatorname{sgn}(q_{e0})$ --- **Shortest-Path Selection**:
    $\boldsymbol{q}_e$ and $-\boldsymbol{q}_e$ represent the same
    attitude. The sign factor always drives the body through the
    rotation of less than $180^\circ$ and prevents unwinding.

-   $-K_d\, \boldsymbol{\omega}$ --- **Derivative (Rate Damping) Term**:
    Dissipates rotational kinetic energy. $K_d$ may be a scalar or a
    $3 \times 3$ matrix.

## Lyapunov Stability Analysis

The candidate function is the rotational kinetic energy plus an attitude
"potential":
$$V(\boldsymbol{q}_e, \boldsymbol{\omega}) = \frac{1}{2}\, \boldsymbol{\omega}^T I\, \boldsymbol{\omega}+ 2 K_p \left( 1 - \lvert q_{e0} \rvert \right)$$
$V \ge 0$, with $V = 0$ only at $\boldsymbol{\omega}= \boldsymbol{0}$
and $q_{e0} = \pm 1$. For a constant target the error kinematics give
$\dot{q}_{e0} = -\tfrac{1}{2}\, \boldsymbol{q}_{ev} \cdot \boldsymbol{\omega}$,
hence $$\begin{aligned}
\dot{V} &= \boldsymbol{\omega}^T \big( \boldsymbol{\tau} - \boldsymbol{\omega}\times I\boldsymbol{\omega}\big) + K_p\, \operatorname{sgn}(q_{e0})\, \boldsymbol{q}_{ev} \cdot \boldsymbol{\omega}\nonumber \\
&= \boldsymbol{\omega}^T \big( \boldsymbol{\tau} + K_p\, \operatorname{sgn}(q_{e0})\, \boldsymbol{q}_{ev} \big)
\qquad (\text{since } \boldsymbol{\omega}\cdot (\boldsymbol{\omega}\times I\boldsymbol{\omega}) = 0)
\end{aligned}$$ Substituting the control law:
$$\dot{V} = -\boldsymbol{\omega}^T K_d\, \boldsymbol{\omega}\;\le\; 0 \quad \text{for } K_d > 0$$
$\dot{V}$ is only negative *semi*-definite (it vanishes whenever
$\boldsymbol{\omega}= \boldsymbol{0}$), so asymptotic convergence
follows from LaSalle's invariance principle:
$\boldsymbol{\omega}\equiv \boldsymbol{0}$ requires
$\dot{\boldsymbol{\omega}} = \boldsymbol{0}$, which forces
$\boldsymbol{q}_{ev} = \boldsymbol{0}$, i.e. the target attitude. Note
that the gyroscopic term does no work and therefore need not be
cancelled for regulation.

## Gain Selection from Settling Time and Damping Ratio

For small errors
$\boldsymbol{q}_{ev} \approx \tfrac{1}{2}\boldsymbol{\theta}$, so each
principal axis behaves as the second-order system
$$I_i \ddot{\theta}_i + K_d\, \dot{\theta}_i + \frac{K_p}{2}\, \theta_i = 0
\quad\Longrightarrow\quad
\omega_n^2 = \frac{K_p}{2 I_i}, \qquad 2 \zeta \omega_n = \frac{K_d}{I_i}$$
Using the 2 % settling-time approximation
$t_s \approx 4 / (\zeta \omega_n)$, and sizing on the largest principal
inertia $I_{\text{ref}} = \max_i I_{ii}$:
$$\omega_n = \frac{4}{\zeta\, t_s},
\qquad
K_p = 2\, I_{\text{ref}}\, \omega_n^2$$ $$K_d =
\begin{cases}
2 \zeta \omega_n I_{\text{ref}} \cdot \mathbf{1}_{3 \times 3}, & \text{isotropic option} \\[1ex]
\operatorname{diag}\!\left( 2 \zeta \omega_n \sqrt{I_{\text{ref}}\, I_{ii}} \right), & \text{per-axis option}
\end{cases}$$ where $\mathbf{1}_{3\times3}$ is the identity matrix. In
the per-axis option the scalar $K_p$ gives axis $i$ a natural frequency
$\omega_n \sqrt{I_{\text{ref}}/I_{ii}}$, and the
$\sqrt{I_{\text{ref}} I_{ii}}$ scaling of $K_d$ makes the damping ratio
equal to $\zeta$ on every axis.

# Attitude Tracking: Feedforward plus PD

## Rate Error

The reference supplies $\boldsymbol{q}_d(t)$, $\boldsymbol{\omega}_d(t)$
and $\dot{\boldsymbol{\omega}}_d(t)$, the latter two expressed in the
desired frame. They are mapped into the body frame with the DCM of the
error quaternion, $C_e = C(\boldsymbol{q}_e)$:
$$\boldsymbol{\omega}_d^B = C_e\, \boldsymbol{\omega}_d,
\qquad
\boldsymbol{\omega}_e = \boldsymbol{\omega}- \boldsymbol{\omega}_d^B$$

## Control Law

$$\boldsymbol{\tau} =
\underbrace{\boldsymbol{\omega}\times (I\boldsymbol{\omega})}_{\text{gyroscopic cancellation}}
+ \underbrace{I \left( C_e\, \dot{\boldsymbol{\omega}}_d - \boldsymbol{\omega}_e \times \boldsymbol{\omega}_d^B \right)}_{\text{reference acceleration feedforward}}
- \underbrace{K_p\, \operatorname{sgn}(q_{e0})\, \boldsymbol{q}_{ev}}_{\text{attitude feedback}}
- \underbrace{K_d\, \boldsymbol{\omega}_e}_{\text{rate-error damping}}$$
The feedforward bracket is the body-frame time derivative of
$\boldsymbol{\omega}_d^B$, because
$\dot{C}_e = -[\boldsymbol{\omega}_e \times]\, C_e$:
$$\frac{d}{dt}\big( C_e\, \boldsymbol{\omega}_d \big) = C_e\, \dot{\boldsymbol{\omega}}_d - \boldsymbol{\omega}_e \times \boldsymbol{\omega}_d^B$$
Substituting into Euler's equation gives the closed-loop error dynamics
$$I\, \dot{\boldsymbol{\omega}}_e = -K_p\, \operatorname{sgn}(q_{e0})\, \boldsymbol{q}_{ev} - K_d\, \boldsymbol{\omega}_e$$
which has the same structure as the regulation problem, so the tracking
error converges with the same Lyapunov argument (with
$\boldsymbol{\omega}$ replaced by $\boldsymbol{\omega}_e$). With a
perfect model and matched initial conditions the tracking error is
identically zero.

# Smooth Reference Generation: Quintic Polynomial Ramp

A step change in the commanded attitude makes $\boldsymbol{\omega}_d$
and $\dot{\boldsymbol{\omega}}_d$ impulsive and produces torque spikes.
Instead, a rotation of total angle $A$ about a fixed unit axis
$\hat{\boldsymbol{e}}$ is commanded between $t_s$ and $t_s + T_r$ using
a $C^2$-smooth profile. With normalised time $\sigma = (t - t_s)/T_r$:
$$\begin{aligned}
s(\sigma) &= 10\sigma^3 - 15\sigma^4 + 6\sigma^5 \\
\dot{s} &= \frac{30\sigma^2 - 60\sigma^3 + 30\sigma^4}{T_r} \\
\ddot{s} &= \frac{60\sigma - 180\sigma^2 + 120\sigma^3}{T_r^2}
\end{aligned}$$ with $(s, \dot{s}, \ddot{s}) = (0, 0, 0)$ for
$t \le t_s$ and $(1, 0, 0)$ for $t \ge t_s + T_r$. This is the unique
fifth-order polynomial satisfying the six boundary conditions
$$s(0) = 0, \quad s(1) = 1, \quad \dot{s}(0) = \dot{s}(1) = 0, \quad \ddot{s}(0) = \ddot{s}(1) = 0$$
so angle, rate and acceleration are all continuous. The reference
signals are
$$\boldsymbol{q}_d(t) = \begin{bmatrix} \cos\!\big( \tfrac{1}{2} A\, s \big) \\ \hat{\boldsymbol{e}}\, \sin\!\big( \tfrac{1}{2} A\, s \big) \end{bmatrix},
\qquad
\boldsymbol{\omega}_d(t) = \hat{\boldsymbol{e}}\, A\, \dot{s},
\qquad
\dot{\boldsymbol{\omega}}_d(t) = \hat{\boldsymbol{e}}\, A\, \ddot{s}$$
For a fixed-axis rotation the components of $\boldsymbol{\omega}_d$ are
identical in the inertial and desired frames. The peak commanded rate is
$\lvert \boldsymbol{\omega}_d \rvert_{\max} = 1.875\, A / T_r$, reached
at $\sigma = \tfrac{1}{2}$.

# Integral Action and Anti-Windup (PID)

## Integral State

A PD law cannot reject a constant disturbance torque: at rest,
$\boldsymbol{\tau}_{\text{dist}} = K_p\, \operatorname{sgn}(q_{e0})\, \boldsymbol{q}_{ev}$
leaves a steady-state attitude error. An integral state
$\boldsymbol{z} \in \mathbb{R}^3$ on the vector part of the error
quaternion is therefore added:
$$\dot{\boldsymbol{z}} = \operatorname{sgn}(q_{e0})\, \boldsymbol{q}_{ev}$$
discretised with forward Euler at the control sample time $\Delta t$,
including both anti-windup mechanisms: $$\boldsymbol{z}_{k+1} =
\begin{cases}
\operatorname{clip}\!\big( \boldsymbol{z}_k + \operatorname{sgn}(q_{e0})\, \boldsymbol{q}_{ev}\, \Delta t,\ -z_{\max},\ +z_{\max} \big), & \text{actuator not saturated} \\[1ex]
\boldsymbol{z}_k, & \text{actuator saturated}
\end{cases}$$

-   **Clamping**: Each component of $\boldsymbol{z}$ is limited to
    $\pm z_{\max}$, bounding the integral torque to $K_i z_{\max}$.

-   **Conditional Integration**: The integrator is frozen while the
    actuator reports saturation, so it does not accumulate error that
    the actuator cannot act upon.

## PID Control Laws

The integral term is subtracted from the PD laws of the previous
sections ($K_p$, $K_i$ may be per-axis vectors applied element-wise,
denoted $\odot$): $$\begin{aligned}
\text{Regulation:} \quad \boldsymbol{\tau} &= -K_p \odot \operatorname{sgn}(q_{e0})\, \boldsymbol{q}_{ev} - K_d\, \boldsymbol{\omega}- K_i \odot \boldsymbol{z} \\
\text{Tracking:} \quad \boldsymbol{\tau} &= \boldsymbol{\omega}\times (I\boldsymbol{\omega}) + I \left( C_e\, \dot{\boldsymbol{\omega}}_d - \boldsymbol{\omega}_e \times \boldsymbol{\omega}_d^B \right) - K_p \odot \operatorname{sgn}(q_{e0})\, \boldsymbol{q}_{ev} - K_d\, \boldsymbol{\omega}_e - K_i \odot \boldsymbol{z}
\end{aligned}$$ Because the controller now has internal state, it is
evaluated exactly once per control sample and *not* inside the RK4
stages.

## Stability Limit on the Integral Gain

Linearising one axis ($\boldsymbol{q}_{ev} \approx \tfrac{1}{2}\theta$)
gives the characteristic polynomial
$$s^3 + \frac{K_d}{I}\, s^2 + \frac{K_p}{2I}\, s + \frac{K_i}{2I} = 0$$
and the Routh--Hurwitz criterion ($a_2 a_1 > a_0$) yields the upper
bound $$0 < K_i < \frac{K_p\, K_d}{I}$$ (neglecting actuator lag, which
reduces the real margin further).

# Actuator Model

## First-Order Lag

The actuator is modelled as a first-order system with time constant
$T_a$:
$$T_a\, \dot{\boldsymbol{\tau}}_{\text{act}} + \boldsymbol{\tau}_{\text{act}} = \boldsymbol{\tau}_{\text{cmd}}$$
Discretising with backward (implicit) Euler, which is unconditionally
stable for any $\Delta t$:
$$\boldsymbol{\tau}_{\text{act}, k+1} = (1 - \alpha_a)\, \boldsymbol{\tau}_{\text{act}, k} + \alpha_a\, \boldsymbol{\tau}_{\text{cmd}, k},
\qquad
\alpha_a = \frac{\Delta t}{T_a + \Delta t}$$

## Magnitude Saturation

The limit is applied to the torque vector norm, so the direction of the
torque is preserved: $$\boldsymbol{\tau}_{\text{act}} \leftarrow
\begin{cases}
\boldsymbol{\tau}_{\text{act}}, & \lVert \boldsymbol{\tau}_{\text{act}} \rVert \le \tau_{\max} \\[1ex]
\tau_{\max}\, \dfrac{\boldsymbol{\tau}_{\text{act}}}{\lVert \boldsymbol{\tau}_{\text{act}} \rVert}, & \lVert \boldsymbol{\tau}_{\text{act}} \rVert > \tau_{\max}
\end{cases}$$ The flag
$\text{saturated} = \big( \lVert \boldsymbol{\tau}_{\text{act}} \rVert > \tau_{\max} \big)$
is passed back to the PID controller on the next sample for conditional
integration.

## Torque Limit from Fin Authority

For the missile, the available control torque is the aerodynamic moment
of a fully deflected fin:
$$\tau_{\max} = \bar{q}\, S\, \bar{c}\, \lvert C_{m\delta} \rvert\, \delta_{\max}$$

# Missile Model: Atmosphere, Inertia and Aerodynamic Moment

## Standard Atmosphere

With $g_0 = 9.80665\ \text{m/s}^2$, $R = 287.05\ \text{J/(kg\,K)}$,
$\gamma = 1.4$ and lapse rate $\lambda = 0.0065\ \text{K/m}$, at
altitude $h$: $$T =
\begin{cases}
288.15 - \lambda h, & h \le 11\,000\ \text{m} \\
216.65\ \text{K}, & h > 11\,000\ \text{m}
\end{cases}$$ $$p =
\begin{cases}
101\,325 \left( \dfrac{T}{288.15} \right)^{g_0 / (\lambda R)}, & h \le 11\,000\ \text{m} \\[2ex]
22\,632.06\, \exp\!\left( -\dfrac{g_0\, (h - 11\,000)}{R\, T} \right), & h > 11\,000\ \text{m}
\end{cases}$$

-   **Air Density** (ideal gas law): $\rho = \dfrac{p}{R\, T}$

-   **Speed of Sound**: $a = \sqrt{\gamma R T}$

-   **Flight Velocity**: $V = M a$

-   **Dynamic Pressure**: $\bar{q} = \tfrac{1}{2} \rho V^2$

## Inertia Matrix

The missile is treated as a uniform solid cylinder of mass $m$, radius
$r$ and length $L$:
$$I_{yy} = I_{zz} = \frac{1}{12} m L^2 + \frac{1}{4} m r^2,
\qquad
I_{xx} = \frac{1}{2} m r^2,
\qquad
I = \operatorname{diag}\!\left( I_{xx},\ I_{yy},\ I_{zz} \right)$$

## Angle of Attack

For the rotation-only demonstration the velocity vector is assumed fixed
along the inertial $x$-axis, so the angle of attack equals the pitch
angle of the yaw--pitch--roll sequence:
$$\alpha = \theta = \arcsin\!\big( 2 (q_0 q_2 - q_1 q_3) \big)$$

## Pitch-Channel Aerodynamic Moment

$$C_m = C_{m\alpha}\, \alpha + C_{m\delta}\, \delta,
\qquad
M_y = \bar{q}\, S\, \bar{c}\, C_m,
\qquad
\boldsymbol{M}_{\text{aero}}(\boldsymbol{q}) = \begin{bmatrix} 0 \\ M_y \\ 0 \end{bmatrix}$$
With $C_{m\alpha} < 0$ the moment is statically stabilising: it opposes
any pitch away from $\alpha = 0$ and therefore acts as a disturbance
that grows with the commanded pitch angle. The fin deflection $\delta$
is set to zero inside this term, because the control moment is supplied
through the actuator torque $\boldsymbol{\tau}_{\text{act}}$.

## Imperfect Model for Feedforward

The controller holds its own copy of the aero model with every
coefficient multiplied by a scale factor $\kappa$:
$$\hat{C}_{m\alpha} = \kappa\, C_{m\alpha},
\qquad
\hat{C}_{m\delta} = \kappa\, C_{m\delta},
\qquad
\hat{\boldsymbol{M}}_{\text{aero}} = \kappa\, \boldsymbol{M}_{\text{aero}}$$
With $\kappa = 0.8$ the controller under-estimates the moment by 20 %.
After cancellation, the residual disturbance
$(1 - \kappa)\, \boldsymbol{M}_{\text{aero}}$ remains, and it is this
residual that the integral term must remove.

# Closed-Loop Missile Pitch Simulation

## Per-Axis Gain Scaling

The roll inertia is roughly $90$ times smaller than the pitch/yaw
inertia. A single scalar $K_p$ sized on $I_{\max}$ would make the roll
loop faster than the actuator bandwidth and unstable. The reference
gains are therefore scaled per axis with $s_i = I_{ii} / I_{\max}$:
$$K_{p,i} = s_i\, K_p^{\text{ref}} = 2\, I_{ii}\, \omega_n^2,
\qquad
K_{d,i} = s_i\, K_d^{\text{ref}} = 2 \zeta \omega_n I_{ii},
\qquad
K_{i,i} = s_i\, K_I$$ so every axis has the same natural frequency
$\omega_n$ and damping ratio $\zeta$.

## Commanded Torque

$$\boldsymbol{\tau}_{\text{cmd}} = \boldsymbol{\tau}_{\text{PID}} - \hat{\boldsymbol{M}}_{\text{aero}}(\boldsymbol{q})$$
where $\boldsymbol{\tau}_{\text{PID}}$ is the tracking PID law and the
second term cancels the (imperfectly known) aerodynamic moment.

## Total Torque on the Plant

$$\boldsymbol{\tau}(t, \boldsymbol{x}) = \underbrace{\boldsymbol{\tau}_{\text{act}, k}}_{\text{held over the sample}} + \underbrace{\boldsymbol{M}_{\text{aero}}\big( \boldsymbol{q}(t) \big)}_{\text{re-evaluated in each RK4 stage}}$$
The actuator output is a zero-order hold over the control sample,
whereas the true aerodynamic moment depends on attitude and is
continuous, so it is evaluated inside the RK4 stages.

## Simulation Algorithm

For $k = 0, 1, \dots, N$ with $N = \operatorname{round}(t_f / \Delta t)$
and $t_k = k\, \Delta t$:

1.  Read the plant state $\boldsymbol{q}_k$, $\boldsymbol{\omega}_k$.

2.  Evaluate the quintic reference: $\boldsymbol{q}_d(t_k)$,
    $\boldsymbol{\omega}_d(t_k)$, $\dot{\boldsymbol{\omega}}_d(t_k)$.

3.  Compute the error quaternion
    $\boldsymbol{q}_e = \boldsymbol{q}_d^{*} \otimes \boldsymbol{q}_k$.

4.  Update the integral state (frozen if the actuator saturated on the
    previous sample) and compute $\boldsymbol{\tau}_{\text{PID}}$.

5.  Subtract the modelled aero moment:
    $\boldsymbol{\tau}_{\text{cmd}} = \boldsymbol{\tau}_{\text{PID}} - \hat{\boldsymbol{M}}_{\text{aero}}(\boldsymbol{q}_k)$.

6.  Pass $\boldsymbol{\tau}_{\text{cmd}}$ through the actuator lag and
    saturation to obtain $\boldsymbol{\tau}_{\text{act}, k}$.

7.  Log $\theta_{\text{err}}$, $\theta$, $\theta_d$,
    $\boldsymbol{\tau}_{\text{cmd}}$, $\boldsymbol{\tau}_{\text{act}}$,
    $\boldsymbol{M}_{\text{aero}}$ and $\boldsymbol{z}$.

8.  If $k < N$: advance the plant one RK4 step of size $\Delta t$ with
    the total torque of the previous subsection.

The simulation is run twice, with $K_I > 0$ (PID) and with $K_I = 0$
(PD), to isolate the effect of integral action on the steady-state
error.

# Simulation Scenarios and Parameters

## Scenario 1: Rest-to-Rest Regulation

A rigid body is regulated from a $70^\circ$ initial rotation about the
body $x$-axis to the identity attitude with the PD law (continuous
feedback). The Lyapunov function $V(t)$ is evaluated along the
trajectory and plotted on a logarithmic scale to confirm monotonic
decay.

## Scenario 2: Smooth Trajectory Tracking

The same rigid body tracks a $40^\circ$ quintic ramp about the $z$-axis
using the feedforward-plus-PD tracking law, starting from rest at the
identity attitude.

::: center
  **Parameter**                                                   **Scenario 1**                             **Scenario 2**
  ----------------------------------------------------- ----------------------------------- ------------------------------------------------
  Inertia $I$ \[kg m$^2$\]                               $\operatorname{diag}(10, 15, 20)$         $\operatorname{diag}(10, 15, 20)$
  Initial attitude                                             $70^\circ$ about $x$                             identity
  Initial rate $\boldsymbol{\omega}_0$                           $\boldsymbol{0}$                           $\boldsymbol{0}$
  Target                                                             identity                $40^\circ$ about $z$, $t_s = 2$ s, $T_r = 6$ s
  Settling time $t_s$ / damping $\zeta$                            $6$ s / $0.8$                             $5$ s / $0.85$
  $\omega_n = 4 / (\zeta t_s)$ \[rad/s\]                              $0.833$                                   $0.941$
  $K_p = 2 I_{\text{ref}}\, \omega_n^2$                               $27.78$                                   $35.43$
  $K_d = 2 \zeta \omega_n I_{\text{ref}}$ (isotropic)                 $26.67$                                   $32.00$
  Step size $h$ / final time                                     $0.01$ s / $15$ s                         $0.01$ s / $15$ s
:::

## Scenario 3: Missile Pitch Tracking with Aero Moment, Actuator and Integral Action

::: center
  **Parameter**                                                                               **Value**
  ------------------------------------------------------------------------ -----------------------------------------------
  Altitude $h$ / Mach number $M$                                                         $10\,000$ m / $2.5$
  Mass $m$, radius $r$, length $L$                                                  $180$ kg, $0.113$ m, $2.65$ m
  Reference area $S$, chord $\bar{c}$                                                  $0.04$ m$^2$, $0.225$ m
  $C_{Z\alpha}$, $C_{m\alpha}$, $C_{m\delta}$, $C_{Z\delta}$ \[per rad\]              $-12$, $-15$, $-20$, $-2$
  Pitch command                                                             $8^\circ$ about $y$, $t_s = 1$ s, $T_r = 3$ s
  Settling time / damping ratio                                                             $1$ s / $0.8$
  Integral gain $K_I$ / integral limit $z_{\max}$                                          $8000$ / $0.5$
  Actuator time constant $T_a$ / fin limit $\delta_{\max}$                              $0.08$ s / $15^\circ$
  Aero model scale $\kappa$                                                                     $0.8$
  Sample time $\Delta t$ / final time                                                     $0.01$ s / $12$ s
:::

Quantities that follow from these parameters through the formulas above:

::: center
  **Derived quantity**                                              **Value**
  --------------------------------------------- --------------------------------------------------
  Temperature $T$                                                   $223.15$ K
  Pressure $p$                                                     $26\,436$ Pa
  Density $\rho$                                                $0.4127$ kg/m$^3$
  Speed of sound $a$ / velocity $V$                        $299.46$ m/s / $748.65$ m/s
  Dynamic pressure $\bar{q}$                                      $115\,657$ Pa
  $I_{yy} = I_{zz}$ / $I_{xx}$                   $105.91$ / $1.149$ kg m$^2$ (ratio $\approx 92$)
  $\omega_n$                                                       $5.0$ rad/s
  $K_p^{\text{ref}}$ / $K_d^{\text{ref}}$                       $5295.6$ / $847.3$
  Integral stability limit $K_p K_d / I_{yy}$      $42\,365$ ($K_I$ is $\approx 19\,\%$ of it)
  Actuator torque limit $\tau_{\max}$                               $5450$ N m
  Aero moment at $\alpha = 8^\circ$                                $-2180$ N m
:::

# Symbol Definitions

-   $\boldsymbol{q}$, $\boldsymbol{q}_d$, $\boldsymbol{q}_e$
    **(Quaternions):** Current, desired and error attitude
    (scalar-first, unit norm).

-   $\boldsymbol{\omega}$, $\boldsymbol{\omega}_d$,
    $\boldsymbol{\omega}_e$ **(Angular Velocities):** Body rate, desired
    rate (desired frame) and rate error (body frame), in rad/s.

-   $I$ **(Inertia Matrix):** Body-frame moment of inertia tensor,
    kg m$^2$.

-   $\boldsymbol{\tau}$, $\boldsymbol{\tau}_{\text{cmd}}$,
    $\boldsymbol{\tau}_{\text{act}}$ **(Torques):** Total, commanded and
    actuator-delivered body torque, N m.

-   $K_p$, $K_d$, $K_i$ **(Controller Gains):** Proportional, derivative
    and integral gains.

-   $\boldsymbol{z}$ **(Integral State):** Accumulated vector part of
    the error quaternion.

-   $c$ **(Constraint Gain):** Quaternion norm stabilisation gain.

-   $\zeta$, $\omega_n$, $t_s$ **(Loop Specifications):** Damping ratio,
    natural frequency and 2 % settling time.

-   $T_a$, $\tau_{\max}$ **(Actuator Parameters):** Time constant and
    torque magnitude limit.

-   $C_e$ **(Error DCM):** Direction cosine matrix from the desired
    frame to the body frame.

-   $\kappa$ **(Model Scale):** Ratio of the controller's aero
    coefficients to the true ones.
