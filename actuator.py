import numpy as np
from kinematics import GuidanceAndControlLaw


# Smooth reference (5th-order polynomial ramp) – removes tracking spikes
def smooth_axis_angle_ramp(axis, angle_deg, t_start, duration, t):
    """
    Returns (q_des, omega_des, omega_des_dot) for a C2-smooth ramp.
    angle goes from 0 → angle_deg in [t_start, t_start+duration].
    omega_des and omega_des_dot are expressed in the desired frame
    (identical to the inertial frame components for a fixed-axis rotation).
    """
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    A = np.radians(angle_deg)

    if t <= t_start:
        s, sdot, sddot = 0.0, 0.0, 0.0
    elif t >= t_start + duration:
        s, sdot, sddot = 1.0, 0.0, 0.0
    else:
        tau = (t - t_start) / duration
        # 5th-order polynomial with zero velocity & acceleration at ends
        s = 10 * tau**3 - 15 * tau**4 + 6 * tau**5
        sdot = (30 * tau**2 - 60 * tau**3 + 30 * tau**4) / duration
        sddot = (60 * tau - 180 * tau**2 + 120 * tau**3) / duration**2

    angle = A * s
    q = np.concatenate([[np.cos(angle / 2)], axis * np.sin(angle / 2)])
    omega = axis * (A * sdot)
    omega_dot = axis * (A * sddot)
    return q, omega, omega_dot


# Controller with integral action (rejects constant disturbances)
class PIDAttitudeController(GuidanceAndControlLaw):
    """
    PD attitude law from GuidanceAndControlLaw plus an integral term on the
    vector part of the error quaternion.
    """

    def __init__(self, Kp, Kd, Ki, inertia_matrix, integral_limit=2.0):
        super().__init__(Kp, Kd, inertia_matrix)
        self.Ki = Ki
        self.integral = np.zeros(3)
        self.integral_limit = integral_limit

    def reset_integral(self):
        self.integral[:] = 0.0

    def _update_integral(self, qe, dt, saturated):
        qev = self._shortest_path_sign(qe[0]) * np.asarray(qe[1:], dtype=float)
        if not saturated:
            self.integral += qev * dt
            self.integral = np.clip(
                self.integral, -self.integral_limit, self.integral_limit
            )
        return self.Ki * self.integral

    def compute_control_torque(self, qe, omega, dt=0.01, saturated=False):
        tau_pd = super().compute_control_torque(qe, omega)
        return tau_pd - self._update_integral(qe, dt, saturated)

    def compute_tracking_torque(
        self,
        qe,
        omega,
        omega_desired,
        omega_desired_dot=None,
        dt=0.01,
        saturated=False,
    ):
        tau_pd, omega_error = super().compute_tracking_torque(
            qe, omega, omega_desired, omega_desired_dot
        )
        return tau_pd - self._update_integral(qe, dt, saturated), omega_error


# First-order actuator + saturation
class ActuatorModel:
    def __init__(self, tau_max=25.0, time_constant=0.05):
        self.tau_max = tau_max
        self.tc = time_constant
        self.tau_actual = np.zeros(3)
        self.saturated = False

    def step(self, tau_cmd, dt):
        alpha = dt / (self.tc + dt)
        self.tau_actual = (1 - alpha) * self.tau_actual + alpha * tau_cmd
        norm = np.linalg.norm(self.tau_actual)
        self.saturated = norm > self.tau_max
        if self.saturated:
            self.tau_actual *= self.tau_max / norm
        return self.tau_actual.copy()


# Aero moment model
class SimpleAeroMoments:
    """
    Pitch-channel aero moment for a rotation-only demonstration.

    The velocity vector is assumed fixed along the inertial x-axis, so the
    angle of attack equals the pitch angle theta of the yaw-pitch-roll
    sequence: sin(theta) = 2 (q0 q2 - q1 q3).
    """

    def __init__(self, missile, model_scale=1.0):
        self.q_bar = missile.q_bar
        self.S = missile.S
        self.c_bar = missile.c_bar
        self.Cm_alpha = model_scale * missile.c_m_alpha
        self.Cm_delta = model_scale * missile.c_m_delta

    @staticmethod
    def angle_of_attack(q):
        return np.arcsin(np.clip(2.0 * (q[0] * q[2] - q[1] * q[3]), -1.0, 1.0))

    def moment(self, q, omega, delta=0.0):
        alpha = self.angle_of_attack(q)
        Cm = self.Cm_alpha * alpha + self.Cm_delta * delta
        My = self.q_bar * self.S * self.c_bar * Cm
        return np.array([0.0, My, 0.0])
