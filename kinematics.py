import scipy.integrate
import numpy as np
import scipy.linalg


class QuaternionSetup:
    def __init__(self, axis, angle_deg):
        self.quaternion = None
        self.quaternion_conjugate = None

        self._axis_angles_to_quaternion(axis, angle_deg)

    def _axis_angles_to_quaternion(self, axis, angle_deg):
        axis = np.array(axis, dtype=np.float64)
        axis_norm = np.linalg.norm(axis)
        if axis_norm < 1e-12:
            raise ValueError("Rotation axis must be non-zero")
        axis = axis / axis_norm
        half_angle = np.radians(angle_deg) / 2.0
        w = np.cos(half_angle)
        xyz = axis * np.sin(half_angle)

        self.quaternion = np.array([w, xyz[0], xyz[1], xyz[2]])
        self.quaternion_conjugate = np.array([w, -xyz[0], -xyz[1], -xyz[2]])

    @staticmethod
    def _quaternion_product(p, q):
        p0, pv = p[0], p[1:]
        q0, qv = q[0], q[1:]

        scalar_part = p0 * q0 - np.dot(pv, qv)
        vector_part = p0 * qv + q0 * pv + np.cross(pv, qv)

        return np.concatenate(([scalar_part], vector_part))

    def calculate_point_rotation(self, point_3d):
        p = np.concatenate(([0.0], point_3d))
        intermediate = self._quaternion_product(self.quaternion, p)
        p_prime = self._quaternion_product(intermediate, self.quaternion_conjugate)
        return p_prime[1:]

    @staticmethod
    def error_quaternion(q, quaternion_desired):
        q = np.asarray(q, dtype=np.float64)
        q = q / np.linalg.norm(q)
        qd = np.asarray(quaternion_desired, dtype=np.float64)
        qd = qd / np.linalg.norm(qd)

        q0, qv = q[0], q[1:]
        qd0, qdv = qd[0], qd[1:]

        q_e0 = q0 * qd0 + np.dot(qv, qdv)
        q_ev = (qd0 * qv) - (q0 * qdv) + np.cross(qv, qdv)

        return np.concatenate(([q_e0], q_ev))

    def calculate_error_quaternion(self, quaternion_desired):
        return self.error_quaternion(self.quaternion, quaternion_desired)


class QuaternionKinematicDifferentialEquation:
    def __init__(self, c=1.0):
        self.c = c

    def calculate_dynamic_time_derivative(self, q, omega):
        wx, wy, wz = omega
        omega_matrix = np.array(
            [
                [0.0, -wx, -wy, -wz],
                [wx, 0.0, wz, -wy],
                [wy, -wz, 0.0, wx],
                [wz, wy, -wx, 0.0],
            ],
            dtype=np.float64,
        )

        kinematic_term = 0.5 * (omega_matrix @ q)
        constraint_term = self.c * (1.0 - np.dot(q, q)) * q
        return kinematic_term + constraint_term


class SecondaryRepresentationsAndConversions:
    def __init__(self, quaternion):
        self.quaternion = np.array(quaternion, dtype=np.float64) / np.linalg.norm(
            quaternion
        )
        self.identity_matrix = np.identity(3)

    def direction_cosine_matrix(self):
        q0, q1, q2, q3 = self.quaternion
        r11 = q0**2 + q1**2 - q2**2 - q3**2
        r12 = 2 * (q1 * q2 + q0 * q3)
        r13 = 2 * (q1 * q3 - q0 * q2)

        r21 = 2 * (q1 * q2 - q0 * q3)
        r22 = q0**2 - q1**2 + q2**2 - q3**2
        r23 = 2 * (q2 * q3 + q0 * q1)

        r31 = 2 * (q1 * q3 + q0 * q2)
        r32 = 2 * (q2 * q3 - q0 * q1)
        r33 = q0**2 - q1**2 - q2**2 + q3**2

        return np.array(
            [[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]], dtype=np.float64
        )

    def euler_yaw_pitch_roll(self):
        q0, q1, q2, q3 = self.quaternion

        theta = -np.arcsin(np.clip(2 * ((q1 * q3) - (q0 * q2)), -1.0, 1.0))
        phi = np.arctan2(2 * ((q2 * q3) + (q0 * q1)), q0**2 - q1**2 - q2**2 + q3**2)
        psi = np.arctan2(2 * ((q1 * q2) + (q0 * q3)), q0**2 + q1**2 - q2**2 - q3**2)

        return phi, theta, psi

    def modified_rodrigues_parameters(self):
        q0, qv = self.quaternion[0], self.quaternion[1:]
        if q0 >= 0:
            sigma = qv / (1.0 + q0)
        else:
            sigma = -qv / (1.0 - q0)

        return sigma


class RungeKuttaDynamicNumericalIntegrator:
    def __init__(self, step_size, initial_quaternion, omega_func, constraint_gain=1.0):
        self.q = np.array(initial_quaternion, dtype=np.float64) / np.linalg.norm(
            initial_quaternion
        )
        self.h = step_size
        self.omega_func = (
            omega_func if callable(omega_func) else (lambda t: np.array(omega_func))
        )
        self.kinematics = QuaternionKinematicDifferentialEquation(c=constraint_gain)

    def step(self, t):
        omega1 = self.omega_func(t)
        k1 = self.h * self.kinematics.calculate_dynamic_time_derivative(self.q, omega1)

        omega2 = self.omega_func(t + 0.5 * self.h)
        k2 = self.h * self.kinematics.calculate_dynamic_time_derivative(
            self.q + 0.5 * k1, omega2
        )
        k3 = self.h * self.kinematics.calculate_dynamic_time_derivative(
            self.q + 0.5 * k2, omega2
        )

        omega4 = self.omega_func(t + self.h)
        k4 = self.h * self.kinematics.calculate_dynamic_time_derivative(
            self.q + k3, omega4
        )

        q_new = self.q + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        self.q = q_new / np.linalg.norm(q_new)
        return self.q.copy()


class RigidBodyDynamics:
    def __init__(self, inertia_matrix, omega, constraint_gain=1.0, use_lu=False):
        self.omega = omega
        self.I = np.array(inertia_matrix, dtype=np.float64)
        self.c = constraint_gain
        self.tau = np.zeros(3)
        self.use_lu = use_lu

        if self.use_lu:
            self.lu, self.piv = scipy.linalg.lu_factor(self.I)
            self.I_inv = None
        else:
            self.I_inv = np.linalg.inv(self.I)
            self.lu, self.piv = None, None

    def compute_angular_acceleration(self):
        gyroscopic_torque = np.cross(self.omega, (self.I @ self.omega))
        b = self.tau - gyroscopic_torque

        if self.use_lu:
            angular_accel = scipy.linalg.lu_solve((self.lu, self.piv), b)
        else:
            angular_accel = self.I_inv @ b

        return angular_accel


class NonLinearStateSpaceFormulation:
    def __init__(self, inertia_matrix, constraint_gain=1.0, use_lu=False):
        dummy_omega = np.zeros(3)
        self.c = constraint_gain
        self.rbd = RigidBodyDynamics(
            inertia_matrix, dummy_omega, constraint_gain=constraint_gain, use_lu=use_lu
        )
        self.qkde = QuaternionKinematicDifferentialEquation(c=constraint_gain)

    def compute_coupled_state_space(self, x, tau):
        q = x[0:4]
        omega = x[4:7]

        self.rbd.omega = omega
        self.rbd.tau = tau

        angular_accel = self.rbd.compute_angular_acceleration()

        state_derivative = np.concatenate(
            [self.qkde.calculate_dynamic_time_derivative(q, omega), angular_accel]
        )
        return state_derivative


class CoupledStateRungeKuttaIntegrator:
    def __init__(
        self,
        step_size,
        initial_quaternion,
        initial_omega,
        inertia_matrix,
        torque_func=None,
        constraint_gain=1.0,
        use_lu=False,
        hold_torque=False,
    ):

        q0 = np.array(initial_quaternion, dtype=np.float64)
        q0 = q0 / np.linalg.norm(q0)
        omega0 = np.array(initial_omega, dtype=np.float64)

        self.x = np.concatenate([q0, omega0])
        self.t = 0.0
        self.h = step_size
        self.hold_torque = hold_torque
        self.torque_func = (
            torque_func if torque_func is not None else (lambda t, x: np.zeros(3))
        )
        self.plant = NonLinearStateSpaceFormulation(
            inertia_matrix, constraint_gain=constraint_gain, use_lu=use_lu
        )

    def _derivative(self, t, x, tau=None):
        if tau is None:
            tau = np.asarray(self.torque_func(t, x), dtype=np.float64)
        return self.plant.compute_coupled_state_space(x, tau)

    def step(self, h=None):
        h = self.h if h is None else h
        t, x = self.t, self.x
        tau = (
            np.asarray(self.torque_func(t, x), dtype=np.float64)
            if self.hold_torque
            else None
        )

        k1 = self._derivative(t, x, tau)
        k2 = self._derivative(t + 0.5 * h, x + 0.5 * h * k1, tau)
        k3 = self._derivative(t + 0.5 * h, x + 0.5 * h * k2, tau)
        k4 = self._derivative(t + h, x + h * k3, tau)

        x_new = x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        x_new[0:4] = x_new[0:4] / np.linalg.norm(x_new[0:4])

        self.x = x_new
        self.t = t + h
        return self.t, self.x.copy()

    def propagate(self, t_final, tol=1e-9):
        t0 = self.t
        duration = t_final - t0
        if duration < 0.0:
            raise ValueError(f"t_final ({t_final}) is earlier than current time ({t0})")

        n_full = int(np.floor(duration / self.h + tol))
        remainder = duration - n_full * self.h
        has_partial = remainder > tol * self.h
        n_steps = n_full + int(has_partial)

        times = np.empty(n_steps + 1)
        states = np.empty((n_steps + 1, 7))
        times[0], states[0] = t0, self.x

        for i in range(1, n_full + 1):
            _, states[i] = self.step()
            self.t = times[i] = t0 + i * self.h

        if has_partial:
            _, states[n_steps] = self.step(remainder)
            self.t = times[n_steps] = t_final

        return times, states


class ConservationLawVerifier:
    def __init__(self, inertia_matrix):
        self.I = np.array(inertia_matrix, dtype=np.float64)

    def compute_rotational_kinetic_energy(self, omega):
        return 0.5 * np.dot(omega, self.I @ omega)

    def compute_angular_momentum(self, q, omega):
        H_B = self.I @ omega
        R = SecondaryRepresentationsAndConversions(q).direction_cosine_matrix()
        H_I = R.T @ H_B
        return H_B, H_I


class GuidanceAndControlLaw:
    def __init__(self, Kp, Kd, inertia_matrix, use_lu=False):
        self.I = np.array(inertia_matrix, dtype=np.float64)
        self.Kp = Kp
        self.Kd = np.array(Kd, dtype=np.float64) if np.ndim(Kd) > 0 else Kd
        self.use_lu = use_lu
        self.clv = ConservationLawVerifier(self.I)

    @staticmethod
    def _shortest_path_sign(qe0):
        return 1.0 if qe0 >= 0.0 else -1.0

    def lyapunov_function(self, qe, omega):
        qe0 = qe[0]
        rate_err = self.clv.compute_rotational_kinetic_energy(omega)
        attitude_err = 2 * self.Kp * (1.0 - abs(qe0))
        return rate_err + attitude_err

    def compute_lyapunov_derivative(self, omega):
        omega = np.asarray(omega, dtype=np.float64)
        if isinstance(self.Kd, np.ndarray) and self.Kd.ndim == 2:
            return -omega @ self.Kd @ omega
        return -np.dot(omega, self.Kd * omega)

    def compute_control_torque(self, qe, omega):
        qe = np.asarray(qe, dtype=np.float64)
        omega = np.asarray(omega, dtype=np.float64)
        qev = self._shortest_path_sign(qe[0]) * qe[1:]
        if isinstance(self.Kd, np.ndarray) and self.Kd.ndim == 2:
            damping_term = self.Kd @ omega
        else:
            damping_term = self.Kd * omega
        return -self.Kp * qev - damping_term

    def compute_tracking_torque(self, qe, omega, omega_desired, omega_desired_dot=None):

        qe = np.asarray(qe, dtype=np.float64)
        omega = np.asarray(omega, dtype=np.float64)
        omega_desired = np.asarray(omega_desired, dtype=np.float64)
        if omega_desired_dot is None:
            omega_desired_dot = np.zeros(3)
        omega_desired_dot = np.asarray(omega_desired_dot, dtype=np.float64)

        C = SecondaryRepresentationsAndConversions(qe).direction_cosine_matrix()
        omega_d_body = C @ omega_desired
        omega_error = omega - omega_d_body

        feedforward = np.cross(omega, self.I @ omega) + self.I @ (
            C @ omega_desired_dot - np.cross(omega_error, omega_d_body)
        )

        qev = self._shortest_path_sign(qe[0]) * qe[1:]
        if isinstance(self.Kd, np.ndarray) and self.Kd.ndim == 2:
            damping_term = self.Kd @ omega_error
        else:
            damping_term = self.Kd * omega_error

        return feedforward - self.Kp * qev - damping_term, omega_error


class GainEstimation:
    def __init__(self, settling_time, zeta=0.707, inertia_matrix=None, isotropic=True):
        self.settling_time = settling_time
        self.zeta = zeta
        self.inertia_matrix = inertia_matrix
        self.isotropic = isotropic

    def estimate_gains(self):
        I_diag = (
            np.diag(self.inertia_matrix)
            if self.inertia_matrix is not None
            else np.array([1.0, 1.0, 1.0])
        )

        omega_n = 4.0 / (self.zeta * self.settling_time)
        if self.isotropic:
            I_ref = np.max(I_diag)
            Kp = 2.0 * I_ref * (omega_n**2)
            Kd = 2.0 * self.zeta * omega_n * I_ref
            Kd = Kd * np.eye(3)
        else:
            I_ref = np.max(I_diag)
            Kp = 2.0 * I_ref * (omega_n**2)
            Kd = np.diag(2.0 * self.zeta * omega_n * np.sqrt(I_ref * I_diag))
        return Kp, Kd


class SingleAxisMissilePitch:
    def __init__(
        self,
        altitude=10000.0,
        reference_area=0.04,
        chord_length=0.225,
        mass=180.0,
        radius=0.113,
        length=2.65,
        mach_number=2.5,
        normal_force_alpha_derivative=-12.0,
        pitching_moment_alpha_derivative=-15.0,
        pitch_control_authority=-20.0,
        pitch_control_force_derivative=-2.0,
    ):
        self.mach_number = mach_number
        self.altitude = altitude
        self.S = reference_area
        self.c_bar = chord_length
        self.m = mass
        self.radius = radius
        self.length = length
        self.c_z_alpha = normal_force_alpha_derivative
        self.c_m_alpha = pitching_moment_alpha_derivative
        self.c_m_delta = pitch_control_authority
        self.c_z_delta = pitch_control_force_derivative
        g0, R_air, gamma = 9.80665, 287.05, 1.4
        if self.altitude <= 11000.0:
            self.temperature = 288.15 - 0.0065 * self.altitude
            self.pressure = 101325.0 * (self.temperature / 288.15) ** (
                g0 / (0.0065 * R_air)
            )
        else:
            self.temperature = 216.65
            self.pressure = 22632.06 * np.exp(
                -g0 * (self.altitude - 11000.0) / (R_air * self.temperature)
            )
        self.density = self.pressure / (R_air * self.temperature)
        self.speed_of_sound = np.sqrt(gamma * R_air * self.temperature)
        self.velocity = self.mach_number * self.speed_of_sound
        self.q_bar = 0.5 * (self.density * self.velocity**2)
        self.pitch_inertia = (1.0 / 12.0) * self.m * (self.length**2) + (
            1.0 / 4.0
        ) * self.m * (self.radius**2)

        self.A = None
        self.B = None
        self.Q = None
        self.R = None

    def _initialise_kinematic_differential_equations(self):
        pass

    def _main(self):
        pass


if __name__ == "__main__":
    samp = SingleAxisMissilePitch()
    samp._main()
