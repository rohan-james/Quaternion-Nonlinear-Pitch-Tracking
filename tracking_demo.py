import numpy as np
import matplotlib.pyplot as plt

from kinematics import (
    QuaternionSetup,
    CoupledStateRungeKuttaIntegrator,
    GuidanceAndControlLaw,
    GainEstimation,
)
from actuator import smooth_axis_angle_ramp


def main():
    I = np.diag([10.0, 15.0, 20.0])
    ge = GainEstimation(settling_time=5.0, zeta=0.85, inertia_matrix=I)
    Kp, Kd = ge.estimate_gains()
    controller = GuidanceAndControlLaw(Kp, Kd, I)

    # Desired trajectory: C2-smooth 40° rotation about Z, 6 s long, from t=2 s
    def reference(t):
        return smooth_axis_angle_ramp(
            axis=[0, 0, 1], angle_deg=40.0, t_start=2.0, duration=6.0, t=t
        )

    def q_desired(t):
        return reference(t)[0]

    def torque_func(t, x):
        q = x[0:4]
        omega = x[4:7]
        qd, wd, wd_dot = reference(t)
        qe = QuaternionSetup.error_quaternion(q, qd)
        tau, _ = controller.compute_tracking_torque(qe, omega, wd, wd_dot)
        return tau

    q0 = np.array([1.0, 0.0, 0.0, 0.0])
    omega0 = np.zeros(3)

    integ = CoupledStateRungeKuttaIntegrator(
        step_size=0.01,
        initial_quaternion=q0,
        initial_omega=omega0,
        inertia_matrix=I,
        torque_func=torque_func,
    )
    times, states = integ.propagate(15.0)

    # Post-process
    q_hist = states[:, 0:4]
    omega_hist = states[:, 4:7]
    qe_hist = np.array(
        [
            QuaternionSetup.error_quaternion(q, q_desired(t))
            for t, q in zip(times, q_hist)
        ]
    )
    att_err = 2.0 * np.arccos(np.clip(np.abs(qe_hist[:, 0]), -1.0, 1.0))

    print(f"Peak attitude error : {np.degrees(att_err.max()):.4f} deg")
    print(f"Final attitude error: {np.degrees(att_err[-1]):.4f} deg")


if __name__ == "__main__":
    main()
