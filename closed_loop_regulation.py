import numpy as np
import matplotlib.pyplot as plt

from kinematics import (
    QuaternionSetup,
    CoupledStateRungeKuttaIntegrator,
    GuidanceAndControlLaw,
    GainEstimation,
    ConservationLawVerifier,
)


def main():
    I = np.diag([10.0, 15.0, 20.0])

    q0 = QuaternionSetup(axis=[1.0, 0.0, 0.0], angle_deg=70.0).quaternion
    omega0 = np.array([0.0, 0.0, 0.0])
    q_desired = np.array([1.0, 0.0, 0.0, 0.0])

    ge = GainEstimation(settling_time=6.0, zeta=0.8, inertia_matrix=I, isotropic=True)
    Kp, Kd = ge.estimate_gains()

    controller = GuidanceAndControlLaw(Kp=Kp, Kd=Kd, inertia_matrix=I)
    print(f"Kp = {Kp}")
    print(f"Kd =\n{Kd}")

    def torque_func(t, x):
        q = x[0:4]
        omega = x[4:7]
        qe = QuaternionSetup.error_quaternion(q, q_desired)
        tau = controller.compute_control_torque(qe, omega)
        return tau

    dt = 0.01
    t_final = 15.0

    integrator = CoupledStateRungeKuttaIntegrator(
        step_size=dt,
        initial_quaternion=q0,
        initial_omega=omega0,
        inertia_matrix=I,
        torque_func=torque_func,
        constraint_gain=1.0,
        hold_torque=False,
    )

    times, states = integrator.propagate(t_final)
    q_hist = states[:, 0:4]
    omega_hist = states[:, 4:7]

    qe_hist = np.array([QuaternionSetup.error_quaternion(q, q_desired) for q in q_hist])
    V_hist = np.array(
        [
            controller.lyapunov_function(qe, omega)
            for qe, omega in zip(qe_hist, omega_hist)
        ]
    )

    tau_hist = np.array([torque_func(t, x) for t, x in zip(times, states)])
    att_error_angle = 2.0 * np.arccos(np.clip(np.abs(qe_hist[:, 0]), -1.0, 1.0))

    print(f"Final quaternion : {np.round(q_hist[-1], 6)}")
    print(f"Final error quat : {np.round(qe_hist[-1], 6)}")
    print(f"Final rates      : {np.round(omega_hist[-1], 6)} rad/s")
    print(f"Final att. error : {np.degrees(att_error_angle[-1]):.4f} deg")
    print(f"Final V          : {V_hist[-1]:.6e}")


if __name__ == "__main__":
    main()
