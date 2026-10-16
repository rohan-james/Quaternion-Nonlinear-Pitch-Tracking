import numpy as np
import matplotlib.pyplot as plt

from kinematics import (
    QuaternionSetup,
    CoupledStateRungeKuttaIntegrator,
    GuidanceAndControlLaw,
    GainEstimation,
)


def run_case(name, I_true, I_ctrl, disturbance_func, t_final=20.0, dt=0.01):
    """I_true = plant, I_ctrl = matrix the controller believes."""
    ge = GainEstimation(settling_time=6.0, zeta=0.8, inertia_matrix=I_ctrl)
    Kp, Kd = ge.estimate_gains()
    controller = GuidanceAndControlLaw(Kp, Kd, I_ctrl)  # controller uses wrong I

    q0 = QuaternionSetup([1, 0.5, 0.2], 80.0).quaternion
    omega0 = np.array([0.1, -0.05, 0.08])
    q_des = np.array([1.0, 0.0, 0.0, 0.0])

    def torque_func(t, x):
        qe = QuaternionSetup.error_quaternion(x[0:4], q_des)
        tau_ctrl = controller.compute_control_torque(qe, x[4:7])
        return tau_ctrl + disturbance_func(t)

    integ = CoupledStateRungeKuttaIntegrator(
        step_size=dt,
        initial_quaternion=q0,
        initial_omega=omega0,
        inertia_matrix=I_true,  # true plant
        torque_func=torque_func,
    )
    times, states = integ.propagate(t_final)

    qe = np.array([QuaternionSetup.error_quaternion(s[0:4], q_des) for s in states])
    att_err = 2.0 * np.arccos(np.clip(np.abs(qe[:, 0]), -1.0, 1.0))
    tau_hist = np.array([torque_func(t, x) for t, x in zip(times, states)])

    return times, att_err, states[:, 4:7], tau_hist


def main():
    I_true = np.diag([10.0, 15.0, 20.0])

    # Case A – nominal
    times_n, err_n, _, tau_n = run_case(
        "Nominal", I_true, I_true, lambda t: np.zeros(3)
    )

    # Case B – inertia mismatch (controller thinks inertias are 30 % higher)
    I_wrong = 1.3 * I_true
    times_m, err_m, _, tau_m = run_case(
        "Inertia +30 %", I_true, I_wrong, lambda t: np.zeros(3)
    )

    # Case C – constant disturbance torque
    dist = np.array([0.8, -0.5, 0.3])  # N·m
    times_d, err_d, _, tau_d = run_case(
        "Const. disturbance", I_true, I_true, lambda t: dist
    )

    # Case D – inertia mismatch + sinusoidal disturbance
    def sino_dist(t):
        return np.array(
            [0.4 * np.sin(0.7 * t), 0.3 * np.cos(1.1 * t), 0.2 * np.sin(0.5 * t)]
        )

    times_b, err_b, _, tau_b = run_case(
        "Mismatch + sinusoid", I_true, I_wrong, sino_dist
    )

    print("Final attitude errors [deg]:")
    print(f"  Nominal              : {np.degrees(err_n[-1]):.4f}")
    print(f"  Inertia +30 %        : {np.degrees(err_m[-1]):.4f}")
    print(f"  Const. disturbance   : {np.degrees(err_d[-1]):.4f}")
    print(f"  Mismatch + sinusoid  : {np.degrees(err_b[-1]):.4f}")


if __name__ == "__main__":
    main()
