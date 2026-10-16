import numpy as np
import matplotlib.pyplot as plt

from kinematics import (
    QuaternionSetup,
    CoupledStateRungeKuttaIntegrator,
    GainEstimation,
    SingleAxisMissilePitch,
)
from actuator import (
    smooth_axis_angle_ramp,
    PIDAttitudeController,
    ActuatorModel,
    SimpleAeroMoments,
)

# Scenario
PITCH_COMMAND_DEG = 8.0  # pitch (= angle of attack)
RAMP_START = 1.0
RAMP_DURATION = 3.0
T_FINAL = 12.0
DT = 0.01
MAX_FIN_DEFLECTION_DEG = 15.0
AERO_MODEL_SCALE = 0.8
SETTLING_TIME = 1.0
KI = 8000.0


def reference(t):
    return smooth_axis_angle_ramp(
        axis=[0, 1, 0],
        angle_deg=PITCH_COMMAND_DEG,
        t_start=RAMP_START,
        duration=RAMP_DURATION,
        t=t,
    )


def simulate(use_integral=True):

    missile = SingleAxisMissilePitch()
    I_roll = 0.5 * missile.m * missile.radius**2
    I_true = np.diag([I_roll, missile.pitch_inertia, missile.pitch_inertia])
    aero = SimpleAeroMoments(missile)
    aero_model = SimpleAeroMoments(missile, model_scale=AERO_MODEL_SCALE)

    ge = GainEstimation(settling_time=SETTLING_TIME, zeta=0.8, inertia_matrix=I_true)
    Kp_ref, Kd_ref = ge.estimate_gains()
    axis_scale = np.diag(I_true) / np.max(np.diag(I_true))
    Kp = Kp_ref * axis_scale
    Kd = Kd_ref * np.diag(axis_scale)
    Ki = (KI if use_integral else 0.0) * axis_scale
    controller = PIDAttitudeController(Kp, Kd, Ki, I_true, integral_limit=0.5)

    tau_max = (
        missile.q_bar
        * missile.S
        * missile.c_bar
        * abs(missile.c_m_delta)
        * np.radians(MAX_FIN_DEFLECTION_DEG)
    )
    actuator = ActuatorModel(tau_max=tau_max, time_constant=0.08)

    q0 = np.array([1.0, 0.0, 0.0, 0.0])
    omega0 = np.zeros(3)
    plant = CoupledStateRungeKuttaIntegrator(
        step_size=DT,
        initial_quaternion=q0,
        initial_omega=omega0,
        inertia_matrix=I_true,
    )

    log = {
        k: [] for k in ("t", "err", "theta", "theta_des", "cmd", "act", "aero", "int")
    }
    n_steps = int(round(T_FINAL / DT))

    for k in range(n_steps + 1):
        t = k * DT
        x = plant.x
        q, omega = x[0:4], x[4:7]

        q_des, w_des, w_des_dot = reference(t)
        qe = QuaternionSetup.error_quaternion(q, q_des)

        tau_pid, _ = controller.compute_tracking_torque(
            qe, omega, w_des, w_des_dot, dt=DT, saturated=actuator.saturated
        )
        tau_cmd = tau_pid - aero_model.moment(q, omega)

        tau_act = actuator.step(tau_cmd, DT)

        log["t"].append(t)
        log["err"].append(2.0 * np.arccos(np.clip(abs(qe[0]), -1.0, 1.0)))
        log["theta"].append(aero.angle_of_attack(q))
        log["theta_des"].append(aero.angle_of_attack(q_des))
        log["cmd"].append(tau_cmd.copy())
        log["act"].append(tau_act.copy())
        log["aero"].append(aero.moment(q, omega))
        log["int"].append(controller.integral.copy())

        if k == n_steps:
            break

        plant.torque_func = lambda tt, xx, ta=tau_act: ta + aero.moment(
            xx[0:4], xx[4:7]
        )
        plant.t = t
        plant.step(DT)

    out = {k: np.array(v) for k, v in log.items()}
    out["tau_max"] = tau_max
    return out


def main():
    pid = simulate(use_integral=True)
    pd = simulate(use_integral=False)

    for name, r in (("With integral   ", pid), ("Without integral", pd)):
        print(
            f"{name}: peak error {np.degrees(r['err'].max()):.3f} deg, "
            f"final error {np.degrees(r['err'][-1]):.4f} deg"
        )
    print(f"Final integral state : {np.round(pid['int'][-1], 4)}")
    print(f"Final aero moment    : {pid['aero'][-1, 1]:.1f} Nm")
    print(
        f"Final actuator torque: {pid['act'][-1, 1]:.1f} Nm"
        f"(limit {pid['tau_max']:.0f} Nm)"
    )


if __name__ == "__main__":
    main()
