import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch
import os

L1, L2 = 0.5, 0.5
m1, m2 = 1.0, 0.5
I1 = (1/12) * m1 * L1**2
I2 = (1/12) * m2 * L2**2
lc1 = L1 / 2
lc2 = L2 / 2
g = 9.81

def mass_matrix(q):
    q2 = q[1]
    c2 = np.cos(q2)
    M11 = I1 + I2 + m1*lc1**2 + m2*(L1**2 + lc2**2 + 2*L1*lc2*c2)
    M12 = I2 + m2*(lc2**2 + L1*lc2*c2)
    M22 = I2 + m2*lc2**2
    return np.array([[M11, M12], [M12, M22]])

def coriolis_matrix(q, dq):
    q2  = q[1]
    dq1 = dq[0]
    dq2 = dq[1]
    h = m2 * L1 * lc2 * np.sin(q2)
    C = np.array([[-h*dq2, -h*(dq1 + dq2)], [h*dq1, 0.0]])
    return C

def gravity_vector(q):
    q1, q2 = q[0], q[1]
    G1 = (m1*lc1 + m2*L1)*g*np.cos(q1) + m2*lc2*g*np.cos(q1 + q2)
    G2 = m2*lc2*g*np.cos(q1 + q2)
    return np.array([G1, G2])

Kp = np.diag([150.0, 150.0])
Kv = np.diag([30.0,  30.0])

def computed_torque_control(q, dq, q_d, dq_d, ddq_d):
    M = mass_matrix(q)
    C = coriolis_matrix(q, dq)
    G = gravity_vector(q)
    e   = q_d  - q
    de  = dq_d - dq
    v   = ddq_d + Kv @ de + Kp @ e
    tau = M @ v + C @ dq + G
    return tau

def inverse_kinematics(x, y):
    r2 = x**2 + y**2
    cos_q2 = (r2 - L1**2 - L2**2) / (2*L1*L2)
    cos_q2 = np.clip(cos_q2, -1.0, 1.0)
    sin_q2 = np.sqrt(1.0 - cos_q2**2)
    q2 = np.arctan2(sin_q2, cos_q2)
    q1 = np.arctan2(y, x) - np.arctan2(L2*np.sin(q2), L1 + L2*np.cos(q2))
    return np.array([q1, q2])

def generate_circular_trajectory(t_arr):
    cx, cy = 0.5, 0.0
    r  = 0.15
    T  = t_arr[-1]
    omega = 2*np.pi / T
    N = len(t_arr)
    q_d   = np.zeros((N, 2))
    dq_d  = np.zeros((N, 2))
    ddq_d = np.zeros((N, 2))

    for i, t in enumerate(t_arr):
        xd  =  cx + r * np.cos(omega * t)
        yd  =  cy + r * np.sin(omega * t)
        dxd = -r * omega * np.sin(omega * t)
        dyd =  r * omega * np.cos(omega * t)
        ddxd = -r * omega**2 * np.cos(omega * t)
        ddyd = -r * omega**2 * np.sin(omega * t)
        q_d[i] = inverse_kinematics(xd, yd)

    dt = np.diff(t_arr)
    for j in range(2):
        dq_d[:-1, j]  = np.diff(q_d[:, j]) / dt
        dq_d[-1,  j]  = dq_d[-2,  j]
        ddq_d[:-1, j] = np.diff(dq_d[:, j]) / dt
        ddq_d[-1,  j] = ddq_d[-2, j]
    return q_d, dq_d, ddq_d

T_end = 5.0
dt    = 0.005
t_arr = np.arange(0, T_end + dt, dt)
N     = len(t_arr)

q_d, dq_d, ddq_d = generate_circular_trajectory(t_arr)

q_hist   = np.zeros((N, 2))
dq_hist  = np.zeros((N, 2))
tau_hist = np.zeros((N, 2))

q_hist[0]  = q_d[0]
dq_hist[0] = dq_d[0]

def forward_dynamics(q, dq, tau):
    M = mass_matrix(q)
    C = coriolis_matrix(q, dq)
    G = gravity_vector(q)
    ddq = np.linalg.solve(M, tau - C @ dq - G)
    return ddq

for i in range(N - 1):
    q  = q_hist[i]
    dq = dq_hist[i]
    tau = computed_torque_control(q, dq, q_d[i], dq_d[i], ddq_d[i])
    tau_hist[i] = tau
    ddq = forward_dynamics(q, dq, tau)
    dq_hist[i+1] = dq + ddq * dt
    q_hist[i+1]  = q  + dq  * dt + 0.5 * ddq * dt**2

tau_hist[-1] = tau_hist[-2]

def forward_kinematics(q_arr):
    x = L1 * np.cos(q_arr[:, 0]) + L2 * np.cos(q_arr[:, 0] + q_arr[:, 1])
    y = L1 * np.sin(q_arr[:, 0]) + L2 * np.sin(q_arr[:, 0] + q_arr[:, 1])
    return x, y

x_actual, y_actual = forward_kinematics(q_hist)
x_desired, y_desired = forward_kinematics(q_d)

error = q_d - q_hist

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

plt.rcParams.update({'font.size': 11, 'figure.dpi': 150})
colors = ['#1f77b4', '#ff7f0e']

fig1, ax = plt.subplots(figsize=(6, 6))
ax.plot(x_desired, y_desired, '--', color='gray', lw=2)
ax.plot(x_actual,  y_actual,  '-',  color=colors[0], lw=1.5)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.axis('equal'); ax.grid(True, alpha=0.4)
fig1.tight_layout()
fig1.savefig(os.path.join(OUTPUT_DIR, 'fig_trajectory.png'), dpi=150)
plt.close(fig1)

fig2, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
for j in range(2):
    axes[j].plot(t_arr, np.rad2deg(q_d[:, j]), '--', color='gray')
    axes[j].plot(t_arr, np.rad2deg(q_hist[:, j]), '-', color=colors[j])
    axes[j].grid(True, alpha=0.4)
fig2.tight_layout()
fig2.savefig(os.path.join(OUTPUT_DIR, 'fig_joint_angles.png'), dpi=150)
plt.close(fig2)

fig3, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
for j in range(2):
    axes[j].plot(t_arr, np.rad2deg(error[:, j]), color=colors[j])
    axes[j].grid(True, alpha=0.4)
    axes[j].axhline(0, color='k', lw=0.8)
fig3.tight_layout()
fig3.savefig(os.path.join(OUTPUT_DIR, 'fig_errors.png'), dpi=150)
plt.close(fig3)

fig4, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
for j in range(2):
    axes[j].plot(t_arr, tau_hist[:, j], color=colors[j])
    axes[j].grid(True, alpha=0.4)
fig4.tight_layout()
fig4.savefig(os.path.join(OUTPUT_DIR, 'fig_torques.png'), dpi=150)
plt.close(fig4)

rms1 = np.sqrt(np.mean(error[:, 0]**2))
rms2 = np.sqrt(np.mean(error[:, 1]**2))
fig5, ax = plt.subplots(figsize=(5, 4))
ax.bar(['Joint 1', 'Joint 2'], [np.rad2deg(rms1), np.rad2deg(rms2)], color=colors, edgecolor='k', alpha=0.85)
ax.grid(True, axis='y', alpha=0.4)
fig5.tight_layout()
fig5.savefig(os.path.join(OUTPUT_DIR, 'fig_rms.png'), dpi=150)
plt.close(fig5)

kp_values = [50, 100, 150, 200]
rms_kp = []

for kp_val in kp_values:
    Kp_test = np.diag([kp_val, kp_val])
    Kv_test = np.diag([2*np.sqrt(kp_val), 2*np.sqrt(kp_val)])
    q_t   = np.zeros((N, 2)); q_t[0]  = q_d[0]
    dq_t  = np.zeros((N, 2)); dq_t[0] = dq_d[0]
    for i in range(N - 1):
        tau_t = (mass_matrix(q_t[i]) @ (ddq_d[i] + Kv_test @ (dq_d[i] - dq_t[i]) + Kp_test @ (q_d[i] - q_t[i])) + coriolis_matrix(q_t[i], dq_t[i]) @ dq_t[i] + gravity_vector(q_t[i]))
        ddq_t = forward_dynamics(q_t[i], dq_t[i], tau_t)
        dq_t[i+1] = dq_t[i] + ddq_t * dt
        q_t[i+1]  = q_t[i] + dq_t[i] * dt + 0.5 * ddq_t * dt**2
    err_t = np.sqrt(np.mean((q_d - q_t)**2, axis=0))
    rms_kp.append(np.rad2deg(np.mean(err_t)))

fig6, ax = plt.subplots(figsize=(6, 4))
ax.plot(kp_values, rms_kp, 'o-', color='#2ca02c', lw=2, markersize=8)
ax.grid(True, alpha=0.4)
fig6.tight_layout()
fig6.savefig(os.path.join(OUTPUT_DIR, 'fig_gain_analysis.png'), dpi=150)
plt.close(fig6)

print(L1, L2, m1, m2)
print(Kp[0,0], Kv[0,0])
print(T_end)
print(np.rad2deg(rms1))
print(np.rad2deg(rms2))
print(np.max(np.abs(tau_hist[:,0])))
print(np.max(np.abs(tau_hist[:,1])))
