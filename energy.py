import numpy as np
import matplotlib.pyplot as plt
from dynamics import M, m, l, g, I_S, f
from integrate import rk4_step, euler_step, rollout


def kinetic_energy(state):
    """Return the kinetic energy of the cart-pole system."""
    p, p_dot, th, th_dot = state[0], state[1], state[2], state[3]
    energy_cart = 0.5 * M * p_dot**2
    #using Königs Theorem to compute the kinetic energy of the pole
    energy_pole = 0.5 * m * ((l * np.sin(th)*th_dot)**2 + (-th_dot * l * np.cos(th) +p_dot)**2) + 0.5 * I_S * th_dot**2
    return energy_cart + energy_pole

def potential_energy(state):
    """Return the potential energy of the cart-pole system."""
    th = state[2]
    return m * g * l * (1 + np.cos(th)) 

def total_energy(state):
    return kinetic_energy(state) + potential_energy(state)

#state 0: only angle is at 0.6, everything else is 0
s0 = np.array([0.0, 0.0, 0.6, 0.0])
#dt = 0,01, n = 3000, u = 0.0
dt1 = 0.01
dt2 = 0.005
T = 30.0
n1 = int(T / dt1)
n2 = int(T / dt2)
#important: the control force is 0.0, so the system is not actuated, and we can see how the energy evolves over time
traj_rk4_dt1 = rollout(f, rk4_step, s0, lambda t, s: 0.0, dt1, n1)
traj_euler_dt1 = rollout(f, euler_step, s0, lambda t, s: 0.0, dt1, n1)
traj_rk4_dt2 = rollout(f, rk4_step, s0, lambda t, s: 0.0, dt2, n2)
traj_euler_dt2 = rollout(f, euler_step, s0, lambda t, s: 0.0, dt2, n2)

E_rk4_dt1 = np.array([total_energy(s) for s in traj_rk4_dt1])
E_euler_dt1 = np.array([total_energy(s) for s in traj_euler_dt1])
E_rk4_dt2 = np.array([total_energy(s) for s in traj_rk4_dt2])
E_euler_dt2 = np.array([total_energy(s) for s in traj_euler_dt2])

print("start:", E_rk4_dt1[0], "  end:", E_rk4_dt1[-1], "  drift:", E_rk4_dt1[-1] - E_rk4_dt1[0])
print("start:", E_euler_dt1[0], "  end:", E_euler_dt1[-1], "  drift:", E_euler_dt1[-1] - E_euler_dt1[0])
print("start:", E_rk4_dt2[0], "  end:", E_rk4_dt2[-1], "  drift:", E_rk4_dt2[-1] - E_rk4_dt2[0])
print("start:", E_euler_dt2[0], "  end:", E_euler_dt2[-1], "  drift:", E_euler_dt2[-1] - E_euler_dt2[0])

#Plotting the energy over time for both integrators and both time steps
t1 = np.arange(len(E_rk4_dt1)) * dt1
t2 = np.arange(len(E_rk4_dt2)) * dt2
plt.figure(figsize=(12, 6))
plt.semilogy(t1, np.abs(E_rk4_dt1 - E_rk4_dt1[0]), label='RK4 dt=0.01', color='blue')
plt.semilogy(t2, np.abs(E_rk4_dt2 - E_rk4_dt2[0]), label='RK4 dt=0.005', color='darkblue')
plt.semilogy(t1, np.abs(E_euler_dt1 - E_euler_dt1[0]), label='Euler dt=0.01', color='green')
plt.semilogy(t2, np.abs(E_euler_dt2 - E_euler_dt2[0]), label='Euler dt=0.005', color='darkgreen')
plt.xlabel("time, t [s]")
plt.ylabel("|E(t) - E(0)|  [J]")
plt.title("Total Energy drift of the Cart-Pole System Over Time")
plt.legend()
plt.grid(True, alpha=0.25)
#plt.savefig("docs/energy_drift.png", dpi=150, bbox_inches="tight")
plt.show()