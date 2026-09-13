from LQR import u_max, k
from dynamics import f, m, M, l, g, I_S
from integrate import rk4_step, rollout
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.animation import FuncAnimation


def swing_up(t, state): # needs two arguments, time and state, to be compatible with rollout
    x, x_dot, theta, theta_dot = state
    theta_wrapped = (theta + np.pi) % (2 * np.pi) - np.pi
    k_energy = 2.0
    
    # check if inside of LQR region of attraction:
    
    if np.abs(4.37*theta_wrapped + theta_dot) < 1.3 and np.abs(theta_dot) < 4: # linear region of attraction with upper cap, to prevent imprecisions 
        # Hand off to linear LQR
        x_linear = np.array([x, x_dot, theta_wrapped, theta_dot])
        return float(np.clip(-k @ x_linear, -u_max, u_max)) 

    # Swing up:
    # E = E_kin + E_pot
    E = 0.5 * (I_S + m * l**2) * (theta_dot**2) + m * g * l * (np.cos(theta_wrapped) + 1)
    E_upright = 2 * m * g * l
    E_err = E - E_upright

    u = -k_energy * E_err * theta_dot * np.cos(theta)
    return float(np.clip(u, -u_max, u_max))


if __name__ == "__main__":
    os.makedirs("docs", exist_ok=True)

    dt   = 0.005
    n    = 3000                                    # 15 s
    s0   = np.array([0.0, 0.0, np.pi - 0.05, 0.0]) # hängend, angestupst
    traj = rollout(f, rk4_step, s0, swing_up, dt, n)
    t    = np.arange(len(traj)) * dt

    wrap = lambda th: (th + np.pi) % (2*np.pi) - np.pi
    J    = I_S + m*l**2
    E_UP = 2*m*g*l
    E_t  = 0.5*J*traj[:,3]**2 + m*g*l*(np.cos(traj[:,2]) + 1)

    # ---------- Energieplot ----------
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
    ax1.plot(t, E_t, color="#1B4D3E", label="Pendelenergie")
    ax1.axhline(E_UP, ls="--", color="0.4", label=f"E_up = {E_UP:.3f} J (Separatrix)")
    ax1.set_ylabel("Energie [J]"); ax1.legend(fontsize=9)
    ax1.set_title("Energy shaping: Aufschwingen und Übergabe an den LQR")

    ax2.plot(t, wrap(traj[:,2]), label="θ (gewickelt) [rad]")
    ax2.plot(t, traj[:,0], label="p [m]")
    ax2.axhline(0, color="0.7", lw=0.8)
    ax2.set_xlabel("Zeit [s]"); ax2.legend(fontsize=9)
    fig.tight_layout(); fig.savefig("docs/swingup_energy.png", dpi=150); plt.close(fig)
    t_switch =  # first time when LQR takes over
    for ax in (ax1, ax2):
        ax.axvline(t_switch, color="#b0452d", lw=1, ls=":")
        ax1.text(t_switch, 0.15, " LQR übernimmt", color="#b0452d", fontsize=9)

    # ---------- GIF ----------
    frames = traj[::10]                            # 0.05 s Abstand -> fps=20 
    CART_W, CART_H = 0.3, 0.15

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.set_xlim(-2.5, 2.5)                         # Wagen läuft bis ca 1 m, Spitze bis ca 2 m
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect("equal")
    ax.axhline(0, color="0.7", lw=1)

    cart = Rectangle((0, -CART_H/2), CART_W, CART_H, facecolor="0.4")
    ax.add_patch(cart)
    pole, = ax.plot([], [], lw=4, color="#1B4D3E")

    def update(i):
        p, _, th, _ = frames[i]
        cart.set_xy((p - CART_W/2, -CART_H/2))
        pole.set_data([p, p - 2*l*np.sin(th)], [0, 2*l*np.cos(th)])
        return cart, pole

    anim = FuncAnimation(fig, update, frames=len(frames), interval=50, blit=True)
    anim.save("docs/swingup.gif", writer="pillow", fps=20)
    plt.close(fig)
    print("docs/swingup_energy.png und docs/swingup.gif geschrieben")