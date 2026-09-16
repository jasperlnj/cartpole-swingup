from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap
from matplotlib.animation import FuncAnimation
from dynamics import M, m, l, g, b, I_S, f
from integrate import rk4_step, rollout
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import os

x0 = np.zeros(4)  # Initial state
u0 = 0  # Initial force input

assert np.allclose(f(x0, u0), 0.0), "linearisation point is not an equilibrium!"

inv_Delta = 1/((M+m)*(I_S+m*l**2)-m**2*l**2)
A_analytic = np.array([[0, 1, 0, 0],
              [0, -b*(I_S+m*l**2)*inv_Delta, m**2*l**2*g*inv_Delta, 0],
              [0, 0, 0, 1],
              [0, -m*l*b*inv_Delta, m*g*l*(M+m)*inv_Delta, 0]])
B_analytic = np.array([0, (I_S+m*l**2)*inv_Delta, 0, m*l*inv_Delta])


def numerical_jacobian(f, x0, u0, eps=1e-6):
    """Compute the Jacobian of f at (x0, u0) numerically using finite differences."""
    A_num = np.zeros((4, 4))
    B_num = (f(x0, u0 + eps) - f(x0, u0 - eps)) / (2 * eps)

    for i in range(4):
        dx = np.zeros(4)
        dx[i] = eps
        A_num[:, i] = (f(x0 + dx, u0) - f(x0 - dx, u0)) / (2 * eps)
        #Central differences cancel out the first-order error term, leaving a truncation error
    return A_num, B_num
        
        

# using Brysons rule to set the state cost matrix Q
p_max, pdot_max = 0.5, 2.0      # m, m/s
th_max, thdot_max = 0.17, 1.0    # rad, rad/s
u_max = 10.0                     # N
Q = np.diag([1/p_max**2, 1/pdot_max**2, 1/th_max**2, 1/thdot_max**2])
R = np.array([[1/u_max**2]])

# solving Continuous-Time Algebraic Riccati Equation (CARE):
# Note: B has to be reshaped to be a 2D array for the matrix multiplication to work correctly
A = A_analytic
B = B_analytic.reshape(4, 1)
P = la.solve_continuous_are(A, B, Q, R)
# K = la.inv(R) @ B_analytic.reshape(4, 1).T @ P 
K = np.linalg.solve(R, B.T @ P) # shape is (1, 4), [[-20.0, -22.55, 129.38, 31.46]]
k = K[0] # convert to row


def make_controller(k, u_max=None):
    """Build a state-feedback controller u = -k·s, optionally force-limited."""
    def u_fn(t, s):
        u = -k @ s
        return u if u_max is None else np.clip(u, -u_max, u_max)
    return u_fn



#-------------------- EXPERIMENTS --------------------#

if __name__ == "__main__":
    A_num, B_num = numerical_jacobian(f, x0, u0)
    # check if the numerical Jacobian matches the analytical Jacobian
    print("max |A diff| =", np.max(np.abs(A_analytic - A_num)))
    print("max |B diff| =", np.max(np.abs(B_analytic - B_num)))

    lqr_free = make_controller(K[0])
    lqr_sat  = make_controller(K[0], u_max=u_max)


    print("open-loop eig  =", np.linalg.eigvals(A))
    #^^gives positive eigenvalues, so the system is unstable
    print("closed-loop eig =", np.linalg.eigvals(A - B @ K))
    #^^gives negative eigenvalues, so the system is stable
    print("K =", K.ravel())
    print("peak |u| =", np.abs(u).max(), "N")
    print("max |p|  =", np.abs(traj[:,0]).max(), "m")


    s0 = np.array([0.0, 0.0, 1.0, 0.0])   # 1.0 rad ≈ 57.3° off vertical, at rest
    dt = 0.005

    fallen = lambda s: np.abs(s[2]) > np.pi/2  # stop if the pole falls too far
    traj = rollout(f, rk4_step, s0, lqr_free, dt, 1200, stop_fn=fallen)    # 6 s

    #-------------------ANIMATION------------------------#
    frames = traj[::6]              
    Cart_W, Cart_H = 0.3, 0.15
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.set_xlim(-3.0, 1.5)               # wide enough for cart travel + pole
    ax.set_ylim(-1.2, 1.2)
    ax.set_aspect('equal')          # or the pole squashes as it rotates
    ax.axhline(0, color='0.7', lw=1)

    cart = Rectangle((0, -Cart_H/2), Cart_W, Cart_H, facecolor='0.4')
    ax.add_patch(cart)

    pole, = ax.plot([], [], lw=4, color='tab:blue')

    def update(i):
        p, p_dot, th, th_dot = frames[i]

        cart.set_xy((p - Cart_W/2, -Cart_H/2))    

        tip_x = p - 2*l * np.sin(th)         
        tip_y = 2*l * np.cos(th)
        pole.set_data([p, tip_x], [0, tip_y])

        return cart, pole

    anim = FuncAnimation(fig, update, frames=len(frames), interval=30, blit=True)
    os.makedirs("docs", exist_ok=True)
    print("writing to:", os.path.abspath("docs/1.0rad.gif"))
    anim.save("docs/1.0rad.gif", writer="pillow", fps=30)
    print("Animation saved to docs/1.0rad.gif")


    #------------------PLOT-----------------#
    t = np.arange(len(traj)) * dt
    u = np.array([lqr_free(0.0, s) for s in traj])     # recompute what the controller asked for

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))

    ax1.plot(t, traj[:, 0], label="p [m]")
    ax1.plot(t, traj[:, 1], label="ṗ [m/s]")
    ax1.plot(t, traj[:, 2], label="θ [rad]")
    ax1.plot(t, traj[:, 3], label="θ̇ [rad/s]")
    ax1.axhline(0, color="0.7", lw=0.8)
    ax1.set_ylabel("state")
    ax1.legend(ncol=4, fontsize=9)
    ax1.set_title(f"LQR from θ₀ = {s0[2]} rad, unsaturated")

    ax2.plot(t, u, color="tab:red", label="u [N]")
    ax2.axhline( u_max, color="0.5", ls="--", lw=1)
    ax2.axhline(-u_max, color="0.5", ls="--", lw=1, label=f"±{u_max:.0f} N budget")
    ax2.set_ylabel("force [N]")
    ax2.set_xlabel("time [s]")
    ax2.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig("docs/lqr_states.png", dpi=150)


    #----------------plot the region of attraction----------------#
    CAUGHT_COLOR, FELL_COLOR = "#1B4D3E", "#ece7e1"


    def scan_roa(controller, theta0_grid, thetadot0_grid, n=800, cache=None):
        """Grid-scan initial (θ, θ̇). True where the controller recovers the pole."""
        if cache and os.path.exists(cache):
            print(f"loaded cached scan: {cache}")
            return np.load(cache)

        caught = np.zeros((len(thetadot0_grid), len(theta0_grid)), dtype=bool)
        for j, thd in enumerate(thetadot0_grid):
            for i, th in enumerate(theta0_grid):
                tr = rollout(f, rk4_step, np.array([0.0, 0.0, th, thd]), controller, dt, n)
                caught[j, i] = (np.isfinite(tr[-1]).all()
                                and abs(tr[-1, 2]) < 0.02      # ended upright
                                and abs(tr[-1, 3]) < 0.05)     # and stopped there
            print(f"\r  row {j+1}/{len(thetadot0_grid)}", end="", flush=True)
        print()

        if cache:
            np.save(cache, caught)
        return caught


    def boundary_angle(caught, theta0_grid, thetadot0_grid):
        """Largest initial angle recovered from rest (θ̇₀ = 0), in radians."""
        row = caught[np.argmin(np.abs(thetadot0_grid))]
        pos = theta0_grid[row & (theta0_grid > 0)]
        return pos.max() if pos.size else 0.0


    def plot_roa(caught, theta0_grid, thetadot0_grid, title, path):
        """Binary region-of-attraction map in the θ–θ̇ plane."""
        th_b = boundary_angle(caught, theta0_grid, thetadot0_grid)

        fig, ax = plt.subplots(figsize=(6.5, 5))
        ax.pcolormesh(theta0_grid, thetadot0_grid, caught,
                    cmap=ListedColormap([FELL_COLOR, CAUGHT_COLOR]), shading="nearest")
        ax.axhline(0, color="0.45", lw=0.8)
        ax.axvline(0, color="0.45", lw=0.8)

        if th_b > 0:
            ax.plot(th_b, 0, "o", color="#161616", ms=5, zorder=3)
            ax.annotate(f"{np.degrees(th_b):.0f}°", (th_b, 0), xytext=(7, 7),
                        textcoords="offset points", color="#161616", fontsize=10)

        ax.set_xlabel("initial angle  θ₀  [rad]")
        ax.set_ylabel("initial angular velocity  θ̇₀  [rad/s]")
        ax.set_title(title)
        ax.legend(handles=[Patch(facecolor=CAUGHT_COLOR, label="recovered"),
                        Patch(facecolor=FELL_COLOR,  label="fell")],
                loc="upper right", framealpha=0.95, fontsize=9)

        fig.tight_layout()
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"wrote {path}")

    os.makedirs("docs", exist_ok=True)
        
    theta0_grid    = np.linspace(-1.2, 1.2, 41)     # rad
    thetadot0_grid = np.linspace(-4.0, 4.0, 41)     # rad/s

    print("scanning saturated…")
    roa_sat  = scan_roa(lqr_sat,  theta0_grid, thetadot0_grid, cache="docs/roa_sat.npy")
    print("scanning unsaturated…")
    roa_free = scan_roa(lqr_free, theta0_grid, thetadot0_grid, cache="docs/roa_free.npy")

    plot_roa(roa_sat,  theta0_grid, thetadot0_grid,
             f"Region of attraction — LQR,  |u| ≤ {u_max:.0f} N", "docs/roa_saturated.png")
    plot_roa(roa_free, theta0_grid, thetadot0_grid,
             "Region of attraction — LQR,  unlimited force", "docs/roa_free.png")

    print(f"boundary from rest —  saturated: "
          f"{np.degrees(boundary_angle(roa_sat,  theta0_grid, thetadot0_grid)):.1f}°   "
          f"unlimited: {np.degrees(boundary_angle(roa_free, theta0_grid, thetadot0_grid)):.1f}°")