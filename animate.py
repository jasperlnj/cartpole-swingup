import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.animation import FuncAnimation
from dynamics import f, l
from integrate import rk4_step, rollout

dt = 0.01
s0 = np.array([0.0, 0.0, 0.6, 0.0])
traj = rollout(f, rk4_step, s0, lambda t, s: 0.0, dt, 800)
frames = traj[::3]              
Cart_W, Cart_H = 0.3, 0.15
fig, ax = plt.subplots(figsize=(7, 4))
ax.set_xlim(-1.3, 1.3)               # wide enough for cart travel + pole
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
import os
os.makedirs("docs", exist_ok=True)
print("writing to:", os.path.abspath("docs/passive.gif"))
anim.save("docs/passive.gif", writer="pillow", fps=30)
print("Animation saved to docs/passive.gif")