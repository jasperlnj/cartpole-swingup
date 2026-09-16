"""Fixed-step integrators."""
import numpy as np


# Here we use the Runge Katta 4th Order method to numerically determine
# the future states of the system
def rk4_step(f, state, u, dt, **kw):
    k1 = f(state,             u, **kw)
    k2 = f(state + dt/2 * k1, u, **kw)
    k3 = f(state + dt/2 * k2, u, **kw)
    k4 = f(state + dt   * k3, u, **kw)
    return state + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

# this is the basic Euler Step which is way less accurate than RK4 (creates energgy)
def euler_step(f, state, u, dt, **kw):
    return state + dt * f(state, u, **kw)


def rollout(f, step, state, u_fn, dt, n, stop_fn = None):
    """Integrate n steps; u_fn(t, state) returns the control force."""
    traj = np.empty((n + 1, len(state))) #n+1 rows, 4 columns
    traj[0] = state
    for i in range(n):
        state = step(f, state, u_fn(i * dt, state), dt)
        traj[i + 1] = state
        if not np.isfinite(state).all():  # if the pole falls too far, stop the simulation
            return traj[:i + 2]
        if stop_fn is not None and stop_fn(state):
            return traj[:i + 2]
    return traj
