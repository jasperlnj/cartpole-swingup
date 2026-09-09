"""Fixed-step integrators. Euler is kept for the Day 2 comparison."""
import numpy as np


def rk4_step(f, state, u, dt, **kw):
    k1 = f(state,             u, **kw)
    k2 = f(state + dt/2 * k1, u, **kw)
    k3 = f(state + dt/2 * k2, u, **kw)
    k4 = f(state + dt   * k3, u, **kw)
    return state + dt/6 * (k1 + 2*k2 + 2*k3 + k4)


def euler_step(f, state, u, dt, **kw):
    return state + dt * f(state, u, **kw)


def rollout(f, step, state, u_fn, dt, n):
    """Integrate n steps; u_fn(t, state) returns the control force."""
    traj = np.empty((n + 1, len(state)))
    traj[0] = state
    for i in range(n):
        state = step(f, state, u_fn(i * dt, state), dt)
        traj[i + 1] = state
    return traj
