"""Cart-pole dynamics, derived by hand (see derivation write-up).

Conventions:
    p      cart position in the ground frame, +x to the right   [m]
    theta  pole angle from upright, counter-clockwise positive  [rad]
    l      pivot to the pole's centre of mass (half the rod)    [m]
"""
import numpy as np

M = 1.0     # cart mass                      [kg]
m = 0.1     # pole mass                      [kg]
l = 0.5     # pivot -> pole centre of mass   [m]
g = 9.81    # gravitational acceleration     [m/s^2]
b = 0.0     # cart viscous damping           [N/(m/s)]

I_S = m * l**2 / 3.0        # uniform slender rod, expressed in terms of l


def f(state, u, dF=0.0):
    """Return d(state)/dt for state = [p, p_dot, theta, theta_dot]."""
    p, p_dot, th, th_dot = state
    s, c = np.sin(th), np.cos(th)

    # A(theta) @ [p_ddot, th_ddot] = rhs
    A = np.array([
        [M + m,       -m * l * c     ],
        [-m * l * c,   I_S + m * l**2],
    ])
    rhs = np.array([
        u - b * p_dot - dF - m * l * th_dot**2 * s,
        m * g * l * s + dF * l * c,
    ])

    p_ddot, th_ddot = np.linalg.solve(A, rhs)
    return np.array([p_dot, p_ddot, th_dot, th_ddot])
