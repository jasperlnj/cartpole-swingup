### Energy-shaping swing-up

With `E` the pendulum's own energy (pivot treated as fixed) and `J = I_S + ml²`:

$$E = \tfrac12 J\dot\theta^2 + mgl(1+\cos\theta), \qquad E_{up} = 2mgl$$

Differentiating with respect to time:

$$\dot E = J\dot\theta\ddot\theta - mgl\sin\theta\,\dot\theta = \dot\theta\,(J\ddot\theta - mgl\sin\theta)$$

Substituting the pole equation (row two of the mass-matrix system), `J θ̈ − mgl sin θ = ml cos θ · p̈`:

$$\dot E = m\,l\,\ddot p\,\dot\theta\cos\theta$$

With `Ẽ = E − E_up`, the choice `p̈ = −k Ẽ θ̇ cos θ` gives

$$\dot{\tilde E} = -k\,m\,l\,\tilde E\,(\dot\theta\cos\theta)^2$$

and with `V = ½Ẽ²`,

$$\dot V = \tilde E\,\dot{\tilde E} = -k\,m\,l\,\tilde E^2(\dot\theta\cos\theta)^2 \le 0$$

so `Ẽ → 0` wherever `θ̇ cos θ ≠ 0`. The exception is `θ = π, θ̇ = 0` exactly where the law produces zero force forever, which is why simulations start from `θ₀ = π − 0.05`.

