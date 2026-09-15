swingup.gif




In this project I built a cart-pole, based on a pole, jointed to a cart that has freedom of movement in one dimension. The goal was to bring the pole in an upright position from any starting position by using the cart's ability to accelerate in one dimension. I achieved this by deriving the equations of motion by hand first from free-body diagrams and then using those to first swing the pole up (if necessary) and then stabilizing it using a linear quadratic regressor (LQR). 




Results:

ROA free vs. ROA sat

I started off entirely ignoring the physical limitations of the force output an electric motor might have. This led to a much wider region of attraction in which the LQR was able to stabilize. In ROA sat the region of attraction is much smaller, however more realistic, as it assumes that the force produced by the motor can be 10 N at most.




## 1: equations deriving - pdf link

I derive the following equations based on this free-body diagram.

Further details can be found in -pdf-. 




## 2: Simulating 

energy.py calculates the system's energy independently of integrate.py. Since there is no damping and no applied force, it must be constant, which makes it a check that couldn't be satisfied by accident. Using this check, I was able to find multiple mistakes later on, including a sign error on the cart's contribution to the pole velocity.




[docs/energy_drift.png]




Further on, I was able to compare the Euler's step to the Runge-Kutta 4 step (RK4) and show that the latter one is much more precise.




## 3: linearizing

I need to linearize the system around the poles' upright equilibrium to use a linear quadratic regressor (LQR) later on. To do so, I first derived the state matrix A and the input matrix by hand and then checked it against their numerical Jacobian counterpart, which I calculated using finite differences in LQR.py.




//mention Eigenvalues




## 4: LQR

I derived the state weight matrix Q and the effort weight R by using Bryson's Rule to then figure out the cost functional, or rather the cost-to-go matrix, by letting scipy.linalg solve the continuous algebraic Riccati equation (CARE) for me. Finally, this allowed me to figure out the gain matrix K, which I used as a controller for the real nonlinear system. 

//mention closed-loop matrix?




## 5: Where LQR stops working

When I started implementing the LQR i was immensely surprised to see that it still worked at angles up to 55°. I then realized that the force exerted by the motor was unrealisticly high. This led me to limit the possibly exerted force to 10 N, to get more realistic results.

I then plotted the region of attraction for both the saturated and unsaturated/free systems to see when my linearized model would stop working.

It is important to note that this only takes into account the region of attraction when the cart's acceleration is zero. 

But even at angles up to +-28° (with no initial angular speed), the LQR controller was able to stabilize the pole, which shows that even though it's a linearized approach, it works quite well for the real, not linear system. 




## 6: Swing up

To solve the issue of the LQR only latching at a limited range of angles, I had to distinguish between the states in which I was able to use the LQR and the ones where I had to swing the pole upwards. I started off just testing whether the angle and angular velocity of the current state were within a small part of the ROA. I used the rectangle set approximately by $-0.25 rad <= theta <= 0.25 rad$ and $-1 rad/s <= theta_dot <= 1 rad/s$ to be exact. However, even though this hardly left any margin of error, it was still a small area, which meant a lot of points that the LQR was able to solve were ignored. I considered using the cost to go xTPx, which sees all four states and is invariant, but ended up using a fitted strip instead, which is tighter on the slice I measured and works sufficiently while being simpler to implement. It is important to note, though, that this approach only takes into account two of the four states. 

Initially I also made the mistake of continuously checking whether the current state was within the boundaries of the ROA. This only needed checking once, though, since it was possible for the angular speed to leave the ROA during the process of the LQR. I fixed this by handing the control over to the LQR entirely once it was in the ROA just once.




To get the pole to land in that ROA when it wasn't initially there, I used the so-called separatrix. You can divide the energy level of the system into three categories: either it has too much energy, e.g., the pendulum spins continuously clockwise around the cart, or it has too little energy, e.g., it swings from left to right slightly. The separatrix sits right in between and is the energy level of the pendulum when it sits upright and still. 

I calculated the difference in energy $E_err$ bewteen the energy of the pole in the current state and the Seperatrix. To figure out how to get closer towards E_upright, I looked at the derivative of E, E_dot = θ˙(Jθ¨−mglsinθ). I then substituted using the pole equation (row two, the mass matrix): Jθ̈ − mgl sin θ = ml cos θ · p̈

which left me with $E_dot = m*l*p_ddot theta_dotcos(theta)$

Since $m*l*theta_dot*cos(theta)$ is a constant known at any time, this controller is linear in $p_ddot$, which means that i can influence the energy of the pendulum using the carts acceleration.

dE_err/dt should have the opposite sign to E_err. 

By using $p_ddot = -k*E_err*theta_dot*cos(theta)$ $p_ddot$ always has the correct sign, because the only possibly non positive variables outside of E_err are always squared:

$E_err_dot = -k_energy*m*l*E_Err(theta_dot*cos(theta))^2$

using this i was able to produce the control law $u = −k·E_err·θ̇·cos θ$.

It is possible to tweak the control law by tuning k_energy.




## Whats next:

i want to continue to work on this project. The next thing i want to do is build a Reinforcement Learning environment, train it and then compare it to the LQR approach. 






## other:

energy.py derivation: 

Using the law of conservation of energy, I was able to predict the state of the system over time based on the initial state. By state I mean the position of the cart $p$, the acceleration of the cart $p_dot$, the angle $theta$ [rad], and its rate of change $theta_dot$ [rad/s]. 

The energy of the system consists of the kinetic energy of the cart, which is $frac{1/2}  M  p_dot^2$ (the cart does not carry potential energy as its height is fixed).

As well as the energy of the pole, which consists of its kinetic energy, which I was able to figure out using König's theorem: $0.5  m  ((l  np.sin(th)th_dot)**2 + (-th_dot  l  np.cos(th) +p_dot)**2) + 0.5  I_S  th_dot**2$


