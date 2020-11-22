# Outer Loop Controls

This component defines controllers for both active and reactive power. The joint design is based
on the fact that many novel control techniques can be based on joint control of active and reactive
power.

## Virtual Inertia and Q-droop ```[OuterControl]```

The following model represent a virtual synchronous machine model to represent how active
power is going to be deployed. The constructor is ```OuterControl{VirtualInertia, ReactivePowerDroop}```.
It defines a new SRF denoted as ``\theta_{\text{olc}}`` for the active power controller and
uses a simple voltage droop for dispatching reactive power:

```math
\begin{align}
    \dot{\delta\omega}_{\text{olc}} &= \frac{p_{\text{ref}}}{T_a} - \frac{p_e}{T_a} - \frac{k_d(\omega_{\text{olc}} - \omega_{\text{pll}})}{T_a} - \frac{k_\omega(\omega_{\text{olc}} - \omega_{\text{ref}})}{T_a} \tag{2a} \\
    \dot{\delta\theta}_{\text{olc}} &= \Omega_b \delta\omega_{\text{olc}} \tag{2b} \\
    \dot{q}_m &= \omega_f (q_e - q_m) \tag{2c}
\end{align}
```

with

```math
\begin{align}
    p_e &= v_ri_r + v_ii_i \tag{2d} \\
    q_e &= v_ii_r - v_ri_i \tag{2e} \\
    v_{\text{olc}}^{\text{ref}} &= v_{\text{ref}} + k_q(q_{\text{ref}} - q_m) \tag{2f}
\end{align}
```

In this case, the measurement of power are being done in the capacitor of the LCL filter.
However, depending on the model, this measurements could be different depending on where
is the point of common coupling.
