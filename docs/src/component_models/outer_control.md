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
    \dot{\delta\omega}_{\text{olc}} &= \frac{p_{\text{ref}}}{T_a} - \frac{p_e}{T_a} - \frac{k_d(\omega_{\text{olc}} - \omega_{\text{pll}})}{T_a} - \frac{k_\omega(\omega_{\text{olc}} - \omega_{\text{ref}})}{T_a} \tag{1a} \\
    \dot{\theta}_{\text{olc}} &= \Omega_b \delta\omega_{\text{olc}} \tag{1b} \\
    \dot{q}_m &= \omega_f (q_e - q_m) \tag{1c}
\end{align}
```

with

```math
\begin{align}
    p_e &= v_ri_r + v_ii_i \tag{1d} \\
    q_e &= v_ii_r - v_ri_i \tag{1e} \\
    v_{\text{olc}}^{\text{ref}} &= v_{\text{ref}} + k_q(q_{\text{ref}} - q_m) \tag{1f}
\end{align}
```

In this case, the measurement of power are being done in the capacitor of the LCL filter.
However, depending on the model, this measurements could be different depending on where
is the point of common coupling.


## Active Power Droop (P-droop) and Q-droop ```[OuterControl]```

The following model represent a ``P\text{-}f`` droop model to represent how active
power is going to be deployed. The constructor is ```OuterControl{ActivePowerControl, ReactivePowerDroop}```.
It defines a new SRF denoted as ``\theta_{\text{olc}}`` for the active power controller and
uses a simple voltage droop for dispatching reactive power. Both active and reactive power are measured via a low-pass filter:

```math
\begin{align}
    \dot{\theta}_{\text{olc}} &= \Omega_b (\omega_{\text{olc}} - \omega_{\text{sys}}) \tag{2a} \\
    \dot{p}_m &= \omega_z (p_e - p_m) \tag{2b} \\
    \dot{q}_m &= \omega_f (q_e - q_m) \tag{2c}
\end{align}
```

with

```math
\begin{align}
    p_e &= v_ri_r + v_ii_i \tag{2d} \\
    q_e &= v_ii_r - v_ri_i \tag{2e} \\
    \omega_{\text{olc}} &= \omega_{\text{ref}} + R_p (p_{\text{ref}} - p_e) \tag{2f} \\
    v_{\text{olc}}^{\text{ref}} &= v_{\text{ref}} + k_q(q_{\text{ref}} - q_m) \tag{2g}
\end{align}
```

In this case, the measurement of power are being done in the capacitor of the LCL filter.
However, depending on the model, this measurements could be different depending on where
is the point of common coupling.


## Active and Reactive Power PI Controllers (Grid Following) ```[OuterControl]```

The following model represents a PI controller for both active and reactive power to generate
the current references that will be used in the Current Controller of the inner control
```CurrentModeControl```. The constructor is ```OuterControl{ActivePowerPI, ReactivePowerPI}```.
The equations are:

```math
\begin{align}
    \dot{\sigma}_{p} &= p_\text{ref} - p_m \tag{3a} \\
    \dot{p}_m &= \omega_z (p_e - p_m) \tag{3b} \\
    \dot{\sigma}_{q} &= q_\text{ref} - q_m \tag{3c} \\
    \dot{q}_m &= \omega_f (q_e - p_m) \tag{3d} \\
\end{align}
```

with

```math
\begin{align}
    p_e &= v_ri_r + v_ii_i \tag{3e} \\
    q_e &= v_ii_r - v_ri_i \tag{3f} \\
    \omega_{\text{olc}} &= \omega_{\text{pll}} \tag{3g} \\
    \theta_{\text{olc}} &= \theta_{\text{pll}} \tag{3h} \\
    i_\text{d,cv}^\text{ref} &= k_p^q (q_\text{ref} - q_m) + k_i^q \sigma_q \tag{3i} \\
    i_\text{q,cv}^\text{ref} &= k_p^p (p_\text{ref} - p_m) + k_i^p \sigma_p \tag{3j} \\
\end{align}
```

This models requires a PLL to have a SRF for an internal ``dq`` reference frame. Contrary
to the Grid-Forming model, it cannot work without a PLL. Since this Outer Control outputs
a current reference, it can only be used with a current mode inner control (i.e. that receives 
a current reference instead of a voltage reference).