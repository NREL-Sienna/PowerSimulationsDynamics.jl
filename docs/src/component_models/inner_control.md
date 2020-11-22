# Inner Loop Controls

This component defines voltage and current controllers to generate the reference signal
for the converter. Although in many controls the current and voltage control are separate blocks
we propose a more general control approach that considers them as a joint control logic.

## Integrated Virtual Impedance, Voltage and Current Controller ```[CurrentControl]```

The following model receives both the outer loop control frequency and reference voltage
signal to generate the reference signal for the converters. The virtual impedance plays a
similar role of the impedance of a synchronous generator. A PI voltage controller is used
to generate the current signal that is used in the PI current controller to finally generate
the voltage reference signal for the converters.

```math
\begin{align*}
    \dot{\xi}_d &= v_{d,\text{vi}}^{\text{ref}} - v_d \tag{3a} \\
    \dot{\xi}_q &= v_{q,\text{vi}}^{\text{ref}} - v_q \tag{3b} \\
    \dot{\gamma}_d &= i_{d,\text{cv}}^{\text{ref}} - i_{d,\text{cv}} \tag{3c} \\
    \dot{\gamma}_q &= i_{q,\text{cv}}^{\text{ref}} - i_{q,\text{cv}} \tag{3d} \\
    \dot{\phi}_d &= \omega_{\text{ad}}(v_d - \phi_d) \tag{3e} \\
    \dot{\phi}_q &= \omega_{\text{ad}}(v_q - \phi_q) \tag{3f}
\end{align*}
```

with

```math
\begin{align}
    v_{d,\text{vi}}^{\text{ref}} &= v_{\text{olc}}^{\text{ref}} - r_v i_d + \omega_{\text{olc}} l_v i_q \tag{3g} \\
    v_{q,\text{vi}}^{\text{ref}} &= - r_v i_q - \omega_{\text{olc}} l_v i_d \tag{3h} \\
    i_{d,\text{cv}}^{\text{ref}} &= k_{pv}\left(v_{d,\text{vi}}^{\text{ref}} - v_q\right) + k_{iv} \xi_d - c_f \omega_{\text{olc}} v_q + k_{\text{ffi}}i_d \tag{3i} \\
    i_{q,\text{cv}}^{\text{ref}} &= k_{pv}\left(v_{q,\text{vi}}^{\text{ref}} - v_q\right) + k_{iv} \xi_q + c_f \omega_{\text{olc}} v_d + k_{\text{ffi}}i_q \tag{3j} \\
    v_d^{\text{ref-signal}} &= k_{pc} \left(i_{d,\text{cv}}^{\text{ref}} - i_{d,\text{cv}}\right) + k_{ic} \gamma_d - \omega_{\text{olc}} l_f i_{q,\text{cv}} + k_{\text{ffv}}v_d - k_{\text{ad}}(v_d - \phi_d) \tag{3k} \\
    v_q^{\text{ref-signal}} &= k_{pc} \left(i_{q,\text{cv}}^{\text{ref}} - i_{q,\text{cv}}\right) + k_{ic} \gamma_q + \omega_{\text{olc}} l_f i_{d,\text{cv}} + k_{\text{ffv}}v_q - k_{\text{ad}}(v_q - \phi_q) \tag{3l}
\end{align}
```

In here the transformation to the ``dq`` reference frame is using the outer-loop reference angle as:

```math
\begin{align*}
v_d + jv_q = (v_r + jv_i)e^{-j\delta\theta_{olc}} \\
i_d + ji_q = (i_r + ji_i)e^{-j\delta\theta_{olc}}
\end{align*}
```

that again ``v_r + jv_i`` could be in the capacitor or the last branch of the filter (i.e.
the point of common coupling). For LCL filters it is considered in the capacitor. In the case
of the converter, the transformation is directly

```math
\begin{align*}
i_{d,\text{cv}} + ji_{q,\text{cv}} = (i_{r,\text{cv}} + ji_{i,\text{cv}})e^{-j\delta\theta_{olc}}
\end{align*}
```
