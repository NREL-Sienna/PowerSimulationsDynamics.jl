
## Prime Movers and Turbine Governors (TG)

This section describes how mechanical power is modified to provide primary frequency control with synchronous generators. It is assumed that ``\tau_{\text{ref}} = P_{\text{ref}}`` since they are decided at nominal frequency ``\omega = 1``.

### Fixed TG ```[TGFixed]```

This a simple model that set the mechanical torque to be equal to a proportion of the desired reference ``\tau_m = \eta P_{\text{ref}}``. To set the mechanical torque to be equal to the desired power, the value of ``\eta`` is set to 1.

### TG Type I ```[TGTypeI]```

This turbine governor is described by a droop controller and a low-pass filter to model the governor and two lead-lag blocks to model the servo and reheat of the turbine governor.

```math
\begin{align}
\dot{x}_{g1} &= \frac{1}{T_s}(p_{\text{in}} - x_{g1}) \tag{13a} \\
\dot{x}_{g2} &= \frac{1}{T_c} \left[ \left(1- \frac{t_3}{T_c}\right)x_{g1} - x_{g2} \right] \tag{13b} \\
\dot{x}_{g3} &= \frac{1}{T_5} \left[\left(1 - \frac{T_4}{T_5}\right)\left(x_{g2} + \frac{T_3}{T_c}x_{g1}\right) - x_{g3}  \right] \tag{13c} \\
\tau_m &= x_{g3} + \frac{T_4}{T_5}\left(x_{g2} + \frac{T_3}{T_c}x_{g1}\right) \tag{13d}
\end{align}
```

with

```math
\begin{align*}
p_{\text{in}} = P_{\text{ref}} + \frac{1}{R}(\omega_s - 1)
\end{align*}
```

### TG Type II ```[TGTypeII]```

This turbine governor is a simplified model of the Type I.

```math
\begin{align}
\dot{x}_g &= \frac{1}{T_2}\left[\frac{1}{R}\left(1 - \frac{T_1}{T_2}\right) (\omega_s - \omega) - x_g\right] \tag{14a} \\
\tau_m &= P_{\text{ref}} + \frac{1}{R}\frac{T_1}{T_2}(\omega_s - \omega) \tag{14b}
\end{align}
```
