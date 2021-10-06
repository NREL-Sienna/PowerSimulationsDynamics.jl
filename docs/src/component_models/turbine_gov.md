
# Prime Movers and Turbine Governors (TG)

This section describes how mechanical power is modified to provide primary frequency control with synchronous generators. It is assumed that ``\tau_{\text{ref}} = P_{\text{ref}}`` since they are decided at nominal frequency ``\omega = 1``.

## Fixed TG ```[TGFixed]```

This a simple model that set the mechanical torque to be equal to a proportion of the desired reference ``\tau_m = \eta P_{\text{ref}}``. To set the mechanical torque to be equal to the desired power, the value of ``\eta`` is set to 1.

## TG Type I ```[TGTypeI]```

This turbine governor is described by a droop controller and a low-pass filter to model the governor and two lead-lag blocks to model the servo and reheat of the turbine governor.

```math
\begin{align}
\dot{x}_{g1} &= \frac{1}{T_s}(p_{\text{in}} - x_{g1}) \tag{1a} \\
\dot{x}_{g2} &= \frac{1}{T_c} \left[ \left(1- \frac{t_3}{T_c}\right)x_{g1} - x_{g2} \right] \tag{1b} \\
\dot{x}_{g3} &= \frac{1}{T_5} \left[\left(1 - \frac{T_4}{T_5}\right)\left(x_{g2} + \frac{T_3}{T_c}x_{g1}\right) - x_{g3}  \right] \tag{1c} \\
\tau_m &= x_{g3} + \frac{T_4}{T_5}\left(x_{g2} + \frac{T_3}{T_c}x_{g1}\right) \tag{1d}
\end{align}
```

with

```math
\begin{align*}
p_{\text{in}} = P_{\text{ref}} + \frac{1}{R}(\omega_s - 1.0)
\end{align*}
```

## TG Type II ```[TGTypeII]```

This turbine governor is a simplified model of the Type I.

```math
\begin{align}
\dot{x}_g &= \frac{1}{T_2}\left[\frac{1}{R}\left(1 - \frac{T_1}{T_2}\right) (\omega_s - \omega) - x_g\right] \tag{2a} \\
\tau_m &= P_{\text{ref}} + \frac{1}{R}\frac{T_1}{T_2}(\omega_s - \omega) \tag{2b}
\end{align}
```

## TGOV1 ```[SteamTurbineGov1]```

This represents a classical Steam-Turbine Governor, known as TGOV1.

```math
\begin{align}
\dot{x}_{g1} &= \frac{1}{T_1} (\text{ref}_{in} - x_{g1}) \tag{3a} \\
\dot{x}_{g2} &= \frac{1}{T_3} \left(x_{g1}^\text{sat} \left(1 - \frac{T_2}{T_3}\right) - x_{g2}\right) \tag{3b}
\end{align}
```

with

```math
\begin{align}
\text{ref}_{in} &= \frac{1}{R} (P_{ref} - (\omega - 1.0)) \tag{3c} \\
x_{g1}^\text{sat} &= \left\{ \begin{array}{cl}
                        x_{g1} & \text{ if } V_{min} \le x_{g1} \le V_{max}\\
                        V_{max} & \text{ if } x_{g1} > V_{max} \\
                        V_{min} & \text{ if } x_{g1} < V_{min}
                    \end{array} \right. \tag{3d} \\
\tau_m &= x_{g2} + \frac{T_2}{T_3} x_{g1} - D_T(\omega - 1.0) \tag{3e}
\end{align}
```

## GAST ```[GasTG]```

This turbine governor represents the Gas Turbine representation, known as GAST.

```math
\begin{align}
\dot{x}_{g1} &= \frac{1}{T_1} (x_{in} - x_{g1}) \tag{4a} \\
\dot{x}_{g2} &= \frac{1}{T_2} \left(x_{g1}^\text{sat} - x_{g2}\right) \tag{4b} \\
\dot{x}_{g3} &= \frac{1}{T_3} (x_{g2} - x_{g3}) \tag{4c}
\end{align}
```

with

```math
\begin{align}
x_{in} &= \min\left\{P_{ref} - \frac{1}{R}(\omega - 1.0), A_T + K_T (A_T - x_{g3}) \right\} \tag{4d} \\
x_{g1}^\text{sat} &= \left\{ \begin{array}{cl}
                        x_{g1} & \text{ if } V_{min} \le x_{g1} \le V_{max}\\
                        V_{max} & \text{ if } x_{g1} > V_{max} \\
                        V_{min} & \text{ if } x_{g1} < V_{min}
                    \end{array} \right. \tag{4e} \\
\tau_m &= x_{g2}  - D_T(\omega - 1.0) \tag{4f}
\end{align}
```

## HYGOV ```[HydroTurbineGov]```

This represents a classical hydro governor, known as HYGOV.

```math
\begin{align}
T_f\dot{x}_{g1} &= P_{in} - x_{g1} \tag{5a} \\
\dot{x}_{g2} &= x_{g1} \tag{5b}\\
T_g \dot{x}_{g3} &= c - x_{g3} \tag{5c}\\
\dot{x}_{g4} &= \frac{1 - h}{T_w} \tag{5d}
\end{align}
```

with

```math
\begin{align}
P_{in} &= P_{ref} - \Delta \omega - R x_{g2} \tag{5e} \\
c &= \frac{x_{g1}}{r} + \frac{x_{g2}}{rT_r} \tag{5f} \\
h &= \left(\frac{x_{g4}}{x_{g3}}\right)^2 \tag{5g}\\
\tau_m &= h\cdot A_t(x_{g4} - q_{NL}) - D_{turb} \Delta\omega \cdot x_{g3} \tag{5h}
\end{align}
``` 
