# Automatic Voltage Regulators (AVR)

AVR are used to determine the voltage in the field winding ``v_f`` in the model.

## Fixed AVR ```[AVRFixed]```

This is a simple model that set the field voltage to be equal to a desired constant value ``v_f = v_{\text{fix}}``.

## Simple AVR ```[AVRSimple]```

This depicts the most basic AVR, on which the field voltage is an integrator over the difference of the measured voltage and a reference:

```math
\begin{align}
\dot{v}_f = K_v(v_{\text{ref}} - v_h) \tag{10a}
\end{align}
```

## AVR Type I ```[AVRTypeI]```

This AVR is a simplified version of the IEEE DC1 AVR model:

```math
\begin{align}
\dot{v}_f &= -\frac{1}{T_e} \left[ V_f(K_e + S_e(v_f))-v_{r1} \right] \tag{11a} \\
\dot{v}_{r1} &= \frac{1}{T_a} \left[ K_a\left(v_{\text{ref}} - v_m - v_{r2} - \frac{K_f}{T_f}v_f\right) - v_{r1} \right]   \tag{11b} \\
\dot{v}_{r2} &=  -\frac{1}{T_f} \left[ \frac{K_f}{T_f}v_f + v_{r2} \right]  \tag{11c} \\
\dot{v}_m &= \frac{1}{T_r} (v_h - v_m) \tag{11d}
\end{align}
```

with the ceiling function:

```math
\begin{align*}
S_e(v_f) = A_e \exp(B_e|v_f|)
\end{align*}
```

## AVR Type II ```[AVRTypeII]```

This model represents a static exciter with higher gains and faster response than the Type I:

```math
\begin{align}
\dot{v}_f &= -\frac{1}{T_e} \left[ V_f(1 + S_e(v_f))-v_{r} \right] \tag{12a} \\
\dot{v}_{r1} &= \frac{1}{T_1} \left[ K_0\left(1 - \frac{T_2}{T_1} \right)(v_{\text{ref}} - v_m) - v_{r1}  \right] \tag{12b} \\
\dot{v}_{r2} &=  \frac{1}{K_0 T_3} \left[ \left( 1 - \frac{T_4}{T_3} \right) \left( v_{r1} + K_0\frac{T_2}{T_1}(v_{\text{ref}} - v_m)\right) - K_0 v_{r2} \right]  \tag{12c} \\
\dot{v}_m &= \frac{1}{T_r} (v_h - v_m) \tag{12d}
\end{align}
```

with

```math
\begin{align*}
v_r &= K_0v_{r2} + \frac{T_4}{T_3} \left( v_{r1} + K_0\frac{T_2}{T_1}(v_{\text{ref}} - v_m)\right) \\
S_e(v_f) &= A_e \exp(B_e|v_f|)
\end{align*}
```
