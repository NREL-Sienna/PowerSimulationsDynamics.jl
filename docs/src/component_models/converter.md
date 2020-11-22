# Converter

This component can be used to model the dynamics of the switching process.

## Average Model ```[AverageConverter]```

The average model outputs the desired reference signal since:

```math
\begin{align}
v_{d}^{\text{cv}} \approx m_{d} v_{\text{dc}} \approx \frac{v_{d}^{\text{ref-signal}}}{v_{\text{dc}}} v_{\text{dc}} \approx v_{d}^{\text{ref-signal}} \tag{4a} \\
v_{q}^{\text{cv}} \approx m_{q} v_{\text{dc}} \approx \frac{v_{q}^{\text{ref-signal}}}{v_{\text{dc}}} v_{\text{dc}} \approx v_{q}^{\text{ref-signal}} \tag{4b}
\end{align}
```

where ``m_{dq}`` is the modulation signal, and ``v_{dq}^{\text{ref-signal}}`` is the voltage reference signal from the inner loop control.
