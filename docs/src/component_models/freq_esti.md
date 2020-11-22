# Frequency Estimators

This component is used to estimate the frequency of the grid based on the voltage at the bus.

## Phase-Locked Loop (PLL) for VSM ```[KauraPLL]```

The following equations present a PLL used to estimate the frequency and PLL angle of
the grid. There are two reference frames considered in this inverter. Those are the VSM
of the outer-loop control ``\delta\theta_{\text{olc}}`` and the PLL one
``\delta\theta_{\text{pll}}``. The notation used a ``\delta\theta`` to refer as the variation
of the respective angle ``\theta`` with respect to the grid SRF (instead of the fixed
``\alpha`` component of the ``\alpha\beta`` transformation):

```math
\begin{align}

\dot{v}_{d,\text{pll}} &= \omega_{\text{lp}} \left [v_{d,\text{out}} - v_{d,\text{pll}} \right] \tag{1a} \\
\dot{v}_{q,\text{pll}} &= \omega_{\text{lp}} \left [v_{q,\text{out}} - v_{q,\text{pll}} \right] \tag{1b} \\
\dot{\varepsilon}_{\text{pll}} &= \tan^{-1}\left(\frac{v_{q,\text{pll}}}{v_{d,\text{pll}}} \right) \tag{1c} \\
\dot{\delta\theta}_{\text{pll}} &= \Omega_b \delta \omega_{\text{pll}} \tag{1d}
\end{align}
```

with

```math
\begin{align}
\delta\omega_{\text{pll}} &= k_{p,\text{pll}} \tan^{-1} \left(\frac{v_{q,\text{pll}}}{v_{d,\text{pll}}} \right) + k_{i,\text{pll}} \varepsilon_{\text{pll}} \tag{1e} \\
v_{d,\text{out}} + jv_{q,\text{out}} &= (v_r + jv_i)e^{-\delta\theta_\text{pll}}  \tag{1f}
\end{align}
```

on which ``v_r + jv_i`` is the voltage in the grid reference frame on which the PLL is
measuring (i.e. point of common coupling), that could be in the capacitor of an LCL filter
or the last branch of such filter.
