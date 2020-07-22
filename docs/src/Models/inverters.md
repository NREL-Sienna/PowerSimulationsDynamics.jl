# Inverter Models

Here we discuss the structure and models used to model inverters in PowerSimulationsDynamics.jl. Each inverter is a data structure that is defined by the following components:

- DC Source: Defines the dynamics of the DC side of the converter.
- Frequency Estimator: That describes how the frequency of the grid can be estimated using the grid voltages. Typically a phase-locked loop (PLL).
- Outer Loop Control: That describes the active and reactive power control dynamics.
- Inner Loop Control: That can describe virtual impedance, voltage control and current control dynamics.
- Converter: That describes the dynamics of the pulse width modulation (PWM) or space vector modulation (SVM).
- Filter: Used to connect the converter output to the grid.

The following figure summarizes the components of a inverter and which variables they share:

```@raw html
<img src="../../assets/inv_metamodel.png" width="75%"/>
``` ⠀

Contrary to the generator, there are many control structures that can be used to model inverter controllers (e.g. grid-following, grid feeding or virtual synchronous machine). For this purpose, more variables are shared among the components in order to cover all these posibilities.

Models are based from the paper: "A Virtual Synchronous Machine implementation for distributed control of power converters in SmartGrids" from S. D'Arco, J.A. Suul and O.B. Fosso, and structures are defined in ```PowerSystems.jl``` abbreviated as ```PSY```.

## DC Source

This component can be used to model the dynamics of the DC side of the converter.

### Fixed DC Source ```[FixedDCSource]```

This is a model that set the DC voltage to a fixed value ``v_{\text{dc}} = v_{\text{dc}}^{\text{fix}}``.


## Frequency Estimators

This component is used to estimate the frequency of the grid based on the voltage at the bus.

### Phase-Locked Loop (PLL) for VSM ```[KauraPLL]```

The following equations present a PLL used to estimate the frequency and PLL angle of the grid. There are two reference frames considered in this inverter. Those are the VSM of the outer-loop control ``\delta\theta_{\text{olc}}`` and the PLL one ``\delta\theta_{\text{pll}}``. The notation used a ``\delta\theta`` to refer as the variation of the respective angle ``\theta`` with respect to the grid SRF (instead of the fixed ``\alpha`` component of the ``\alpha\beta`` transformation):

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
on which ``v_r + jv_i`` is the voltage in the grid reference frame on which the PLL is measuring (i.e. point of common coupling), that could be in the capacitor of an LCL filter or the last branch of such filter.

## Outer Loop Controls

This component defines controllers for both active and reactive power

### Virtual Inertia and Q-droop ```[OuterControl]```

The following model represent a virtual synchronous machine model to represent how active power is going to be deployed. The constructor is ```OuterControl{VirtualInertia, ReactivePowerDroop}```. It defines a new SRF denoted as ``\theta_{\text{olc}}`` for the active power controller and uses a simple voltage droop for dispatching reactive power:

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
In this case, the measurement of power are being done in the capacitor of the LCL filter. However, depending on the model, this measurements could be different depending on where is the point of common coupling.

## Inner Loop Controls

This component defines voltage and current controllers to generate the reference signal for the converter.

### Integrated Virtual Impedance, Voltage and Current Controller ```[CurrentControl]```

The following model receives both the outer loop control frequency and reference voltage signal to generate the reference signal for the converters. The virtual impedance plays a similar role of the impedance of a synchronous generator. A PI voltage controller is used to generate the current signal that is used in the PI current controller to finally generate the voltage reference signal for the converters.

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
that again ``v_r + jv_i`` could be in the capacitor or the last branch of the filter (i.e. the point of common coupling). For LCL filters it is considered in the capacitor. In the case of the converter, the transformation is directly
```math
\begin{align*}
i_{d,\text{cv}} + ji_{q,\text{cv}} = (i_{r,\text{cv}} + ji_{i,\text{cv}})e^{-j\delta\theta_{olc}}
\end{align*}
```

## Converter

This component can be used to model the dynamics of the switching process.

### Average Model ```[AverageConverter]```

The average model simply output the desired reference signal since:

```math
\begin{align}
v_{d}^{\text{cv}} \approx m_{d} v_{\text{dc}} \approx \frac{v_{d}^{\text{ref-signal}}}{v_{\text{dc}}} v_{\text{dc}} \approx v_{d}^{\text{ref-signal}} \tag{4a} \\
v_{q}^{\text{cv}} \approx m_{q} v_{\text{dc}} \approx \frac{v_{q}^{\text{ref-signal}}}{v_{\text{dc}}} v_{\text{dc}} \approx v_{q}^{\text{ref-signal}} \tag{4b}
\end{align}
```

where ``m_{dq}`` is the modulation signal, and ``v_{dq}^{\text{ref-signal}}`` is the voltage reference signal from the inner loop control.

## Filters

### LCL Filter ```[LCLFilter]```

A standard LCL filter is proposed to connect the output of the converter to the grid. In this case, ``v_r`` and ``v_i`` are voltages in the capacitor, while ``v_r^{\text{grid}}`` and ``v_i^{\text{grid}}`` represent the voltage at the bus. The L filter after the capacitor can also include a step-up transformer to increase the voltage, that is model as an extra impedance.


```math
\begin{align}
    \dot{i}_{r,\text{cv}} &= \frac{\Omega_b}{l_f}\left( v_r^{\text{cv}} - v_r  - r_f i_{r,\text{cv}} + \omega_{\text{grid}} l_f i_{i,\text{cv}} \right) \tag{5a} \\
    \dot{i}_{i,\text{cv}} &= \frac{\Omega_b}{l_f}\left( v_i^{\text{cv}} - v_i  - r_f i_{i,\text{cv}} - \omega_{\text{grid}} l_f i_{r,\text{cv}} \right) \tag{5b} \\
    \dot{v}_{r} &=  \frac{\Omega_b}{c_f}\left( i_r^{\text{cv}} - i_r + \omega_{\text{grid}} c_f v_i \right) \tag{5c} \\
    \dot{v}_{i} &=  \frac{\Omega_b}{c_f}\left( i_i^{\text{cv}} - i_i - \omega_{\text{grid}} c_f v_r \right) \tag{5d} \\
    \dot{i}_{r} &= \frac{\Omega_b}{l_g}\left( v_r^{\text{cv}} - v_r^{\text{grid}} - r_g i_{r} + \omega_{\text{grid}} l_g i_{i,\text{cv}} \right) \tag{5e} \\
    \dot{i}_{i} &= \frac{\Omega_b}{l_g}\left( v_i^{\text{cv}} - v_i^{\text{grid}} - r_g i_{i} - \omega_{\text{grid}} l_g i_{r,\text{cv}} \right) \tag{5f}
\end{align}
```
on which 
```math
\begin{align*}
v_r^\text{cv} + jv_i^\text{cv} = (v_d^\text{cv} + jv_q^\text{cv})e^{j\delta\theta_{olc}}
\end{align*}
```
that comes from the converter model.

## Reference

For constructors check the API on [PowerSystems.jl documentation](https://nrel.github.io/PowerSystems.jl/latest/api/PowerSystems/)
