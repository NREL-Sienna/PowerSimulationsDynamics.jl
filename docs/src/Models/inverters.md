# Inverter Models

Here we discuss the structure and models used to model inverters in LITS.jl. Each inverter is a data structure that is defined by the following components:

- DC Source: Defines the dynamics of the DC side of the converter.
- Converter: That describes the dynamics of the pulse width modulation (PWM) or space vector modulation (SVM).
- Frequency Estimator: That describes how the frequency of the grid can be estimated using the grid voltages. Typically a phase-locked loop (PLL).
- Outer Loop Control: That describes the active and reactive power control dynamics.
- Inner Loop Control: That can describe virtual impedance, voltage control and current control dynamics.
- Filter: Used to connect the converter output to the grid.

The following figure summarizes the components of a inverter and which variables they share:

```@raw html
<img src="../../assets/inv_metamodel.png" width="75%"/>
``` ⠀

Contrary to the generator, there are many control structures that can be used to model inverter controllers (e.g. grid-following, grid feeding or virtual synchronous machine). For this purpose, more variables are shared among the components in order to cover all these posibilities.

## DC Source

This component can be used to model the dynamics of the DC side of the converter.

### Fixed DC Source

This is a model that set the DC voltage to a fixed value ``v_{\text{dc}} = v_{\text{dc}}^{\text{fix}}``.

## Converter

This component can be used to model the dynamics of the switching process.

### Average Model

The average model simply output the desired reference signal since:

```math
\begin{align}
v_{dq}^{\text{cnv}} \approx m_{dq} v_{\text{dc}} \approx \frac{v_{dq}^{\text{ref-signal}}}{v_{\text{dc}}} v_{\text{dc}} \approx v_{dq}^{\text{ref-signal}}
\end{align}
```

where ``m_{dq}`` is the modulation signal, and ``v_{dq}^{\text{ref-signal}}`` is the voltage reference signal from the inner loop control.

## Frequency Estimator

a
