# Network model

Here we discuss the models used to describe the network in `PowerSimulationsDynamics.jl`.
This is based on a standard current injection model as defined in [Power System Modelling and Scripting](https://www.springer.com/gp/book/9783642136689).
The numerical advantages of current injection models outweigh the complexities of
implementing constant power loads for longer-term transient stability analysis.
The network is defined in a synchronous reference frame (SRF), named the RI (real-imaginary)
reference frame, rotating at the constant base frequency ``\Omega_b``.

In simple terms, `PowerSimulationsDynamics.jl` internally tracks the current-injection
balances at the nodal level from all the devices on the system. Based on the buses and
branches information, the system constructor computes the admittance matrix ``\boldsymbol{Y}``
assuming nominal frequency and this is used for static branch modeling. The algebraic equations
for the static portions of the network are as follows:

```math
 \begin{align}
 0 = \boldsymbol{i}(\boldsymbol{x}, \boldsymbol{v}) - \boldsymbol{Y}\boldsymbol{v}
 \end{align}
```

where ``\boldsymbol{i} = i_r + ji_i`` is the vector of the sum of complex current injections from devices
, ``\boldsymbol{x}`` is the vector of states and ``\boldsymbol{v} = v_r + jv_i`` is the vector of complex
bus voltages. Equations (1) connect all the port variables, i.e., currents, defined for each
injection device. Components that contribute to (1) by modifying the current ``\boldsymbol{i}``
are (i) static injection devices, (ii) dynamic injection devices, and (iii) dynamic network
branches. Components that contribute to modify the admittance matrix ``\boldsymbol{Y}``
are static branches.

## Static Branches (or Algebraic Branches)

### Lines

Each line is defined using a ``\pi`` model connecting two buses ``(n,m)``, with a series
resistance ``r`` and reactance ``x``, and a shunt capacitance at both ends ``(c_n, c_m)``.
The values are already in system per unit. Then each branch contributes to the admittance matrix as follows:

```math
\begin{align}
Y_{nn} &+\!= \frac{1}{r+jx} + jc_n \\
Y_{nm} &+\!= \frac{-1}{r+jx} \\
Y_{mm} &+\!= \frac{1}{r+jx} + jc_m \\
Y_{mn} &+\!= \frac{-1}{r+jx} \\
\end{align}
```

### Two-Windings Transformers

Similarly to lines these are defined by a series reactance and impedance. The equations are
equivalently of the lines without the shunt capacitance.

## Dynamic Branches

Dynamic network branches contribute directly to (1) by modifying the vector of complex currents.
Their parameters are also the series resistance ``r`` and reactance ``x``, and a shunt
capacitance at both ends ``(c_n, c_m)`` for a line ``\ell``. In addition, they define 3
new additional differential equations per line (6 in total for real and imaginary part):

```math
\begin{align}
    \frac{l}{\Omega_b} \frac{d\boldsymbol{i}_{\ell}}{dt} &= (\boldsymbol{v}_n - \boldsymbol{v}_m) - (r+jl) \boldsymbol{i}_{\ell} \\
     \frac{c_n}{\Omega_b} \frac{d\boldsymbol{v}_n}{dt} &=  \boldsymbol{i}_n^{\text{cap}} - jc_n\boldsymbol{v}_n   \\
      \frac{c_m}{\Omega_b} \frac{d\boldsymbol{v}_m}{dt} &= \boldsymbol{i}_m^{\text{cap}} - jc_m\boldsymbol{v}_m
\end{align}
```

Since all the values are in per unit, the reactance is equal to the inductance.

A detail discussion about the effects of different line models in the modeling of inverters
is presented in [Grid Forming Inverter Small Signal Stability: Examining Role of Line and Voltage Dynamics](https://ieeexplore.ieee.org/document/9255030)
