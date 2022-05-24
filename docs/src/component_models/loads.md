Here we discuss the models used to describe load modeling in `PowerSimulationsDynamics.jl`. 
In a similar fashion of other devices, loads will withdraw power (i.e. current) from the current-injection balances at the nodal level. Based on the specified parameters and model chosen, the equations for computing such withdrawal will change.

## Static Loads (or Algebraic Loads)

### ZIP + Exponential Load Model

`PowerSimulationsDynamics.jl` uses all the static ZIP and exponential loads at each bus to obtain a single structure that creates an aggregate ZIP load model and a collection of all exponential loads.

The ZIP load model given by:

```math
\begin{align}
P_\text{zip} &= P_\text{power} + P_\text{current} \cdot \frac{V}{V_0} + P_\text{impedance} \cdot \left(\frac{V}{V_0}\right)^2 \\
Q_\text{zip} &= Q_\text{power} + Q_\text{current}  \cdot \frac{V}{V_0} + Q_\text{impedance} \cdot \left(\frac{V}{V_0}\right)^2
\end{align}
```

with ``V = \sqrt{V_r^2 + V_i^2}`` and ``V_0`` the voltage magnitude from the power flow solution.

The current taken for the load is computed as:

```math
\begin{align}
I_\text{zip} &= \frac{(P_\text{zip} + j Q_\text{zip})^*}{(V_r + j V_i)^*} \\
I_\text{zip} &= \frac{P_\text{zip} - j Q_\text{zip}}{V_r - j V_i}
\end{align}
```

For constant impedance load, the current obtained is:

```math
\begin{align}
I_\text{re}^z &= \frac{1}{V_0^2} \cdot (V_r \cdot P_\text{impedance} + V_i \cdot Q_\text{impedance}) \\
I_\text{im}^z &= \frac{1}{V_0^2} \cdot (V_i \cdot P_\text{impedance} - V_r \cdot Q_\text{impedance})
\end{align}
```

For constant current load, the current obtained is:

```math
\begin{align}
I_\text{re}^i  &= \frac{1}{V_0} \cdot \frac{V_r * P_\text{current} + V_i * Q_\text{current}}{V} \\
I_\text{im}^i  &= \frac{1}{V_0} \cdot \frac{V_i * P_\text{current} - V_r * Q_\text{current}}{V}
\end{align}
```

For constant power load, the current obtained is:

```math
\begin{align}
I_\text{re}^p  =  \frac{V_r \cdot P_\text{power} + V_i \cdot Q_\text{power}}{V^2} \\
I_\text{im}^p =  \frac{V_i \cdot P_\text{power} - V_r \cdot Q_\text{power}}{V^2}
\end{align}
```

Then the total current withdrawed from the ZIP load model is simply
```math
\begin{align}
I_\text{zip}^\text{re}  &=  I_\text{re}^z + I_\text{re}^i + I_\text{re}^p \\
I_\text{zip}^\text{im}  &=  I_\text{im}^z + I_\text{im}^i + I_\text{im}^p
\end{align}
```

On the case of Exponential Loads, the model is given by:

```math
\begin{align}
P_\text{exp} = P_0 \cdot \left(\frac{V}{V_0}\right)^\alpha \\
Q_\text{exp} = Q_0 \cdot \left(\frac{V}{V_0}\right)^\beta
\end{align}
```

The current taken for the load is computed as:
```math
\begin{align}
I_\text{exp} &= \frac{(P_\text{exp} + j Q_\text{exp})^*}{(V_r + j V_i)^*} \\
I_\text{exp} &= \frac{P_\text{exp} - j Q_\text{exp}}{V_r - j V_i}
\end{align}
```

that results:
```math
\begin{align}
I_\text{exp}^\text{re}  &= V_r \cdot P_0 \cdot \frac{V^{\alpha - 2}}{V_0^\alpha} + V_i \cdot Q_0 \cdot \frac{V^{\beta - 2}}{V_0^\beta} \\
I_\text{exp}^\text{im}  &= V_i \cdot P_0 \cdot \frac{V^{\alpha - 2}}{V_0^\alpha} - V_r \cdot Q_0 \cdot \frac{V^{\beta - 2}}{V_0^\beta}
\end{align}
```

## Dynamic loads

### 5th-order Single Cage Induction Machine ```[SingleCageInductionMachine]```

The following model is used to model a 5th-order induction machine with a quadratic relationship speed-torque.

