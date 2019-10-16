# Generator Models

Here we discuss the structure and models used to model generators in LITS.jl. Each generator is a data structure that is defined by the following components:

- Machine: That defines the stator electro-magnetic dynamics.
- Shaft: That describes the rotor electro-mechanical dynamics.
- Automatic Voltage Regulator: Electromotive dynamics to model an AVR controller.
- Power System Stabilizer: Control dynamics to define an stabilization signal for the AVR.
- Prime Mover (Governor): Thermo-mechanical dynamics and associated controllers.

The following figure summarizes the components of a generator and which variables they share:

```@raw html
<img src="../../assets/gen_metamodel.png" width="75%"/>
``` ⠀

Each generator is defined in its own ``dq`` reference frame. Let ``\delta`` be the rotor angle of the generator. If ``v_r + jv_i = v_h\angle \theta`` defines the voltage in the bus in the network reference frame ``RI`` rotating at nominal frequency ``\Omega_b``, then the following equations (both are equivalent) can be used to convert between reference frames:
```math
\begin{align}
\left[ \begin{array}{c} v_d \\ v_q \end{array} \right] &=  \left[ \begin{array}{c} v_h \sin(\delta - \theta) \\ v_h \cos(\delta - \theta) \end{array} \right]  \tag{0a} \\
\left[ \begin{array}{c} v_d \\ v_q \end{array} \right] &= \left[ \begin{array}{cc} \sin(\delta) & -\cos(\delta) \\ \cos(\delta) & \sin(\delta) \end{array} \right] \left[ \begin{array}{c} v_r \\ v_i \end{array} \right] \tag{0b}
\end{align}
```

## Machines

### Classical Model
This is the classical order model that does not have differential equations (``\delta`` and ``\omega`` are defined in the shaft):

```math
\begin{align}
  \left[ \begin{array}{c} i_d \\ i_q \end{array} \right] &= \left[ \begin{array}{cc} r_a & -x_d' \\ x_d' & r_a \end{array} \right]  \left[ \begin{array}{c} -v_d \\ e_q' - v_q \end{array} \right] \tag{1a}\\
p_e &= (v_q + r_a i_q)i_q + (v_d + r_ai_d)i_d \tag{1b}
\end{align}
```



## Shafts



### Rotor Mass

### Multi-Mass



### 4th Order

bla
