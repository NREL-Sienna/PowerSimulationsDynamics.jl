# Dynamic Injection Models

## Generator Models

Here we discuss the structure and models used to model generators in `PowerSimulationsDynamics.jl`.
Each generator is a data structure that is defined by the following components:

- Machine: That defines the stator electro-magnetic dynamics.
- Shaft: That describes the rotor electro-mechanical dynamics.
- Automatic Voltage Regulator: Electromotive dynamics to model an AVR controller.
- Power System Stabilizer: Control dynamics to define an stabilization signal for the AVR.
- Prime Mover and Turbine Governor: Thermo-mechanical dynamics and associated controllers.

The following figure summarizes the components of a generator and which variables they share:

![gen_meta](https://github.com/NREL-SIIP/PowerSystems.jl/blob/master/docs/src/assets/gen_metamodel.png)

Models are based from Federico Milano's book: "Power System Modelling and Scripting" and
Prabha Kundur's book: "Power System's Stability and Control" and structures are defined
in ```PowerSystems.jl```.

## Inverter Models

Here we discuss the structure and models used to model inverters in `PowerSimulationsDynamics.jl`. Each inverter is a data structure that is defined by the following components:

- DC Source: Defines the dynamics of the DC side of the converter.
- Frequency Estimator: That describes how the frequency of the grid can be estimated using the grid voltages. Typically a phase-locked loop (PLL).
- Outer Loop Control: That describes the active and reactive power control dynamics.
- Inner Loop Control: That can describe virtual impedance, voltage control and current control dynamics.
- Converter: That describes the dynamics of the pulse width modulation (PWM) or space vector modulation (SVM).
- Filter: Used to connect the converter output to the grid.

The following figure summarizes the components of a inverter and which variables they share:

![inv_meta](https://github.com/NREL-SIIP/PowerSystems.jl/blob/master/docs/src/assets/inv_metamodel.png)

Contrary to the generator, there are many control structures that can be used to model inverter controllers (e.g. grid-following, grid feeding or virtual synchronous machine). For this purpose, more variables are shared among the components in order to cover all these posibilities.

Models are based from the paper: "A Virtual Synchronous Machine implementation for distributed control of power converters in SmartGrids" from S. D'Arco, J.A. Suul and O.B. Fosso, and structures are defined in ```PowerSystems.jl``` abbreviated as ```PSY```.

## Reference

For constructors check the API on [PowerSystems.jl documentation](https://nrel.github.io/PowerSystems.jl/latest/api/PowerSystems/)
