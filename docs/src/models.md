
# Models

## Simulation Models

PowerSimulations dynamics supports two formulations for the simulation model and define different methods for each simulation model. You can pass `ImplicitModel` or `MassMatrixModel` to a call to Simulation to define the preferred formulation.

In this way, we provide a common set of development requirements for contributors of new models that maintains the same flexibility in choosing the solving algorithm.

- *MassMatrixModel*: Defines models that can be solved using [Mass-Matrix Solvers](https://diffeq.sciml.ai/stable/solvers/dae_solve/#OrdinaryDiffEq.jl-(Mass-Matrix)). The model is formulated as follows:

```math
\begin{align}
M\frac{dx(t)}{dt} = x(t)
\end{align}
```

At this stage we have not conducted extensive tests with all the solvers in [DifferentialEquations](https://diffeq.sciml.ai/) most of our tests use `Rodas5()`.


- *ImplicitModel*: Define models that can be solved using [Implicit ODE solvers](https://diffeq.sciml.ai/stable/solvers/dae_solve/#OrdinaryDiffEq.jl-(Implicit-ODE)) and also the solver IDA from [Sundials](https://diffeq.sciml.ai/stable/solvers/dae_solve/#Sundials.jl). The model is formulated to solved the following problem:

```math
\begin{align}
r(t) = \frac{dx(t)}{dt} - x(t)
\end{align}
```

At this stage we have not conducted extensive tests with all the solvers in [DifferentialEquations](https://diffeq.sciml.ai/) if you are solving a larger system use `IDA()`.

### The dynamic systen model in PowerSimulationsDynamics

In order to support both formulations, the default implementation of the ImplicitModel solves the following problem:

```math
\begin{align}
r(t) = M\frac{dx(t)}{dt} - x(t)
\end{align}
```

## Generator Models

Here we discuss the structure and models used to model generators in `PowerSimulationsDynamics.jl`. See [`PowerSystems.jl` dynamic devices](https://nrel-siip.github.io/PowerSystems.jl/stable/modeler_guide/example_dynamic_data/)
for details.

Each generator is a data structure composed of the following components defined in `PowerSystems.jl`:

- [`Machine`](https://nrel-siip.github.io/PowerSystems.jl/stable/model_library/generated_Machine/#Machine): That defines the stator electro-magnetic dynamics.
- [`Shaft`](https://nrel-siip.github.io/PowerSystems.jl/stable/model_library/generated_Shaft/#Shaft): That describes the rotor electro-mechanical dynamics.
- [`Automatic Voltage Regulator`](https://nrel-siip.github.io/PowerSystems.jl/stable/model_library/generated_AVR/#AVR): Electromotive dynamics to model an AVR controller.
- [`Power System Stabilizer`](https://nrel-siip.github.io/PowerSystems.jl/stable/model_library/generated_PSS/#PSS): Control dynamics to define an stabilization signal for the AVR.
- [`Prime Mover and Turbine Governor`](https://nrel-siip.github.io/PowerSystems.jl/stable/model_library/generated_TurbineGov/#TurbineGov): Thermo-mechanical dynamics and associated controllers.

The implementation of Synchronous generators as components uses the following structure to
share values across components.

```@raw html
<img src="https://github.com/NREL-SIIP/PowerSystems.jl/blob/master/docs/src/assets/gen_metamodel.png?raw=true" width="75%">
```

## Inverter Models

Here we discuss the structure and models used to model inverters in `PowerSimulationsDynamics.jl`. See [`PowerSystems.jl` dynamic devices](https://nrel-siip.github.io/PowerSystems.jl/stable/modeler_guide/example_dynamic_data/)
for details. One of the key contributions in this software package is a separation of the
components in a way that resembles current practices for synchronoues machine modeling.

- [`DC Source`](https://nrel-siip.github.io/PowerSystems.jl/stable/model_library/generated_DCSource/#DCSource): Defines the dynamics of the DC side of the converter.
- [`Frequency Estimator`](https://nrel-siip.github.io/PowerSystems.jl/stable/model_library/generated_FrequencyEstimator/#FrequencyEstimator): That describes how the frequency of the grid can be estimated using the grid voltages. Typically a phase-locked loop (PLL).
- [`Outer Loop Control`](https://nrel-siip.github.io/PowerSystems.jl/stable/model_library/outer_control/#OuterControl): That describes the active and reactive power control dynamics.
- [`Inner Loop Control`](https://nrel-siip.github.io/PowerSystems.jl/stable/model_library/generated_InnerControl/#InnerControl): That can describe virtual impedance, voltage control and current control dynamics.
- [`Converter`](https://nrel-siip.github.io/PowerSystems.jl/stable/model_library/generated_Converter/#Converter): That describes the dynamics of the pulse width modulation (PWM) or space vector modulation (SVM).
- [`Filter`](https://nrel-siip.github.io/PowerSystems.jl/stable/model_library/generated_Filter/): Used to connect the converter output to the grid.

The following figure summarizes the components of a inverter and which variables they share:

```@raw html
<img src="https://github.com/NREL-SIIP/PowerSystems.jl/blob/master/docs/src/assets/inv_metamodel.png?raw=true" width="75%">
```

Contrary to the generator, there are many control structures that can be used to model
inverter controllers (e.g. grid-following, grid feeding or virtual synchronous machine).
For this purpose, more variables are shared among the components in order to cover all
these posibilities.

## Reference

For models, check the library in [PowerSystems.jl](https://nrel-siip.github.io/PowerSystems.jl/stable/)
