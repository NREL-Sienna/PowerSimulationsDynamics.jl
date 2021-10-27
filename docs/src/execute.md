# Executing a Simulation

Once a Simulation is constructed and properly initialized, the `execute!` command is used to run the Simulation. If no perturbation was included, then a steady state simulation will be run over the time span defined. The constructor is the following:

```raw
execute!(
    Simulation Structure,
    Solver,
    Optional: Keyword Arguments
)
```

## Solvers

Solvers must be chosen accordingly depending on the type of model used in the Simulation. For example, a Residual model can be executed using Sundials IDA solver:

```julia
using Sundials
sim = Simulation(
    ResidualModel,
    sys,
    pwd(),
    (0.0, 20.0),
    perturbation,
)
execute!(sim, IDA())
```

Results can be explored using:
```julia
results = read_results(sim)
```

Similarly, a Mass Matrix model can be executed using `Rodas4` solver.

```julia
using OrdinaryDiffEq
sim2 = Simulation(
    ResidualModel,
    sys,
    pwd(),
    (0.0, 20.0),
    perturbation,
)
execute!(sim2, Rodas4())
```

## Exploring the Solution

Once a Simulation is executed and the results are stored via `results = read_results(sim)`, the following functions can be used to explore the Simulation solution:

### Show initial conditions

The function `show_states_initial_value(results)` can be used to display the initial condition of the voltages and dynamic states of each dynamic component.

### Explore bus voltages

The function `get_voltage_magnitude_series(results, BusNumber)` can be used to obtain the voltage magnitude time series of the specified bus. Similarly, `get_voltage_angle_series(results, BusNumber)` can be used to obtain the voltage angle time series of the specified bus.

### Explore output currents

The functions `get_real_current_series(results, "DeviceName")` and `get_imaginary_current_series(results, "DeviceName")` can be used to obtain the output current time series of the specified device.

### Explore output power

The functions `get_activepower_series(results, "DeviceName")` and `get_reactivepower_series(results, "DeviceName")` can be used to obtain the output power time series of the specified device.

### Explore dynamic states

The function `get_state_series(results, ("DeviceName", :StateSymbol)` can be used to obtain the specified state time series of the specified device.

### Explore Reference Setpoints

The function `get_setpoints(sim)` can be used to obtain the reference setpoints of each dynamic device. **Note:** If a setpoint was changed via a perturbation, this function will return the modified setpoint.


## Keyword Arguments

Any solver option available in `DifferentialEquations.jl` can be passed as keyword arguments in the `execute!` function. Please see the [Common Solver Options](https://diffeq.sciml.ai/stable/basics/common_solver_opts/) in the `DifferentialEquations.jl` documentation for more details.

Most common solver options used are `dtmax` to control the maximum dt for adaptive timestepping. `abstol` and `reltol` are also commonly used to control the tolerance in the adaptive timestepping. `saveat` is also used to store the results at a specified time stamps. For example, the following code is valid to further specify your solver options:

```julia
execute!(sim, IDA(), dtmax = 0.01, abstol = 1e-9, reltol = 1e-6, saveat = 0.01)
```

In addition, the keyword argument `enable_progress_bar = false` can be used to disable the progress bar.