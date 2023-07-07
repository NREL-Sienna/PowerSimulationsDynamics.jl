# Quick Start Guide

The data for these tutorials is provided in [PowerSystemCaseBuilder](https://github.com/nrel-sienna/PowerSystemCaseBuilder.jl). If you want to build your own case, take a look at the tutorial [Creating and Handling Data for Dynamic Simulations](@ref)

For more details about loading data and adding more dynamic components check the
[Creating a System with Dynamic devices](https://nrel-sienna.github.io/PowerSystems.jl/stable/modeler_guide/system_dynamic_data/)
section of the documentation in `PowerSystems.jl`.

For a detailed tutorial about this case visit [One Machine against Infinite Bus (OMIB) Simulation](@ref)

## Loading data

Data can be loaded from a pss/e raw file and a pss/e dyr file.

```@repl quick_start_guide
using PowerSystems
using PowerSimulationsDynamics
using PowerSystemCaseBuilder
using Sundials
using Plots

omib_sys = build_system(PSIDSystems, "OMIB System")
```

## Define the Simulation

```@repl quick_start_guide
time_span = (0.0, 30.0)
perturbation_trip = BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1")
sim = Simulation!(ResidualModel, omib_sys, pwd(), time_span, perturbation_trip)
```

## Explore initial conditions for the simulation

```@repl quick_start_guide
x0_init = read_initial_conditions(sim)
```

```@repl quick_start_guide
show_states_initial_value(sim)
```

## Obtain small signal results for initial conditions

Show eigenvalues for operating point
```@repl quick_start_guide
    small_sig = small_signal_analysis(sim)
    summary_eigenvalues(small_sig)
```

## Execute the simulation

```@repl quick_start_guide
    execute!(sim, IDA(), dtmax = 0.02, saveat = 0.02, enable_progress_bar = false)
```

## Make a plot of the results

```@repl quick_start_guide
results = read_results(sim)
angle = get_state_series(results, ("generator-102-1", :δ));
plot(angle, xlabel = "time", ylabel = "rotor angle [rad]", label = "gen-102-1")
```

![plot](assets/f-plot.svg)

If you miss PSS/e's plotting aesthetics and want something that resembles that, you can use [`UnicodePlots`](https://github.com/Evizero/UnicodePlots.jl).

```@repl quick_start_guide
using UnicodePlots
unicodeplots()
plot(angle, xlabel = "time", ylabel = "rotor angle [rad]", label = "gen-102-1");
```

![plot](assets/unicode.png)
