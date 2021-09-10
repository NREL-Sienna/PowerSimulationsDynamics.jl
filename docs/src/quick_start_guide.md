# Quick Start Guide

You can access example data in the [Power Systems Test Data Repository](https://github.com/NREL-SIIP/PowerSystemsTestData),
the data can be downloaded with the `PowerSystems.jl` submodule `UtilsData`.

## Loading data

Data can be loaded from a pss/e raw file and a pss/e dyr file.

```@repl quick_start_guide
using PowerSystems, PowerSimulationsDynamics, Sundials, Plots, Logging
logger = configure_logging(console_level = Logging.Error, file_level = Logging.Info)
DATA_DIR = download(PowerSystems.UtilsData.TestData, folder = pwd())
omib_sys = System(joinpath(DATA_DIR, "psse_raw/OMIB.raw"),
                     joinpath(DATA_DIR, "psse_dyr/OMIB.dyr");
                     time_series_in_memory = true)
```

For more details about loading data and adding more dynamic components check the
[Creating a System with Dynamic devices](https://nrel-siip.github.io/PowerSystems.jl/stable/modeler_guide/system_dynamic_data/)
section of the documentation in `PowerSystems.jl`.

## Define the Simulation

```@repl quick_start_guide
time_span = (0.0, 30.0)
perturbation_trip = BranchTrip(1.0, "BUS 1-BUS 2-i_1")
sim = Simulation!(ResidualModel, omib_sys, pwd(), time_span, perturbation_trip)
```

## Explore initial conditions for the simulation

```@repl quick_start_guide
x0_init = read_initial_conditions(sim)
```

## Obtain small signal results for initial conditions

```@repl quick_start_guide
    small_sig = small_signal_analysis(sim)
```

## Execute the simulation

```@repl quick_start_guide
    execute!(sim, IDA())
```

## Make a plot of the results

```@repl quick_start_guide
angle = get_state_series(res, ("generator-102-1", :δ));
plot(angle, xlabel = "time", ylabel = "rotor angle [rad]", label = "gen-102-1");
```

![plot](assets/f-plot.svg)

If you miss PSS/e's plotting aesthetics and want something that resembles that, you can use [`UnicodePlots`](https://github.com/Evizero/UnicodePlots.jl).

```@repl quick_start_guide
using UnicodePlots
unicodeplots()
plot(angle, xlabel = "time", ylabel = "rotor angle [rad]", label = "gen-102-1");
```

![plot](assets/unicode.png)
