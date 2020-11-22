# Quick Start Guide

You can access example data in the [Power Systems Test Data Repository](https://github.com/NREL-SIIP/PowerSystemsTestData),
the data can be downloaded with the `PowerSystems.jl` submodule `UtilsData`

```julia
using PowerSystems
DATA_DIR = download(PowerSystems.UtilsData.TestData, folder = pwd())
```

## Loading data

Data can be loaded from several file formats and return a summary of the system's components and
time-series.

```@repl generated_quick_start_guide
using PowerSystems
DATA_DIR = "../../data" #hide
system_data = System(joinpath(DATA_DIR, "matpower/RTS_GMLC.m"))
```
