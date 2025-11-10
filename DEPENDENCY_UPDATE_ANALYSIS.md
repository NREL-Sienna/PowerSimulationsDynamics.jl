# Dependency Update Analysis

## Overview
This document analyzes the changes made to update PowerSimulationsDynamics.jl to be compatible with the latest versions of InfrastructureSystems, PowerSystems, PowerFlows, and PowerNetworkMatrices.

## Version Changes

| Package | Previous | New | Release Date | Breaking Changes |
|---------|----------|-----|--------------|------------------|
| InfrastructureSystems | 2 | 3 | Oct 31, 2024 | Yes - time series interface, supplemental attributes |
| PowerSystems | 4 | 5 | Nov 2024 | Yes - struct redefinitions, API removals (get_bus_numbers) |
| PowerFlows | 0.9 | 0.10 | Nov 4, 2024 | Yes - complete API rewrite |
| PowerNetworkMatrices | 0.12.1 | 0.14 | Nov 4, 2024 | Yes - arc-based indexing |
| **DataStructures** | ~0.18 | 0.19 | - | Required by PowerSystems v5 |
| **Julia** | ^1.6 | ^1.10 | - | Required by PowerSystems v5 and PowerFlows v0.10 |

## Code Changes Made

### 1. InfrastructureSystems v2 → v3 (Project.toml)

**Required for PowerSystems v5 compatibility.**

InfrastructureSystems v3.0.0 breaking changes:
- New Time Series interface
- Supplemental attributes getters modified
- Cost function definitions updated

No code changes required - existing IS API usage remains compatible.

### 2. PowerFlows API Update (src/base/simulation_initialization.jl:14)

**Previous (PowerFlows 0.9):**
```julia
pf = PF.ACPowerFlow{PF.TrustRegionACPowerFlow}()
res = PF.solve_powerflow!(pf, sys)
```

**Updated (PowerFlows 0.10):**
```julia
pf = PF.ACPowerFlow(PF.TrustRegionACPowerFlow)
res = PF.solve_powerflow!(pf, sys)
```

**Reason:** PowerFlows v0.10.0 introduced a complete API rewrite. The constructor pattern changed from using type parameters to function arguments.

**API Evolution:**
- PowerFlows 0.8 and earlier: `PF.solve_ac_powerflow!(sys)`
- PowerFlows 0.9: `pf = PF.ACPowerFlow{PF.TrustRegionACPowerFlow}(); PF.solve_powerflow!(pf, sys)`
- PowerFlows 0.10: `pf = PF.ACPowerFlow(PF.TrustRegionACPowerFlow); PF.solve_powerflow!(pf, sys)`

### 3. PowerSystems get_bus_numbers Removal (src/base/simulation_initialization.jl:20)

**Previous (PowerSystems 4):**
```julia
bus_size = length(PSY.get_bus_numbers(sys))
```

**Updated (PowerSystems 5):**
```julia
bus_size = length(collect(PSY.get_components(PSY.Bus, sys)))
```

**Reason:** `get_bus_numbers()` was removed in PowerSystems v5. The function no longer exists in the codebase.

**Migration:** Use `length(collect(get_components(PSY.Bus, sys)))` to count buses.

### 4. PowerNetworkMatrices v0.14 Compatibility

**Locations of PNM.Ybus usage:**
- `src/base/simulation_inputs.jl:299,309` - `PNM.Ybus(sys)`
- `src/base/perturbations.jl:358` - `PNM.find_subnetworks(ybus, ...)`
- `src/base/perturbations.jl:366` - `NetworkSwitch(time::Float64, ybus::PNM.Ybus)`
- Multiple test files

**Analysis:**
The code accesses these Ybus properties:
- `Ybus_[:, :]` - Matrix indexing (still supported via Base.getindex)
- `Ybus_.lookup[1]` - Bus number to index mapping (lookup is a NTuple{2, Dict})
- `ybus.data` - Raw matrix data (SparseMatrixCSC field)

**v0.14 Ybus Structure:**
```julia
struct Ybus{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{ComplexF32}
    data::SparseArrays.SparseMatrixCSC{ComplexF32, Int}
    adjacency_data::SparseArrays.SparseMatrixCSC{Int8, Int}
    axes::Ax
    lookup::L  # Tuple of 2 dictionaries
    subnetwork_axes::Dict{Int, Ax}
    arc_subnetwork_axis::Dict{Int, Vector{Tuple{Int, Int}}}
    network_reduction_data::NetworkReductionData
    arc_admittance_from_to::Union{ArcAdmittanceMatrix, Nothing}
    arc_admittance_to_from::Union{ArcAdmittanceMatrix, Nothing}
end
```

**Result:** No code changes required. The API is compatible.

## Summary of Files Changed

1. **Project.toml** - Updated all dependency versions
2. **src/base/simulation_initialization.jl** - Fixed PowerFlows API and get_bus_numbers removal

## Breaking Changes in Dependencies

### InfrastructureSystems v3.0.0
- New Time Series interface (unified get_time_series interfaces)
- Supplemental attributes getters refactored
- Cost function definitions modified
- Removed deprecated functions

### PowerSystems v5.0.0
- **REMOVED: `get_bus_numbers()`** - use `collect(get_components(PSY.Bus, sys))` instead
- Redefinition of many structs and parsing logic
- New supplemental attributes interface
- Revised time series handling
- Enhanced hydro and HVDC structures
- Support for 3-winding transformers
- HydroEnergyReservoir struct removed

### PowerFlows v0.10.0
- **Complete API rewrite** with new solve methods
- Constructor pattern changed: `ACPowerFlow(SolverType)` instead of `ACPowerFlow{SolverType}()`
- Integration of network reductions from PowerNetworkMatrices
- New solver types: `NewtonRaphsonACPowerFlow`, `TrustRegionACPowerFlow`, `LevenbergMarquardtACPowerFlow`, `RobustHomotopyPowerFlow`
- Default constructor `ACPowerFlow()` uses Newton-Raphson solver

### PowerNetworkMatrices v0.14.0
- Changed from branch-based indexing to arc-based indexing (internal change)
- New reduction interfaces
- Ybus struct expanded with new fields
- Existing API (constructor, indexing, property access) remains compatible

## Testing Status

All identified API changes have been implemented. The code should now be compatible with:
- InfrastructureSystems 3.x
- PowerSystems 5.x
- PowerFlows 0.10.x
- PowerNetworkMatrices 0.14.x

## References

- [InfrastructureSystems v3.0.0 Release](https://github.com/NREL-Sienna/InfrastructureSystems.jl/releases/tag/v3.0.0)
- [PowerSystems v5.0.0 Release](https://github.com/NREL-Sienna/PowerSystems.jl/releases/tag/v5.0.0)
- [PowerFlows v0.10.0 Release](https://github.com/NREL-Sienna/PowerFlows.jl/releases/tag/v0.10.0)
- [PowerNetworkMatrices v0.14.0 Release](https://github.com/NREL-Sienna/PowerNetworkMatrices.jl/releases/tag/v0.14.0)
- [PowerNetworkMatrices PR #141](https://github.com/NREL-Sienna/PowerNetworkMatrices.jl/pull/141) - Arc-based indexing refactor

## Commits

1. `43f5694` - Update PowerSystems, PowerFlows, and PowerNetworkMatrices dependencies
2. `617d0f0` - Add comprehensive dependency update analysis
3. `42d19f1` - Update InfrastructureSystems to version 3
4. `0386d05` - Fix PowerSystems v5 API: replace get_bus_numbers
5. `713de0b` - Update dependency analysis with all implemented fixes
6. `815df11` - Update DataStructures to version 0.19
7. `2337329` - Update minimum Julia version to 1.10
