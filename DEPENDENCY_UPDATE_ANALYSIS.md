# Dependency Update Analysis

## Overview
This document analyzes the changes made to update PowerSimulationsDynamics.jl to be compatible with the latest versions of PowerSystems, PowerFlows, and PowerNetworkMatrices.

## Version Changes

| Package | Previous | New | Release Date | Breaking Changes |
|---------|----------|-----|--------------|------------------|
| PowerSystems | 4 | 5 | Nov 2024 | Yes - struct redefinitions, time series changes |
| PowerFlows | 0.9 | 0.10 | Nov 4, 2024 | Yes - complete API rewrite |
| PowerNetworkMatrices | 0.12.1 | 0.14 | Nov 4, 2024 | Yes - arc-based indexing |

## Code Changes Made

### 1. PowerFlows API Update (src/base/simulation_initialization.jl:14)

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

### 2. PowerNetworkMatrices Usage

**Locations of PNM.Ybus usage:**
- `src/base/simulation_inputs.jl:299,309` - `PNM.Ybus(sys)`
- `src/base/perturbations.jl:358` - `PNM.find_subnetworks(ybus, ...)`
- `src/base/perturbations.jl:366` - `NetworkSwitch(time::Float64, ybus::PNM.Ybus)`
- Multiple test files

**Analysis:**
The code accesses these Ybus properties:
- `Ybus_[:, :]` - Matrix indexing
- `Ybus_.lookup[1]` - Bus number to index mapping
- `ybus.data` - Raw matrix data

**v0.14 Changes:**
- Eliminated branch-based indexing in favor of arc-based indexing
- This change is primarily internal to PowerNetworkMatrices
- The Ybus constructor API appears to remain stable
- Property access patterns (`.data`, `.lookup`) need verification via testing

**Risk Assessment:** MEDIUM - The basic API appears unchanged, but property access patterns may have been affected.

## Testing Requirements

Since Julia cannot be installed in this environment, testing will occur via GitHub Actions CI when the PR is opened.

### Critical Test Areas:

1. **Power Flow Initialization**
   - File: `src/base/simulation_initialization.jl:9-30`
   - Tests: All test cases that call `Simulation()` constructor
   - Verify: Power flow convergence with new ACPowerFlow API

2. **Ybus Construction and Usage**
   - Files: `src/base/simulation_inputs.jl`, `src/base/perturbations.jl`
   - Tests: Network switch perturbations, system initialization
   - Verify: `.data`, `.lookup[1]`, and `[:, :]` access patterns still work

3. **System Building**
   - Tests: All data tests in `test/data_tests/`
   - Verify: PowerSystems v5 struct changes don't break system construction

4. **Dynamic Simulations**
   - Tests: Full simulation tests (OMIB, three-bus systems, etc.)
   - Verify: End-to-end functionality with all three updated dependencies

## Expected CI Workflow

When the PR is created, GitHub Actions will:

1. Install Julia 1.x (latest stable)
2. Run `julia-buildpkg` to install all dependencies including:
   - PowerSystems 5.x
   - PowerFlows 0.10.x
   - PowerNetworkMatrices 0.14.x
3. Run `julia-runtest` to execute the full test suite
4. Report any failures

## Potential Issues to Watch For

### High Priority:
1. **Ybus property access** - If `.data` or `.lookup` structure changed in v0.14
2. **PowerSystems struct changes** - v5.0 redefined many structs

### Medium Priority:
3. **Time series handling** - PowerSystems v5 changed time series interfaces
4. **Solver convergence** - Different default settings in PowerFlows v0.10

### Low Priority:
5. **Hydro/HVDC components** - If test systems use these features
6. **3-winding transformers** - New handling in PowerSystems v5

## Verification Steps Post-CI

If CI tests fail:

1. Check error messages for specific API mismatches
2. Review PowerNetworkMatrices v0.14 struct definitions for Ybus
3. Check PowerSystems v5 migration guide for struct changes
4. Update property access patterns as needed

## References

- [PowerSystems v5.0.0 Release](https://github.com/NREL-Sienna/PowerSystems.jl/releases/tag/v5.0.0)
- [PowerFlows v0.10.0 Release](https://github.com/NREL-Sienna/PowerFlows.jl/releases/tag/v0.10.0)
- [PowerNetworkMatrices v0.14.0 Release](https://github.com/NREL-Sienna/PowerNetworkMatrices.jl/releases/tag/v0.14.0)
- [PowerNetworkMatrices PR #141](https://github.com/NREL-Sienna/PowerNetworkMatrices.jl/pull/141) - Arc-based indexing refactor
