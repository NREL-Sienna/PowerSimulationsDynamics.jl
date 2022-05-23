using PowerSimulationsDynamics
using Sundials

"""
Case 25:
This case study a three-bus system, with two Marconato generators (in buses 1 and 2), and a constant impedance load in bus 3.
The perturbation increase the reference of mechanical power of generator-2 from 0.8 to 0.9 at t=1.0s.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test25.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

# time span
tspan = (0.0, 40.0);

# PSCAD benchmark data
csv_file = joinpath(TEST_FILES_DIR, "benchmarks/pscad/Test25/Test25_v102.csv")
t_offset = 49.0

# Define Fault using Callbacks
gen2 = get_dynamic_injector(get_component(Generator, sys, "generator-102-1"));
Pref_change = ControlReferenceChange(1.0, gen2, :P_ref, 0.9);

@testset "Test 25 Marconato with Dynamic Lines ResidualModel" begin
    path = (joinpath(pwd(), "test-25"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation!(ResidualModel, sys, path, tspan, Pref_change)

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test25_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        @test small_sig.stable

        # Solve problem in equilibrium
        @test execute!(sim, Sundials.IDA(), dtmax = 0.01, saveat = 0.01) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain voltage magnitude data
        series = get_voltage_magnitude_series(results, 102)
        t = series[1]
        v = series[2]

        # Obtain benchmark data from PSCAD
        M = get_csv_data(csv_file)
        t_pscad = M[:, 1] .- t_offset
        v_pscad = M[:, 2]

        # Relaxed constraint to account for mismatch in damping
        @test LinearAlgebra.norm(v - v_pscad) <= 0.1
        @test LinearAlgebra.norm(t - round.(t_pscad, digits = 3)) == 0.0
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 25 Marconato with Dynamic Lines MassMatrixModel" begin
    path = (joinpath(pwd(), "test-25"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation!(MassMatrixModel, sys, path, tspan, Pref_change)

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test25_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        @test small_sig.stable

        # Solve problem in equilibrium
        @test execute!(sim, Rodas4(), dtmax = 0.01, saveat = 0.01) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain voltage magnitude data
        series = get_voltage_magnitude_series(results, 102)
        t = series[1]
        v = series[2]

        # Obtain benchmark data from PSCAD
        M = get_csv_data(csv_file)
        t_pscad = M[:, 1] .- t_offset
        v_pscad = M[:, 2]

        # Relaxed constraint to account for mismatch in damping
        @test LinearAlgebra.norm(v - v_pscad) <= 0.1
        @test LinearAlgebra.norm(t - round.(t_pscad, digits = 3)) == 0.0
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
