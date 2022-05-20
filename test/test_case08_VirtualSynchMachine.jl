using PowerSimulationsDynamics
using Sundials

"""
Case 8:
This case study a 19-state virtual synchronous machine against an infinite bus located at bus 1, with VSM located at bus 2.
The perturbation increase the reference power (analogy for mechanical power) from 0.5 to 0.7.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test08.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#PSCAD benchmark data
csv_file = joinpath(TEST_FILES_DIR, "benchmarks/pscad/Test08/Test08_omega.csv")
t_offset = 9.0

#Define Fault using Callbacks
Pref_change = ControlReferenceChange(1.0, case_inv, :P_ref, 0.7)

@testset "Test 08 VSM Inverter Infinite Bus ResidualModel" begin
    path = (joinpath(pwd(), "test-08"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation(
            ResidualModel,
            omib_sys, # system
            path,
            (0.0, 4.0),
            Pref_change,
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test08_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test08_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain frequency data
        series = get_state_series(results, ("generator-102-1", :ω_oc))
        t = series[1]
        ω = series[2]

        # Should return zeros and a warning
        series3 = get_field_current_series(results, "generator-102-1")
        series4 = get_field_voltage_series(results, "generator-102-1")

        # Obtain PSCAD benchmark data
        M = get_csv_data(csv_file)
        t_pscad = M[:, 1] .- t_offset
        ω_pscad = M[:, 2]

        power = PSID.get_activepower_series(results, "generator-102-1")
        rpower = PSID.get_reactivepower_series(results, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
        @test LinearAlgebra.norm(ω - ω_pscad) <= 1e-4
        @test LinearAlgebra.norm(t - round.(t_pscad, digits = 3)) == 0.0

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 08 VSM Inverter Infinite Bus MassMatrixModel" begin
    path = (joinpath(pwd(), "test-08"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation(
            MassMatrixModel,
            omib_sys, # system
            path,
            (0.0, 4.0),
            Pref_change,
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test08_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test08_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, Rodas5(), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain frequency data
        series = get_state_series(results, ("generator-102-1", :ω_oc))
        t = series[1]
        ω = series[2]

        # Obtain PSCAD benchmark data
        M = get_csv_data(csv_file)
        t_pscad = M[:, 1] .- t_offset
        ω_pscad = M[:, 2]

        power = PSID.get_activepower_series(results, "generator-102-1")
        rpower = PSID.get_reactivepower_series(results, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
        @test LinearAlgebra.norm(ω - ω_pscad) <= 1e-4
        @test LinearAlgebra.norm(t - round.(t_pscad, digits = 3)) == 0.0

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
