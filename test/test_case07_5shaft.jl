"""
Case 7:
This case study a three bus system with 2 machine located at bus 2 and 3.
The generator at bus 3 uses the model of a one d- one q- machine, and has a 5-mass shaft and a turbine governor.
The generator at bus 2 uses the model of a one d- one q- machine, and single mass shaft.
The fault disconnects a circuit between buses 1 and 2, doubling its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test07.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

;

# Define Fault: Change of YBus
Ybus_change = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

@testset "Test 07 5-Mass-shaft model ResidualModel" begin
    path = (joinpath(pwd(), "test-07"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            threebus_sys, #system,
            path,
            Ybus_change, #Type of Fault
        ) #initial guess

        # Test Initial Condition
        diffvals = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test07_x0_init
            diffvals[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diffvals[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test07_eigvals) < 1e-3

        #Solve problem
        @test execute!(sim, IDA(), (0.0, 20.0), dtmax = 0.001) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for angles
        series = get_state_series(results, ("generator-103-1", :δ))
        series2 = get_state_series(results, ("generator-103-1", :δ_hp))
        series3 = get_state_series(results, ("generator-103-1", :δ_ip))
        series4 = get_state_series(results, ("generator-103-1", :δ_ex))
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 07 5-Mass-shaft model MassMatrixModel" begin
    path = (joinpath(pwd(), "test-07"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            threebus_sys, #system,
            path,
            Ybus_change, #Type of Fault
        ) #initial guess

        # Test Initial Condition
        diffvals = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test07_x0_init
            diffvals[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diffvals[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test07_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, Rodas4(), (0.0, 20.0), dtmax = 0.001) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for angles
        series = get_state_series(results, ("generator-103-1", :δ))
        series2 = get_state_series(results, ("generator-103-1", :δ_hp))
        series3 = get_state_series(results, ("generator-103-1", :δ_ip))
        series4 = get_state_series(results, ("generator-103-1", :δ_ex))
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
