"""
Case 12:
This case study a three bus system with 2 machines (Classic Model - Single Shaft: 2 State model) without loads.
The machine at bus 1 is used as a reference machine, while machine at bus 2 has a simplified droop governor (TGTypeII).
The perturbation trips four (out of 5) circuits of line between buses 1 and 2, multiplying by 4 its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test12.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

# Time span
tspan = (0.0, 5.0)

# Define Fault: Change of YBus
Ybus_change = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

@testset "Test 12 Multi Machine ResidualModel" begin
    path = (joinpath(pwd(), "test-12"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation(
            ResidualModel,
            threebus_sys, #system,
            path,
            Ybus_change, #Type of Fault
        )

        # Test Initial Condition
        diffvals = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test12_x0_init
            diffvals[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diffvals[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test12_eigvals) < 1e-3

        # Run simulation
        @test execute!(
            sim, #simulation structure
            IDA(),#Sundials DAE Solver
            tspan,
            dtmax = 0.02, #keywords arguments
        ) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)
        series = get_state_series(results, ("generator-102-1", :ω))
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 12 Multi Machine MassMatrixModel" begin
    path = (joinpath(pwd(), "test-12"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation(
            MassMatrixModel,
            threebus_sys, #system,
            path,
            Ybus_change, #Type of Fault
        )

        # Test Initial Condition
        diffvals = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test12_x0_init
            diffvals[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diffvals[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test12_eigvals) < 1e-3

        # Run simulation
        @test execute!(
            sim, #simulation structure
            Rodas4(),#Sundials DAE Solver
            tspan,
            dtmax = 0.02, #keywords arguments
        ) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)
        series = get_state_series(results, ("generator-102-1", :ω))
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
