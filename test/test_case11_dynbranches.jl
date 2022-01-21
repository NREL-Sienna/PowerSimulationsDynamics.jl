"""
Case 11:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source. The connection between buses 2 and 3 is modeled using dynamic lines.
The perturbation trips two of the three circuits of line between buses 1 and 2, triplicating its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test11.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

# time span
tspan = (0.0, 40.0)

# Define Fault: Change of YBus
Ybus_change = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

@testset "Test 11 Dynamic Branches ResidualModel" begin
    path = (joinpath(pwd(), "test-11"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            threebus_sys, #system,
            path,
            Ybus_change, #Type of Fault
        )
        # Get dictionary of initial conditions
        init_conds = read_initial_conditions(sim)

        # Test Initial Condition
        diffvals = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test10_x0_init
            diffvals[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diffvals[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test11_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, IDA(), tspan) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for voltages
        series = get_voltage_magnitude_series(results, 102)
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 11 Dynamic Branches MassMatrixModel" begin
    path = (joinpath(pwd(), "test-11"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            threebus_sys, #system,
            path,
            Ybus_change, #Type of Fault
        )

        # Test Initial Condition
        diffvals = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test10_x0_init
            diffvals[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diffvals[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test11_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, Rodas5(), tspan) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for voltages
        series = get_voltage_magnitude_series(results, 102)
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
