"""
Case 13:
This case study a three bus system with 2 machines (One d- One q-: 4th order model) and an infinite source.
The case is similar to case 04, with different AVR and TG models.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test13.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

# Time span


#Define Fault: Change of YBus
Ybus_change = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

@testset "Test 13 AVR ResidualModel" begin
    path = (joinpath(pwd(), "test-13"))
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
        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test13_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test13_eigvals) < 1e-3

        #Solve problem
        @test execute!(sim, IDA(), (0.0, 20.0), dtmax = 0.02) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        #Obtain data for angles
        series = get_state_series(results, ("generator-102-1", :δ))
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 13 AVR MassMarixcModel" begin
    path = (joinpath(pwd(), "test-13"))
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
        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test13_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test13_eigvals) < 1e-3

        #Solve problem
        @test execute!(sim, Rodas4(), (0.0, 20.0), dtmax = 0.02) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        #Obtain data for angles
        series = get_state_series(results, ("generator-102-1", :δ))
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
