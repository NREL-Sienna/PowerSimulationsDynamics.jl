"""
Case GENROU + AVR Type I:
This case study defines a three bus system with an infinite bus, GENROU (with AVR TypeI) and a load.
The fault drop the line connecting the infinite bus and GENROU.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test17.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

@testset "Test 17 GENROU AVR ResidualModel" begin
    path = (joinpath(pwd(), "test-psse-genrou-avr"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            sys, #system
            path,
            (0.0, 30.0), #time span
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault,
        ) #Type of Fault

        # Test Initial Condition
        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test_psse_genrou_avr_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test17_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, IDA(), (0.0, 20.0), dtmax = 0.01) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for angles
        series = get_state_series(results, ("generator-102-1", :δ))
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 17 GENROU AVR MassMatrixModel" begin
    path = (joinpath(pwd(), "test-psse-genrou-avr"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            sys, #system
            path,
            (0.0, 30.0), #time span
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault,
        ) #Type of Fault

        # Test Initial Condition
        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test_psse_genrou_avr_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test17_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, Rodas4(), (0.0, 20.0), dtmax = 0.01) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for angles
        series = get_state_series(results, ("generator-102-1", :δ))
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
