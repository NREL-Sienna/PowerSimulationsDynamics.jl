"""
Case 5:
This case study a three bus system with 2 machines (Simple Anderson-Fouad: 4th order model) and an infinite source.
The fault drop the connection between buses 1 and 3, eliminating the direct connection between the infinite source
and the generator located in bus 3.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test05.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

# Define Fault: Change of YBus
Ybus_change = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

@testset "Test 05 Simple Anderson ResidualModel" begin
    path = (joinpath(pwd(), "test-05"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation(
            ResidualModel,
            threebus_sys, #system
            path,
            Ybus_change, #Type of Fault
        ) #initial guess

        # Test Initial Condition
        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test05_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test05_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, IDA(), (0.0, 200.0), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        power = PSID.get_activepower_series(results, "generator-102-1")
        rpower = PSID.get_reactivepower_series(results, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 05 Simple Anderson MassMatrixModel" begin
    path = (joinpath(pwd(), "test-05"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation(
            MassMatrixModel,
            threebus_sys, #system
            path,
            Ybus_change, #Type of Fault
        ) #initial guess

        # Test Initial Condition
        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test05_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test05_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, Rodas4(), (0.0, 200.0), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        power = PSID.get_activepower_series(results, "generator-102-1")
        rpower = PSID.get_reactivepower_series(results, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
