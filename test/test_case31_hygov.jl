"""
Validation PSSE/HYGOV:
This case study defines a three bus system with an infinite bus, GENROU+SEXS+HYGOV and a load.
The fault drop the line connecting the infinite bus and GENROU
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/HYGOV/ThreeBusMulti.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/HYGOV/ThreeBus_HYGOV.dyr")
csv_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/HYGOV/HYGOV_RESULTS.csv")

@testset "Test 31 HYGOV ResidualModel" begin
    path = (joinpath(pwd(), "test-psse-gast"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)

        # Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            sys, #system
            path,
            (0.0, 20.0), #time span
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault
        )

        # Test Initial Condition
        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test31_hygov_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test31_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for angles
        t, δ = get_state_series(results, ("generator-102-1", :δ))
        _, Vt = get_voltage_magnitude_series(results, 102)
        _, ω = get_state_series(results, ("generator-102-1", :ω))

        # Obtain PSSE results
        M = get_csv_data(csv_file)
        t_psse = M[:, 1]
        Vt_psse = M[:, 2]
        δ_psse = M[:, 3]
        ω_psse = M[:, 4] .+ 1.0

        # Test Transient Simulation Results
        # PSSE results are in Degrees
        @test LinearAlgebra.norm(δ - (δ_psse .* pi / 180), Inf) <= 1e-2
        @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0
        @test LinearAlgebra.norm(Vt - Vt_psse, Inf) <= 1e-3
        @test LinearAlgebra.norm(ω - ω_psse, Inf) <= 1e-3
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 31 HYGOV MassMatrixModel" begin
    path = (joinpath(pwd(), "test-psse-gast"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)

        # Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            sys, #system
            path,
            (0.0, 20.0), #time span
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault
        )

        # Test Initial Condition
        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test31_hygov_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test31_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, Rodas4(), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for angles
        series = get_state_series(results, ("generator-102-1", :δ))
        t = series[1]
        δ = series[2]

        # Obtain PSSE results
        #t_psse, δ_psse = get_csv_delta(csv_file)

        # Test Transient Simulation Results
        # PSSE results are in Degrees
        #@test LinearAlgebra.norm(δ - (δ_psse .* pi / 180), Inf) <= 1e-2
        #@test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
