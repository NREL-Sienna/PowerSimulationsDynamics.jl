"""
Case 1:
This case study defines a classical machine against an infinite bus. The fault
drop a circuit on the (double circuit) line connecting the two buses, doubling its impedance
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test01.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################
# Define Fault: Change of YBus
Ybus_change = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

@testset "Test 01 OMIB ResidualModel" begin
    path = (joinpath(pwd(), "test-01"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            omib_sys, #system
            path,
            Ybus_change,
        )

        # Test Initial Condition
        diffvals = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test01_x0_init
            diffvals[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diffvals[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test01_eigvals) < 1e-3
        @test LinearAlgebra.norm(eigs - test01_eigvals_psat, Inf) < 5.0

        # Solve problem
        @test execute!(sim, IDA(), (0.0, 20.0), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for angles
        series = get_state_series(results, ("generator-102-1", :δ))
        t = series[1]
        δ = series[2]

        series2 = get_voltage_magnitude_series(results, 102)
        series3 = get_voltage_angle_series(results, 102)

        #Obtain PSAT and PSS/e benchmark data
        psat_csv = joinpath(TEST_FILES_DIR, "benchmarks/psat/Test01/Test01_delta.csv")
        psse_csv = joinpath(TEST_FILES_DIR, "benchmarks/psse/Test01/Test01_delta.csv")
        t_psat, δ_psat = get_csv_delta(psat_csv)
        t_psse, δ_psse = get_csv_delta(psse_csv)

        # Test Transient Simulation Results
        # We test that the time series have the same number of items for convenience.
        @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0
        # PSSE results are in Degrees
        @test LinearAlgebra.norm(δ - (δ_psse .* pi / 180), Inf) <= 2e-3
        @test LinearAlgebra.norm(t - t_psat) == 0.0
        @test LinearAlgebra.norm(δ - δ_psat, Inf) <= 1e-3

        power = get_activepower_series(results, "generator-102-1")
        rpower = get_reactivepower_series(results, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 01 OMIB MassMatrixModel" begin
    path = (joinpath(pwd(), "test-01"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            omib_sys, #system
            path,
            Ybus_change,
        )

        # Test Initial Condition
        diffvals = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test01_x0_init
            diffvals[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diffvals[1] < 1e-3)

        # Test Initial Condition
        diffvals = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test01_x0_init
            diffvals[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diffvals[1] < 1e-3)

        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        #Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test01_eigvals) < 1e-3
        @test LinearAlgebra.norm(eigs - test01_eigvals_psat, Inf) < 5.0
        #Test Small Signal
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test01_eigvals) < 1e-3
        @test LinearAlgebra.norm(eigs - test01_eigvals_psat, Inf) < 5.0

        # Solve problem
        @test execute!(sim, Rodas4(), (0.0, 20.0), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for angles
        series = get_state_series(results, ("generator-102-1", :δ))
        t = series[1]
        δ = series[2]

        series2 = get_voltage_magnitude_series(results, 102)
        series3 = get_voltage_angle_series(results, 102)

        # Obtain PSAT and PSS/e benchmark data
        psat_csv = joinpath(TEST_FILES_DIR, "benchmarks/psat/Test01/Test01_delta.csv")
        psse_csv = joinpath(TEST_FILES_DIR, "benchmarks/psse/Test01/Test01_delta.csv")
        t_psat, δ_psat = get_csv_delta(psat_csv)
        t_psse, δ_psse = get_csv_delta(psse_csv)

        # Test Transient Simulation Results
        @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0
        # PSSE results are in Degrees
        @test LinearAlgebra.norm(δ - (δ_psse .* pi / 180), Inf) <= 2e-3
        @test LinearAlgebra.norm(t - t_psat) == 0.0
        @test LinearAlgebra.norm(δ - δ_psat, Inf) <= 1e-3

        power = get_activepower_series(results, "generator-102-1")
        rpower = get_reactivepower_series(results, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
