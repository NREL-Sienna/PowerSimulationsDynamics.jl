"""
Validation PSSE/TGOV1:
This case study defines a three bus system with an infinite bus, GENROU+AC1A+GAST and a load.
The fault drop the line connecting the infinite bus and GENROU
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

raw_file = joinpath(dirname(@__FILE__), "benchmarks/psse/GAST/ThreeBusMulti.raw")
dyr_file = joinpath(dirname(@__FILE__), "benchmarks/psse/GAST/ThreeBus_GAST.dyr")
csv_file = joinpath(dirname(@__FILE__), "benchmarks/psse/GAST/GAST_TEST.csv")

@testset "Test 21 GAST ImplicitModel" begin
    path = (joinpath(pwd(), "test-psse-gast"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)
        sim = Simulation!(
            ImplicitModel,
            sys, #system
            path,
            (0.0, 20.0), #time span
            BranchTrip(1.0, "BUS 1-BUS 2-i_1"), #Type of Fault
        )

        #Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        #Solve problem
        execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005)

        series = get_state_series(sim, ("generator-102-1", :δ))
        t = series[1]
        δ = series[2]

        #Obtain PSSE results
        t_psse, δ_psse = get_csv_delta(csv_file)

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test_psse_gast_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        #Test Initial Condition
        @test (diff[1] < 1e-3)
        #Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test21_eigvals) < 1e-3
        #Test Solution DiffEq
        @test sim.solution.retcode == :Success
        #Test Transient Simulation Results
        # PSSE results are in Degrees
        @test LinearAlgebra.norm(δ - (δ_psse .* pi / 180), Inf) <= 1e-2
        @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 21 GAST MassMatrixModel" begin
    path = (joinpath(pwd(), "test-psse-gast"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)
        sim = Simulation!(
            MassMatrixModel,
            sys, #system
            path,
            (0.0, 20.0), #time span
            BranchTrip(1.0, "BUS 1-BUS 2-i_1"), #Type of Fault
        )

        #Solve problem
        execute!(sim, Rodas5(autodiff = true), dtmax = 0.005, saveat = 0.005)

        small_sig = small_signal_analysis(sim; reset_simulation = true)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        series = get_state_series(sim, ("generator-102-1", :δ))
        t = series[1]
        δ = series[2]

        #Obtain PSSE results
        t_psse, δ_psse = get_csv_delta(csv_file)

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test_psse_gast_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        #Test Initial Condition
        @test (diff[1] < 1e-3)
        #Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test21_eigvals) < 1e-3
        #Test Solution DiffEq
        @test sim.solution.retcode == :Success
        #Test Transient Simulation Results
        # PSSE results are in Degrees
        @test LinearAlgebra.norm(δ - (δ_psse .* pi / 180), Inf) <= 1e-2
        @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
