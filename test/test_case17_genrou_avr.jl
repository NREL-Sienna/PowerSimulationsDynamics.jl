"""
Case GENROU + AVR Type I:
This case study defines a three bus system with an infinite bus, GENROU (with AVR TypeI) and a load.
The fault drop the line connecting the infinite bus and GENROU.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test17.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

@testset "Test 17 GENROU AVR ImplicitModel" begin
    path = (joinpath(pwd(), "test-psse-genrou-avr"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            ImplicitModel,
            sys, #system
            path,
            (0.0, 30.0), #time span
            BranchTrip(1.0, "BUS 1-BUS 2-i_1"), #Type of Fault,
        ) #Type of Fault

        #Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        @test small_sig.stable

        #Solve problem in equilibrium
        execute!(sim, IDA(), dtmax = 0.01)

        #Obtain data for angles
        series = get_state_series(sim, ("generator-102-1", :δ))

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test_psse_genrou_avr_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        #Test Initial Condition
        @test (diff[1] < 1e-3)
        #Test Solution DiffEq
        @test sim.solution.retcode == :Success

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 17 GENROU AVR MassMatrixModel" begin
    path = (joinpath(pwd(), "test-psse-genrou-avr"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            sys, #system
            path,
            (0.0, 30.0), #time span
            BranchTrip(1.0, "BUS 1-BUS 2-i_1"), #Type of Fault,
        ) #Type of Fault

        #Obtain small signal results for initial conditions
        # small_sig = small_signal_analysis(sim)
        # @test small_sig.stable

        #Solve problem in equilibrium
        execute!(sim, Rodas5(), dtmax = 0.01)

        #Obtain data for angles
        series = get_state_series(sim, ("generator-102-1", :δ))

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test_psse_genrou_avr_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        #Test Initial Condition
        @test (diff[1] < 1e-3)
        #Test Solution DiffEq
        @test sim.solution.retcode == :Success

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
