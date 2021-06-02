"""
Case 6:
This case study a three bus system with 2 machines (Anderson-Fouad: 6th order model) and an infinite source.
The fault drop the connection between buses 1 and 3, eliminating the direct connection between the infinite source
and the generator located in bus 3.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test06.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 200.0);
#Define Fault: Change of YBus
Ybus_change = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

@testset "Test 06 Anderson ImplicitModel" begin
    path = (joinpath(pwd(), "test-06"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            ImplicitModel,
            threebus_sys, #system,
            path,
            tspan, #time span
            Ybus_change, #Type of Fault
        )

        small_sig = small_signal_analysis(sim)
        @test small_sig.stable

        #Solve problem in equilibrium
        execute!(sim, IDA(), dtmax = 0.02)

        #Obtain data for angles
        series = get_state_series(sim, ("generator-102-1", :δ))

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test06_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test sim.solution.retcode == :Success

        power = PSID.get_activepower_series(sim, "generator-102-1")
        rpower = PSID.get_reactivepower_series(sim, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 06 Anderson MassMatrixModel" begin
    path = (joinpath(pwd(), "test-06"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            threebus_sys, #system,
            path,
            tspan, #time span
            Ybus_change, #Type of Fault
        )

        # small_sig = small_signal_analysis(sim)
        # @test small_sig.stable

        #Solve problem in equilibrium
        execute!(sim, Rodas5(), dtmax = 0.02)

        #Obtain data for angles
        series = get_state_series(sim, ("generator-102-1", :δ))

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test06_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test sim.solution.retcode == :Success

        power = PSID.get_activepower_series(sim, "generator-102-1")
        rpower = PSID.get_reactivepower_series(sim, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
