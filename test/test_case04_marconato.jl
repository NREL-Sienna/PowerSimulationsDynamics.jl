"""
Case 4:
This case study a three bus system with 2 machines (Marconato: 8th order model) and an infinite source.
The fault drop the connection between buses 1 and 3, eliminating the direct connection between the infinite source
and the generator located in bus 3.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test04.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#Define Fault: Change of YBus
Ybus_change = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

@testset "Test 04 Marconato ImplicitModel" begin
    path = (joinpath(pwd(), "test-04"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            ImplicitModel,
            threebus_sys, #system,
            path,
            (0.0, 20.0), #time span
            Ybus_change, #Type of Fault
        ) #initial guess

        small_sig = small_signal_analysis(sim)
        @test small_sig.stable

        #Solve problem in equilibrium
        execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005)

        #Obtain data for angles
        series = get_state_series(sim, ("generator-102-1", :δ))
        t = series[1]
        δ = series[2]

        #Obtain PSAT benchmark data
        psat_csv = joinpath(dirname(@__FILE__), "benchmarks/psat/Test03/Test03_delta.csv")
        t_psat, δ_psat = get_csv_delta(psat_csv)

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test04_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test sim.solution.retcode == :Success
        @test LinearAlgebra.norm(t - t_psat) == 0.0
        @test LinearAlgebra.norm(δ - δ_psat, Inf) <= 1e-3

        power = PSID.get_activepower_series(sim, "generator-102-1")
        rpower = PSID.get_reactivepower_series(sim, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 04 Marconato MassMatrixModel" begin
    path = (joinpath(pwd(), "test-04"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            threebus_sys, #system,
            path,
            (0.0, 20.0), #time span
            Ybus_change, #Type of Fault
        ) #initial guess

        #small_sig = small_signal_analysis(sim)
        #@test small_sig.stable

        #Solve problem in equilibrium
        execute!(sim, Rodas5(), dtmax = 0.005, saveat = 0.005)

        #Obtain data for angles
        series = get_state_series(sim, ("generator-102-1", :δ))
        t = series[1]
        δ = series[2]

        #Obtain PSAT benchmark data
        psat_csv = joinpath(dirname(@__FILE__), "benchmarks/psat/Test03/Test03_delta.csv")
        t_psat, δ_psat = get_csv_delta(psat_csv)

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test04_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test sim.solution.retcode == :Success
        @test LinearAlgebra.norm(t - t_psat) == 0.0
        @test LinearAlgebra.norm(δ - δ_psat, Inf) <= 1e-3

        power = PSID.get_activepower_series(sim, "generator-102-1")
        rpower = PSID.get_reactivepower_series(sim, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
