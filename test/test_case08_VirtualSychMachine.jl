using PowerSimulationsDynamics
using Sundials

"""
Case 8:
This case study a 19-state virtual synchronous machine against an infinite bus located at bus 1, with VSM located at bus 2.
The perturbation increase the reference power (analogy for mechanical power) from 0.5 to 0.7.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test08.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 4.0);

#PSCAD benchmark data
csv_file = joinpath(dirname(@__FILE__), "benchmarks/pscad/Test08/Test08_omega.csv")
t_offset = 9.0

#Define Fault using Callbacks
Pref_change = ControlReferenceChange(1.0, case_inv, PSID.P_ref_index, 0.7)

@testset "Test 08 VSM Inverter Infinite Bus ImplicitModel" begin
    path = (joinpath(pwd(), "test-08"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            ImplicitModel,
            omib_sys, # system
            path,
            tspan,
            Pref_change,
        )

        small_sig = small_signal_analysis(sim)
        @test small_sig.stable

        #Solve problem
        execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005)
        test_case08_VirtualInertia
        #Obtain frequency data
        series = get_state_series(sim, ("generator-102-1", :ω_oc))
        t = series[1]
        ω = series[2]

        #Obtain PSCAD benchmark data
        M = get_csv_data(csv_file)
        t_pscad = M[:, 1] .- t_offset
        ω_pscad = M[:, 2]

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test08_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test sim.solution.retcode == :Success

        power = PSID.get_activepower_series(sim, "generator-102-1")
        rpower = PSID.get_reactivepower_series(sim, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
        @test LinearAlgebra.norm(ω - ω_pscad) <= 1e-4
        @test LinearAlgebra.norm(t - round.(t_pscad, digits = 3)) == 0.0

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 08 VSM Inverter Infinite Bus MassMatrixModel" begin
    path = (joinpath(pwd(), "test-08"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            omib_sys, # system
            path,
            tspan,
            Pref_change,
        )

        # small_sig = small_signal_analysis(sim)
        # @test small_sig.stable

        #Solve problem
        execute!(sim, Rodas5(), dtmax = 0.005, saveat = 0.005)

        #Obtain frequency data
        series = get_state_series(sim, ("generator-102-1", :ω_oc))
        t = series[1]
        ω = series[2]

        #Obtain PSCAD benchmark data
        M = get_csv_data(csv_file)
        t_pscad = M[:, 1] .- t_offset
        ω_pscad = M[:, 2]

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test08_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test sim.solution.retcode == :Success

        power = PSID.get_activepower_series(sim, "generator-102-1")
        rpower = PSID.get_reactivepower_series(sim, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
        @test LinearAlgebra.norm(ω - ω_pscad) <= 1e-4
        @test LinearAlgebra.norm(t - round.(t_pscad, digits = 3)) == 0.0

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
