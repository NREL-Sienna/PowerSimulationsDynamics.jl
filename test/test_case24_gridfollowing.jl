using PowerSimulationsDynamics
using Sundials

"""
Case 24:
This case study a 15-state grid following inverter against an infinite bus located at bus 1, with the inverter located at bus 2.
The perturbation increase the reference power (analogy for mechanical power) from 0.5 to 0.7.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test24.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 2.0);

#Define Fault using Callbacks
Pref_change = ControlReferenceChange(1.0, case_inv, PSID.P_ref_index, 0.7)

@testset "Test 24 Grid Following Inverter ImplicitModel" begin
    path = (joinpath(pwd(), "test-24"))
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
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        #Solve problem in equilibrium
        execute!(sim, Sundials.IDA(), dtmax = 0.001)

        #Obtain frequency data
        series = get_state_series(sim, ("generator-102-1", :p_oc))
        t = series[1]
        p = series[2]

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test24_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test LinearAlgebra.norm(eigs - test24_eigvals) < 1e-3
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


@testset "Test 24 Grid Following Inverter ImplicitModel" begin
    path = (joinpath(pwd(), "test-24"))
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

        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        #Solve problem in equilibrium
        execute!(sim, Rodas5(), dtmax = 0.001)

        #Obtain frequency data
        series = get_state_series(sim, ("generator-102-1", :p_oc))
        t = series[1]
        p = series[2]

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test24_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test LinearAlgebra.norm(eigs - test24_eigvals) < 1e-3
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
