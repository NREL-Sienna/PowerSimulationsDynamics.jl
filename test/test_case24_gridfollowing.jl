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

#PSCAD benchmark data
csv_file = joinpath(dirname(@__FILE__), "benchmarks/pscad/Test24/Test24_p.csv")
t_offset = 9.0

#Define Fault using Callbacks
Pref_change = ControlReferenceChange(1.0, case_inv, PSID.P_ref_index, 0.7)

@testset "Test 24 Grid Following Inverter ResidualModel" begin
    path = (joinpath(pwd(), "test-24"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            omib_sys, # system
            path,
            tspan,
            Pref_change,
        )

        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        #Solve problem in equilibrium
        execute!(sim, Sundials.IDA(), dtmax = 0.001, saveat = 0.005)

        #Obtain frequency data
        series = get_state_series(res, ("generator-102-1", :p_oc))
        t = series[1]
        p = series[2]

        #Obtain PSCAD benchmark data
        M = get_csv_data(csv_file)
        t_pscad = M[:, 1] .- t_offset
        p_pscad = M[:, 2]

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test24_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test LinearAlgebra.norm(eigs - test24_eigvals) < 1e-3
        @test res.solution.retcode == :Success

        power = PSID.get_activepower_series(res, "generator-102-1")
        rpower = PSID.get_reactivepower_series(res, "generator-102-1")
        ir = PSID.get_real_current_series(res, "generator-102-1")
        ii = PSID.get_imaginary_current_series(res, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(ir, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(ii, Tuple{Vector{Float64}, Vector{Float64}})
        @test LinearAlgebra.norm(p - p_pscad) <= 5e-3
        @test LinearAlgebra.norm(t - round.(t_pscad, digits = 3)) == 0.0

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 24 Grid Following Inverter MassMatrixModel" begin
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
        execute!(sim, Rodas5(), dtmax = 0.001, saveat = 0.005)

        #Obtain frequency data
        series = get_state_series(res, ("generator-102-1", :p_oc))
        t = series[1]
        p = series[2]

        #Obtain PSCAD benchmark data
        M = get_csv_data(csv_file)
        t_pscad = M[:, 1] .- t_offset
        p_pscad = M[:, 2]

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test24_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test LinearAlgebra.norm(eigs - test24_eigvals) < 1e-3
        @test res.solution.retcode == :Success

        power = PSID.get_activepower_series(res, "generator-102-1")
        rpower = PSID.get_reactivepower_series(res, "generator-102-1")
        ir = PSID.get_real_current_series(res, "generator-102-1")
        ii = PSID.get_imaginary_current_series(res, "generator-102-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(ir, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(ii, Tuple{Vector{Float64}, Vector{Float64}})
        @test LinearAlgebra.norm(p - p_pscad) <= 5e-3
        @test LinearAlgebra.norm(t - round.(t_pscad, digits = 3)) == 0.0

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
