using PowerSimulationsDynamics
using Sundials

"""
Case 23:
This case study a 15-state droop grid forming inverter against an infinite bus located at bus 1, with inverter located at bus 2.
The perturbation increase the reference power (analogy for mechanical power) from 0.5 to 0.7.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test23.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 4.0);

#PSCAD benchmark data
csv_file = joinpath(TEST_FILES_DIR, "benchmarks/pscad/Test23/Test23_theta.csv")
t_offset = 9.0

#Define Fault using Callbacks
Pref_change = ControlReferenceChange(1.0, case_inv, PSID.P_ref_index, 0.7)

@testset "Test 23 Droop Inverter ResidualModel" begin
    path = (joinpath(pwd(), "test-23"))
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

        #Solve problem
        execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005)

        #Obtain data for angles
        series = get_state_series(res, ("generator-102-1", :θ_oc))
        t = series[1]
        θ = series[2]

        #Obtain PSCAD benchmark data
        M = get_csv_data(csv_file)
        t_pscad = M[:, 1] .- t_offset
        θ_pscad = M[:, 2]

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test23_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test LinearAlgebra.norm(eigs - test23_eigvals) < 1e-3
        @test res.solution.retcode == :Success
        @test LinearAlgebra.norm(θ - θ_pscad) <= 3e-2
        @test LinearAlgebra.norm(t - round.(t_pscad, digits = 3)) == 0.0

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 23 Droop Inverter MassMatrixModel" begin
    path = (joinpath(pwd(), "test-23"))
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

        #Solve problem
        execute!(sim, Rodas5(), dtmax = 0.005, saveat = 0.005)

        #Obtain data for angles
        series = get_state_series(res, ("generator-102-1", :θ_oc))
        t = series[1]
        θ = series[2]

        #Obtain PSCAD benchmark data
        M = get_csv_data(csv_file)
        t_pscad = M[:, 1] .- t_offset
        θ_pscad = M[:, 2]

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test23_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test LinearAlgebra.norm(eigs - test23_eigvals) < 1e-3
        @test res.solution.retcode == :Success
        @test LinearAlgebra.norm(θ - θ_pscad) <= 3e-2
        @test LinearAlgebra.norm(t - round.(t_pscad, digits = 3)) == 0.0

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
