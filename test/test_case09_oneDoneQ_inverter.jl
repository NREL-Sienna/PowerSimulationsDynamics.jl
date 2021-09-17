"""
Case 9:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source.
The perturbation increase the reference power (analogy for mechanical power) of the inverter from 1.0 to 1.2.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test09.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#time span
tspan = (0.0, 20.0);
case_inv = collect(PSY.get_components(PSY.DynamicInverter, threebus_sys))[1]

#Define Fault using Callbacks
Pref_change = ControlReferenceChange(1.0, case_inv, PSID.P_ref_index, 1.2)

@testset "Test 09 VSM Inverter and OneDoneQ ResidualModel" begin
    path = (joinpath(pwd(), "test-09"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            threebus_sys, # system
            path,
            tspan,
            Pref_change,
        )

        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        #Solve problem
        execute!(sim, IDA(), dtmax = 0.02)

        #Obtain data for angles
        series = get_state_series(res, ("generator-103-1", :θ_oc))

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test09_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test LinearAlgebra.norm(eigs - test09_eigvals) < 1e-3
        @test res.solution.retcode == :Success
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 09 VSM Inverter and OneDoneQ MassMatrixModel" begin
    path = (joinpath(pwd(), "test-09"))
    !isdir(path) && mkdir(path)
    try
        #Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            threebus_sys, # system
            path,
            tspan,
            Pref_change,
        )

        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        #Solve problem
        execute!(sim, Rodas5(), dtmax = 0.02)

        #Obtain data for angles
        series = get_state_series(res, ("generator-103-1", :θ_oc))

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test09_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff[1] < 1e-3)
        @test LinearAlgebra.norm(eigs - test09_eigvals) < 1e-3
        @test res.solution.retcode == :Success

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
