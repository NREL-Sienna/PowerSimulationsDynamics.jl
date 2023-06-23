"""
Case 9:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source.
The perturbation increase the reference power (analogy for mechanical power) of the inverter from 1.0 to 1.2.
"""

##################################################
############### LOAD DATA ########################
##################################################

threebus_sys = build_system(PSIDTestSystems, "psid_test_threebus_machine_vsm")

##################################################
############### SOLVE PROBLEM ####################
##################################################

@testset "Test 09 VSM Inverter and OneDoneQ ResidualModel" begin
    path = mktempdir()
    #time span
    tspan = (0.0, 20.0)
    case_inv = collect(PSY.get_components(PSY.DynamicInverter, threebus_sys))[1]
    case_gen = collect(PSY.get_components(PSY.DynamicGenerator, threebus_sys))[1]

    #Define Fault using Callbacks
    Pref_change = ControlReferenceChange(1.0, case_inv, :P_ref, 1.2)
    gen_trip = GeneratorTrip(1.5, case_gen)
    try
        #Define Simulation Problem
        sim = Simulation(
            ResidualModel,
            threebus_sys, # system
            path,
            tspan,
            [Pref_change, gen_trip],
        )

        # Test Initial Conditions
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test09_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test09_eigvals) < 1e-3

        #Solve problem
        @test execute!(sim, IDA(); dtmax = 0.02) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        #Obtain data for angles
        t, θ_oc = get_state_series(results, ("generator-103-1", :θ_oc))
        @test length(t) == length(θ_oc)
        t, ω_oc = get_state_series(results, ("generator-103-1", :ω_oc))
        @test length(t) == length(ω_oc)
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test 09 VSM Inverter and OneDoneQ MassMatrixModel" begin
    path = (joinpath(pwd(), "test-09"))
    !isdir(path) && mkdir(path)
    #time span
    tspan = (0.0, 20.0)
    case_inv = collect(PSY.get_components(PSY.DynamicInverter, threebus_sys))[1]
    case_gen = collect(PSY.get_components(PSY.DynamicGenerator, threebus_sys))[1]

    #Define Fault using Callbacks
    Pref_change = ControlReferenceChange(1.0, case_inv, :P_ref, 1.2)
    gen_trip = GeneratorTrip(1.5, case_gen)
    try
        #Define Simulation Problem
        sim = Simulation(
            MassMatrixModel,
            threebus_sys, # system
            path,
            tspan,
            [Pref_change, gen_trip],
        )

        # Test Initial Conditions
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test09_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test09_eigvals) < 1e-3

        #Solve problem
        @test execute!(sim, Rodas4(); dtmax = 0.02) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        #Obtain data for angles
        t, θ_oc = get_state_series(results, ("generator-103-1", :θ_oc))
        @test length(t) == length(θ_oc)
        t, ω_oc = get_state_series(results, ("generator-103-1", :ω_oc))
        @test length(t) == length(ω_oc)

        power = PSID.get_activepower_series(results, "generator-103-1")
        rpower = PSID.get_reactivepower_series(results, "generator-103-1")
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})

    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end
