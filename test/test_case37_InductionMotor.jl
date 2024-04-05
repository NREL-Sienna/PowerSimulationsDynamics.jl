"""
Test InductionMotor model:
This case study defines a four bus system with an infinite bus in 1,
a GENSAL in bus 2, and the (5th order model) induction motor in bus 3
The GENSAL machine has a SEXS and a HYGOV.
The disturbance is the outage of one line between buses 1 and 4
"""
##################################################
############### LOAD DATA ########################
##################################################

raw_file = joinpath(TEST_FILES_DIR, "data_tests/TVC_System_motor.raw")
dyr_file = joinpath(TEST_FILES_DIR, "data_tests/TVC_System_motor.dyr")

##################################################
############### SOLVE PROBLEM ####################
##################################################

@testset "Test 37 5th_order Ind. Motor" begin
    path = (joinpath(pwd(), "test-ind-mot_5th"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)
        time_span = (0.0, 20.0)
        perturbation_trip = BranchTrip(1.0, Line, "BUS1-BUS4-i_1")
        # Create Simulation with Constant Power
        sim_P = Simulation(ResidualModel, sys, path, time_span, perturbation_trip)

        # Motor parameters
        load = first(get_components(PSY.StandardLoad, sys))
        # Include the induction motor
        dynamic_injector = Ind_Motor(load)
        set_dynamic_injector!(load, dynamic_injector)
        sim = Simulation(ResidualModel, sys, path, time_span, perturbation_trip)

        # Test initial voltages between two systems are equivalent
        voltages_P = sim_P.x0_init[1:8]
        voltages_motor = sim.x0_init[1:8]
        @test LinearAlgebra.norm(voltages_P - voltages_motor) < 1e-4

        # Test Initial Condition
        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test37_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        #Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test37_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, IDA(); dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # No comparison since model is not available in other tools
        power = get_activepower_series(results, PSY.get_name(load))
        rpower = get_reactivepower_series(results, PSY.get_name(load))
    finally
        CRC.@ignore_derivatives @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test 37 5th_order Ind. Motor" begin
    path = (joinpath(pwd(), "test-ind-mot_5th"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)
        time_span = (0.0, 20.0)
        perturbation_trip = BranchTrip(1.0, Line, "BUS1-BUS4-i_1")
        # Create Simulation with Constant Power
        sim_P = Simulation(MassMatrixModel, sys, path, time_span, perturbation_trip)

        # Motor parameters
        load = first(get_components(PSY.StandardLoad, sys))
        # Include the induction motor
        dynamic_injector = Ind_Motor(load)
        set_dynamic_injector!(load, dynamic_injector)
        sim = Simulation(MassMatrixModel, sys, path, time_span, perturbation_trip)

        # Test initial voltages between two systems are equivalent
        voltages_P = sim_P.x0_init[1:8]
        voltages_motor = sim.x0_init[1:8]
        @test LinearAlgebra.norm(voltages_P - voltages_motor) < 1e-4

        # Test Initial Condition
        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test37_x0_init
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        #Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test37_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, Rodas5(); dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # No comparison since model is not available in other tools
        power = get_activepower_series(results, PSY.get_name(load))
        rpower = get_reactivepower_series(results, PSY.get_name(load))
    finally
        CRC.@ignore_derivatives @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end
