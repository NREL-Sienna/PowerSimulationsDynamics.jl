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

@testset "Test 5th_order Ind. Motor" begin
    path = (joinpath(pwd(), "test-ind-mot_5th"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)
        # correct X_thevenin
        source = first(get_components(Source, sys))
        PSY.set_X_th!(source, 0.01)
        # Motor parameters 
 
        Ind_Motor() = PSY.SingleCageInductionMachine(
            name="motor",
            Rs=0.013,
            Rr=0.009,
            Xs=0.067,
            Xr=0.17,
            Xm=3.8,
            H=1.5,
            A=1.0,
            B=0.0,
            base_power=1000
        )

        # Include the induction motor
        load = first(get_components(PSY.PowerLoad, sys))
        dynamic_injector = Ind_Motor()
        set_dynamic_injector!(load, dynamic_injector)

        time_span = (0.0, 20.0)
        perturbation_trip = BranchTrip(1.0, Line, "BUS1-BUS4-i_1")

        sim = Simulation!(ResidualModel, sys, pwd(), time_span, perturbation_trip)

        # Test Initial Condition
        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in init_cond
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
 
        @test (diff[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        #Test Eigenvalues
        @test LinearAlgebra.norm(eigs - eigs_value) < 1e-3

        # Solve problem
        @test execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for bus voltages
        V_1 = get_voltage_magnitude_series(results,1)
        V_2 = get_voltage_magnitude_series(results,2)
        V_3 = get_voltage_magnitude_series(results,3)
        V_4 = get_voltage_magnitude_series(results,4)

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
