"""
Case 14:
This case study a three bus system with 1 machines (Classic Model - Single Shaft: 2 State model) and 1 Inverter without loads.
The inverter at bus 1 is used as a reference device, while machine at bus 2 has a simplified droop governor (TGTypeII).
The perturbation trips four (out of 5) circuits of line between buses 1 and 2, multiplying by 4 its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test14.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

#Time span
tspan = (0.0, 10.0)

#Define Fault: Change of YBus
Ybus_change = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

@testset "Test 14 Inverter Ref ResidualModel" begin
    path = (joinpath(pwd(), "test-14"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation(
            ResidualModel,
            threebus_sys, #system,
            path,
            tspan, #time span
            Ybus_change, #Type of Fault
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test14_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @show eigs
        @test LinearAlgebra.norm(eigs - test14_eigvals) < 1e-3

        #Run simulation
        @test execute!(
            sim, #simulation structure
            IDA(),#Sundials DAE Solver
            dtmax = 0.001, #keywords arguments
        ) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        series = get_state_series(results, ("generator-101-1", :ω_oc))
        series2 = get_mechanical_torque_series(results, "generator-101-1")
        series3 = get_mechanical_torque_series(results, "generator-102-1")
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 14 Inverter Ref MassMatrixModel" begin
    path = (joinpath(pwd(), "test-14"))
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation(
            MassMatrixModel,
            threebus_sys, #system,
            path,
            tspan, #time span
            Ybus_change, #Type of Fault
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test14_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test14_eigvals) < 1e-3

        #Run simulation
        @test execute!(
            sim, #simulation structure
            Rodas4(),#Sundials DAE Solver
            dtmax = 0.001, #keywords arguments
        ) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        series = get_state_series(results, ("generator-101-1", :ω_oc))
        series2 = get_mechanical_torque_series(results, "generator-101-1")
        series3 = get_mechanical_torque_series(results, "generator-102-1")
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
