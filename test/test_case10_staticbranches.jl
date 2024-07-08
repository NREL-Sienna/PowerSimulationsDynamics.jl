"""
Case 10:
This case study a three bus system with 1 machine (One d- One q-: 4th order model), a VSM of 19 states and an infinite source. All lines are modeled as a static lines.
The perturbation trips two of the three circuits of line between buses 1 and 2, triplicating its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/test10.jl"))

##################################################
############### SOLVE PROBLEM ####################
##################################################

# time span
tspan = (0.0, 40.0)

# Define Fault: Change of YBus
Ybus_change = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

@testset "Test 10 Static Branches ResidualModel" begin
    path = mktempdir()
    try
        # Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            threebus_sys, #system,
            path,
            tspan, #time span
            Ybus_change, #Type of Fault
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test10_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        try
            small_sig = small_signal_analysis(sim)
        catch e
            T = ResidualModel
            inputs = PSID.get_simulation_inputs(sim)
            x_eval = PSID.get_initial_conditions(sim)
            jacwrapper = PSID.get_jacobian(T, inputs, x_eval, 0)
            jacobian = jacwrapper.Jv
            println("jacobian:")
            for row in eachrow(jacobian)
                println(row)
            end
            mass_matrix = PSID.get_mass_matrix(inputs)
            println("mass matrix:")
            for row in eachrow(mass_matrix)
                println(row)
            end
            diff_states = PSID.get_DAE_vector(inputs)
            println("diff states:")
            @error(diff_states)
            global_index = PSID.make_global_state_map(inputs)
            jac_index = PSID._make_reduced_jacobian_index(global_index, diff_states)
            println("jac index:")
            @error(jac_index)
            reduced_jacobian =
                PSID._reduce_jacobian(jacobian, diff_states, mass_matrix, global_index)
            println("reduced jacobian:")
            for row in eachrow(reduced_jacobian)
                println(row)
            end
            eigen_vals, R_eigen_vect = LinearAlgebra.eigen(Matrix(reduced_jacobian))
            println("eigvals:")
            for eig in eigen_vals
                println(eig)
            end
            println("eigvecs")
            for row in eachrow(R_eigen_vect)
                println(row)
            end
        end
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test10_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, IDA()) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for voltages
        series = get_voltage_magnitude_series(results, 102)

        # Create Zoom for comparison with dynamic lines
        zoom = [
            (series[1][ix], series[2][ix]) for
            (ix, s) in enumerate(series[1]) if (s > 0.90 && s < 1.6)
        ]
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test 10 Static Branches MassMatrixModel" begin
    path = mktempdir()
    try
        # Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            threebus_sys, #system,
            path,
            tspan, #time span
            Ybus_change, #Type of Fault
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test10_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test10_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, Rodas4()) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for voltages
        series = get_voltage_magnitude_series(results, 102)

        # Create Zoom for comparison with dynamic lines
        zoom = [
            (series[1][ix], series[2][ix]) for
            (ix, s) in enumerate(series[1]) if (s > 0.90 && s < 1.6)
        ]
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end
