"""
Validation PSSE/ST6B:
This case study defines a three bus system with an infinite bus, GENROU and a load.
The GENROU machine has connected an ST8C Excitation System.
The fault drop the line connecting the infinite bus and GENROU.
"""
##################################################
############### LOAD DATA ########################
##################################################

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/ST8C/ThreeBusMulti.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/ST8C/ThreeBus_ST8C.dyr")
#csv_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/ST8C/results_PSSe.csv")

@testset "Test 60 ST8C ResidualModel" begin
    path = joinpath(pwd(), "test-psse-ST8C")
    !isdir(path) && mkdir(path)
    try
        # Define system
        sys = System(raw_file, dyr_file)
        for l in get_components(PSY.StandardLoad, sys)
            transform_load_to_constant_impedance(l)
        end

        # Define Simulation Problem
        sim = Simulation(
            ResidualModel,
            sys, #system
            path,
            (0.0, 20.0), #time span
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test60_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test diff_val[1] < 1e-3

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test60_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, IDA(); dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain results

        t_psid, v2_psid = get_voltage_magnitude_series(results, 102)
        _, v3_psid = get_voltage_magnitude_series(results, 103)
        _, ω_psid = get_state_series(results, ("generator-102-1", :ω))
        _, Vf = get_field_voltage_series(results, "generator-102-1")

        #=
        # Obtain PSSE results
        M = get_csv_data(csv_file)
        t_psse = M[:, 1]
        v1_psse = M[:, 2]
        v2_psse = M[:, 3]
        v3_psse = M[:, 4]
        v4_psse = M[:, 5]
        ω_psse = M[:, 6] .+ 1.0
        efd_psse = M[:, 7]

        # Test Transient Simulation Results

        @test LinearAlgebra.norm(t_psid - round.(t_psse, digits = 3)) == 0.0
        @test LinearAlgebra.norm(v2_psid - v2_psse, Inf) <= 1e-3
        @test LinearAlgebra.norm(v3_psid - v3_psse, Inf) <= 1e-3
        @test LinearAlgebra.norm(ω_psid - ω_psse, Inf) <= 1e-3
        =#
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test 60 ST8C MassMatrixModel" begin
    path = joinpath(pwd(), "test-psse-ST6B")
    !isdir(path) && mkdir(path)
    try
        # Define system
        sys = System(raw_file, dyr_file)
        for l in get_components(PSY.StandardLoad, sys)
            transform_load_to_constant_impedance(l)
        end

        # Define Simulation Problem
        sim = Simulation(
            MassMatrixModel,
            sys, #system
            path,
            (0.0, 20.0), #time span
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test60_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test diff_val[1] < 1e-3

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test60_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, Rodas4(); dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain results

        t_psid, v2_psid = get_voltage_magnitude_series(results, 102)
        _, v3_psid = get_voltage_magnitude_series(results, 103)
        _, ω_psid = get_state_series(results, ("generator-102-1", :ω))
        _, Vf = get_field_voltage_series(results, "generator-102-1")

        #=
        # Obtain PSSE results
        M = get_csv_data(csv_file)
        t_psse = M[:, 1]
        v1_psse = M[:, 2]
        v2_psse = M[:, 3]
        v3_psse = M[:, 4]
        v4_psse = M[:, 5]
        ω_psse = M[:, 6] .+ 1.0
        efd_psse = M[:, 7]

        # Test Transient Simulation Results

        @test LinearAlgebra.norm(t_psid - round.(t_psse, digits = 3)) == 0.0
        @test LinearAlgebra.norm(v2_psid - v2_psse, Inf) <= 1e-3
        @test LinearAlgebra.norm(v3_psid - v3_psse, Inf) <= 1e-3
        @test LinearAlgebra.norm(ω_psid - ω_psse, Inf) <= 1e-3
        =#
    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end
