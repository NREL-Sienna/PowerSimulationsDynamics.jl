"""
Case 42:
This case study a three bus system with one AggregateDistributedGenerationA model, one load, and one infinite source. 
The fault drops the line connecting the infinite bus and AggregateDistributedGenerationA.

"""
##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(TEST_FILES_DIR, "data_tests/data_utils.jl"))
raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/DERA/ThreeBusMulti.raw")

##################################################
############### SOLVE PROBLEM ####################
##################################################

names = ["DERA: FreqFlag=0", "DERA: FreqFlag=1"]

#TODO - include set of dyr values once parser includes DERA. 
FreqFlag_values = [0, 1]
csv_files = (
    joinpath(TEST_FILES_DIR, "benchmarks/psse/DERA/dera_freqflag0.csv"),
    joinpath(TEST_FILES_DIR, "benchmarks/psse/DERA/dera_freqflag1.csv"),
)
init_conditions = [test_psse_dera_freqflag0_init, test_psse_dera_freqflag1_init]
eigs_values = [test42_eigvals_freqflag0, test42_eigvals_freqflag1]

tspan = (0.0, 4.0)
csv_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/DERA/TEST_DERA.csv")

function test_dera_residual(freqflag_value, csv_file, init_cond, eigs_value)
    path = (joinpath(pwd(), "test-psse-dera"))
    !isdir(path) && mkdir(path)
    try
        threebus_sys = System(raw_file, runchecks = false)
        for g in get_components(ThermalStandard, threebus_sys)
            g.bus.bustype == BusTypes.REF && remove_component!(threebus_sys, g)
        end
        add_source_to_ref(threebus_sys)
        for g in get_components(ThermalStandard, threebus_sys)
            case_dera = dera(g, freqflag_value)
            add_component!(threebus_sys, case_dera, g)
        end
        for l in get_components(PSY.StandardLoad, threebus_sys)
            PSID.transform_load_to_constant_impedance(l)
        end

        sim = Simulation(
            ResidualModel,
            threebus_sys,
            path,
            tspan,
            BranchTrip(2.0, Line, "BUS 1-BUS 2-i_1"),
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in init_cond
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff_val[1] < 1e-3)

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

        # Obtain data for angles
        t_psid = get_voltage_angle_series(results, 102)[1]
        θ_psid = get_voltage_angle_series(results, 102)[2]
        V_psid = get_voltage_magnitude_series(results, 102)[2]
        power = get_activepower_series(results, "generator-102-1")

        M = get_csv_data(csv_file)
        t_psse, V_psse = clean_extra_timestep!(M[:, 1], M[:, 2])
        _, θ_psse = clean_extra_timestep!(M[:, 1], M[:, 3])

        @test LinearAlgebra.norm(θ_psid - (θ_psse .* pi / 180), Inf) <= 1e-1
        @test LinearAlgebra.norm(V_psid - V_psse, Inf) <= 1e-1
        @test LinearAlgebra.norm(t_psid - round.(t_psse, digits = 3)) == 0.0
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

function test_dera_massmatrix(freqflag_value, csv_file, init_cond, eigs_value)
    path = (joinpath(pwd(), "test-psse-dera"))
    !isdir(path) && mkdir(path)
    try
        threebus_sys = System(raw_file, runchecks = false)
        for g in get_components(ThermalStandard, threebus_sys)
            g.bus.bustype == BusTypes.REF && remove_component!(threebus_sys, g)
        end
        add_source_to_ref(threebus_sys)
        for g in get_components(ThermalStandard, threebus_sys)
            case_dera = dera(g, freqflag_value)
            add_component!(threebus_sys, case_dera, g)
        end
        for l in get_components(PSY.StandardLoad, threebus_sys)
            PSID.transform_load_to_constant_impedance(l)
        end

        sim = Simulation(
            MassMatrixModel,
            threebus_sys,
            path,
            tspan,
            BranchTrip(2.0, Line, "BUS 1-BUS 2-i_1"),
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in init_cond
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end
        @test (diff_val[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        #Test Eigenvalues
        @test LinearAlgebra.norm(eigs - eigs_value) < 1e-3

        # Solve problem
        @test execute!(sim, Rodas4(), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for angles
        t_psid = get_voltage_angle_series(results, 102)[1]
        θ_psid = get_voltage_angle_series(results, 102)[2]
        V_psid = get_voltage_magnitude_series(results, 102)[2]
        power = get_activepower_series(results, "generator-102-1")

        M = get_csv_data(csv_file)
        t_psse, V_psse = clean_extra_timestep!(M[:, 1], M[:, 2])
        _, θ_psse = clean_extra_timestep!(M[:, 1], M[:, 3])

        @test LinearAlgebra.norm(θ_psid - (θ_psse .* pi / 180), Inf) <= 1e-1
        @test LinearAlgebra.norm(V_psid - V_psse, Inf) <= 1e-1
        @test LinearAlgebra.norm(t_psid - round.(t_psse, digits = 3)) == 0.0

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 42 DERA ResidualModel" begin
    for (ix, name) in enumerate(names)
        @testset "$(name)" begin
            freqflag_value = FreqFlag_values[ix]
            csv_file = csv_files[ix]
            init_cond = init_conditions[ix]
            eigs_value = eigs_values[ix]
            test_dera_residual(freqflag_value, csv_file, init_cond, eigs_value)
        end
    end
end

@testset "Test 42 DERA MassMatrixModel" begin
    for (ix, name) in enumerate(names)
        @testset "$(name)" begin
            freqflag_value = FreqFlag_values[ix]
            csv_file = csv_files[ix]
            init_cond = init_conditions[ix]
            eigs_value = eigs_values[ix]
            test_dera_massmatrix(freqflag_value, csv_file, init_cond, eigs_value)
        end
    end
end
