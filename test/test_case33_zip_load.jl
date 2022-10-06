"""
Validation PSSE/ZIP:
This case study defines a three bus system with an infinite bus, GENROU and a load.
The fault drop the line connecting the infinite bus and GENROU.
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

names = ["Constant Power", "Constant Current"]
load_models = [PSY.LoadModels.ConstantPower, PSY.LoadModels.ConstantCurrent]
raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/LOAD/ThreeBusMulti.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/LOAD/ThreeBus_GENROU.dyr")

initial_conditions = test33_zipload_x0_init
eigs_values = [test33_eigvals_constantP, test33_eigvals_constantI]
csv_files = (
    joinpath(TEST_FILES_DIR, "benchmarks/psse/LOAD/constantpower.csv"),
    joinpath(TEST_FILES_DIR, "benchmarks/psse/LOAD/constantcurrent.csv"),
)

tspan = (0.0, 20.0)

function test_zipload_implicit(csv_file, eigs_value, load_model)
    path = (joinpath(pwd(), "test-psse-zipload"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)
        for l in get_components(PSY.PowerLoad, sys)
            PSY.set_model!(l, load_model)
        end

        # Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            sys,
            path,
            tspan,
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in initial_conditions
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
        @test execute!(sim, IDA(), abstol = 1e-9, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for voltages
        series = get_voltage_magnitude_series(results, 102)
        t_psid = series[1]
        v2_psid = series[2]
        _, v3_psid = get_voltage_magnitude_series(results, 103)

        M = get_csv_data(csv_file)
        t_psse = M[:, 1]
        v2_psse = M[:, 2]
        v3_psse = M[:, 4]

        #TODO: Test for LoadPower
        p = get_activepower_series(results, "load1031")

        # Test Transient Simulation Results
        @test LinearAlgebra.norm(v2_psid - v2_psse, Inf) <= 5e-2
        @test LinearAlgebra.norm(v3_psid - v3_psse, Inf) <= 5e-2
        @test LinearAlgebra.norm(t_psid - round.(t_psse, digits = 3)) == 0.0
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

function test_zipload_mass_matrix(csv_file, eigs_value, load_model)
    path = (joinpath(pwd(), "test-psse-zipload"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)
        for l in get_components(PSY.PowerLoad, sys)
            PSY.set_model!(l, load_model)
        end

        # Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            sys,
            path,
            tspan,
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in initial_conditions
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
        @test execute!(sim, Rodas4(), abstol = 1e-9, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for voltages
        series = get_voltage_magnitude_series(results, 102)
        t_psid = series[1]
        v2_psid = series[2]
        _, v3_psid = get_voltage_magnitude_series(results, 103)

        M = get_csv_data(csv_file)
        t_psse = M[:, 1]
        v2_psse = M[:, 2]
        v3_psse = M[:, 4]

        # Test Transient Simulation Results
        @test LinearAlgebra.norm(v2_psid - v2_psse, Inf) <= 5e-2
        @test LinearAlgebra.norm(v3_psid - v3_psse, Inf) <= 5e-2
        @test LinearAlgebra.norm(t_psid - round.(t_psse, digits = 3)) == 0.0
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 33 ZIPLoad ResidualModel" begin
    for (ix, name) in enumerate(names)
        @testset "$(name)" begin
            csv_file = csv_files[ix]
            eigs_value = eigs_values[ix]
            load_model = load_models[ix]
            test_zipload_implicit(csv_file, eigs_value, load_model)
        end
    end
end

@testset "Test 33 ZIPLoad ResidualModel" begin
    for (ix, name) in enumerate(names)
        @testset "$(name)" begin
            csv_file = csv_files[ix]
            eigs_value = eigs_values[ix]
            load_model = load_models[ix]
            test_zipload_mass_matrix(csv_file, eigs_value, load_model)
        end
    end
end
