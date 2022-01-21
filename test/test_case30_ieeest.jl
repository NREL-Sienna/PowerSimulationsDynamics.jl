"""
Validation PSSE/IEEEST
This case study defines a three bus system with an infinite bus, GENROU+SEXS+IEEEST and a load.
The fault drop the line connecting the infinite bus and GENROU
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

names = ["IEEEST no Filter", "IEEEST with Filter"]

raw_file_dir = joinpath(TEST_FILES_DIR, "benchmarks/psse/IEEEST/ThreeBusMulti.raw")
dyr_files = [
    joinpath(TEST_FILES_DIR, "benchmarks/psse/IEEEST/ThreeBus_IEEEST.dyr")
    joinpath(TEST_FILES_DIR, "benchmarks/psse/IEEEST/ThreeBus_IEEEST_with_filter.dyr")
]

csv_files = [
    joinpath(TEST_FILES_DIR, "benchmarks/psse/IEEEST/IEEEST_SEXS_RESULTS_NOFILT.csv")
    joinpath(TEST_FILES_DIR, "benchmarks/psse/IEEEST/IEEEST_SEXS_RESULTS_FILT.csv")
]

init_conditions = [test_psse_ieeest_no_filt_init, test_psse_ieeest_with_filt_init]

eigs_values = [test30_eigvals_no_filt, test30_eigvals_with_filt]

function test_ieeest_implicit(dyr_file, csv_file, init_cond, eigs_value)
    path = (joinpath(pwd(), "test-psse-ieeest"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file_dir, dyr_file)

        # Define Simulation Problem
        sim = Simulation(
            ResidualModel,
            sys, #system
            path,
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault
        ) #Type of Fault

        # Test Initial Condition
        diffvals = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in init_cond
            diffvals[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diffvals[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        #Test Eigenvalues
        @test LinearAlgebra.norm(eigs - eigs_value) < 1e-3

        # Solve problem
        @test execute!(sim, IDA(), (0.0, 20.0), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for angles
        series = get_state_series(results, ("generator-102-1", :Vf))
        t = series[1]
        Efd = series[2]
        series2 = get_voltage_magnitude_series(results, 102)

        # Obtain PSS/E data
        M = get_csv_data(csv_file)
        t_psse = M[:, 1]
        Efd_psse = M[:, 2]
        V2_psse = M[:, 3]
        t_psse, Efd_psse = clean_extra_timestep!(t_psse, Efd_psse)
        t_psse, V2_psse = clean_extra_timestep!(t_psse, V2_psse)

        # Test Transient Simulation Results
        @test LinearAlgebra.norm(Efd - Efd_psse, Inf) <= 1e-2
        @test LinearAlgebra.norm(series2[2] - V2_psse, Inf) <= 1e-2
        @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0

        power = PSID.get_activepower_series(results, "generator-102-1")
        rpower = PSID.get_reactivepower_series(results, "generator-102-1")
        @test isa(series2, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
        sim = nothing
    end
    return
end

function test_ieeest_mass_matrix(dyr_file, csv_file, init_cond, eigs_value)
    path = (joinpath(pwd(), "test-psse-genrou"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file_dir, dyr_file)

        # Define Simulation Problem
        sim = Simulation(
            MassMatrixModel,
            sys, #system
            path,
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault
        ) #Type of Fault

        # Test Initial Condition
        diffvals = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in init_cond
            diffvals[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test (diffvals[1] < 1e-3)

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        #Test Eigenvalues
        @test LinearAlgebra.norm(eigs - eigs_value) < 1e-3

        # Solve problem
        @test execute!(sim, Rodas4(), (0.0, 20.0), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # Obtain data for angles
        series = get_state_series(results, ("generator-102-1", :Vf))
        t = series[1]
        Efd = series[2]
        series2 = get_voltage_magnitude_series(results, 102)

        # Obtain PSS/E data
        M = get_csv_data(csv_file)
        t_psse = M[:, 1]
        Efd_psse = M[:, 2]
        V2_psse = M[:, 3]
        t_psse, Efd_psse = clean_extra_timestep!(t_psse, Efd_psse)
        t_psse, V2_psse = clean_extra_timestep!(t_psse, V2_psse)

        # Test Transient Simulation Results
        @test LinearAlgebra.norm(Efd - Efd_psse, Inf) <= 1e-2
        @test LinearAlgebra.norm(series2[2] - V2_psse, Inf) <= 1e-2
        @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0

        power = PSID.get_activepower_series(results, "generator-102-1")
        rpower = PSID.get_reactivepower_series(results, "generator-102-1")
        @test isa(series2, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(power, Tuple{Vector{Float64}, Vector{Float64}})
        @test isa(rpower, Tuple{Vector{Float64}, Vector{Float64}})
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
        sim = nothing
    end
    return
end

@testset "Test 30 IEEEST ResidualModel" begin
    for (ix, name) in enumerate(names)
        @testset "$(name)" begin
            dyr_file = dyr_files[ix]
            csv_file = csv_files[ix]
            init_cond = init_conditions[ix]
            eigs_value = eigs_values[ix]
            test_ieeest_implicit(dyr_file, csv_file, init_cond, eigs_value)
        end
    end
end

@testset "Test 30 IEEEST MassMatrixModel" begin
    for (ix, name) in enumerate(names)
        @testset "$(name)" begin
            dyr_file = dyr_files[ix]
            csv_file = csv_files[ix]
            init_cond = init_conditions[ix]
            eigs_value = eigs_values[ix]
            test_ieeest_mass_matrix(dyr_file, csv_file, init_cond, eigs_value)
        end
    end
end
