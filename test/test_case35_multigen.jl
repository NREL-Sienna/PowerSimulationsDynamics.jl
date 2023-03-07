"""
Validation PSSE/MultiGen:
This case study defines a four bus system multimachine with several units in one bus.
The fault trips a generator
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/MultiGen/FourBusMulti.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/MultiGen/FourBus_multigen.dyr")
gen_trip_csv_file =
    joinpath(TEST_FILES_DIR, "benchmarks/psse/MultiGen/gen_trip_results.csv")
line_trip_csv_file =
    joinpath(TEST_FILES_DIR, "benchmarks/psse/MultiGen/line_trip_results.csv")

@testset "Test 35 MultiGen ResidualModel" begin
    path = (joinpath(pwd(), "test-multigen"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)
        for l in get_components(PSY.StandardLoad, sys)
            PSID.transform_load_to_constant_impedance(l)
        end

        dc = get_dynamic_injector(get_component(ThermalStandard, sys, "generator-102-SW"))
        # Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            sys, #system
            path,
            (0.0, 20.0), #time span
            GeneratorTrip(1.0, dc),
        )

        # Solve problem
        @test execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # # Obtain results
        M = get_csv_data(gen_trip_csv_file)
        t_psse = M[:, 1]
        for (ix, b) in enumerate(101:104)
            t, Vt = get_voltage_magnitude_series(results, b)
            Vt_psse = M[:, 47 + 2 * (ix - 1)]
            @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0
            @test LinearAlgebra.norm(Vt - Vt_psse, Inf) <= 1e-3
        end

        for (ix, g) in enumerate(get_components(ThermalStandard, sys))
            if isa(get_dynamic_injector(g), PSY.DynamicInverter)
                continue
            end
            gen_name = get_name(g)
            t, ω = get_state_series(results, (gen_name, :ω))
            ω_psse = M[:, 5 + 5 * (ix - 1)] .+ 1.0
            @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0
            @test LinearAlgebra.norm(ω - ω_psse, Inf) <= 0.1 # relaxed due to inconsistent PSSE outputs
        end

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 35 MultiGen MassMatrixModel" begin
    path = (joinpath(pwd(), "test-psse-gast"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file, dyr_file)
        for l in get_components(PSY.StandardLoad, sys)
            PSID.transform_load_to_constant_impedance(l)
        end

        dc = get_dynamic_injector(get_component(ThermalStandard, sys, "generator-102-SW"))
        # Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            sys, #system
            path,
            (0.0, 20.0), #time span
            GeneratorTrip(1.0, dc),
        )

        # Solve problem
        @test execute!(sim, Rodas4(), dtmax = 0.005, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        M = get_csv_data(gen_trip_csv_file)
        t_psse = M[:, 1]
        for (ix, b) in enumerate(101:104)
            t, Vt = get_voltage_magnitude_series(results, b)
            Vt_psse = M[:, 47 + 2 * (ix - 1)]
            @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0
            @test LinearAlgebra.norm(Vt - Vt_psse, Inf) <= 1e-3
        end

        for (ix, g) in enumerate(get_components(ThermalStandard, sys))
            if isa(get_dynamic_injector(g), PSY.DynamicInverter)
                continue
            end
            gen_name = get_name(g)
            t, ω = get_state_series(results, (gen_name, :ω))
            ω_psse = M[:, 5 + 5 * (ix - 1)] .+ 1.0
            @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0
            @test LinearAlgebra.norm(ω - ω_psse, Inf) <= 0.1 # relaxed due to inconsistent PSSE outputs
        end
    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end
