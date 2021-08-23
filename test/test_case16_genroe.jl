"""
Validation PSSE/GENROE:
This case study defines a three bus system with an infinite bus, GENROE and a load.
The fault drop the line connecting the infinite bus and GENROE.
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

#Define dyr files

names = ["GENROE: Normal Saturation", "GENROE: High Saturation"]

dyr_files = [
    joinpath(dirname(@__FILE__), "benchmarks/psse/GENROE/ThreeBus_GENROE.dyr"),
    joinpath(dirname(@__FILE__), "benchmarks/psse/GENROE/ThreeBus_GENROE_HIGH_SAT.dyr"),
]

csv_files = (
    joinpath(dirname(@__FILE__), "benchmarks/psse/GENROE/TEST_GENROE.csv"),
    joinpath(dirname(@__FILE__), "benchmarks/psse/GENROE/TEST_GENROE_HIGH_SAT.csv"),
)

init_conditions = [test_psse_genroe_init, test_psse_genroe_high_sat_init]

eigs_values = [test16_eigvals, test16_eigvals_high_sat]

raw_file_dir = joinpath(dirname(@__FILE__), "benchmarks/psse/GENROE/ThreeBusMulti.raw")
tspan = (0.0, 20.0)

function test_genroe_implicit(dyr_file, csv_file, init_cond, eigs_value)
    path = (joinpath(pwd(), "test-psse-genrou"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file_dir, dyr_file)

        #Compute Ybus_fault
        sys2 = System(raw_file_dir)
        #Remove line connecting bus 1 to 2.
        remove_component!(Line, sys2, "BUS 1-BUS 2-i_1")
        Ybus_fault = Ybus(sys2).data

        #Define Fault: Change of YBus
        Ybus_change = NetworkSwitch(
            1.0, #change at t = 1.0
            Ybus_fault, #New YBus
        )

        #Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            sys, #system
            path,
            tspan, #time span
            Ybus_change,
        ) #Type of Fault

        #Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        #Solve problem
        execute!(sim, IDA(), dtmax = 0.005, saveat = 0.005)

        #Obtain data for angles
        series = get_state_series(sim, ("generator-102-1", :δ))
        t = series[1]
        δ = series[2]

        series2 = get_voltage_magnitude_series(sim, 102)

        t_psse, δ_psse = get_csv_delta(csv_file)

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in init_cond
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        #Test Initial Condition
        @test (diff[1] < 1e-3)
        #Test Eigenvalues
        @test LinearAlgebra.norm(eigs - eigs_value) < 1e-3
        #Test Solution DiffEq
        @test sim.solution.retcode == :Success

        #Test Transient Simulation Results
        # PSSE results are in Degrees
        @test LinearAlgebra.norm(δ - (δ_psse .* pi / 180), Inf) <= 1e-1
        @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

function test_genroe_mass_matrix(dyr_file, csv_file, init_cond, eigs_value)
    path = (joinpath(pwd(), "test-psse-genrou"))
    !isdir(path) && mkdir(path)
    try
        sys = System(raw_file_dir, dyr_file)

        #Compute Ybus_fault
        sys2 = System(raw_file_dir)
        #Remove line connecting bus 1 to 2.
        remove_component!(Line, sys2, "BUS 1-BUS 2-i_1")
        Ybus_fault = Ybus(sys2).data

        #Define Fault: Change of YBus
        Ybus_change = NetworkSwitch(
            1.0, #change at t = 1.0
            Ybus_fault, #New YBus
        )

        #Define Simulation Problem
        sim = Simulation!(
            MassMatrixModel,
            sys, #system
            path,
            tspan, #time span
            Ybus_change,
        ) #Type of Fault

        #Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        #Solve problem
        execute!(sim, Rodas5(), dtmax = 0.005, saveat = 0.005)

        #Obtain data for angles
        series = get_state_series(sim, ("generator-102-1", :δ))
        t = series[1]
        δ = series[2]

        series2 = get_voltage_magnitude_series(sim, 102)

        t_psse, δ_psse = get_csv_delta(csv_file)

        diff = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in init_cond
            diff[1] += LinearAlgebra.norm(res[k] - v)
        end
        #Test Initial Condition
        @test (diff[1] < 1e-3)
        #Test Eigenvalues
        @test LinearAlgebra.norm(eigs - eigs_value) < 1e-3
        #Test Solution DiffEq
        @test sim.solution.retcode == :Success
        #Test Transient Simulation Results
        # PSSE results are in Degrees
        @test LinearAlgebra.norm(δ - (δ_psse .* pi / 180), Inf) <= 1e-1
        @test LinearAlgebra.norm(t - round.(t_psse, digits = 3)) == 0.0

    finally
        @info("removing test files")
        rm(path, force = true, recursive = true)
    end
end

@testset "Test 16 GENROE ResidualModel" begin
    for (ix, name) in enumerate(names)
        @testset "$(name)" begin
            dyr_file = dyr_files[ix]
            csv_file = csv_files[ix]
            init_cond = init_conditions[ix]
            eigs_value = eigs_values[ix]
            test_genroe_implicit(dyr_file, csv_file, init_cond, eigs_value)
        end
    end
end

@testset "Test 16 GENROE MassMatrixModel" begin
    for (ix, name) in enumerate(names)
        @testset "$(name)" begin
            dyr_file = dyr_files[ix]
            csv_file = csv_files[ix]
            init_cond = init_conditions[ix]
            eigs_value = eigs_values[ix]
            test_genroe_mass_matrix(dyr_file, csv_file, init_cond, eigs_value)
        end
    end
end
