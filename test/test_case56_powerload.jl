"""
Validation Power Load
This case study defines a three bus system with an infinite bus, GENROU and a load.
The fault drop the line connecting the infinite bus and GENROU. The test validates
that the PowerLoad Model matches the StandardLoad model with only constant power component.
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/LOAD/ThreeBusMulti.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/LOAD/ThreeBus_GENROU.dyr")

# Create PowerLoad system
sys_power = System(raw_file, dyr_file)

# Create StandardLoad system
sys_standard = System(raw_file, dyr_file)

# Replace StandardLoad with PowerLoad
for l in collect(get_components(PSY.StandardLoad, sys_power))
    power_load = PSY.PowerLoad(;
        name = PSY.get_name(l),
        available = PSY.get_available(l),
        bus = PSY.get_bus(l),
        active_power = PSY.get_constant_active_power(l),
        reactive_power = PSY.get_constant_reactive_power(l),
        #α = 0.0, # Constant Power
        #β = 0.0, # Constant Power
        base_power = PSY.get_base_power(l),
        max_active_power = PSY.get_max_constant_active_power(l),
        max_reactive_power = PSY.get_max_constant_reactive_power(l),
    )
    PSY.remove_component!(sys_power, l)
    PSY.add_component!(sys_power, power_load)
end

@testset "Test 56 PowerLoad ResidualModel" begin
    path = (joinpath(pwd(), "test-56"))
    !isdir(path) && mkdir(path)
    try
        # Instantiate Simulations
        sim_power = Simulation(
            ResidualModel,
            sys_power,
            path,
            tspan,
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault
        )
        sim_standard = Simulation(
            ResidualModel,
            sys_standard,
            path,
            tspan,
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault
        )

        # Test Initial Conditions
        @test LinearAlgebra.norm(sim_power.x0_init - sim_standard.x0_init) < 1e-4

        # Test Small Signal
        ss_power = small_signal_analysis(sim_power)
        @test ss_power.stable
        ss_standard = small_signal_analysis(sim_standard)
        @test ss_standard.stable
        # Compare Eigenvalues
        @test LinearAlgebra.norm(ss_power.eigenvalues - ss_standard.eigenvalues) < 1e-4

        # Solve Problems
        @test execute!(sim_power, IDA(); abstol = 1e-9, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results_power = read_results(sim_power)
        @test execute!(sim_standard, IDA(); abstol = 1e-9, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results_standard = read_results(sim_standard)

        # Store results
        _, v102_power = get_voltage_magnitude_series(results_power, 102)
        _, v102_standard = get_voltage_magnitude_series(results_standard, 102)
        _, v103_power = get_voltage_magnitude_series(results_power, 103)
        _, v103_standard = get_voltage_magnitude_series(results_standard, 103)

        _, p_standard = get_activepower_series(results_standard, "load1031")
        _, p_power = get_activepower_series(results_power, "load1031")

        # Test Transient Simulation Results
        @test LinearAlgebra.norm(v102_power - v102_standard, Inf) <= 1e-3
        @test LinearAlgebra.norm(v103_power - v103_standard, Inf) <= 1e-3
        @test LinearAlgebra.norm(p_power - p_standard, Inf) <= 1e-3

    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test 56 PowerLoad MassMatrixModel" begin
    path = (joinpath(pwd(), "test-56"))
    !isdir(path) && mkdir(path)
    try
        # Instantiate Simulations
        sim_power = Simulation(
            MassMatrixModel,
            sys_power,
            path,
            tspan,
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault
        )
        sim_standard = Simulation(
            MassMatrixModel,
            sys_standard,
            path,
            tspan,
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault
        )

        # Test Initial Conditions
        @test LinearAlgebra.norm(sim_power.x0_init - sim_standard.x0_init) < 1e-4

        # Test Small Signal
        ss_power = small_signal_analysis(sim_power)
        @test ss_power.stable
        ss_standard = small_signal_analysis(sim_standard)
        @test ss_standard.stable
        # Compare Eigenvalues
        @test LinearAlgebra.norm(ss_power.eigenvalues - ss_standard.eigenvalues) < 1e-4

        # Solve Problems
        @test execute!(sim_power, Rodas5P(); abstol = 1e-9, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results_power = read_results(sim_power)
        @test execute!(sim_standard, Rodas5P(); abstol = 1e-9, saveat = 0.005) ==
              PSID.SIMULATION_FINALIZED
        results_standard = read_results(sim_standard)

        # Store results
        _, v102_power = get_voltage_magnitude_series(results_power, 102)
        _, v102_standard = get_voltage_magnitude_series(results_standard, 102)
        _, v103_power = get_voltage_magnitude_series(results_power, 103)
        _, v103_standard = get_voltage_magnitude_series(results_standard, 103)

        _, p_standard = get_activepower_series(results_standard, "load1031")
        _, p_power = get_activepower_series(results_power, "load1031")

        # Test Transient Simulation Results
        @test LinearAlgebra.norm(v102_power - v102_standard, Inf) <= 1e-3
        @test LinearAlgebra.norm(v103_power - v103_standard, Inf) <= 1e-3
        @test LinearAlgebra.norm(p_power - p_standard, Inf) <= 1e-3

    finally
        @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end
