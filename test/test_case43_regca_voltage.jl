"""
Test for REGCA Voltage model.
This case study defines a three bus system with an infinite bus in 1,
and a REGCA with a voltage model in bus 3.
The disturbance trips the line connecting buses 1 and 2.
There is no equivalent model for this one in PSS/e.
Note: The model will not work with a capacitor value different than zero.
"""
##################################################
############### LOAD DATA ########################
##################################################

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/RENA/ThreeBusRenewable.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/RENA/ThreeBus_REN_A_DEFAULT_FLAG.dyr")

sys = System(raw_file, dyr_file)
for l in get_components(PSY.StandardLoad, sys)
    transform_load_to_constant_impedance(l)
end

# Modify system to use a REGCA voltage model
inverter = get_component(ThermalStandard, sys, "generator-103-1")
dynamic_injector = get_dynamic_injector(inverter)

filt_rl = deepcopy(get_filter(dynamic_injector))
converter = deepcopy(get_converter(dynamic_injector))
outer = deepcopy(get_outer_control(dynamic_injector))
inner = deepcopy(get_inner_control(dynamic_injector))
dc_sorc = deepcopy(get_dc_source(dynamic_injector))
freq_est = deepcopy(get_freq_estimator(dynamic_injector))

conv_voltage() = PSY.RenewableEnergyVoltageConverterTypeA(
    0.02,
    10.0,
    0.9,
    0.4,
    1.22,
    1.2,
    (min = 0.5, max = 0.9),
    -1.3,
    0.2,
    0.0,
    (min = -100.0, max = 100.0),
    0.0,
    0,
)

filt_noc() = LCLFilter(; lf = 0.08, rf = 0.003, cf = 0.0, lg = 0.2, rg = 0.01)

new_dynamic_injector = PSY.DynamicInverter(;
    name = PSY.get_name(inverter),
    ω_ref = 1.0,
    converter = conv_voltage(),
    outer_control = outer,
    inner_control = inner,
    dc_source = dc_sorc,
    freq_estimator = freq_est,
    filter = filt_noc(),
    base_power = 100.0,
)

PSY.set_dynamic_injector!(inverter, nothing)
PSY.set_dynamic_injector!(inverter, new_dynamic_injector)

##################################################
############### SOLVE PROBLEM ####################
##################################################

@testset "Test 43 REGCA Voltage ResidualModel" begin
    path = joinpath(pwd(), "test-regca-voltage")
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation(
            ResidualModel,
            sys, #system
            path,
            (0.0, 2.0), #time span
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault,
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test43_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test diff_val[1] < 1e-3

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test43_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, IDA(); abstol = 1e-9) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # No comparison
    finally
        CRC.@ignore_derivatives @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end

@testset "Test 43 REGCA Voltage MassMatrixModel" begin
    path = joinpath(pwd(), "test-regca-voltage")
    !isdir(path) && mkdir(path)
    try
        # Define Simulation Problem
        sim = Simulation(
            MassMatrixModel,
            sys, #system
            path,
            (0.0, 2.0), #time span
            BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1"), #Type of Fault,
        )

        # Test Initial Condition
        diff_val = [0.0]
        res = get_init_values_for_comparison(sim)
        for (k, v) in test43_x0_init
            diff_val[1] += LinearAlgebra.norm(res[k] - v)
        end

        @test diff_val[1] < 1e-3

        # Obtain small signal results for initial conditions
        small_sig = small_signal_analysis(sim)
        eigs = small_sig.eigenvalues
        @test small_sig.stable

        # Test Eigenvalues
        @test LinearAlgebra.norm(eigs - test43_eigvals) < 1e-3

        # Solve problem
        @test execute!(sim, Rodas4(); abstol = 1e-9) == PSID.SIMULATION_FINALIZED
        results = read_results(sim)

        # No comparison
    finally
        CRC.@ignore_derivatives @info("removing test files")
        rm(path; force = true, recursive = true)
    end
end
