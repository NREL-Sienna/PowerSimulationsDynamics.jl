##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test01.jl"))
omib_sys_file = System(PowerModelsData(omib_file_dir), runchecks = false)

##################################################
############### SOLVE PROBLEM ####################
##################################################
#Define Fault: Change of YBus
Ybus_change = NetworkSwitch(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

@testset "Make Simulation" begin
    path1 = (joinpath(pwd(), "test-Base-1"))
    !isdir(path1) && mkdir(path1)
    try
        #Define Simulation Problem
        sim = Simulation(
            ResidualModel,
            omib_sys, #system
            path1,
            (0.0, 30.0), #time span
            Ybus_change;
            system_to_file = true,
        )
        @test sim.status == PSID.BUILT
        o_system = System(joinpath(path1, "input_system.json"))
        for b in get_components(Bus, o_system)
            b_sys = get_component(Bus, omib_sys, get_name(b))
            b_file = get_component(Bus, omib_sys_file, get_name(b))
            @test get_angle(b) == get_angle(b_sys)
            @test get_angle(b) == get_angle(b_file)
        end

    finally
        @info("removing test files")
        rm(path1, force = true, recursive = true)
    end

    path2 = (joinpath(pwd(), "test-Base-2"))
    !isdir(path2) && mkdir(path2)
    try
        #Define Simulation Problem
        sim = Simulation!(
            ResidualModel,
            omib_sys, #system
            path2,
            (0.0, 30.0), #time span
            Ybus_change;
            system_to_file = true,
        )
        @test sim.status == PSID.BUILT
        m_system = System(joinpath(path2, "initialized_system.json"))
        for b in get_components(Bus, m_system)
            b_sys = get_component(Bus, omib_sys, get_name(b))
            b_file = get_component(Bus, omib_sys_file, get_name(b))
            @test get_angle(b) == get_angle(b_sys)
            if get_bustype(b) == PSY.BusTypes.REF
                @test get_angle(b) == get_angle(b_file)
            else
                @test get_angle(b) != get_angle(b_file)
            end
        end

    finally
        @info("removing test files")
        rm(path2, force = true, recursive = true)
    end
end

@testset "Hybrid Line Indexing" begin
    ## Create threebus system with more dyn lines ##
    three_bus_file_dir = joinpath(dirname(@__FILE__), "data_tests/ThreeBusInverter.raw")
    threebus_sys_dyns = System(three_bus_file_dir, runchecks = false)
    add_source_to_ref(threebus_sys_dyns)

    dyn_branch12 =
        DynamicBranch(get_component(Branch, threebus_sys_dyns, "BUS 1-BUS 2-i_2"))
    dyn_branch23 =
        DynamicBranch(get_component(Branch, threebus_sys_dyns, "BUS 2-BUS 3-i_3"))

    add_component!(threebus_sys_dyns, dyn_branch12)
    add_component!(threebus_sys_dyns, dyn_branch23)

    # Attach dyn devices
    for g in get_components(Generator, threebus_sys_dyns)
        if get_number(get_bus(g)) == 102
            case_gen = dyn_gen_second_order(g)
            add_component!(threebus_sys_dyns, case_gen, g)
        elseif get_number(get_bus(g)) == 103
            case_inv = inv_case78(g)
            add_component!(threebus_sys_dyns, case_inv, g)
        end
    end

    # Tests for all Dynamic Lines
    sim = Simulation(ResidualModel, threebus_sys_dyns, pwd(), (0.0, 10.0))
    @test sim.status == PSID.BUILT
    sim_inputs = sim.simulation_inputs
    DAE_vector = PSID.get_DAE_vector(sim_inputs)
    @test all(DAE_vector)
    total_shunts = PSID.get_total_shunts(sim_inputs)
    for v in LinearAlgebra.diag(total_shunts)
        @test imag(v) > 0
    end

    for entry in LinearAlgebra.diag(sim_inputs.mass_matrix)
        @test entry > 0
    end
    voltage_buses_ix = PSID.get_voltage_buses_ix(sim_inputs)
    @test length(voltage_buses_ix) == 3

    # Tests for dynamic lines with b = 0
    set_b!(dyn_branch12, (from = 0.0, to = 0.0))
    set_b!(dyn_branch23, (from = 0.0, to = 0.0))
    sim = Simulation(ResidualModel, threebus_sys_dyns, pwd(), (0.0, 10.0))
    @test sim.status == PSID.BUILT
    sim_inputs = sim.simulation_inputs
    DAE_vector = PSID.get_DAE_vector(sim_inputs)
    @test sum(.!DAE_vector) == 6
    for (ix, entry) in enumerate(DAE_vector)
        if !entry
            @test LinearAlgebra.diag(sim_inputs.mass_matrix)[ix] == 0
        elseif entry
            @test LinearAlgebra.diag(sim_inputs.mass_matrix)[ix] == 1
        else
            @test false
        end
    end

    total_shunts = PSID.get_total_shunts(sim_inputs)
    @test sum(LinearAlgebra.diag(total_shunts)) == 0.0
    voltage_buses_ix = PSID.get_voltage_buses_ix(sim_inputs)
    @test isempty(voltage_buses_ix)
end

@testset "Test initial conditions externally" begin
    #Create system with 2 GENROU generators
    threebus_file_dir =
        joinpath(dirname(@__FILE__), "benchmarks/psse/GENROU/ThreeBusMulti_LessLoad.raw")
    dyr_file_dir =
        joinpath(dirname(@__FILE__), "benchmarks/psse/GENROU/ThreeBus_GENROU_SEXS.dyr")
    sys = System(threebus_file_dir, dyr_file_dir)
    x0_test = zeros(23)
    x0_test[1:6] .= 1.0
    #Initialize System from given initial condition
    sim = Simulation(
        ImplicitModel,
        sys,
        pwd(),
        (0.0, 20.0),
        BranchTrip(1.0, "BUS 1-BUS 2-i_1");
        initialize_simulation = false,
        initial_conditions = x0_test,
    )
    @test LinearAlgebra.norm(sim.x0_init - x0_test) <= 1e-6
    #Initialize System normally
    sim_normal = Simulation(
        ImplicitModel,
        sys,
        pwd(),
        (0.0, 20.0),
        BranchTrip(1.0, "BUS 1-BUS 2-i_1"),
    )
    #Save states without generator at bus 2
    x0 = sim_normal.x0_init
    x0_no_gen = vcat(x0[1:6], x0[13:end])
    #Make Generator 2 unavailable and transform bus into PQ bus
    gen = PSY.get_component(ThermalStandard, sys, "generator-102-1")
    PSY.set_available!(gen, false)
    b = PSY.get_component(Bus, sys, "BUS 2")
    PSY.set_bustype!(b, PSY.BusTypes.PQ)
    #Create Simulation without Gen 2 starting from steady-state with Gen 2
    sim_trip_gen = Simulation(
        ImplicitModel,
        sys,
        pwd(),
        (0.0, 20.0);
        initialize_simulation = false,
        initial_conditions = x0_no_gen,
    )
    @test LinearAlgebra.norm(sim_trip_gen.x0_init - x0_no_gen) <= 1e-6
    #Create Simulation without Gen 2 at steady state
    sim_normal_no_gen = Simulation(
        ImplicitModel,
        sys,
        pwd(),
        (0.0, 20.0),
        BranchTrip(1.0, "BUS 1-BUS 2-i_1"),
    )
    @test length(sim_normal_no_gen.x0_init) == 17
    #Ignore Initial Conditions without passing initialize_simulation = false
    sim_ignore_init =
        Simulation(ImplicitModel, sys, pwd(), (0.0, 20.0); initial_conditions = x0_no_gen)
    @test LinearAlgebra.norm(sim_ignore_init.x0_init - sim_normal_no_gen.x0_init) <= 1e-6
    #Pass wrong vector size
    x0_wrong = zeros(20)
    @test_throws IS.ConflictingInputsError Simulation(
        ImplicitModel,
        sys,
        pwd(),
        (0.0, 20.0);
        initialize_simulation = false,
        initial_conditions = x0_wrong,
    )
    #Flat start initialization
    sim_flat =
        Simulation(ImplicitModel, sys, pwd(), (0.0, 20.0); initialize_simulation = false)
    x0_flat = zeros(17)
    x0_flat[1:3] .= 1.0
    @test LinearAlgebra.norm(sim_flat.x0_init - x0_flat) <= 1e-6
end

@testset "Test Network Kirchoff Calculation" begin
    voltages = [
        1.05
        1.019861574381897
        0.0
        -0.016803841801155062
    ]
    V_r = voltages[1:2]
    V_i = voltages[3:end]
    ybus_ = PSY.Ybus(omib_sys).data
    I_balance_ybus = -1 * ybus_ * (V_r + V_i .* 1im)
    inputs = PSID.SimulationInputs(ResidualModel, omib_sys)
    I_balance_sim = zeros(4)
    PSID.network_model(inputs, I_balance_sim, voltages)
    for i in 1:2
        @test isapprox(real(I_balance_ybus[i]), I_balance_sim[i])
        @test isapprox(imag(I_balance_ybus[i]), I_balance_sim[i + 2])
    end
end
