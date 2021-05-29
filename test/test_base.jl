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
            ImplicitModel,
            omib_sys, #system
            path1,
            (0.0, 30.0), #time span
            Ybus_change;
            system_to_file = true,
        )

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
            ImplicitModel,
            omib_sys, #system
            path2,
            (0.0, 30.0), #time span
            Ybus_change;
            system_to_file = true,
        )

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
    sim = Simulation(ImplicitModel, threebus_sys_dyns, pwd(), (0.0, 10.0))
    sim_inputs = sim.simulation_inputs
    dae_vector = PSID.get_DAE_vector(sim_inputs)
    @test all(dae_vector)
    total_shunts = PSID.get_total_shunts(sim_inputs)
    for (k, v) in total_shunts
        @test v > 0
    end
    voltage_buses_ix = PSID.get_voltage_buses_ix(sim_inputs)
    @test length(voltage_buses_ix) == 3

    # Tests for dynamic lines with b = 0
    set_b!(dyn_branch12, (from = 0.0, to = 0.0))
    set_b!(dyn_branch23, (from = 0.0, to = 0.0))
    sim = Simulation(ImplicitModel, threebus_sys_dyns, pwd(), (0.0, 10.0))
    sim_inputs = sim.simulation_inputs
    dae_vector = PSID.get_DAE_vector(sim_inputs)
    @test sum(.!dae_vector) == 6
    total_shunts = PSID.get_total_shunts(sim_inputs)
    @test isempty(total_shunts)
    voltage_buses_ix = PSID.get_voltage_buses_ix(sim_inputs)
    @test isempty(voltage_buses_ix)
end
