##################################################
############### LOAD DATA ########################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/test01.jl"))
omib_sys_file = System(PowerModelsData(omib_file_dir), runchecks = false)

##################################################
############### SOLVE PROBLEM ####################
##################################################
#Define Fault: Change of YBus
Ybus_change = ThreePhaseFault(
    1.0, #change at t = 1.0
    Ybus_fault,
) #New YBus

path1 = (joinpath(pwd(), "test-Base-1"))
!isdir(path1) && mkdir(path1)
try
    #Define Simulation Problem
    sim = Simulation(
        path1,
        omib_sys, #system
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
        path2,
        omib_sys, #system
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
