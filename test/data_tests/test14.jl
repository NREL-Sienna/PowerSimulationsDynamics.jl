using PowerSystems
using NLsolve
const PSY = PowerSystems

############### Data Network ########################
include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
threebus_file_dir = joinpath(dirname(@__FILE__), "ThreeBusMulti.raw")
threebus_sys = System(PowerModelsData(threebus_file_dir), runchecks = false)

function inv_case78(static_device)
    return PSY.DynamicInverter(
        static_device,
        1.0, # ω_ref,
        converter_high_power(), #converter
        outer_control(), #outer control
        inner_control(), #inner control voltage source
        dc_source_lv(), #dc source
        pll(), #pll
        filt(), #filter
    )
end

function dyn_gen_multi_tg(generator)
    return PSY.DynamicGenerator(
        generator,
        1.0, # ω_ref,
        machine_classic(), #machine
        shaft_damping(), #shaft
        avr_none(), #avr
        tg_type2(), #tg
        pss_none(),
    )
end

# Add dynamic generators to the system (each gen is linked through a static one)
for g in get_components(Generator, threebus_sys)
    if get_number(get_bus(g)) == 101
        case_gen = inv_case78(g)
        add_component!(threebus_sys, case_gen)
    elseif get_number(get_bus(g)) == 102
        case_gen = dyn_gen_multi_tg(g)
        add_component!(threebus_sys, case_gen)
    end
end

#Create Ybus_Fault
fault_branches = deepcopy(collect(get_components(Branch, threebus_sys)))
for br in fault_branches
    if get_name(br) == "1"
        br.r = 4 * br.r
        br.x = 4 * br.x
    end
end
Ybus_fault = PSY.Ybus(fault_branches, get_components(Bus, threebus_sys))[:, :];
