using PowerSystems
using NLsolve
const PSY = PowerSystems

############### Data Network ########################
include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
threebus_file_dir = joinpath(dirname(@__FILE__), "ThreeBusNetwork.raw")
threebus_sys = System(
    PowerModelsData(threebus_file_dir),
    runchecks = false,
    unit_system = "device_base",
)
add_source_to_ref(threebus_sys)

function dyn_gen_avr_type2(generator)
    return PSY.DynamicGenerator(
        generator, #static generator
        1.0, # ω_ref
        machine_oneDoneQ(), #machine
        shaft_no_damping(), #shaft
        avr_type2(), #avr
        tg_type1(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen_simple_avr(generator)
    return PSY.DynamicGenerator(
        generator, #static generator
        1.0, # ω_ref
        machine_oneDoneQ(), #machine
        shaft_no_damping(), #shaft
        avr_propr(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

# Add dynamic generators to the system (each gen is linked through a static one)
for g in get_components(Generator, threebus_sys)
    if get_number(get_bus(g)) == 102
        case_gen = dyn_gen_avr_type2(g)
        add_component!(threebus_sys, case_gen)
    elseif get_number(get_bus(g)) == 103
        case_gen = dyn_gen_simple_avr(g)
        add_component!(threebus_sys, case_gen)
    end
end

#Compute Y_bus after fault
fault_branches = deepcopy(collect(get_components(Branch, threebus_sys))[2:end])
Ybus_fault = PSY.Ybus(fault_branches, get_components(Bus, threebus_sys))[:, :]
