using PowerSystems
using NLsolve
const PSY = PowerSystems

############### Data Network ########################
include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
threebus_file_dir = joinpath(dirname(@__FILE__), "ThreeBusNetwork.raw")
threebus_sys = System(PowerModelsData(threebus_file_dir), runchecks = false)
add_source_to_ref(threebus_sys)

### Case 2 Generators ###

function dyn_gen_simple_marconato(generator)
    return PSY.DynamicGenerator(
        generator, #static generator
        1.0, # ω_ref
        machine_simple_marconato(), #machine
        shaft_no_damping(), #shaft
        avr_type1(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

for g in get_components(Generator, threebus_sys)
    case_gen = dyn_gen_simple_marconato(g)
    add_component!(threebus_sys, case_gen)
end

#Compute Y_bus after fault
fault_branches = deepcopy(collect(get_components(Branch, threebus_sys))[2:end])
Ybus_fault = PSY.Ybus(fault_branches, get_components(Bus, threebus_sys))[:, :]
