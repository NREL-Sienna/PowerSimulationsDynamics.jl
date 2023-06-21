using PowerSystems
using NLsolve
const PSY = PowerSystems

############### Data Network ########################
include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
threebus_file_dir = joinpath(dirname(@__FILE__), "ThreeBusNetwork.raw")
threebus_sys = System(threebus_file_dir; runchecks = false)
add_source_to_ref(threebus_sys)

### Case 2 Generators ###

function dyn_gen_simple_anderson(generator)
    return PSY.DynamicGenerator(;
        name = get_name(generator), #static generator
        ω_ref = 1.0, # ω_ref
        machine = machine_simple_anderson(), #machine
        shaft = shaft_no_damping(), #shaft
        avr = avr_type1(), #avr
        prime_mover = tg_none(), #tg
        pss = pss_none(),
    ) #pss
end

for g in get_components(Generator, threebus_sys)
    case_gen = dyn_gen_simple_anderson(g)
    add_component!(threebus_sys, case_gen, g)
end

for l in get_components(PSY.StandardLoad, threebus_sys)
    transform_load_to_constant_impedance(l)
end

#Compute Y_bus after fault
fault_branches = deepcopy(collect(get_components(Branch, threebus_sys))[2:end])
Ybus_fault = PNM.Ybus(fault_branches, collect(get_components(Bus, threebus_sys)))[:, :]
