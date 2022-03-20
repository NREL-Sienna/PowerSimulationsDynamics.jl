using PowerSystems
using NLsolve
const PSY = PowerSystems

############### Data Network ########################
include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
threebus_file_dir = joinpath(dirname(@__FILE__), "ThreeBusNetwork.raw")
threebus_sys = System(threebus_file_dir, runchecks = false)
add_source_to_ref(threebus_sys)
res = run_powerflow!(threebus_sys)

### Case 2 Generators ###

function dyn_gen_oneDoneQ(generator)
    return PSY.DynamicGenerator(
        name = get_name(generator),
        ω_ref = 1.0,
        machine = machine_oneDoneQ(), #machine
        shaft = shaft_no_damping(), #shaft
        avr = avr_type1(), #avr
        prime_mover = tg_none(), #tg
        pss = pss_none(),
    ) #pss
end

# Add dynamic generators to the system (each gen is linked through a static one)
for g in get_components(Generator, threebus_sys)
    case_gen = dyn_gen_oneDoneQ(g)
    add_component!(threebus_sys, case_gen, g)
end

#Compute Y_bus after fault
fault_branches = deepcopy(collect(get_components(Branch, threebus_sys))[2:end])
Ybus_fault = PSY.Ybus(fault_branches, get_components(Bus, threebus_sys))[:, :]
