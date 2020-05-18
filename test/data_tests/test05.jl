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
#Reduce generator output
for g in get_components(Generator, threebus_sys)
    g.activepower = 0.75
end
res = solve_powerflow!(threebus_sys, nlsolve)

function dyn_gen_five_mass_shaft_order(generator)
    return PSY.DynamicGenerator(
        generator,
        1.0, # ω_ref,
        machine_4th(), #machine
        shaft_fivemass(), #shaft
        avr_type1(), #avr
        tg_none(),
        pss_none(),
    ) #pss
end

function dyn_gen_first_order(generator)
    return PSY.DynamicGenerator(
        generator,
        1.0, # ω_ref,
        machine_4th(), #machine
        shaft_damping(), #shaft
        avr_type1(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

gens = collect(get_components(Generator, threebus_sys))
#At node 3
case5_gen = dyn_gen_five_mass_shaft_order(gens[1])
add_component!(threebus_sys, case5_gen)
#At node 2
case1_gen = dyn_gen_first_order(gens[2])
add_component!(threebus_sys, case1_gen)

#Compute Y_bus after fault
fault_branches = deepcopy(collect(get_components(Branch, threebus_sys))[2:end])
Ybus_fault = PSY.Ybus(fault_branches, get_components(Bus, threebus_sys))[:, :]
