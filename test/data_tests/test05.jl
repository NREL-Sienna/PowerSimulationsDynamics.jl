using PowerSystems
const PSY = PowerSystems

############### Data Network ########################
include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
threebus_file_dir= joinpath(dirname(@__FILE__), "ThreeBusNetwork.raw")
threebus_sys = System(PowerModelsData(threebus_file_dir), runchecks=false)
add_source_to_ref(threebus_sys)
#Reduce generator output
for g in get_components(Generator, threebus_sys)
    g.activepower = 0.5
end
res = solve_powerflow!(threebus_sys, nlsolve)

function dyn_gen_five_mass_shaft_order(generator)
    return PSY.DynamicGenerator(
        1, #Number
        "Case5_$(get_name(generator))",
        get_bus(generator), #bus
        1.0, # ω_ref,
        get_voltage(get_bus(generator)), #V_ref
        get_activepower(generator), #P_ref
        get_reactivepower(generator), #Q_ref
        machine_4th(), #machine
        shaft_fivemass(), #shaft
        avr_type1(), #avr
        tg_type2(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen_first_order(generator)
    return PSY.DynamicGenerator(
        1, #Number
        "Case1Gen",
        get_bus(generator), #bus
        1.0, # ω_ref,
        get_voltage(get_bus(generator)), #V_ref
        get_activepower(generator), #P_ref
        get_reactivepower(generator), #Q_ref
        machine_OMIB(), #machine
        shaft_damping(), #shaft
        avr_none(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

gens = collect(get_components(Generator, threebus_sys))
case5_gen = dyn_gen_five_mass_shaft_order(gens[2])
add_component!(threebus_sys, case5_gen)
case1_gen = dyn_gen_first_order(gens[1])
add_component!(threebus_sys, case1_gen)

#Compute Y_bus after fault
fault_branches = deepcopy(collect(get_components(Branch, threebus_sys))[2:end])
Ybus_fault = PSY.Ybus(fault_branches, get_components(Bus, threebus_sys))[:, :]
