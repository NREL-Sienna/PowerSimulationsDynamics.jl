using PowerSystems
const PSY = PowerSystems

############### Data Network ########################
include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
threebus_file_dir = joinpath(dirname(@__FILE__), "ThreeBusNetwork.raw")
threebus_sys = System(PowerModelsData(threebus_file_dir), runchecks = false)
add_source_to_ref(threebus_sys)
res = solve_powerflow!(threebus_sys, nlsolve)
#Make line 3 the dynamic line

make_dynamic_branch!(branches[3], sys)

######## Machine Data #########
function dyn_gen_second_order(generator)
    return PSY.DynamicGenerator(
        1, #Number
        "Case9_$(get_name(generator))",
        get_bus(generator), #bus
        1.0, # ω_ref,
        1.0, #V_ref
        get_activepower(generator), #P_ref
        get_reactivepower(generator), #Q_ref
        machine_4th(), #machine
        shaft_no_damping(), #shaft
        avr_type1(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

############### Inverter Data ########################
function inv_case9(buses)
    return PSY.DynamicInverter(
        1, #Number
        "DARCO", #name
        buses[3], #bus
        1.0, # ω_ref,
        0.8, #V_ref
        0.5, #P_ref
        -0.3, #Q_ref
        100.0, #MVABase
        converter_case78(), #converter
        outer_control_test(), #outer control
        vsc_test(), #inner control voltage source
        dc_source_case78(), #dc source
        pll_test(), #pll
        filter_test(),
    ) #filter
end

for g in get_components(Generator, threebus_sys)
    case_gen = dyn_gen_second_order(g)
    add_component!(threebus_sys, case_gen)
end

#Compute Y_bus after fault
fault_branches = deepcopy(collect(get_components(Branch, threebus_sys))[2:end])
Ybus_fault = PSY.Ybus(fault_branches, get_components(Bus, threebus_sys))[:, :]
