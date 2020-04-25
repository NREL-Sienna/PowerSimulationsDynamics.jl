using PowerSystems
using NLsolve
const PSY = PowerSystems

include("test/data_tests/dynamic_test_data.jl")
include("test/data_tests/data_utils.jl")
############### Data Network ########################
omib_file_dir = joinpath(pwd(), "test/data_tests/OMIB.raw")
omib_sys = System(PowerModelsData(omib_file_dir), runchecks=false)
add_source_to_ref(omib_sys)
res = solve_powerflow!(omib_sys, nlsolve)
############### Data Dynamic devices ########################
function dyn_gen_first_order(generator)
    return PSY.DynamicGenerator(
        1, #Number
        "Case1Gen",
        get_bus(generator), #bus
        1.0, # ω_ref,
        1.0, #V_ref
        get_activepower(generator), #P_ref
        get_reactivepower(generator), #Q_ref
        machine_OMIB(), #machine
        shaft_damping(), #shaft
        avr_none(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

#Attach dynamic generator. Currently use PSS/e format based on bus #.
gen = [g for g in get_components(Generator, omib_sys)][1]
case1_gen = dyn_gen_first_order(gen)
add_component!(omib_sys, case1_gen)
#Compute Y_bus after fault
fault_branch = deepcopy(collect(get_components(Branch, omib_sys))[1])
fault_branch.r = 0.02; fault_branch.x = 0.1
Ybus_fault = PSY.Ybus([fault_branch], get_components(Bus, omib_sys))[:, :]
######################### Dynamical System ########################
