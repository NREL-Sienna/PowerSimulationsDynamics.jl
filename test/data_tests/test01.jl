using PowerSystems
using NLsolve
const PSY = PowerSystems

include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
omib_file_dir = joinpath(dirname(@__FILE__), "OMIB.raw")
omib_sys = System(PowerModelsData(omib_file_dir), runchecks = false)
add_source_to_ref(omib_sys)
res = solve_powerflow!(omib_sys)
############### Data Dynamic devices ########################
function dyn_gen_classic(generator)
    return PSY.DynamicGenerator(
        generator,
        1.0, #ω_ref
        machine_classic(), #machine
        shaft_damping(), #shaft
        avr_none(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

#Attach dynamic generator. Currently use PSS/e format based on bus #.
gen = [g for g in get_components(Generator, omib_sys)][1]
case_gen = dyn_gen_classic(gen)
add_component!(omib_sys, case_gen)

#Compute Y_bus after fault
fault_branch = deepcopy(collect(get_components(Branch, omib_sys))[1])
fault_branch.r = 0.02;
fault_branch.x = 0.1;
Ybus_fault = PSY.Ybus([fault_branch], get_components(Bus, omib_sys))[:, :]
