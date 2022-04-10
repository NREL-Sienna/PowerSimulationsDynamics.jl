using PowerSystems
using NLsolve

include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
omib_file_dir = joinpath(dirname(@__FILE__), "OMIB.raw")
omib_sys = System(omib_file_dir, runchecks = false)
add_source_to_ref(omib_sys)
############### Data Dynamic devices ########################
function dyn_gen_classic(generator)
    return DynamicGenerator(
        name = get_name(generator),
        ω_ref = 1.0,
        machine = machine_classic(),
        shaft = shaft_damping(),
        avr = avr_none(),
        prime_mover = tg_none(),
        pss = pss_none(),
    )
end

#Attach dynamic generator. Currently use PSS/e format based on bus #.
gen = [g for g in get_components(Generator, omib_sys)][1]
case_gen = dyn_gen_classic(gen)
add_component!(omib_sys, case_gen, gen)

for l in get_components(PSY.PowerLoad, omib_sys)
    PSY.set_model!(l, PSY.LoadModels.ConstantImpedance)
end

#Compute Y_bus after fault
fault_branch = deepcopy(collect(get_components(Branch, omib_sys))[1])
fault_branch.r = 0.00;
fault_branch.x = 0.1;
Ybus_fault = PSY.Ybus([fault_branch], get_components(Bus, omib_sys))[:, :]
