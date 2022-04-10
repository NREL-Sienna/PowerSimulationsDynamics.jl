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
#Reduce generator output
for g in get_components(Generator, threebus_sys)
    g.active_power = 0.75
end
res = run_powerflow!(threebus_sys)

function dyn_gen_five_mass_shaft_order(generator)
    return PSY.DynamicGenerator(
        name = get_name(generator),
        ω_ref = 1.0, # ω_ref,
        machine = machine_oneDoneQ(), #machine
        shaft = shaft_fivemass(), #shaft
        avr = avr_type1(), #avr
        prime_mover = tg_none(),
        pss = pss_none(),
    ) #pss
end

function dyn_gen_first_order(generator)
    return PSY.DynamicGenerator(
        name = get_name(generator),
        ω_ref = 1.0, # ω_ref,
        machine = machine_oneDoneQ(), #machine
        shaft = shaft_damping(), #shaft
        avr = avr_type1(), #avr
        prime_mover = tg_none(), #tg
        pss = pss_none(),
    ) #pss
end

for g in get_components(Generator, threebus_sys)
    if get_number(get_bus(g)) == 103
        case_gen = dyn_gen_five_mass_shaft_order(g)
        add_component!(threebus_sys, case_gen, g)
    elseif get_number(get_bus(g)) == 102
        case_inv = dyn_gen_first_order(g)
        add_component!(threebus_sys, case_inv, g)
    end
end

for l in get_components(PSY.PowerLoad, threebus_sys)
    PSY.set_model!(l, PSY.LoadModels.ConstantImpedance)
end

#Compute Y_bus after fault
fault_branches = deepcopy(collect(get_components(Branch, threebus_sys))[2:end])
Ybus_fault = PSY.Ybus(fault_branches, get_components(Bus, threebus_sys))[:, :]
