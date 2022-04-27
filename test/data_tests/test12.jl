using PowerSystems
using NLsolve
const PSY = PowerSystems

############### Data Network ########################
include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
threebus_file_dir = joinpath(dirname(@__FILE__), "ThreeBusMulti.raw")
threebus_sys = System(threebus_file_dir, runchecks = false)

function dyn_gen_multi(generator)
    return PSY.DynamicGenerator(
        name = get_name(generator),
        ω_ref = 1.0, # ω_ref,
        machine = machine_classic(), #machine
        shaft = shaft_damping(), #shaft
        avr = avr_none(), #avr
        prime_mover = tg_none(), #tg
        pss = pss_none(),
    )
end

function dyn_gen_multi_tg(generator)
    return PSY.DynamicGenerator(
        name = get_name(generator),
        ω_ref = 1.0, # ω_ref,
        machine = machine_classic(), #machine
        shaft = shaft_damping(), #shaft
        avr = avr_none(), #avr
        prime_mover = tg_type2(), #tg
        pss = pss_none(),
    )
end

# Add dynamic generators to the system (each gen is linked through a static one)
for g in get_components(Generator, threebus_sys)
    if get_number(get_bus(g)) == 101
        case_gen = dyn_gen_multi(g)
        add_component!(threebus_sys, case_gen, g)
    elseif get_number(get_bus(g)) == 102
        case_gen = dyn_gen_multi_tg(g)
        add_component!(threebus_sys, case_gen, g)
    end
end

for l in get_components(PSY.PowerLoad, threebus_sys)
    PSY.set_model!(l, PSY.LoadModels.ConstantImpedance)
end

#Create Ybus_Fault
fault_branches = deepcopy(collect(get_components(Branch, threebus_sys)))
for br in fault_branches
    if get_name(br) == "BUS 1-BUS 2-i_1"
        br.r = 4 * br.r
        br.x = 4 * br.x
    end
end
Ybus_fault = PSY.Ybus(fault_branches, get_components(Bus, threebus_sys))[:, :];
