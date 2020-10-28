using PowerSystems
using NLsolve
const PSY = PowerSystems

include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################

raw_file_dir = joinpath(dirname(@__FILE__), "../benchmarks/psse/GENROU/ThreeBusMulti.raw")
sys = System(raw_file_dir);

#Replace Gen101 by Source
remove_component!(ThermalStandard, sys, "generator-101-1");
add_source_to_ref(sys)

function dyn_gen_genrou(generator)
    return PSY.DynamicGenerator(
        name = get_name(generator),
        ω_ref = 1.0, #ω_ref
        machine = machine_genrou(), #machine
        shaft = shaft_genrou(), #shaft
        avr = avr_type1(), #avr
        prime_mover = tg_none(), #tg
        pss = pss_none(),
    ) #pss
end

#Add GENROU to System
g = get_component(ThermalStandard, sys, "generator-102-1")
dyn_gen = dyn_gen_genrou(g)
add_component!(sys, dyn_gen, g);

#Compute Ybus_fault
sys2 = System(raw_file_dir)
#Remove line connecting bus 1 to 2.
remove_component!(Line, sys2, "1")
Ybus_fault = Ybus(sys2).data
