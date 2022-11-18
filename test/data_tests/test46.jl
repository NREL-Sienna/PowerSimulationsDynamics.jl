using PowerSystems
using NLsolve
const PSY = PowerSystems

############### Data Network ########################
include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
sys_dir = joinpath(dirname(@__FILE__), "OMIB_ActiveLoad.raw")
sys = System(sys_dir, runchecks = false)

function inv_model(static_device)
    return DynamicInverter(
        name = get_name(static_device),
        ω_ref = 1.0, # ω_ref,
        converter = converter_high_power(), #converter
        outer_control = outer_control(), #outer control
        inner_control = inner_control(), #inner control voltage source
        dc_source = dc_source_lv(), #dc source
        freq_estimator = pll(), #pll
        filter = filt(), #filter
    )
end

for g in get_components(Generator, sys)
    case_gen = inv_model(g)
    add_component!(sys, case_gen, g)
end

for l in get_components(PSY.PowerLoad, sys)
    dyn_load = ActiveLoad(l)
    add_component!(sys, dyn_load, l)
end
