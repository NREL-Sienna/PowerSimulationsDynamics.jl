using PowerSystems
using NLsolve
const PSY = PowerSystems

include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
omib_file_dir= joinpath(dirname(@__FILE__), "OMIB_DARCO_PSR.raw")
omib_sys = System(PowerModelsData(omib_file_dir), runchecks=false)
gen_bus = get_components_by_name(Component, omib_sys, "BUS 2       ")[1]
gen_bus.bustype = BusTypes.PQ
add_source_to_ref(omib_sys)
res = solve_powerflow!(omib_sys, nlsolve)

function inv_DAIB(bus)
    return PSY.DynamicInverter(
        1, #Number
        "DARCO", #name
        bus, #bus
        1.0, # ω_ref,
        get_voltage(bus), #V_ref
        0.5, #P_ref
        0.0, #Q_ref
        2.75, #MVABase
        converter_DAIB(), #converter
        outer_control_test(), #outer control
        vsc_test(), #inner control voltage source
        dc_source_DAIB(), #dc source
        pll_test(), #pll
        filter_test(),
    ) #filter
end

darco_inverter = inv_DAIB(gen_bus)
add_component!(omib_sys, darco_inverter)
