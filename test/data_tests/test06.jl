using PowerSystems
using NLsolve
const PSY = PowerSystems

include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
omib_file_dir = joinpath(dirname(@__FILE__), "OMIB_DARCO_PSR.raw")
omib_sys = System(PowerModelsData(omib_file_dir), runchecks = false)
gen_bus = get_components_by_name(Component, omib_sys, "BUS 2       ")[1]
gen_bus.bustype = BusTypes.PQ
add_source_to_ref(omib_sys)
res = solve_powerflow!(omib_sys, nlsolve)

function inv_DAIB(bus)
    return PSY.DynamicInverter(
        number = 1,
        name = "VSM",
        bus = get_bus(battery),
        ω_ref = 1.0,
        V_ref = get_voltage(get_bus(battery)),
        P_ref = get_activepower(battery),
        Q_ref = get_reactivepower(battery),
        MVABase = 2.75,
        converter = converter(),
        outer_control = outer_control(),
        inner_control = inner_control(),
        dc_source = dc_source(),
        freq_estimator = pll(),
        filter = filt(),
    )
end

darco_inverter = inv_DAIB(gen_bus)
add_component!(omib_sys, darco_inverter)
