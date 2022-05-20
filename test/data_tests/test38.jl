using PowerSystems
using NLsolve
const PSY = PowerSystems

############### Data Network ########################
include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
raw_file = joinpath(dirname(@__FILE__), "ThreeBusMultiLoad.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/DERA/ThreeBus_DERA.dyr")

#Eventually, include DERA in the dyr File (need parser). Include zero inertia generator at bus 1 which is represented as IB. 
threebus_sys = System(raw_file, runchecks = false)
for g in get_components(ThermalStandard, threebus_sys)
    g.bus.bustype == BusTypes.REF && remove_component!(threebus_sys, g)
end
add_source_to_ref(threebus_sys)
