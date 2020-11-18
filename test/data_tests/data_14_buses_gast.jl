using PowerSystems
using NLsolve
const PSY = PowerSystems

include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################

raw_file_dir = joinpath(dirname(@__FILE__), "../benchmarks/psse/GAST/14bus.raw")
dyr_file_dir = joinpath(dirname(@__FILE__), "../benchmarks/psse/GAST/14_bus_dyn_data.dyr")
sys = System(raw_file_dir, dyr_file_dir)

#Compute Ybus_fault
sys2 = System(raw_file_dir)
#Remove line connecting bus 1 to 2.
remove_component!(Line, sys2, "BUS 02-BUS 04-i_4")
Ybus_fault = Ybus(sys2).data
