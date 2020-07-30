using PowerSystems
using NLsolve
const PSY = PowerSystems

raw_file_dir = joinpath(dirname(@__FILE__), "../benchmarks/psse/GENROU/ThreeBusMulti.raw")
dyr_file_dir = joinpath(dirname(@__FILE__), "../benchmarks/psse/GENROU/ThreeBus_GENROU.dyr")

#Create System
sys = System(raw_file_dir, dyr_file_dir)

#Compute Ybus_fault
sys2 = System(raw_file_dir)
#Remove line connecting bus 1 to 2.
remove_component!(Line, sys2, "1")
Ybus_fault = Ybus(sys2).data
