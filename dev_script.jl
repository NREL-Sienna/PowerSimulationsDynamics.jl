using PowerSystems
using NLsolve
const PSY = PowerSystems
omib_file_dir = joinpath(pwd(), "test/data_tests/ThreeBusNetwork.raw")
wscc_file_dir = joinpath(pwd(), "test/data_tests/WSCC 9 bus.raw")
# We Should export this function
omib_sys = System(PowerModelsData(omib_file_dir), runchecks = false)

wscc_sys = System(PSY.PowerModelsData(wscc_file_dir))
solve_powerflow!(omib_sys, nlsolve)
solve_powerflow!(wscc_sys, nlsolve, method = :newton)
