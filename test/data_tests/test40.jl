using PowerSystems
using NLsolve
const PSY = PowerSystems

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/EXAC1/TVC_System_32.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/EXAC1/TVC_System.dyr")

sys = System(raw_file, dyr_file)
source = first(get_components(Source, sys))
PSY.set_X_th!(source, 0.01)

for l in get_components(PSY.PowerLoad, sys)
    PSY.set_model!(l, PSY.LoadModels.ConstantImpedance)
end
