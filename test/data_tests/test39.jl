using PowerSystems
using NLsolve
const PSY = PowerSystems

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/EXST1/TVC_System_32.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/EXST1/TVC_System.dyr")

sys = System(raw_file, dyr_file)
source = first(get_components(Source, sys))
PSY.set_X_th!(source, 0.01)

for l in get_components(PSY.PowerLoad, sys)
    PSY.set_model!(l, PSY.LoadModels.ConstantImpedance)
end

for gen in get_components(ThermalStandard, sys)
    #Find the generator at bus 2
    if get_number(get_bus(gen)) == 2
        dynamic_injector = get_dynamic_injector(gen)
        machine = deepcopy(get_machine(dynamic_injector))
        shaft = deepcopy(get_shaft(dynamic_injector))
        tg = deepcopy(get_prime_mover(dynamic_injector))
        pss = deepcopy(get_pss(dynamic_injector))
        avr = deepcopy(get_avr(dynamic_injector))
        new_dynamic_injector = DynamicGenerator(
            name = get_name(gen),
            machine = machine,
            shaft = shaft,
            avr = avr_exst1(),
            prime_mover = tg,
            ω_ref = 1.0,
            pss = pss,
        )
        set_dynamic_injector!(gen, nothing)
        add_component!(sys, new_dynamic_injector, gen)
    end
end
