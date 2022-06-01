using PowerSystems
using NLsolve
const PSY = PowerSystems

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/STAB1/OMIB_SSS.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/STAB1/OMIB_SSS.dyr")

for gen in get_components(ThermalStandard, sys)
    #Find the generator at bus 1
    if get_number(get_bus(gen)) == 1
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
            avr = avr,
            prime_mover = tg,
            ω_ref = 1.0,
            pss = pss_stab1(),
        )
        set_dynamic_injector!(gen, nothing)
        add_component!(sys, new_dynamic_injector, gen)
    end
end
