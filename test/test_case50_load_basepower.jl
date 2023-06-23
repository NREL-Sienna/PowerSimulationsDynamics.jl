"""
Test basepowers load models:
Builds a OMIB system with a load. Tests that the active power flow
across the single branch does not change when the base power of the
load is changed. Includes a LoadChange and LoadTrip to ensure that
base power scaling is handled properly in the perturbation.
"""

include(joinpath(TEST_FILES_DIR, "data_tests", "dynamic_test_data.jl"))

function standard_load(b)
    return StandardLoad(;
        name = "test-load",
        available = true,
        bus = b,
        base_power = 100.0,
        constant_active_power = 0.0,
        constant_reactive_power = 0.0,
        impedance_active_power = 0.1,
        impedance_reactive_power = 0.01,
        current_active_power = 0.0,
        current_reactive_power = 0.0,
        max_constant_active_power = 0.0,
        max_constant_reactive_power = 0.0,
        max_impedance_active_power = 1.0,
        max_impedance_reactive_power = 1.0,
        max_current_active_power = 0.0,
        max_current_reactive_power = 0.0,
    )
end

function dyn_gen_classic(generator)
    return DynamicGenerator(;
        name = get_name(generator),
        ω_ref = 1.0,
        machine = machine_classic(),
        shaft = shaft_damping(),
        avr = avr_none(),
        prime_mover = tg_none(),
        pss = pss_none(),
    )
end

OMIB_dir = joinpath(TEST_FILES_DIR, "data_tests", "OMIB.raw")
omib_sys = System(OMIB_dir; runchecks = false)
add_source_to_ref(omib_sys)
b = get_component(Bus, omib_sys, "BUS 2")
l = standard_load(b)
add_component!(omib_sys, l)
gen = [g for g in get_components(Generator, omib_sys)][1]
case_gen = dyn_gen_classic(gen)
add_component!(omib_sys, case_gen, gen)

pert = LoadChange(0.1, l, :P_ref_impedance, 0.3)
sim = Simulation!(MassMatrixModel, omib_sys, pwd(), (0.0, 1.0), pert)
execute!(sim, Rodas5())
results = read_results(sim)
_, P_refchange_original = get_activepower_branch_flow(results, "BUS 1-BUS 2-i_1", :to)

pert = LoadTrip(0.1, l)
sim = Simulation!(MassMatrixModel, omib_sys, pwd(), (0.0, 1.0), pert)
execute!(sim, Rodas5())
results = read_results(sim)
_, P_trip_original = get_activepower_branch_flow(results, "BUS 1-BUS 2-i_1", :to)

# Change base power of load, setpoints, and perturbation such that
# the result should be identical to the original system.
OMIB_dir = joinpath(TEST_FILES_DIR, "data_tests", "OMIB.raw")
omib_sys = System(OMIB_dir; runchecks = false)
add_source_to_ref(omib_sys)
b = get_component(Bus, omib_sys, "BUS 2")
l = standard_load(b)
add_component!(omib_sys, l)
gen = [g for g in get_components(Generator, omib_sys)][1]
case_gen = dyn_gen_classic(gen)
add_component!(omib_sys, case_gen, gen)

set_base_power!(l, get_base_power(l) / 2)
set_impedance_active_power!(l, get_impedance_active_power(l) * 2)
set_impedance_reactive_power!(l, get_impedance_reactive_power(l) * 2)
pert = LoadChange(0.1, l, :P_ref_impedance, 0.6)

sim = Simulation!(MassMatrixModel, omib_sys, pwd(), (0.0, 1.0), pert)
execute!(sim, Rodas5())
results = read_results(sim)
_, P_refchange_new = get_activepower_branch_flow(results, "BUS 1-BUS 2-i_1", :to)

pert = LoadTrip(0.1, l)
sim = Simulation!(MassMatrixModel, omib_sys, pwd(), (0.0, 1.0), pert)
execute!(sim, Rodas5())
results = read_results(sim)
_, P_trip_new = get_activepower_branch_flow(results, "BUS 1-BUS 2-i_1", :to)

#Test starting power is the same
@test isapprox(P_refchange_original[1], P_refchange_new[1], atol = 1e-10)
#Test end power is the same after load reference change
@test isapprox(P_refchange_original[end], P_refchange_new[end], atol = 1e-10)

#Test starting power is the same
@test isapprox(P_trip_original[1], P_trip_new[1], atol = 1e-10)
#Test end power is the same after load trip
@test isapprox(P_trip_original[end], P_trip_new[end], atol = 1e-10)
