using PowerSystems, PowerSystemCaseBuilder, Sundials, Plots, Logging, InfrastructureSystems 
const PSY = PowerSystems


using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics


avr_exst1() = PSY.EXST1_PTI(
    Tr = 0.01,
    Vi_lim = (-0.2,0.2) , 
    Tc = 1.0,
    Tb = 10.0,
    Ka = 200.0,
    Ta = 0.01,
    Vr_lim = (0.0,6.2) , 
    Kc = 0.15,
    Kf = 0.0,
    Tf = 1.0
)

avr_ex4vsa1() = PSY.EX4VSA(
    Iflim=1.2,
    d=-0.1,
    f=1.0,
    Spar=0.0,
    K1=1.0,
    K2=-1.0,
    Oel_lim = (-12.0,10.0),
    G=70.0,
    Ta=10.0,
    Tb=20.0,
    Te=0.1,
    E_lim = (0.0,4.0)
)


pss_stab1() = PSY.STAB1(
    Ks = 0,
    Ts = 1.5,
    T1T3 = 13,
    T3 = 0.045,
    T2T4 = 13,
    T4 = 0.045,
    H_lim=0.2
)


  synmachineRR()= PSY.RoundRotorQuadratic(
    R= 0.001,
    Td0_p=5.0,
    Td0_pp=0.04,
    Tq0_p=0.4,
    Tq0_pp=0.04,
    Xd= 1.9,
    Xq= 1.8,
    Xd_p= 0.3,
    Xq_p=0.5,
    Xd_pp=0.2,
    Xl = 0.15,
    Se = (0.08,0.45)
    )




sys = PowerSystemCaseBuilder.build_system(
                       PowerSystemCaseBuilder.PSSETestSystems,
                       "psse_OMIB_sys",
                       )
                       
thermal_gen = first(get_components(ThermalStandard, sys))
dynamic_injector = get_dynamic_injector(thermal_gen)
machine = deepcopy(get_machine(dynamic_injector))
shaft = deepcopy(get_shaft(dynamic_injector))
tg = deepcopy(get_prime_mover(dynamic_injector))
pss = deepcopy(get_pss(dynamic_injector))

new_dynamic_injector = DynamicGenerator(name = get_name(thermal_gen), machine = synmachineRR(), shaft=shaft, avr = avr_ex4vsa1(), prime_mover = tg, ω_ref = 1.0, pss = pss_stab1())
set_dynamic_injector!(thermal_gen, nothing) # delete gen
add_component!(sys, new_dynamic_injector, thermal_gen)

sys

time_span = (0.0, 20.0)	
perturbation_trip = BranchTrip(1.0, Line, "BUS 1-BUS 2-i_1")

sim = Simulation!(ResidualModel, sys, pwd(), time_span, perturbation_trip)

x0_init = read_initial_conditions(sim)
show_states_initial_value(sim)


execute!(sim, IDA(), dtmax = 0.02, saveat = 0.02, enable_progress_bar = false)


results = read_results(sim)

angle = get_state_series(results, ("generator-102-1", :δ));

Voutput = get_state_series(results, ("generator-102-1", :Vex));

Voel = get_state_series(results, ("generator-102-1", :oel));

PSS_out = get_state_series(results, ("generator-102-1", :x_p3));


plot(angle, xlabel = "time", ylabel = "rotor angle [rad]", label = "gen-102-1")

plot(Voutput , xlabel = "time", ylabel = "field voltage [pu]", label = "gen-102-1")

plot(Voel , xlabel = "time", ylabel = "OEL [pu]", label = "gen-102-1")

plot(PSS_out , xlabel = "time", ylabel = "Vpss [pu]", label = "gen-102-1")