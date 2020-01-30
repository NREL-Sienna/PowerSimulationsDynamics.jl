nodes_OMIB = [
    PSY.Bus(
        1, #number
        "Bus 1", #Name
        "REF", #BusType (REF, PV, PQ)
        0, #Angle in radians
        1.06, #Voltage in pu
        (min = 0.94, max = 1.06), #Voltage limits in pu
        69,
    ), #Base voltage in kV
    PSY.Bus(2, "Bus 2", "PV", 0, 1.045, (min = 0.94, max = 1.06), 69),
]

branch_OMIB = [PSY.Line(
    "Line1", #name
    true, #available
    0.0, #active power flow initial condition (from-to)
    0.0, #reactive power flow initial condition (from-to)
    Arc(from = nodes_OMIB[1], to = nodes_OMIB[2]), #Connection between buses
    0.01, #resistance in pu
    0.05, #reactance in pu
    (from = 0.0, to = 0.0), #susceptance in pu
    18.046, #rate in MW
    1.04,
)]  #angle limits (-min and max)

#### Generators #####
Basic = PSY.BaseMachine(
    0.0, #R
    0.2995, #Xd_p
    1.05, #eq_p
    615.0,
)  #MVABase

BaseShaft = PSY.SingleMass(
    5.148, #H
    2.0,
) #D

fixed_avr = PSY.AVRFixed(1.05) #Emf

proportional_avr = PSY.AVRSimple(5000.0) #Kv

fixed_tg = PSY.TGFixed(1.0) #eff

no_pss = PSY.PSSFixed(0.0)

Gen1AVR = PSY.DynamicGenerator(
    1, #Number
    "TestGen",
    nodes_OMIB[2],#bus
    1.0, # ω_ref,
    1.05,
    0.4,
    Basic,
    BaseShaft,
    proportional_avr, #avr
    fixed_tg, #tg
    no_pss,
);

Gen2AVR = PSY.DynamicGenerator(
    1, #Number
    "TestGen2",
    nodes_OMIB[2],#bus
    1.0, # ω_ref,
    1.05,
    0.4,
    Basic,
    BaseShaft,
    proportional_avr, #avr
    fixed_tg, #tg
    no_pss,
);

#### Inverter ####

converter = PSY.AvgCnvFixedDC(690.0, 2750000.0)

dc_source = PSY.FixedDCSource(1500.0) #Not in the original data, guessed.

filt = PSY.LCLFilter(
    0.08, #Series inductance lf in pu
    0.003, #Series resitance rf in pu
    0.074, #Shunt capacitance cf in pu
    0.2, #Series reactance rg to grid connection (#Step up transformer or similar)
    0.01,
) #Series resistance lg to grid connection (#Step up transformer or similar)

pll = PSY.PLL(
    500.0, #ω_lp: Cut-off frequency for LowPass filter of PLL filter.
    0.084, #k_p: PLL proportional gain
    4.69,
) #k_i: PLL integral gain

virtual_H = PSY.VirtualInertia(
    2.0, #Ta:: VSM inertia constant
    400.0, #kd:: VSM damping coefficient
    20.0, #kω:: Frequency droop gain in pu
    2 * pi * 50.0,
) #ωb:: Rated angular frequency

Q_control = ReactivePowerDroop(
    0.2, #kq:: Reactive power droop gain in pu
    1000.0,
) #ωf:: Reactive power cut-off low pass filter frequency

outer_control = VirtualInertiaQdroop(virtual_H, Q_control)

vsc = CombinedVIwithVZ(
    0.59, #kpv:: Voltage controller proportional gain
    736.0, #kiv:: Voltage controller integral gain
    0.0, #kffv:: Binary variable enabling the voltage feed-forward in output of current controllers
    0.0, #rv:: Virtual resistance in pu
    0.2, #lv: Virtual inductance in pu
    1.27, #kpc:: Current controller proportional gain
    14.3, #kiv:: Current controller integral gain
    0.0, #kffi:: Binary variable enabling the current feed-forward in output of current controllers
    50.0, #ωad:: Active damping low pass filter cut-off frequency
    0.2,
) #kad:: Active damping gain

test_inverter = PSY.DynamicInverter(
    2, #number
    "DARCO", #name
    nodes_OMIB[1], #bus location
    1.0, #ω_ref
    1.02, #V_ref
    0.5, #P_ref
    0.0, #Q_ref
    100.0, #MVABase
    converter, #Converter
    outer_control, #OuterControl
    vsc, #Voltage Source Controller
    dc_source, #DC Source
    pll, #Frequency Estimator
    filt,
); #Output Filter


@testset "Dynamic Generator in System" begin

    sys = PSY.System(100)
    for bus in nodes_OMIB
        PSY.add_component!(sys, bus)
    end
    for lines in branch_OMIB
        PSY.add_component!(sys, lines)
    end
    PSY.add_component!(sys, Gen1AVR)
    PSY.add_component!(sys, Gen2AVR)
    PSY.add_component!(sys, test_inverter)

    @test collect(PSY.get_components(DynamicGenerator, sys))[1] == Gen1AVR

end
