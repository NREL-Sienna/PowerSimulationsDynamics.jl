using PowerSystems

##Nodes
nodes_DAIB= [Bus(1 , #number
                 "Bus 1", #Name
                 "REF" , #BusType (REF, PV, PQ)
                 0, #Angle in radians
                 1.04, #Voltage in pu
                 (min=0.94, max=1.06), #Voltage limits in pu
                 0.69),  #Base voltage in kV
                 Bus(2 , "Bus 2"  , "PV" ,  0 , 1.0 , (min=0.94, max=1.06), 0.69)]


branch_DAIB = [Line("Line1", #name
                    true, #available
                    0.0, #active power flow initial condition (from-to)
                    0.0, #reactive power flow initial condition (from-to)
                    Arc(from=nodes_DAIB[1],to=nodes_DAIB[2]), #Connection between buses
                    0.0, #resistance in pu
                    0.075, #reactance in pu
                    (from=0.0, to=0.0), #susceptance in pu
                    5.0, #rate in MW
                    1.04)]  #angle limits (-min and max)

##Inf bus
inf_gen_DAIB = StaticSource(1, #number
                    :InfBus, #name
                    nodes_DAIB[1],#bus
                    1.00, #VR
                    0.0, #VI
                    0.000005) #Xth

## Inverter

converter = AvgCnvFixedDC(690.0, #Rated Voltage
                          2.75) #Rated MVA
dc_source = FixedDCSource(600.0) #Not in the original data, guessed.

filter = LCLFilter(0.08, #Series inductance lf in pu
                   0.003, #Series resitance rf in pu
                   0.074, #Shunt capacitance cf in pu
                   0.2, #Series ractance rg to grid connection (#Step up transformer or similar)
                   0.01) #Series resistance lg to grid connection (#Step up transformer or similar)
pll = PLL(500.0, #ω_lp: Cut-off frequency for LowPass filter of PLL filter.
          0.084, #k_p: PLL proportional gain
          4.69) #k_i: PLL integral gain

virtual_H = VirtualInertia(2.0, #Ta:: VSM inertia constant
                           400.0, #kd:: VSM damping coefficient
                           20.0, #kω:: Frequency droop gain in pu
                           2*pi*50.0) #ωb:: Rated angular frequency

Q_control = ReactivePowerDroop(0.2, #kq:: Reactive power droop gain in pu
                              1000.0) #ωf:: Reactive power cut-off low pass filter frequency

outer_control = VirtualInertiaQdroop(virtual_H, Q_control)

vsc = CombinedVIwithVZ(0.59, #kpv:: Voltage controller proportional gain
                       736.0, #kiv:: Voltage controller integral gain
                       0.0, #kffv:: Binary variable enabling the voltage feed-forward in output of current controllers
                       0.0, #rv:: Virtual resistance in pu
                       0.2, #lv: Virtual inductance in pu
                       1.27, #kpc:: Current controller proportional gain
                       14.3, #kiv:: Current controller integral gain
                       0.0, #kffi:: Binary variable enabling the current feed-forward in output of current controllers
                       50.0, #ωad:: Active damping low pass filter cut-off frequency
                       0.2) #kad:: Active damping gain

Darco_Inverter = DynInverter(1, #number
                             :DARCO, #name
                             nodes_DAIB[2], #bus location
                             1.0, #ω_ref
                             1.02, #V_ref
                             0.5, #P_ref
                             0.0, #Q_ref
                             2.75, #MVABase
                             converter, #Converter
                             outer_control, #OuterControl
                             vsc, #Voltage Source Controller
                             dc_source, #DC Source
                             pll, #Frequency Estimator
                             filter) #Output Filter

loads_DAIB = [PowerLoad("Bus1", true, nodes_DAIB[1], PowerSystems.ConstantPower, 0.0, 0.0, 0.0, 0.0)]

DAIB = DynamicSystem(nodes_DAIB, branch_DAIB, [Darco_Inverter], [inf_gen_DAIB], 100.0, 50.0);

#DAIB = DynamicSystem(nodes_DAIB, branch_DAIB, [Darco_Inverter], loads_DAIB, 100.0, 50.0)
