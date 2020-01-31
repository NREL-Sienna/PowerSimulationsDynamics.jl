"""
Case 3:
This case study a three bus system with 2 machines (Simple Marconato: 6th order model) and an infinite source.
The fault drop the connection between buses 1 and 3, eliminating the direct connection between the infinite source
and the generator located in bus 3.
"""

##################################################
############### LOAD DATA ########################
##################################################

############### Data Network ########################

nodes_case234 = [PSY.Bus(1 , "Bus 1"  , "REF" , 0 , 1.02  , (min=0.94, max=1.06), 138),
                 PSY.Bus(2 , "Bus 2"  , "PV" , 0 , 1.00 , (min=0.94, max=1.06), 138),
                 PSY.Bus(3 , "Bus 3"  , "PV" , 0 , 1.00 , (min=0.94, max=1.06), 138)]

branch_case234 = [PSY.Line("Line1", true, 0.0, 0.0, Arc(from=nodes_case234[1], to=nodes_case234[3]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04),
                  PSY.Line("Line2", true, 0.0, 0.0, Arc(from=nodes_case234[1], to=nodes_case234[2]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04),
                  PSY.Line("Line3", true, 0.0, 0.0, Arc(from=nodes_case234[2], to=nodes_case234[3]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04)]

#Trip of Line 1.
branch_case234_fault = [PSY.Line("Line2", true, 0.0, 0.0, Arc(from=nodes_case234[1], to=nodes_case234[2]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04),
                        PSY.Line("Line3", true, 0.0, 0.0, Arc(from=nodes_case234[2], to=nodes_case234[3]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04)]

loads_case234 = [PSY.PowerLoad("Bus1", true, nodes_case234[1], PowerSystems.ConstantPower, 1.5, 0.8, 1.5, 0.8),
                 PSY.PowerLoad("Bus2", true, nodes_case234[2], PowerSystems.ConstantPower, 1.5, 0.7, 1.5, 0.8),
                 PSY.PowerLoad("Bus3", true, nodes_case234[3], PowerSystems.ConstantPower, 0.5, 0.3, 0.5, 0.3)]


############### Data devices ########################

inf_gen_case234 = PSY.Source(
                "InfBus", #name
                true, #availability
                nodes_case234[1], #bus
                1.02, #VR
                0.0, #VI
                0.000001) #Xth

######## Machine Data #########

### Case 3: 6th Order Model with AVR (3-bus case) ###
case3_machine2 = PSY.SimpleMarconatoMachine(0.0,
                                        1.3125, #Xd
                                        1.2578, #Xq
                                        0.1813, #Xd_p
                                        0.25, #Xq_p
                                        0.14, #Xd_pp
                                        0.18, #Xq_pp
                                        5.89, #Td0_p
                                        0.6, #Tq0_p
                                        0.5, #Td0_pp
                                        0.023, #Tq0_pp
                                        0.0, #T_AA
                                        100.0) #MVABase

case3_machine3 = PSY.SimpleMarconatoMachine(0.0,
                                        1.3125, #Xd
                                        1.2578, #Xq
                                        0.1813, #Xd_p
                                        0.25, #Xq_p
                                        0.14, #Xd_pp
                                        0.18, #Xq_pp
                                        5.89, #Td0_p
                                        0.6, #Tq0_p
                                        0.5, #Td0_pp
                                        0.023, #Tq0_pp
                                        0.0, #T_AA
                                        100.0) #MVABase

######## Shaft Data #########

### Shafts for Gens 2 and 3: Cases 2, 3 and 4 ###
case234_shaft2 = PSY.SingleMass(3.01, #H (M = 6.02 -> H = M/2)
                            0.0) #D


case234_shaft3 = PSY.SingleMass(3.01, #H (M = 6.02 -> H = M/2)
                            0.0) #D

######## PSS Data #########
cases_no_pss = PSY.PSSFixed(0.0)


######## TG Data #########

### No TG for Cases 1, 2, 3, 4 ###
case1234_no_tg = PSY.TGFixed(1.0) #eff


########  AVR Data #########
### AVRs for Case 2, 3, 4 and 5 ###
case2345_avr2 = PSY.AVRTypeI(20.0, #Ka - Gain
                        0.01, #Ke
                        0.063, #Kf
                        0.2, #Ta
                        0.314, #Te
                        0.35, #Tf
                        0.001, #Tr
                        5.0, #Vrmax
                        -5.0, #Vrmin
                        0.0039, #Ae - 1st ceiling coefficient
                        1.555) #Be - 2nd ceiling coefficient

case2345_avr3 = PSY.AVRTypeI(20.0, #Ka - Gain
                        0.01, #Ke
                        0.063, #Kf
                        0.2, #Ta
                        0.314, #Te
                        0.35, #Tf
                        0.001, #Tr
                        5.0, #Vrmax
                        -5.0, #Vrmin
                        0.0039, #Ae - 1st ceiling coefficient
                        1.555) #Be - 2nd ceiling coefficient

### Case 3 Generators ###
case3_gen2 = PSY.DynamicGenerator(1, #Number
                         "Case3Gen2",
                         nodes_case234[2], #bus
                         1.0, # ω_ref,
                         1.0142, #V_ref
                         1.0, #P_ref
                         0.0, #Q_ref
                         case3_machine2, #machine
                         case234_shaft2, #shaft
                         case2345_avr2, #avr
                         case1234_no_tg, #tg
                         cases_no_pss) #pss

case3_gen3 = PSY.DynamicGenerator(2, #Number
                         "Case3Gen3",
                         nodes_case234[3], #bus
                         1.0, # ω_ref,
                         1.0059, #V_ref
                         1.0, #P_ref
                         0.0, #Q_ref
                         case3_machine3, #machine
                         case234_shaft3, #shaft
                         case2345_avr3, #avr
                         case1234_no_tg, #tg
                         cases_no_pss) #pss


######################### Dynamical System ########################

#Create system with BasePower = 100 MVA and nominal frequency 60 Hz.
sys = PSY.System(100.0, frequency = 60.0);

#Add buses
for bus in nodes_case234
    PSY.add_component!(sys,bus)
end

#Add lines
for lines in branch_case234
    PSY.add_component!(sys,lines)
end

#Add loads
for loads in loads_case234
    PSY.add_component!(sys,loads)
end

#Add infinite source
PSY.add_component!(sys,inf_gen_case234)

#Add generators
PSY.add_component!(sys,case3_gen2)
PSY.add_component!(sys,case3_gen3)


##################################################
############### SOLVE PROBLEM ####################
##################################################



#Compute Y_bus after fault
sys2 = PSY.System(100.0, frequency = 60.0);
#Add buses
for bus in nodes_case234
    PSY.add_component!(sys2, bus)
end
#Add lines
for lines in branch_case234_fault
    PSY.add_component!(sys2, lines)
end
Ybus_fault = PSY.Ybus(sys2)[:,:]

#time span
tspan = (0.0, 20.0);

#Initial guess
x0_guess = [1.02, 1.0, 1.0, 0.0, -0.01, -0.01,
      1.0, #eq_p
      0.47, #ed_p
      0.95, #eq_pp
      0.8, #ed_pp
      0.6, #δ
      1.0, #ω
      2.1, #Vf
      0.28, #Vr1
      -0.39, #Vr2,
      1.0, #Vm
      0.81, #eq_p
      0.59, #ed_p
      0.75, #eq_pp
      0.6, #ed_pp
      0.86, #δ
      1.0, #ω
      1.7, #Vf
      0.11, #Vr1
      -0.31, #Vr2,
      1.0] #Vm

#Define Fault: Change of YBus
Ybus_change = ThreePhaseFault(1.0, #change at t = 1.0
                            Ybus_fault) #New YBus

#Define Simulation Problem
sim = Simulation(sys, #system
                 tspan, #time span
                 Ybus_change, #Type of Fault
                 initial_guess = x0_guess) #initial guess

#Solve problem in equilibrium
run_simulation!(sim, IDA());

#Obtain data for angles
series = get_state_series(sim, ("Case3Gen2", :δ));

@test sim.solution.retcode == :Success
