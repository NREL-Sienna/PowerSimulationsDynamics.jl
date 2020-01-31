"""
Case 5:
This case study a three bus system with 1 machine located at bus 2.
The generator uses the model of a one d- one q- machine, and has a 5-mass shaft and a turbine governor.
The fault disconnects a circuit between buses 1 and 2, doubling its impedance.
"""

##################################################
############### LOAD DATA ########################
##################################################

############### Data Network ########################

nodes_case5 = [PSY.Bus(1 , "Bus 1"  , "REF" , 0 , 1.05  , (min=0.94, max=1.06), 138),
               PSY.Bus(2 , "Bus 2"  , "PV" , 0 , 1.00 , (min=0.94, max=1.06), 138),
               PSY.Bus(3 , "Bus 3"  , "PV" , 0 , 1.00 , (min=0.94, max=1.06), 138)]

branch_case5 = [PSY.Line("Line1", true, 0.0, 0.0, Arc(from=nodes_case5[1], to=nodes_case5[2]), 0.01, 0.05, (from=0.0, to=0.0), 100, 1.04),
                PSY.Line("Line2", true, 0.0, 0.0, Arc(from=nodes_case5[2], to=nodes_case5[3]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04),
                PSY.Line("Line3", true, 0.0, 0.0, Arc(from=nodes_case5[1], to=nodes_case5[3]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04)]

#Trip of a single circuit of Line 1 -> Resistance and Reactance doubled.
branch_case5_fault = [Line("Line1", true, 0.0, 0.0, Arc(from=nodes_case5[1], to=nodes_case5[2]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04),
                      Line("Line2", true, 0.0, 0.0, Arc(from=nodes_case5[2], to=nodes_case5[3]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04),
                      Line("Line3", true, 0.0, 0.0, Arc(from=nodes_case5[1], to=nodes_case5[3]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04)]

loads_case5 = [PSY.PowerLoad("Bus2", true, nodes_case5[2], PowerSystems.ConstantPower, 0.3, 0.05, 0.3, 0.05),
               PSY.PowerLoad("Bus3", true, nodes_case5[3], PowerSystems.ConstantPower, 0.3, 0.05, 0.3, 0.05)]


############### Data devices ########################

inf_gen_case5 = PSY.Source(
                "InfBus", #name
                true, #availability
                nodes_case5[1], #bus
                1.00, #VR
                0.0, #VI
                0.000001) #Xth

######## Machine Data #########

### Case 5: 4th Order Model with AVR + TG + Multishaft ###
case5_machine =  PSY.OneDOneQMachine(0.0, #R
                                  1.3125, #Xd
                                  1.2578, #Xq
                                  0.1813, #Xd_p
                                  0.25, #Xq_p
                                  5.89, #Td0_p
                                  0.6, #Tq0_p
                                  100.0)   #MVABase

######## Shaft Data #########

case5_shaft = PSY.FiveMassShaft(3.01, #5.148, #H
                            0.3348, #H_hp
                            0.7306, #H_ip
                            0.8154, #H_lp
                            0.0452, #H_ex,
                            0.0, #2.0, #D
                            0.5180, #D_hp
                            0.2240, #D_ip
                            0.2240, #D_lp
                            0.1450, #D_ex
                            0.0518, #D_12
                            0.0224, #D_23
                            0.0224, #D_34
                            0.0145, #D_45
                            33.07, #K_hp
                            28.59, #K_ip
                            44.68, #K_lp
                            21.984) #K_ex


######## PSS Data #########
cases_no_pss = PSY.PSSFixed(0.0)


######## TG Data #########

### TG for Case 5 ###
case5_tg = PSY.TGTypeII(0.05, #R
                    1.0, #T1
                    2.0, #T2
                    1.5, #τ_max
                    0.1) #τ_min

########  AVR Data #########
### AVRs for Case 2, 3, 4 and 5 ###
case5_avr = PSY.AVRTypeI(20.0, #Ka - Gain
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



### Case 5 Generator ###
case5_gen = PSY.DynamicGenerator(1, #Number
                         "Case5Gen",
                         nodes_case5[2], #bus
                         1.0, # ω_ref,
                         1.0155, #V_ref
                         0.5, #P_ref
                         0.0, #Q_ref
                         case5_machine, #machine
                         case5_shaft, #shaft
                         case5_avr, #avr
                         case5_tg, #tg
                         cases_no_pss) #pss


######################### Dynamical System ########################

#Create system with BasePower = 100 MVA and nominal frequency 60 Hz.
sys = PSY.System(100.0, frequency = 60.0);

#Add buses
for bus in nodes_case5
    PSY.add_component!(sys,bus)
end

#Add lines
for lines in branch_case5
    PSY.add_component!(sys,lines)
end

#Add loads
for loads in loads_case5
    PSY.add_component!(sys,loads)
end

#Add infinite source
PSY.add_component!(sys,inf_gen_case5)

#Add generators
PSY.add_component!(sys,case5_gen)


##################################################
############### SOLVE PROBLEM ####################
##################################################

#Compute Y_bus after fault
sys2 = PSY.System(100.0, frequency = 60.0);
#Add buses
for bus in nodes_case5
    PSY.add_component!(sys2, bus)
end
#Add lines
for lines in branch_case5_fault
    PSY.add_component!(sys2, lines)
end
Ybus_fault = PSY.Ybus(sys2)[:,:]

#time span
tspan = (0.0, 20.0);

#Initial guess
x0_guess = [1.02, 1.0, 1.0, 0.0, -0.01, -0.01,
      1.0, #eq_p
      0.0, #ed_p
      0.05, #δ
      1.0, #ω
      0.05, #δ_HP
      1.0, #ω_HP
      0.05, #δ_IP
      1.0, #ω_LP
      0.05, #δ_LP
      1.0, #ω_LP
      0.05, #δ_Ex
      1.0, #ω_Ex
      1.0, #Vf
      0.01, #Vr1
      -0.1, #Vr2,
      1.0, #Vm
      0.0] #xg

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
series = get_state_series(sim, ("Case5Gen", :δ));

@test sim.solution.retcode == :Success
