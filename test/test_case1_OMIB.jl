"""
Case 1:
This case study defines a classical machine against an infinite bus. The fault
drop a circuit on the (double circuit) line connecting the two buses, doubling its impedance
"""


##################################################
############### LOAD DATA ########################
##################################################

############### Data Network ########################

nodes_case1 = [PSY.Bus(1 , #number
                   "Bus 1", #Name
                   "REF" , #BusType (REF, PV, PQ)
                   0, #Angle in radians
                   1.05, #Voltage in pu
                   (min=0.94, max=1.06), #Voltage limits in pu
                   69), #Base voltage in kV
                   PSY.Bus(2 , "Bus 2"  , "PV" , 0 , 1.0 , (min=0.94, max=1.06), 69)]

branch_case1 = [PSY.Line("Line1", #name
                     true, #available
                     0.0, #active power flow initial condition (from-to)
                     0.0, #reactive power flow initial condition (from-to)
                     Arc(from=nodes_case1[1], to=nodes_case1[2]), #Connection between buses
                     0.01, #resistance in pu
                     0.05, #reactance in pu
                     (from=0.0, to=0.0), #susceptance in pu
                     18.046, #rate in MW
                     1.04)]  #angle limits (-min and max)

#Trip of a single circuit of Line 1 -> Resistance and Reactance doubled.
branch_case1_fault = [PSY.Line("Line1", #name
                           true, #available
                           0.0, #active power flow initial condition (from-to)
                           0.0, #reactive power flow initial condition (from-to)
                           Arc(from=nodes_case1[1], to=nodes_case1[2]), #Connection between buses
                           0.02, #resistance in pu
                           0.1, #reactance in pu
                           (from=0.0, to=0.0), #susceptance in pu
                           18.046, #rate in MW
                           1.04)]  #angle limits (-min and max)

loads_case1 = [PSY.PowerLoad("LBus1", #name
                            true, #availability
                            nodes_case1[2], #bus
                            PowerSystems.ConstantPower, #type
                            0.3, #P
                            0.01, #Q
                            0.3, #P_max
                            0.01)] #Q_max


############### Data devices ########################

inf_gen_case1 = PSY.Source(
                "InfBus", #name
                true, #availability
                nodes_case1[1], #bus
                1.05, #VR
                0.0, #VI
                0.000001) #Xth

######## Machine Data #########

### Case 1: Classical machine against infinite bus ###
case1_machine = PSY.BaseMachine(0.0, #R
                            0.2995, #Xd_p
                            0.7087, #eq_p
                            100.0)  #MVABase

######## Shaft Data #########

### Shaft for Case 1 ###
case1_shaft = PSY.SingleMass(3.148, #H
                     2.0) #D



######## PSS Data #########
cases_no_pss = PSY.PSSFixed(0.0)


######## TG Data #########

### No TG for Cases 1, 2, 3, 4 ###
case1234_no_tg = PSY.TGFixed(1.0) #eff

########  AVR Data #########
case1_avr = PSY.AVRFixed(0.0) #Vf not applicable in Classic Machines

### Case 1 Generator ###
case1_gen = PSY.DynamicGenerator(1, #Number
                      "Case1Gen",
                      nodes_case1[2], #bus
                      1.0, # ω_ref,
                      1.0, #V_ref
                      0.5, #P_ref
                      0.0, #Q_ref
                      case1_machine, #machine
                      case1_shaft, #shaft
                      case1_avr, #avr
                      case1234_no_tg, #tg
                      cases_no_pss) #pss


######################### Dynamical System ########################

#Create system with BasePower = 100 MVA and nominal frequency 60 Hz.
sys = PSY.System(100.0, frequency = 60.0);

#Add buses
for bus in nodes_case1
    PSY.add_component!(sys,bus)
end

#Add lines
for lines in branch_case1
    PSY.add_component!(sys,lines)
end

#Add loads
for loads in loads_case1
    PSY.add_component!(sys,loads)
end

#Add infinite source
PSY.add_component!(sys,inf_gen_case1)

#Add generator
PSY.add_component!(sys,case1_gen)


##################################################
############### SOLVE PROBLEM ####################
##################################################

#Compute Y_bus after fault
sys2 = PSY.System(100.0, frequency = 60.0);
#Add buses
for bus in nodes_case1
    PSY.add_component!(sys2, bus)
end

#Add lines
for lines in branch_case1_fault
    PSY.add_component!(sys2, lines)
end
Ybus_fault = PSY.Ybus(sys2)[:,:]

tspan = (0.0, 30.0);

#Define Fault: Change of YBus
Ybus_change = ThreePhaseFault(1.0, #change at t = 1.0
                            Ybus_fault) #New YBus

#Define Simulation Problem
sim = Simulation(sys, #system
                 tspan, #time span
                 Ybus_change) #Type of Fault

#Solve problem in equilibrium
run_simulation!(sim, IDA(), dtmax=0.02);

#Obtain data for angles
series = get_state_series(sim, ("Case1Gen", :δ));
series2 = get_voltagemag_series(sim, 2)
LITS.print_init_states(sim)

@test sim.solution.retcode == :Success
