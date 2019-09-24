"""
Case 6:
This case study a 19-state virtual synchronous machine against an infinite bus located at bus 1, with VSM located at bus 2.
The perturbation increase the reference power (analogy for mechanical power) from 0.5 to 0.7.
"""

##################################################
############### LOAD DATA ########################
##################################################

############### Data Network ########################

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

############### Data devices ########################

inf_gen_DAIB = StaticSource(1, #number
                    :InfBus, #name
                    nodes_DAIB[1],#bus
                    1.00, #VR
                    0.0, #VI
                    0.000005) #Xth

############### Inverter Data ########################

converter = AvgCnvFixedDC(690.0, #Rated Voltage
                          2.75) #Rated MVA

dc_source = FixedDCSource(600.0) #Not in the original data, guessed.

filt = LCLFilter(0.08, #Series inductance lf in pu
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
                             filt) #Output Filter


######################### Dynamical System ########################

DAIB = DynamicSystem(nodes_DAIB, branch_DAIB, [Darco_Inverter], [inf_gen_DAIB], 100.0, 50.0);


##################################################
############### SOLVE PROBLEM ####################
##################################################

#Initialize variables
dx0 = zeros(LITS.get_total_rows(DAIB))
x0 = [1.00, #V1_R
          1.0648, #V2_R
          0.0, #V1_I
          0.001, #V2_I
          0.0, #δω_vsm
          0.2, #δθ_vsm
          0.025, #qm
          0.0015, #ξ_d
          -0.07, #ξ_q
          0.05, #γ_d
        -0.001, #γ_q
         0.95, #ϕ_d
         -0.10, #ϕ_q
         1.004, #vpll_d
         0.0, #vpll_q
         0.0, #ε_pll
         0.1, #δθ_pll
         0.5, #id_cv
         0.0, #iq_cv
         0.95, #vod
         -0.1, #voq
         0.49, #iod
        -0.1] #ioq
diff_vars = DAIB.DAE_vector
Darco_Inverter.inner_vars[13] = 0.95 #Vd_cnv var
Darco_Inverter.inner_vars[14] = -0.1 #Vq_cnv var
u0 = [0.2] #Increase in P_ref
tspan = (0.0, 4.0);

#Find initial condition
inif! = (out,x) -> system_model!(out, dx0 ,x, (u0,DAIB), 0.0)
sys_solve = nlsolve(inif!, x0, xtol=:1e-8,ftol=:1e-8,method=:trust_region)
x0_init = sys_solve.zero

#Define problem
prob = DiffEqBase.DAEProblem(system_model!, dx0, x0_init, tspan, (u0, DAIB), differential_vars = diff_vars)

#Solve problem in equilibrium
sol = solve(prob, IDA());

#Define data for using callbacks for defining the perturbation
tstop = [1.0]
cb = DiffEqBase.DiscreteCallback(LITS.change_t_one, LITS.step_change3!)

#Solve DAE system
sol2 = solve(prob, IDA(init_all=:false), callback=cb, tstops=tstop);

#Obtain data for virtual rotor angle speed
series = LITS.get_state_series(sol2, DAIB, (:DARCO, :δω_vsm))

@test sol2.retcode == :Success
