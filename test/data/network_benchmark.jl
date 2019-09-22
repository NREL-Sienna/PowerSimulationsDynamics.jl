using PowerSystems

############### OMIB Benchmark (Case 1 and 5) ########################

nodes_case1 = [Bus(1 , #number
                   "Bus 1", #Name
                   "REF" , #BusType (REF, PV, PQ)
                   0, #Angle in radians
                   1.05, #Voltage in pu
                   (min=0.94, max=1.06), #Voltage limits in pu
                   69), #Base voltage in kV
                   Bus(2 , "Bus 2"  , "PV" , 0 , 1.0 , (min=0.94, max=1.06), 69)]

branch_case1 = [Line("Line1", #name
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
branch_case1_fault = [Line("Line1", #name
                           true, #available
                           0.0, #active power flow initial condition (from-to)
                           0.0, #reactive power flow initial condition (from-to)
                           Arc(from=nodes_case1[1], to=nodes_case1[2]), #Connection between buses
                           0.02, #resistance in pu
                           0.1, #reactance in pu
                           (from=0.0, to=0.0), #susceptance in pu
                           18.046, #rate in MW
                           1.04)]  #angle limits (-min and max)

loads_case1 = [PowerLoad("Bus1", true, nodes_case1[2], PowerSystems.ConstantPower, 0.3, 0.01, 0.3, 0.01)]


################ Benchmark 3 Bus case (Cases 2, 3 and 4)########################

nodes_case234   = [Bus(1 , "Bus 1"  , "REF" , 0 , 1.02  , (min=0.94, max=1.06), 138),
                    Bus(2 , "Bus 2"  , "PV" , 0 , 1.00 , (min=0.94, max=1.06), 138),
                    Bus(3 , "Bus 3"  , "PV" , 0 , 1.00 , (min=0.94, max=1.06), 138)]

branch_case234  =  [Line("Line1", true, 0.0, 0.0, Arc(from=nodes_case234[1], to=nodes_case234[3]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04),
                    Line("Line2", true, 0.0, 0.0, Arc(from=nodes_case234[1], to=nodes_case234[2]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04),
                    Line("Line3", true, 0.0, 0.0, Arc(from=nodes_case234[2], to=nodes_case234[3]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04)]

#Trip of Line 1.
branch_case234_fault = [Line("Line2", true, 0.0, 0.0, Arc(from=nodes_case234[1], to=nodes_case234[2]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04),
                        Line("Line3", true, 0.0, 0.0, Arc(from=nodes_case234[2], to=nodes_case234[3]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04)]

loads_case234 =   [PowerLoad("Bus1", true, nodes_case234[1], PowerSystems.ConstantPower, 1.5, 0.8, 1.5, 0.8),
                   PowerLoad("Bus2", true, nodes_case234[2], PowerSystems.ConstantPower, 1.5, 0.7, 1.5, 0.8),
                   PowerLoad("Bus3", true, nodes_case234[3], PowerSystems.ConstantPower, 0.5, 0.3, 0.5, 0.3)]


################ Benchmark 3 Bus case TG (Case 5)########################

nodes_case5   =    [Bus(1 , "Bus 1"  , "REF" , 0 , 1.05  , (min=0.94, max=1.06), 138),
                    Bus(2 , "Bus 2"  , "PV" , 0 , 1.00 , (min=0.94, max=1.06), 138),
                    Bus(3 , "Bus 3"  , "PV" , 0 , 1.00 , (min=0.94, max=1.06), 138)]

branch_case5 =     [Line("Line1", true, 0.0, 0.0, Arc(from=nodes_case5[1], to=nodes_case5[2]), 0.01, 0.05, (from=0.0, to=0.0), 100, 1.04),
                    Line("Line2", true, 0.0, 0.0, Arc(from=nodes_case5[2], to=nodes_case5[3]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04),
                    Line("Line3", true, 0.0, 0.0, Arc(from=nodes_case5[1], to=nodes_case5[3]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04)]

#Trip of a single circuit of Line 1 -> Resistance and Reactance doubled.
branch_case5_fault =     [Line("Line1", true, 0.0, 0.0, Arc(from=nodes_case5[1], to=nodes_case5[2]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04),
                          Line("Line2", true, 0.0, 0.0, Arc(from=nodes_case5[2], to=nodes_case5[3]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04),
                          Line("Line3", true, 0.0, 0.0, Arc(from=nodes_case5[1], to=nodes_case5[3]), 0.02, 0.1, (from=0.0, to=0.0), 100, 1.04)]

loads_case5 = [PowerLoad("Bus2", true, nodes_case5[2], PowerSystems.ConstantPower, 0.3, 0.05, 0.3, 0.05),
               PowerLoad("Bus3", true, nodes_case5[3], PowerSystems.ConstantPower, 0.3, 0.05, 0.3, 0.05)]
