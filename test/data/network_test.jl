nodes_test= [Bus(1 , #number
                 "Bus 1", #Name
                 "REF" , #BusType (REF, PV, PQ)
                 0, #Angle in radians
                 1.02, #Voltage in pu
                 (min=0.94, max=1.06), #Voltage limits in pu
                 69), #Base voltage in kV
             Bus(2 , "Bus 2"  , "PV" , 0 , 1.0 , (min=0.94, max=1.06), 69)]

nodes_test2= [Bus(1 , #number
              "Bus 1", #Name
              "REF" , #BusType (REF, PV, PQ)
              0, #Angle in radians
              1.02, #Voltage in pu
              (min=0.94, max=1.06), #Voltage limits in pu
              69), #Base voltage in kV
          Bus(2 , "Bus 2"  , "PV" , 0 , 1.0 , (min=0.94, max=1.06), 69),
          Bus(3 , "Bus 3"  , "PQ" , 0 , 1.0 , (min=0.94, max=1.06), 69)]

branch_test = [Line("Line1", #name
                    true, #available
                    0.0, #active power flow initial condition (from-to)
                    0.0, #reactive power flow initial condition (from-to)
                    Arc(from=nodes_test[1], to=nodes_test[2]), #Connection between buses
                    0.01, #resistance in pu
                    0.05, #reactance in pu
                    (from=0.0, to=0.0), #susceptance in pu
                    18.046, #rate in MW
                    1.04)]  #angle limits (-min and max)




branch_test_fault = [Line("Line1", #name
                    true, #available
                    0.0, #active power flow initial condition (from-to)
                    0.0, #reactive power flow initial condition (from-to)
                    Arc(from=nodes_test[1], to=nodes_test[2]), #Connection between buses
                    0.02, #resistance in pu
                    0.1, #reactance in pu
                    (from=0.0, to=0.0), #susceptance in pu
                    18.046, #rate in MW
                    1.04)]  #angle limits (-min and max)


####

branch_test2 = [Line("Line1", #name
                    true, #available
                    0.0, #active power flow initial condition (from-to)
                    0.0, #reactive power flow initial condition (from-to)
                    Arc(from=nodes_test2[1], to=nodes_test2[2]), #Connection between buses
                    0.01, #resistance in pu
                    0.05, #reactance in pu
                    (from=0.0, to=0.0), #susceptance in pu
                    18.046, #rate in MW
                    1.04),
                    Line("Line2", #name
                            true, #available
                            0.0, #active power flow initial condition (from-to)
                            0.0, #reactive power flow initial condition (from-to)
                            Arc(from=nodes_test2[2], to=nodes_test2[3]), #Connection between buses
                            0.02, #resistance in pu
                            0.1, #reactance in pu
                            (from=0.0, to=0.0), #susceptance in pu
                            18.046, #rate in MW
                            1.04)]  #angle limits (-min and max)

branch_test2_fault = [Line("Line1", #name
                true, #available
                0.0, #active power flow initial condition (from-to)
                0.0, #reactive power flow initial condition (from-to)
                Arc(from=nodes_test2[1], to=nodes_test2[2]), #Connection between buses
                0.02, #resistance in pu
                0.1, #reactance in pu
                (from=0.0, to=0.0), #susceptance in pu
                18.046, #rate in MW
                1.04),
                Line("Line2", #name
                        true, #available
                        0.0, #active power flow initial condition (from-to)
                        0.0, #reactive power flow initial condition (from-to)
                        Arc(from=nodes_test2[2], to=nodes_test2[3]), #Connection between buses
                        0.02, #resistance in pu
                        0.1, #reactance in pu
                        (from=0.0, to=0.0), #susceptance in pu
                        18.046, #rate in MW
                        1.04)]  #angle limits (-min and max)
