using PowerSystems


################## OMIB data #####################

nodes_OMIB= [Bus(1 , #number
                 "Bus 1", #Name
                 "REF" , #BusType (REF, PV, PQ)
                 0, #Angle in radians
                 1.06, #Voltage in pu
                 (min=0.94, max=1.06), #Voltage limits in pu
                 69), #Base voltage in kV
             Bus(2 , "Bus 2"  , "PV" , 0 , 1.045 , (min=0.94, max=1.06), 69)]

branch_OMIB = [Line("Line1", #name
                    true, #available
                    0.0, #active power flow initial condition (from-to)
                    0.0, #reactive power flow initial condition (from-to)
                    Arc(from=nodes_OMIB[1], to=nodes_OMIB[2]), #Connection between buses
                    0.01938, #resistance in pu
                    0.05917, #reactance in pu
                    (from=0.0264, to=0.0264), #susceptance in pu
                    18.046, #rate in MW
                    1.04)]  #angle limits (-min and max)

Dynbranch_OMIB = [DynLine("Line1", #name
                    true, #available
                    Arc(from=nodes_OMIB[1], to=nodes_OMIB[2]), #Connection between buses
                    0.01938, #resistance in pu
                    0.05917, #reactance in pu
                    (from=0.0264, to=0.0264))]  #susceptance in pu

loads = [PowerLoad("Bus2", #name
                   true, #available
                   nodes_OMIB[1], #bus
                   PowerSystems.ConstantPower, #model
                   0.217, #active power in pu
                   0.127, #reactive power in pu
                   0.217, #max active power in pu
                   0.127)] #max reactive power in pu

################## Kundur machine case data #####################

nodes_kundur= [Bus(1 , "Bus 1"  , "PQ" , 0 , 1.02  , (min=0.94, max=1.06), 24),
                    Bus(2 , "Bus 2"  , "REF" , 0 , 1.01 , (min=0.94, max=1.06), 24)]

branch_kundur = [Line("Line1", true, 0.0, 0.0, Arc(from=nodes_kundur[1], to=nodes_kundur[2]),
                      0.01938, 0.05917, (from=0.0, to=0.0), 18.046, 1.04)]

loads_kundur = [PowerLoad("Bus1", true, nodes_kundur[1], PowerSystems.ConstantPower, 0.0, 0.0, 0.0, 0.0)]

################ Benchmark 3 Bus case ########################

nodes_benchmark = [Bus(1 , "Bus 1"  , "REF" , 0 , 1.02  , (min=0.94, max=1.06), 138),
                    Bus(2 , "Bus 2"  , "PV" , 0 , 1.00 , (min=0.94, max=1.06), 138),
                    Bus(3 , "Bus 3"  , "PV" , 0 , 1.00 , (min=0.94, max=1.06), 138)]

branch_benchmark = [Line("Line1", true, 0.0, 0.0, Arc(from=nodes_benchmark[1], to=nodes_benchmark[3]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04),
                    Line("Line2", true, 0.0, 0.0, Arc(from=nodes_benchmark[1], to=nodes_benchmark[2]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04),
                    Line("Line3", true, 0.0, 0.0, Arc(from=nodes_benchmark[2], to=nodes_benchmark[3]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04)]


branch_benchmark_noline = [Line("Line2", true, 0.0, 0.0, Arc(from=nodes_benchmark[1], to=nodes_benchmark[2]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04),
                           Line("Line3", true, 0.0, 0.0, Arc(from=nodes_benchmark[2], to=nodes_benchmark[3]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04)]

loads_benchmark = [PowerLoad("Bus1", true, nodes_benchmark[1], PowerSystems.ConstantPower, 1.5, 0.8, 1.5, 0.8),
                   PowerLoad("Bus2", true, nodes_benchmark[2], PowerSystems.ConstantPower, 1.5, 0.7, 1.5, 0.8),
                   PowerLoad("Bus3", true, nodes_benchmark[3], PowerSystems.ConstantPower, 0.5, 0.3, 0.5, 0.3)]

################ Benchmark 2 Bus case ########################

nodes_bench = [Bus(1 , "Bus 1"  , "REF" , 0 , 1.02  , (min=0.94, max=1.06), 138),
                    Bus(2 , "Bus 2"  , "PV" , 0 , 1.00 , (min=0.94, max=1.06), 138)]

branch_bench = [Line("Line1", true, 0.0, 0.0, Arc(from=nodes_bench[1], to=nodes_bench[2]), 0.01, 0.12, (from=0.0, to=0.0), 100, 1.04)]

loads_bench = [PowerLoad("Bus1", true, nodes_bench[2], PowerSystems.ConstantPower, 1.0, 0.3, 1.5, 0.8)]

############## Static Kundur Machines for PF #######################

machine_Kundur = ThermalStandard("Bus1", true, nodes_kundur[1], 2.22, 0,
                    TechThermal(5.55, PowerSystems.ST, PowerSystems.COAL, (min=0.0, max=5.55), (min=-1.0, max=1.0), nothing, nothing),
                    ThreePartCost((430.292599, 2000.0), 0.0, 0.0, 0.0)
                    )
infinite_bus = ThermalStandard("Bus2", true, nodes_kundur[2], -2.22, 0,
                    TechThermal(5.55, PowerSystems.ST, PowerSystems.COAL, (min=-500.55, max=500.55), (min=-5.55, max=5.55), nothing, nothing),
                    ThreePartCost((430.292599, 2000.0), 0.0, 0.0, 0.0)
                    )
