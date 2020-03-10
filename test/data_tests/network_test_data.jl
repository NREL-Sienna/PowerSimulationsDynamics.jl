using PowerSystems
const PSY = PowerSystems

############### Buses Data ########################

nodes_OMIB() = [
    PSY.Bus(
        1, #number
        "Bus 1", #Name
        "REF", #BusType (REF, PV, PQ)
        0, #Angle in radians
        1.05, #Voltage in pu
        (min = 0.94, max = 1.06), #Voltage limits in pu
        69, #Base voltage in kV
        nothing,
        nothing,
    ),
    PSY.Bus(2, "Bus 2", "PV", 0, 1.0, (min = 0.94, max = 1.06), 69, nothing, nothing),
]

nodes_3bus() = [
    PSY.Bus(1, "Bus 1", "REF", 0, 1.02, (min = 0.94, max = 1.06), 138, nothing, nothing),
    PSY.Bus(2, "Bus 2", "PV", 0, 1.00, (min = 0.94, max = 1.06), 138, nothing, nothing),
    PSY.Bus(3, "Bus 3", "PV", 0, 1.00, (min = 0.94, max = 1.06), 138, nothing, nothing),
]

nodes_3bus_case5() = [
    PSY.Bus(1, "Bus 1", "REF", 0, 1.05, (min = 0.94, max = 1.06), 138, nothing, nothing)
    PSY.Bus(2, "Bus 2", "PV", 0, 1.00, (min = 0.94, max = 1.06), 138, nothing, nothing)
    PSY.Bus(3, "Bus 3", "PV", 0, 1.00, (min = 0.94, max = 1.06), 138, nothing, nothing)
]

nodes_DArco_IB() = [
    PSY.Bus(
        1, #number
        "Bus 1", #Name
        "REF", #BusType (REF, PV, PQ)
        0, #Angle in radians
        1.04, #Voltage in pu
        (min = 0.94, max = 1.06), #Voltage limits in pu
        0.69, #Base voltage in kV
        nothing,
        nothing,
    ),
    PSY.Bus(2, "Bus 2", "PV", 0, 1.0, (min = 0.94, max = 1.06), 0.69, nothing, nothing),
]

nodes_multimachine() = [
    PSY.Bus(1, "Bus 1", "REF", 0, 1.05, (min = 0.94, max = 1.06), 138, nothing, nothing),
    PSY.Bus(2, "Bus 2", "PQ", 0, 1.00, (min = 0.94, max = 1.06), 138, nothing, nothing),
    PSY.Bus(3, "Bus 3", "PV", 0, 1.00, (min = 0.94, max = 1.06), 138, nothing, nothing),
]

############### Branches Data ########################

branches_OMIB(nodes_OMIB) = [PSY.Line(
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

branches_OMIB_fault(nodes_OMIB) = [PSY.Line(
    "Line1", #name
    true, #available
    0.0, #active power flow initial condition (from-to)
    0.0, #reactive power flow initial condition (from-to)
    Arc(from = nodes_OMIB[1], to = nodes_OMIB[2]), #Connection between buses
    0.02, #resistance in pu
    0.1, #reactance in pu
    (from = 0.0, to = 0.0), #susceptance in pu
    18.046, #rate in MW
    1.04,
)]  #angle limits (-min and max)

branches_3lines(nodes_3bus) = [
    PSY.Line(
        "Line1",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus[1], to = nodes_3bus[3]),
        0.01,
        0.12,
        (from = 0.0, to = 0.0),
        100,
        1.04,
    ),
    PSY.Line(
        "Line2",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus[1], to = nodes_3bus[2]),
        0.01,
        0.12,
        (from = 0.0, to = 0.0),
        100,
        1.04,
    ),
    PSY.Line(
        "Line3",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus[2], to = nodes_3bus[3]),
        0.01,
        0.12,
        (from = 0.0, to = 0.0),
        100,
        1.04,
    ),
]

branches_3lines_fault(nodes_3bus) = [
    PSY.Line(
        "Line2",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus[1], to = nodes_3bus[2]),
        0.01,
        0.12,
        (from = 0.0, to = 0.0),
        100,
        1.04,
    ),
    PSY.Line(
        "Line3",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus[2], to = nodes_3bus[3]),
        0.01,
        0.12,
        (from = 0.0, to = 0.0),
        100,
        1.04,
    ),
]

branches_3lines_case5(nodes_3bus_case5) = [
    PSY.Line(
        "Line1",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus_case5[1], to = nodes_3bus_case5[2]),
        0.01,
        0.05,
        (from = 0.0, to = 0.0),
        100,
        1.04,
    ),
    PSY.Line(
        "Line2",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus_case5[2], to = nodes_3bus_case5[3]),
        0.02,
        0.1,
        (from = 0.0, to = 0.0),
        100,
        1.04,
    ),
    PSY.Line(
        "Line3",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus_case5[1], to = nodes_3bus_case5[3]),
        0.02,
        0.1,
        (from = 0.0, to = 0.0),
        100,
        1.04,
    ),
]

branches_3lines_case5_fault(nodes_3bus_case5) = [
    Line(
        "Line1",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus_case5[1], to = nodes_3bus_case5[2]),
        0.02,
        0.1,
        (from = 0.0, to = 0.0),
        100,
        1.04,
    ),
    Line(
        "Line2",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus_case5[2], to = nodes_3bus_case5[3]),
        0.02,
        0.1,
        (from = 0.0, to = 0.0),
        100,
        1.04,
    ),
    Line(
        "Line3",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus_case5[1], to = nodes_3bus_case5[3]),
        0.02,
        0.1,
        (from = 0.0, to = 0.0),
        100,
        1.04,
    ),
]

branches_DArco_IB(nodes_DArco_IB) = [PSY.Line(
    "Line1", #name
    true, #available
    0.0, #active power flow initial condition (from-to)
    0.0, #reactive power flow initial condition (from-to)
    Arc(from = nodes_DArco_IB[1], to = nodes_DArco_IB[2]), #Connection between buses
    0.0, #resistance in pu
    0.075, #reactance in pu
    (from = 0.0, to = 0.0), #susceptance in pu
    5.0, #rate in MW
    1.04,
)]  #angle limits (-min and max)

branches_3lines_case8(nodes_3bus) = [
    PSY.Line(
        "Line1",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus[1], to = nodes_3bus[3]),
        0.01,
        0.12,
        (from = 0.1, to = 0.1),
        100,
        1.04,
    ),
    PSY.Line(
        "Line2",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus[1], to = nodes_3bus[2]),
        0.01,
        0.12,
        (from = 0.1, to = 0.1),
        100,
        1.04,
    ),
    PSY.Line(
        "Line3",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus[2], to = nodes_3bus[3]),
        0.02,
        0.9,
        (from = 0.5, to = 0.5),
        100,
        1.04,
    ),
]

branches_3lines_case8_fault(nodes_3bus) = [
    PSY.Line(
        "Line1",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus[1], to = nodes_3bus[3]),
        0.03,
        0.36,
        (from = 0.03, to = 0.03),
        100,
        1.04,
    ),
    PSY.Line(
        "Line2",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus[1], to = nodes_3bus[2]),
        0.01,
        0.12,
        (from = 0.1, to = 0.1),
        100,
        1.04,
    ),
    PSY.Line(
        "Line3",
        true,
        0.0,
        0.0,
        Arc(from = nodes_3bus[2], to = nodes_3bus[3]),
        0.02,
        0.9,
        (from = 0.5, to = 0.5),
        100,
        1.04,
    ),
]

branch_3bus_case9(nodes) = [
    PSY.Line(
        "Line1",
        true,
        0.0,
        0.0,
        Arc(from = nodes[1], to = nodes[3]),
        0.01,
        0.12,
        (from = 0.1, to = 0.1),
        100,
        1.04,
    ),
    PSY.Line(
        "Line2",
        true,
        0.0,
        0.0,
        Arc(from = nodes[1], to = nodes[2]),
        0.01,
        0.12,
        (from = 0.1, to = 0.1),
        100,
        1.04,
    ),
    PSY.Line(
        "Line3",
        true,
        0.0,
        0.0,
        Arc(from = nodes[2], to = nodes[3]),
        0.02,
        0.9,
        (from = 0.5, to = 0.5),
        100,
        1.04,
    ),
]

branch_3bus_case9_fault(nodes) = [
    PSY.Line(
        "Line1",
        true,
        0.0,
        0.0,
        Arc(from = nodes[1], to = nodes[3]),
        0.03,
        0.36,
        (from = 0.03, to = 0.03),
        100,
        1.04,
    ),
    PSY.Line(
        "Line2",
        true,
        0.0,
        0.0,
        Arc(from = nodes[1], to = nodes[2]),
        0.01,
        0.12,
        (from = 0.1, to = 0.1),
        100,
        1.04,
    ),
    PSY.Line(
        "Line3",
        true,
        0.0,
        0.0,
        Arc(from = nodes[2], to = nodes[3]),
        0.02,
        0.9,
        (from = 0.5, to = 0.5),
        100,
        1.04,
    ),
]

branch_multimachine(nodes) = [
    Line(
        "Line1", #name
        true, #available
        0.0, #active power flow initial condition (from-to)
        0.0, #reactive power flow initial condition (from-to)
        Arc(from = nodes[1], to = nodes[2]), #Connection between buses
        0.01, #resistance in pu
        0.05, #reactance in pu
        (from = 0.0, to = 0.0), #susceptance in pu
        18.046, #rate in MW
        1.04,
    ),
    Line(
        "Line2", #name
        true, #available
        0.0, #active power flow initial condition (from-to)
        0.0, #reactive power flow initial condition (from-to)
        Arc(from = nodes[2], to = nodes[3]), #Connection between buses
        0.01, #resistance in pu
        0.05, #reactance in pu
        (from = 0.0, to = 0.0), #susceptance in pu
        18.046, #rate in MW
        1.04,
    ),  #angle limits (-min and max)
]

branch_multimachine_fault(nodes) = [
    Line(
        "Line1", #name
        true, #available
        0.0, #active power flow initial condition (from-to)
        0.0, #reactive power flow initial condition (from-to)
        Arc(from = nodes[1], to = nodes[2]), #Connection between buses
        0.01, #resistance in pu
        0.05, #reactance in pu
        (from = 0.0, to = 0.0), #susceptance in pu
        18.046, #rate in MW
        1.04,
    ),
    Line(
        "Line2", #name
        true, #available
        0.0, #active power flow initial condition (from-to)
        0.0, #reactive power flow initial condition (from-to)
        Arc(from = nodes[2], to = nodes[3]), #Connection between buses
        0.02, #resistance in pu
        0.10, #reactance in pu
        (from = 0.0, to = 0.0), #susceptance in pu
        18.046, #rate in MW
        1.04,
    ),  #angle limits (-min and max)
]

############### Load Data ########################

loads_OMIB(nodes_OMIB) = [PSY.PowerLoad(
    "LBus1", #name
    true, #availability
    nodes_OMIB[2], #bus
    PSY.LoadModels.ConstantPower, #type
    0.3, #P
    0.01, #Q
    0.3, #P_max
    0.01,
)] #Q_max

loads_3bus(nodes_3bus) = [
    PSY.PowerLoad(
        "Bus1",
        true,
        nodes_3bus[1],
        PSY.LoadModels.ConstantPower,
        1.5,
        0.8,
        1.5,
        0.8,
    ),
    PSY.PowerLoad(
        "Bus2",
        true,
        nodes_3bus[2],
        PSY.LoadModels.ConstantPower,
        1.5,
        0.7,
        1.5,
        0.8,
    ),
    PSY.PowerLoad(
        "Bus3",
        true,
        nodes_3bus[3],
        PSY.LoadModels.ConstantPower,
        0.5,
        0.3,
        0.5,
        0.3,
    ),
]

loads_3bus_case5(nodes_3bus) = [
    PSY.PowerLoad(
        "Bus2",
        true,
        nodes_3bus[2],
        PSY.LoadModels.ConstantPower,
        0.3,
        0.05,
        0.3,
        0.05,
    ),
    PSY.PowerLoad(
        "Bus3",
        true,
        nodes_3bus[3],
        PSY.LoadModels.ConstantPower,
        0.3,
        0.05,
        0.3,
        0.05,
    ),
]

loads_3bus_case7(nodes_3bus) = [
    PSY.PowerLoad(
        "Bus1",
        true,
        nodes_3bus[1],
        PSY.LoadModels.ConstantPower,
        0.5,
        0.1,
        1.5,
        0.8,
    ),
    PSY.PowerLoad(
        "Bus2",
        true,
        nodes_3bus[2],
        PSY.LoadModels.ConstantPower,
        1.0,
        0.3,
        1.5,
        0.8,
    ),
    PSY.PowerLoad(
        "Bus3",
        true,
        nodes_3bus[3],
        PSY.LoadModels.ConstantPower,
        0.3,
        0.1,
        0.5,
        0.3,
    ),
]

loads_3bus_case8(nodes_3bus) = [
    PSY.PowerLoad(
        "Bus1",
        true,
        nodes_3bus[1],
        PSY.LoadModels.ConstantPower,
        0.5,
        0.1,
        1.5,
        0.8,
    ),
    PSY.PowerLoad(
        "Bus2",
        true,
        nodes_3bus[2],
        PSY.LoadModels.ConstantPower,
        1.0,
        0.3,
        1.5,
        0.8,
    ),
    PSY.PowerLoad(
        "Bus3",
        true,
        nodes_3bus[3],
        PSY.LoadModels.ConstantPower,
        0.3,
        0.1,
        0.5,
        0.3,
    ),
]

loads_3bus_case9(nodes) = [
    PowerLoad("Bus1", true, nodes[1], PSY.LoadModels.ConstantPower, 0.5, 0.1, 1.5, 0.8),
    PowerLoad("Bus2", true, nodes[2], PSY.LoadModels.ConstantPower, 1.0, 0.3, 1.5, 0.8),
    PowerLoad("Bus3", true, nodes[3], PSY.LoadModels.ConstantPower, 0.3, 0.1, 0.5, 0.3),
]

loads_multimachine(nodes) = [PowerLoad(
    "LBus1", #name
    true, #availability
    nodes[2], #bus
    LoadModels.ConstantPower, #type
    0.8, #P
    0.01, #Q
    0.8, #P_max
    0.01, #Q_max
)]

############### Infinite Sources Data ########################

inf_gen_1_pu(nodes) = PSY.Source(
    "InfBus", #name
    true, #availability
    nodes[1], #bus
    1.00, #VR
    0.0, #VI
    0.000005,
) #Xth

inf_gen_102_pu(nodes) = PSY.Source(
    "InfBus", #name
    true, #availability
    nodes[1], #bus
    1.02, #VR
    0.0, #VI
    0.000001,
) #Xth

inf_gen_105_pu(nodes) = PSY.Source(
    "InfBus", #name
    true, #availability
    nodes[1], #bus
    1.05, #VR
    0.0, #VI
    0.000001,
) #Xth

############# Obtain Ybus ###########

function get_admittance_matrix(nodes, branches)
    sys2 = PSY.System(100.0, frequency = 60.0)
    for bus in nodes
        PSY.add_component!(sys2, bus)
    end

    #Add lines
    for lines in branches
        PSY.add_component!(sys2, lines)
    end
    return PSY.Ybus(sys2)[:, :]
end

function get_admittance_matrix_case9(nodes, branches)
    sys2 = PSY.System(100.0, frequency = 60.0)
    for bus in nodes
        PSY.add_component!(sys2, bus)
    end

    #Add lines
    for lines in branches
        PSY.add_component!(sys2, lines)
    end

    #Make line 3 the dynamic line
    make_dynamic_branch!(branches[3], sys2)

    return PSY.Ybus(sys2)[:, :]
end
