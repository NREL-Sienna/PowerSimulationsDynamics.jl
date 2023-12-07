function add_source_to_ref(sys::PSY.System, X_th::Float64)
    for g in PSY.get_components(StaticInjection, sys)
        isa(g, ElectricLoad) && continue
        g.bus.bustype == ACBusTypes.REF &&
            error("A device is already attached to the REF bus")
    end

    slack_bus =
        [b for b in PSY.get_components(ACBus, sys) if b.bustype == ACBusTypes.REF][1]
    inf_source = Source(;
        name = "InfBus", #name
        available = true, #availability
        active_power = 0.0,
        reactive_power = 0.0,
        bus = slack_bus, #bus
        R_th = 0.0,
        X_th = X_th, #Xth
    )
    PSY.add_component!(sys, inf_source)
    return
end

function add_source_to_ref(sys::PSY.System)
    for g in PSY.get_components(StaticInjection, sys)
        isa(g, ElectricLoad) && continue
        g.bus.bustype == ACBusTypes.REF &&
            error("A device is already attached to the REF bus")
    end

    slack_bus =
        [b for b in PSY.get_components(ACBus, sys) if b.bustype == ACBusTypes.REF][1]
    inf_source = Source(;
        name = "InfBus", #name
        available = true, #availability
        active_power = 0.0,
        reactive_power = 0.0,
        bus = slack_bus, #bus
        R_th = 0.0,
        X_th = 5e-6, #Xth
    )
    PSY.add_component!(sys, inf_source)
    return
end

function get_ybus_fault_threebus_sys(sys)
    fault_branch =
        filter!(x -> get_name(x) != "BUS 1-BUS 3-i_1", collect(get_components(Branch, sys)))
    sorted_buses =
        sort!(collect(get_components(ACBus, threebus_sys)); by = x -> get_number(x))
    Ybus_fault = PNM.Ybus(fault_branch, sorted_buses)[:, :]
    return Ybus_fault
end
