function add_source_to_ref(sys::PSY.System)
    for g in PSY.get_components(StaticInjection, sys)
        isa(g, ElectricLoad) && continue
        g.bus.bustype == BusTypes.REF && error("A device is already attached to the REF bus")
    end

    slack_bus = [b for b in PSY.get_components(Bus, sys) if b.bustype == BusTypes.REF][1]
    inf_source = Source(
        name = "InfBus", #name
        available = true, #availability
        activepower = 0.0,
        reactivepower = 0.0,
        bus = slack_bus, #bus
        X_th = 0.000005
    ) #Xth
    PSY.add_component!(sys, inf_source)
    return
end
