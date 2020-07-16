function get_init_values_for_comparison(sim::Simulation)
    bus_size = LITS.get_bus_count(sim.system)
    V_R = Vector{Float64}(undef, bus_size)
    V_I = Vector{Float64}(undef, bus_size)
    Vm = Vector{Float64}(undef, bus_size)
    θ = Vector{Float64}(undef, bus_size)
    for bus in PSY.get_components(PSY.Bus, sim.system)
        bus_n = PSY.get_number(bus)
        bus_ix = PSY.get_ext(sim.system)[LITS.LOOKUP][bus_n]
        V_R[bus_ix] = sim.x0_init[bus_ix]
        V_I[bus_ix] = sim.x0_init[bus_ix + bus_size]
        Vm[bus_ix] = sqrt(V_R[bus_ix]^2 + V_I[bus_ix]^2)
        θ[bus_ix] = angle(V_R[bus_ix] + V_I[bus_ix] * 1im)
    end
    results =
        Dict{String, Vector{Float64}}("V_R" => V_R, "V_I" => V_I, "Vm" => Vm, "θ" => θ)
    for device in PSY.get_components(PSY.DynamicInjection, sim.system)
        states = PSY.get_states(device)
        name = PSY.get_name(device)
        global_index = PSY.get_ext(sim.system)[LITS.GLOBAL_INDEX][name]
        x0_device = Vector{Float64}(undef, length(states))
        for (i, s) in enumerate(states)
            x0_device[i] = sim.x0_init[global_index[s]]
        end
        results[name] = x0_device
    end
    dyn_branches = PSY.get_components(PSY.DynamicBranch, sim.system)
    if !isempty(dyn_branches)
        for br in dyn_branches
            states = PSY.get_states(br)
            name = PSY.get_name(br)
            global_index = PSY.get_ext(sim.system)[LITS.GLOBAL_INDEX][name]
            x0_br = Vector{Float64}(undef, length(states))
            for (i, s) in enumerate(states)
                x0_br[i] = sim.x0_init[global_index[s]]
            end
            printed_name = "Line " * name
            results[printed_name] = x0_br
        end
    end
    return results
end
