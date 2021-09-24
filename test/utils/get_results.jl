function get_init_values_for_comparison(sim::Simulation)
    bus_size = PSID.get_bus_count(sim.inputs)
    system = PSID.get_system(sim)
    V_R = Vector{Float64}(undef, bus_size)
    V_I = Vector{Float64}(undef, bus_size)
    Vm = Vector{Float64}(undef, bus_size)
    θ = Vector{Float64}(undef, bus_size)
    for bus in PSY.get_components(PSY.Bus, system)
        bus_n = PSY.get_number(bus)
        bus_ix = PSID.get_lookup(sim.inputs)[bus_n]
        V_R[bus_ix] = sim.x0_init[bus_ix]
        V_I[bus_ix] = sim.x0_init[bus_ix + bus_size]
        Vm[bus_ix] = sqrt(V_R[bus_ix]^2 + V_I[bus_ix]^2)
        θ[bus_ix] = angle(V_R[bus_ix] + V_I[bus_ix] * 1im)
    end
    results =
        Dict{String, Vector{Float64}}("V_R" => V_R, "V_I" => V_I, "Vm" => Vm, "θ" => θ)
    for device in PSID.get_dynamic_injectors(sim.inputs)
        states = PSY.get_states(device)
        name = PSY.get_name(device)
        global_index = PSID.get_global_index(device)
        x0_device = Vector{Float64}(undef, length(states))
        for (i, s) in enumerate(states)
            x0_device[i] = sim.x0_init[global_index[s]]
        end
        results[name] = x0_device
    end
    for br in PSID.get_dynamic_branches(sim.inputs)
        states = PSY.get_states(br)
        name = PSY.get_name(br)
        global_index = PSID.get_global_index(br)
        x0_br = Vector{Float64}(undef, length(states))
        for (i, s) in enumerate(states)
            x0_br[i] = sim.x0_init[global_index[s]]
        end
        printed_name = "Line " * name
        results[printed_name] = x0_br
    end

    return results
end

function clean_extra_timestep!(t::Vector{Float64}, δ::Vector{Float64})
    idx = unique(i -> t[i], 1:length(t))
    return t[idx], δ[idx]
end

function get_csv_delta(str::AbstractString)
    M = readdlm(str, ',')
    return clean_extra_timestep!(M[:, 1], M[:, 2])
end

function get_csv_data(str::AbstractString)
    M = readdlm(str, ',')
    return M
end
