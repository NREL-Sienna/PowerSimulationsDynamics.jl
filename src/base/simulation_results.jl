struct SimulationResults
    global_index::MAPPING_DICT
    voltage_buses::Vector{Int}
    current_buses::Vector{Int}
    solution::SciMLBase.AbstractODESolution
end

function make_global_state(wrapped_devices)
    dict = MAPPING_DICT()
    global_state_index[PSY.get_name(device)] = Dict{Symbol, Int}()
    for s in PSY.get_states(device)
        state_space_ix[1] += 1
        global_state_index[PSY.get_name(device)][s] = state_space_ix[1]
    end
    #This V_ix should be V_number.
    # global_state_index["V_$(ix)"] = Dict(:R => ix, :I => ix + n_buses)
    return
end
