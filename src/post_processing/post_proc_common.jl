function make_global_state_map(inputs::SimulationInputs)
    dic = inputs.global_state_map
    if !isempty(dic)
        return dic
    end
    device_wrappers = get_dynamic_injectors(inputs)
    branches_wrappers = get_dynamic_branches(inputs)
    buses_diffs = get_voltage_buses_ix(inputs)
    n_buses = get_bus_count(inputs)
    for d in device_wrappers
        dic[PSY.get_name(d)] = OrderedDict{Symbol, Int}(get_global_index(d)...)
    end
    if !isempty(branches_wrappers)
        for br in branches_wrappers
            dic[PSY.get_name(br)] = OrderedDict{Symbol, Int}(get_global_index(br)...)
        end
    end
    if !isempty(buses_diffs)
        for ix in buses_diffs
            dic["V_$(ix)"] = OrderedDict{Symbol, Int}(:R => ix, :I => ix + n_buses)
        end
    end
    return dic
end

function get_state_from_ix(global_index::MAPPING_DICT, idx::Int)
    for (name, device_ix) in global_index
        if idx ∈ values(device_ix)
            state = [state for (state, number) in device_ix if number == idx]
            IS.@assert_op length(state) == 1
            return name, state[1]
        end
    end
    error("State with index $(idx) not found in the global index")
end

function _obtain_shunt_current(
    branch::B,
    V_R_from,
    V_I_from,
    V_R_to,
    V_I_to,
) where {B <: PSY.Branch}
    error(
        "Current branch $(PSY.get_name(branch)) of type $(typeof(branch)) is not supported for this calculation",
    )
end

function _obtain_shunt_current(branch::PSY.Line, V_R_from, V_I_from, V_R_to, V_I_to)
    b_from, b_to = PSY.get_b(branch)
    z_from = 1.0 / (b_from * 1im)
    z_to = 1.0 / (b_to * 1im)
    I_shunt_from = (0.0 - (V_R_from + V_I_from * 1im)) / z_from
    I_shunt_to = (0.0 - (V_R_to + V_I_to * 1im)) / z_to
    return real.(I_shunt_from), imag.(I_shunt_from), real.(I_shunt_to), imag.(I_shunt_to)
end
