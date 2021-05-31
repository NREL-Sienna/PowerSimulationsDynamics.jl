function _attach_ports!(component::PSY.DynamicComponent)
    PSY.get_ext(component)[PORTS] = Ports(component)
    return
end

function index_local_states(component::PSY.DynamicComponent, device_states::Vector{Symbol})
    component_state_index = Vector{Int}(undef, PSY.get_n_states(component))
    component_states = PSY.get_states(component)
    for (ix, s) in enumerate(component_states)
        component_state_index[ix] = findfirst(x -> x == s, device_states)
    end
    return component_state_index
end

function index_port_mapping!(component::PSY.DynamicComponent, device_states::Vector{Symbol})
    _attach_ports!(component)
    index_component_inputs = Vector{Int}()
    for i in PSY.get_ext(component)[PORTS][:states]
        tmp = [(ix, var) for (ix, var) in enumerate(device_states) if var == i]
        isempty(tmp) && continue
        push!(index_component_inputs, tmp[1][1])
    end
    return index_component_inputs
end

function _add_dynamic_bus_states!(
    DAE_vector::Vector{Bool},
    voltage_buses_ix::Vector{Int},
    bus_ix::Int,
    n_buses::Int,
)
    push!(voltage_buses_ix, bus_ix)
    DAE_vector[bus_ix] = DAE_vector[bus_ix + n_buses] = true
    return
end

function _add_to_total_shunts!(total_shunts::Dict{Int, Float64}, pairs...)
    merge!(+, total_shunts, Dict(pairs...))
    return
end

function index_dynamic_lines!(
    inputs::SimulationInputs,
    branch::PSY.DynamicBranch,
    n_buses::Int,
)
    DAE_vector = get_DAE_vector(inputs)
    voltage_buses_ix = get_voltage_buses_ix(inputs)
    arc = PSY.get_arc(branch)
    from_bus_number = PSY.get_number(arc.from)
    to_bus_number = PSY.get_number(arc.to)
    bus_ix_from = get_lookup(inputs)[from_bus_number]
    bus_ix_to = get_lookup(inputs)[to_bus_number]
    b_from = PSY.get_b(branch).from
    b_to = PSY.get_b(branch).to
    total_shunts = get_total_shunts(inputs)
    b_from > 0.0 && _add_to_total_shunts!(total_shunts, bus_ix_from => b_from)
    b_to > 0.0 && _add_to_total_shunts!(total_shunts, bus_ix_to => b_to)
    b_from > 0.0 &&
        _add_dynamic_bus_states!(DAE_vector, voltage_buses_ix, bus_ix_from, n_buses)
    b_to > 0.0 && _add_dynamic_bus_states!(DAE_vector, voltage_buses_ix, bus_ix_to, n_buses)
    n_states = PSY.get_n_states(branch)
    DAE_vector = push!(DAE_vector, collect(trues(n_states))...)
    return
end

function index_static_injection(inputs::SimulationInputs, ::PSY.StaticInjection)
    return
end

function index_dynamic_injection(
    inputs::SimulationInputs,
    dynamic_device::PSY.DynamicInjection,
    state_space_ix::Vector{Int},
)
    state_types = make_device_index!(dynamic_device)
    add_states_to_global!(inputs.global_index, state_space_ix, dynamic_device)
    DAE_vector = get_DAE_vector(inputs)
    push!(DAE_vector, state_types...)
    return
end

"""
Indexes the devices states and maps to the ports
"""
function make_device_index!(dynamic_device::PSY.DynamicInjection)
    device_states = PSY.get_states(dynamic_device)
    device_state_mapping = DEVICE_INTERNAL_MAPPING()
    input_port_mapping = DEVICE_INTERNAL_MAPPING()
    _attach_inner_vars!(dynamic_device)
    dae_vector = collect(trues(PSY.get_n_states(dynamic_device)))
    for c in PSY.get_dynamic_components(dynamic_device)
        device_state_mapping[typeof(c)] = index_local_states(c, device_states)
        input_port_mapping[typeof(c)] = index_port_mapping!(c, device_states)
    end
    dynamic_device.ext[LOCAL_STATE_MAPPING] = device_state_mapping
    dynamic_device.ext[INPUT_PORT_MAPPING] = input_port_mapping
    return dae_vector
end
