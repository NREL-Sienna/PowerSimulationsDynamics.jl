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
    total_shunts[bus_ix_from, bus_ix_from] += 1im * b_from
    total_shunts[bus_ix_to, bus_ix_to] += 1im * b_to
    b_from > 0.0 &&
        _add_dynamic_bus_states!(DAE_vector, voltage_buses_ix, bus_ix_from, n_buses)
    b_to > 0.0 && _add_dynamic_bus_states!(DAE_vector, voltage_buses_ix, bus_ix_to, n_buses)
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
    DAE_vector = get_DAE_vector(inputs)
    make_device_index!(dynamic_device, DAE_vector)
    add_states_to_global!(inputs.global_index, state_space_ix, dynamic_device)
    return
end

function _attach_inner_vars!(
    dynamic_device::PSY.DynamicGenerator,
    ::Type{T} = Real,
) where {T <: Real}
    dynamic_device.ext[INNER_VARS] = zeros(T, 9)
    return
end

function _attach_inner_vars!(
    dynamic_device::PSY.DynamicInverter,
    ::Type{T} = Real,
) where {T <: Real}
    dynamic_device.ext[INNER_VARS] = zeros(T, 14)
    return
end

function _attach_control_refs!(device::PSY.StaticInjection)
    dynamic_device = PSY.get_dynamic_injector(device)
    dynamic_device.ext[CONTROL_REFS] = [
        PSY.get_V_ref(dynamic_device),
        PSY.get_ω_ref(dynamic_device),
        PSY.get_P_ref(dynamic_device),
        PSY.get_reactive_power(device),
    ]
    return
end

function _get_internal_mapping(
    dynamic_device::PSY.DynamicInjection,
    key::AbstractString,
    ty::Type{T},
) where {T <: PSY.DynamicComponent}
    device_index = PSY.get_ext(dynamic_device)[key]
    val = get(device_index, ty, nothing)
    @assert !isnothing(val)
    return val
end

function get_local_state_ix(
    dynamic_device::PSY.DynamicInjection,
    ty::Type{T},
) where {T <: PSY.DynamicComponent}
    return _get_internal_mapping(dynamic_device, LOCAL_STATE_MAPPING, ty)
end

function get_input_port_ix(
    dynamic_device::PSY.DynamicInjection,
    ty::Type{T},
) where {T <: PSY.DynamicComponent}
    return _get_internal_mapping(dynamic_device, INPUT_PORT_MAPPING, ty)
end

"""
Default implementation to index the devices states and maps to the ports
"""
function make_device_index!(dynamic_device::PSY.DynamicInjection, ::Vector{Bool})
    device_states = PSY.get_states(dynamic_device)
    device_state_mapping = DEVICE_INTERNAL_MAPPING()
    input_port_mapping = DEVICE_INTERNAL_MAPPING()
    _attach_inner_vars!(dynamic_device)
    for c in PSY.get_dynamic_components(dynamic_device)
        device_state_mapping[typeof(c)] = index_local_states(c, device_states)
        input_port_mapping[typeof(c)] = index_port_mapping!(c, device_states)
    end
    dynamic_device.ext[LOCAL_STATE_MAPPING] = device_state_mapping
    dynamic_device.ext[INPUT_PORT_MAPPING] = input_port_mapping
    return
end
