function _attach_ports!(component::PSY.DynamicComponent)
    PSY.get_ext(component)[PORTS] = Ports(component)
    return
end

function index_local_states(component::PSY.DynamicComponent, device_states::Vector{Symbol})
    component_state_index = Vector{Int64}(undef, PSY.get_n_states(component))
    component_states = PSY.get_states(component)
    for (ix, s) in enumerate(component_states)
        component_state_index[ix] = findfirst(x -> x == s, device_states)
    end
    return component_state_index
end

function index_port_mapping!(component::PSY.DynamicComponent, device_states::Vector{Symbol})
    _attach_ports!(component)
    index_component_inputs = Vector{Int64}()
    for i in PSY.get_ext(component)[PORTS][:states]
        tmp = [(ix, var) for (ix, var) in enumerate(device_states) if var == i]
        isempty(tmp) && continue
        push!(index_component_inputs, tmp[1][1])
    end
    return index_component_inputs
end

"""
Indexes the devices states and maps to the ports
"""
function make_device_index!(device::PSY.StaticInjection)
    dynamic_device = PSY.get_dynamic_injector(device)
    device_states = PSY.get_states(dynamic_device)
    device_state_mapping = Dict{Type{<:PSY.DynamicComponent}, Vector{Int64}}()
    input_port_mapping = Dict{Type{<:PSY.DynamicComponent}, Vector{Int64}}()
    _attach_inner_vars!(dynamic_device)
    _attach_control_refs!(device)
    dae_vector = collect(trues(PSY.get_n_states(dynamic_device)))
    for c in PSY.get_dynamic_components(dynamic_device)
        device_state_mapping[typeof(c)] = index_local_states(c, device_states)
        input_port_mapping[typeof(c)] = index_port_mapping!(c, device_states)
    end
    dynamic_device.ext[LOCAL_STATE_MAPPING] = device_state_mapping
    dynamic_device.ext[INPUT_PORT_MAPPING] = input_port_mapping
    return dae_vector
end
