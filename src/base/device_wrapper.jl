get_inner_vars_count(::PSY.DynamicGenerator) = GEN_INNER_VARS_SIZE
get_inner_vars_count(::PSY.DynamicInverter) = INV_INNER_VARS_SIZE

index(::PSY.TurbineGov) = 1
index(::PSY.PSS) = 2
index(::PSY.Machine) = 3
index(::PSY.Shaft) = 4
index(::PSY.AVR) = 5

index(::PSY.DCSource) = 1
index(::PSY.FrequencyEstimator) = 2
index(::PSY.OuterControl) = 3
index(::PSY.InnerControl) = 4
index(::PSY.Converter) = 5
index(::PSY.Filter) = 6

"""
Wraps DynamicInjection devices from PowerSystems to handle changes in controls and connection
status, and allocate the required indexes of the state space.
"""
struct DeviceWrapper{T <: PSY.DynamicInjection}
    device::T
    connection_status::Base.RefValue{Float64}
    v_ref::Base.RefValue{Float64}
    ω_ref::Base.RefValue{Float64}
    P_ref::Base.RefValue{Float64}
    Q_ref::Base.RefValue{Float64}
    inner_vars_index::Vector{Int}
    ix_range::Vector{Int}
    ode_range::Vector{Int}
    bus_ix::Int
    global_index::Base.ImmutableDict{Symbol, Int}
    component_state_mapping::Base.ImmutableDict{Int, Vector{Int}}
    input_port_mapping::Base.ImmutableDict{Int, Vector{Int}}
    function DeviceWrapper(
        device::T,
        bus_ix::Int,
        ix_range,
        ode_range,
        inner_var_range,
    ) where {T <: PSY.Device}
        dynamic_device = PSY.get_dynamic_injector(device)
        device_states = PSY.get_states(dynamic_device)
        component_state_mapping = Dict{Int, Vector{Int}}()
        input_port_mapping = Dict{Int, Vector{Int}}()

        for c in PSY.get_dynamic_components(dynamic_device)
            component_state_mapping[PSID.index(c)] =
                PSID._index_local_states(c, device_states)
            input_port_mapping[PSID.index(c)] = PSID._index_port_mapping!(c, device_states)
        end

        new{typeof(dynamic_device)}(
            dynamic_device,
            Base.Ref(1.0),
            Base.Ref(PSY.get_V_ref(dynamic_device)),
            Base.Ref(PSY.get_ω_ref(dynamic_device)),
            Base.Ref(PSY.get_P_ref(dynamic_device)),
            Base.Ref(PSY.get_reactive_power(device)),
            inner_var_range,
            ix_range,
            ode_range,
            bus_ix,
            Base.ImmutableDict(Dict(device_states .=> ix_range)...),
            Base.ImmutableDict(component_state_mapping...),
            Base.ImmutableDict(input_port_mapping...),
        )
    end
end

function _index_local_states(component::PSY.DynamicComponent, device_states::Vector{Symbol})
    component_state_index = Vector{Int}(undef, PSY.get_n_states(component))
    component_states = PSY.get_states(component)
    for (ix, s) in enumerate(component_states)
        component_state_index[ix] = findfirst(x -> x == s, device_states)
    end
    return component_state_index
end

function _index_port_mapping!(
    component::PSY.DynamicComponent,
    device_states::Vector{Symbol},
)
    ports = Ports(component)
    index_component_inputs = Vector{Int}()
    for i in ports[:states]
        tmp = [(ix, var) for (ix, var) in enumerate(device_states) if var == i]
        isempty(tmp) && continue
        push!(index_component_inputs, tmp[1][1])
    end
    return index_component_inputs
end

get_inner_vars_index(wrapper::DeviceWrapper) = wrapper.inner_vars_index
get_ix_range(wrapper::DeviceWrapper) = wrapper.ix_range
get_ode_range(wrapper::DeviceWrapper) = wrapper.ode_range
get_bus_ix(wrapper::DeviceWrapper) = wrapper.bus_ix
get_global_index(wrapper::DeviceWrapper) = wrapper.global_index
get_component_state_mapping(wrapper::DeviceWrapper) = wrapper.component_state_mapping
get_input_port_mapping(wrapper::DeviceWrapper) = wrapper.input_port_mapping

get_P_ref(wrapper::DeviceWrapper) = wrapper.P_ref[]
get_Q_ref(wrapper::DeviceWrapper) = wrapper.P_ref[]
get_V_ref(wrapper::DeviceWrapper) = wrapper.V_ref[]
get_ω_ref(wrapper::DeviceWrapper) = wrapper.ω_ref[]

set_P_ref(wrapper::DeviceWrapper, val) = wrapper.P_ref[] = val
set_Q_ref(wrapper::DeviceWrapper, val) = wrapper.P_ref[] = val
set_V_ref(wrapper::DeviceWrapper, val) = wrapper.V_ref[] = val
set_ω_ref(wrapper::DeviceWrapper, val) = wrapper.ω_ref[] = val

# PSY overloads for the wrapper
PSY.get_name(wrapper::DeviceWrapper) = wrapper.device.name
PSY.get_ext(wrapper::DeviceWrapper) = wrapper.device.ext
PSY.get_states(wrapper::DeviceWrapper) = wrapper.device.states
PSY.get_n_states(wrapper::DeviceWrapper) = wrapper.device.n_states
PSY.get_base_power(wrapper::DeviceWrapper) = wrapper.device.base_power

PSY.get_machine(wrapper::DeviceWrapper{T}) where {T <: PSY.DynamicGenerator} =
    wrapper.device.machine
PSY.get_shaft(wrapper::DeviceWrapper{T}) where {T <: PSY.DynamicGenerator} =
    wrapper.device.shaft
PSY.get_avr(wrapper::DeviceWrapper{T}) where {T <: PSY.DynamicGenerator} =
    wrapper.device.avr
PSY.get_prime_mover(wrapper::DeviceWrapper{T}) where {T <: PSY.DynamicGenerator} =
    wrapper.device.prime_mover
PSY.get_pss(wrapper::DeviceWrapper{T}) where {T <: PSY.DynamicGenerator} =
    wrapper.device.pss

PSY.get_converter(wrapper::DeviceWrapper{T}) where {T <: PSY.DynamicInverter} =
    wrapper.device.converter
PSY.get_outer_control(wrapper::DeviceWrapper{T}) where {T <: PSY.DynamicInverter} =
    wrapper.device.outer_control
PSY.get_inner_control(wrapper::DeviceWrapper{T}) where {T <: PSY.DynamicInverter} =
    wrapper.device.inner_control
PSY.get_dc_source(wrapper::DeviceWrapper{T}) where {T <: PSY.DynamicInverter} =
    wrapper.device.dc_source
PSY.get_freq_estimator(wrapper::DeviceWrapper{T}) where {T <: PSY.DynamicInverter} =
    wrapper.device.freq_estimator
PSY.get_filter(wrapper::DeviceWrapper{T}) where {T <: PSY.DynamicInverter} =
    wrapper.device.filter

function get_local_state_ix(
    wrapper::DeviceWrapper,
    component::T,
) where {T <: PSY.DynamicComponent}
    return wrapper.device_state_mapping[index(component)]
end

function get_input_port_ix(
    wrapper::DeviceWrapper,
    component::T,
) where {T <: PSY.DynamicComponent}
    return wrapper.input_port_mapping[index(component)]
end
