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
    connection_status::Base.RefValue{Bool}
    v_ref::Base.RefValue{Float64}
    ω_ref::Base.RefValue{Float64}
    P_ref::Base.RefValue{Float64}
    Q_ref::Base.RefValue{Float64}
    inner_vars_index::Vector{Int}
    ix_range::Vector{Int}
    ode_range::Vector{Int}
    bus_ix::Int
    global_index::Base.ImmutableDict{Symbol, Int}
    device_state_mapping::Base.ImmutableDict{Int, Vector{Int}}
    input_port_mapping::Base.ImmutableDict{Int, Vector{Int}}
    function DeviceWrapper(
        device::T,
        connection_status::Bool
    ) where T <: PSY.Device
        dynamic_device = PSY.get_dynamic_injector(device)
        device_states = PSY.get_states(dynamic_device)
        device_state_mapping = Dict{Int, Vector{Int}}()
        input_port_mapping = Dict{Int, Vector{Int}}()

        for c in PSY.get_dynamic_components(dynamic_device)
            device_state_mapping[index(c)] = PSID.index_local_states(c, device_states)
            input_port_mapping[index(c)] = PSID.index_port_mapping!(c, device_states)
        end

        new{typeof(dynamic_device)}(
            dynamic_device,
            Base.Ref(connection_status),
            Base.Ref(PSY.get_V_ref(dynamic_device)),
            Base.Ref(PSY.get_ω_ref(dynamic_device)),
            Base.Ref(PSY.get_P_ref(dynamic_device)),
            Base.Ref(PSY.get_reactive_power(device)),
            zeros(Int, get_inner_vars_count(dynamic_device)),
            Base.ImmutableDict{Symbol, Int}(),
            Base.ImmutableDict(device_state_mapping...),
            Base.ImmutableDict(input_port_mapping...)
        )
    end
end

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
