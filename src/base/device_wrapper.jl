get_inner_vars_count(::PSY.DynamicGenerator) = GEN_INNER_VARS_SIZE
get_inner_vars_count(::PSY.DynamicInverter) = INV_INNER_VARS_SIZE
get_inner_vars_count(::PSY.PeriodicVariableSource) = 0

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
struct DynamicWrapper{T <: PSY.DynamicInjection}
    device::T
    bus_category::Type{<:BusCategory}
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
    function DynamicWrapper(
        device::T,
        bus_ix::Int,
        ix_range,
        ode_range,
        inner_var_range,
    ) where {T <: PSY.Device}
        dynamic_device = PSY.get_dynamic_injector(device)
        @assert dynamic_device !== nothing
        device_states = PSY.get_states(dynamic_device)
        component_state_mapping = Dict{Int, Vector{Int}}()
        input_port_mapping = Dict{Int, Vector{Int}}()

        for c in PSY.get_dynamic_components(dynamic_device)
            component_state_mapping[index(c)] =
                _index_local_states(c, device_states)
            input_port_mapping[index(c)] = _index_port_mapping!(c, device_states)
        end

        new{typeof(dynamic_device)}(
            dynamic_device,
            BUS_MAP[PSY.get_bustype(PSY.get_bus(device))],
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

get_bus_category(wrapper::DynamicWrapper) = wrapper.bus_category
get_inner_vars_index(wrapper::DynamicWrapper) = wrapper.inner_vars_index
get_ix_range(wrapper::DynamicWrapper) = wrapper.ix_range
get_ode_range(wrapper::DynamicWrapper) = wrapper.ode_range
get_bus_ix(wrapper::DynamicWrapper) = wrapper.bus_ix
get_global_index(wrapper::DynamicWrapper) = wrapper.global_index
get_component_state_mapping(wrapper::DynamicWrapper) = wrapper.component_state_mapping
get_input_port_mapping(wrapper::DynamicWrapper) = wrapper.input_port_mapping


get_P_ref(wrapper::DynamicWrapper) = wrapper.P_ref[]
get_Q_ref(wrapper::DynamicWrapper) = wrapper.P_ref[]
get_V_ref(wrapper::DynamicWrapper) = wrapper.V_ref[]
get_ω_ref(wrapper::DynamicWrapper) = wrapper.ω_ref[]

set_P_ref(wrapper::DynamicWrapper, val::Float64) = wrapper.P_ref[] = val
set_Q_ref(wrapper::DynamicWrapper, val::Float64) = wrapper.P_ref[] = val
set_V_ref(wrapper::DynamicWrapper, val::Float64) = wrapper.V_ref[] = val
set_ω_ref(wrapper::DynamicWrapper, val::Float64) = wrapper.ω_ref[] = val

# PSY overloads for the wrapper
PSY.get_name(wrapper::DynamicWrapper) = wrapper.device.name
PSY.get_ext(wrapper::DynamicWrapper) = wrapper.device.ext
PSY.get_states(wrapper::DynamicWrapper) = wrapper.device.states
PSY.get_n_states(wrapper::DynamicWrapper) = wrapper.device.n_states
PSY.get_base_power(wrapper::DynamicWrapper) = wrapper.device.base_power

PSY.get_machine(wrapper::DynamicWrapper{T}) where {T <: PSY.DynamicGenerator} =
    wrapper.device.machine
PSY.get_shaft(wrapper::DynamicWrapper{T}) where {T <: PSY.DynamicGenerator} =
    wrapper.device.shaft
PSY.get_avr(wrapper::DynamicWrapper{T}) where {T <: PSY.DynamicGenerator} =
    wrapper.device.avr
PSY.get_prime_mover(wrapper::DynamicWrapper{T}) where {T <: PSY.DynamicGenerator} =
    wrapper.device.prime_mover
PSY.get_pss(wrapper::DynamicWrapper{T}) where {T <: PSY.DynamicGenerator} =
    wrapper.device.pss

PSY.get_converter(wrapper::DynamicWrapper{T}) where {T <: PSY.DynamicInverter} =
    wrapper.device.converter
PSY.get_outer_control(wrapper::DynamicWrapper{T}) where {T <: PSY.DynamicInverter} =
    wrapper.device.outer_control
PSY.get_inner_control(wrapper::DynamicWrapper{T}) where {T <: PSY.DynamicInverter} =
    wrapper.device.inner_control
PSY.get_dc_source(wrapper::DynamicWrapper{T}) where {T <: PSY.DynamicInverter} =
    wrapper.device.dc_source
PSY.get_freq_estimator(wrapper::DynamicWrapper{T}) where {T <: PSY.DynamicInverter} =
    wrapper.device.freq_estimator
PSY.get_filter(wrapper::DynamicWrapper{T}) where {T <: PSY.DynamicInverter} =
    wrapper.device.filter

function get_local_state_ix(wrapper::DynamicWrapper, component::PSY.DynamicComponent)
    return wrapper.device_state_mapping[index(component)]
end

function get_input_port_ix(wrapper::DynamicWrapper, component::PSY.DynamicComponent)
    return wrapper.input_port_mapping[index(component)]
end

struct StaticWrapper{T <: PSY.StaticInjection, V}
    device::T
    connection_status::Base.RefValue{Float64}
    V_ref::Base.RefValue{Float64}
    θ_ref::Base.RefValue{Float64}
    P_ref::Base.RefValue{Float64}
    Q_ref::Base.RefValue{Float64}
    bus_ix::Int
end

function DynamicWrapper(device::T, bus_ix::Int) where {T <: PSY.Device}
    bus = PSY.get_bus(device)
    StaticWrapper{T, BUS_MAP[PSY.get_bustype(bus)]}(
        device,
        Base.Ref(1.0),
        Base.Ref(PSY.get_magnitude(bus)),
        Base.Ref(PSY.get_angle(bus)),
        Base.Ref(PSY.get_active_power(device)),
        Base.Ref(PSY.get_reactive_power(device)),
        bus_ix,
    )
end

function StaticWrapper(device::T, bus_ix::Int) where {T <: PSY.Source}
    bus = PSY.get_bus(device)
    return StaticWrapper{T, BUS_MAP[PSY.get_bustype(bus)]}(
        device,
        Base.Ref(1.0),
        Base.Ref(PSY.get_internal_voltage(device)),
        Base.Ref(PSY.get_internal_angle(device)),
        Base.Ref(PSY.get_active_power(device)),
        Base.Ref(PSY.get_reactive_power(device)),
        bus_ix,
    )
end

# get_bustype is already exported in PSY. So this is named this way to avoid name collisions
get_bus_category(::StaticWrapper{<: PSY.StaticInjection, U}) where {U} = U

function StaticWrapper(device::T, bus_ix::Int) where {T <: PSY.ElectricLoad}
    bus = PSY.get_bus(device)
    return StaticWrapper{T, LOAD_MAP[PSY.get_model(device)]}(
        device,
        Base.Ref(1.0),
        Base.Ref(PSY.get_magnitude(bus)),
        Base.Ref(PSY.get_angle(bus)),
        Base.Ref(PSY.get_active_power(device)),
        Base.Ref(PSY.get_reactive_power(device)),
        bus_ix,
    )
end

# get_bustype is already exported in PSY. So this is named this way to avoid name collisions
get_bus_category(::StaticWrapper{<: PSY.ElectricLoad, U}) where {U} = PQBus
get_load_category(::StaticWrapper{<: PSY.ElectricLoad, U}) where {U} = U

get_P_ref(wrapper::StaticWrapper) = wrapper.P_ref[]
get_Q_ref(wrapper::StaticWrapper) = wrapper.P_ref[]
get_V_ref(wrapper::StaticWrapper) = wrapper.V_ref[]
get_θ_ref(wrapper::StaticWrapper) = wrapper.θ_ref[]

set_P_ref(wrapper::StaticWrapper, val::Float64) = wrapper.P_ref[] = val
set_Q_ref(wrapper::StaticWrapper, val::Float64) = wrapper.P_ref[] = val
set_V_ref(wrapper::StaticWrapper, val::Float64) = wrapper.V_ref[] = val
set_θ_ref(wrapper::StaticWrapper, val::Float64) = wrapper.θ_ref[] = val
