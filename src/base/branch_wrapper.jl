"""
Wraps DynamicBranch devices from PowerSystems to handle changes in controls and connection
status, and allocate the required indexes of the state space.
"""
struct BranchWrapper
    branch::PSY.DynamicBranch
    system_base_power::Float64
    system_base_frequency::Float64
    connection_status::Base.RefValue{Float64}
    bus_ix_from::Int
    bus_ix_to::Int
    ix_range::Vector{Int}
    ode_range::Vector{Int}
    global_index::Base.ImmutableDict{Symbol, Int}
    function BranchWrapper(
        branch::PSY.DynamicBranch,
        bus_ix_from::Int,
        bus_ix_to::Int,
        ix_range,
        ode_range,
        sys_base_power,
        sys_base_freq,
    )
        branch_states = PSY.get_states(branch)
        new(
            branch,
            sys_base_power,
            sys_base_freq,
            Base.Ref(1.0),
            bus_ix_from,
            bus_ix_to,
            ix_range,
            ode_range,
            Base.ImmutableDict(Dict(branch_states .=> ix_range)...),
        )
    end
end

get_connection_status(wrapper::BranchWrapper) = wrapper.connection_status[]
get_bus_ix_from(wrapper::BranchWrapper) = wrapper.bus_ix_from
get_bus_ix_to(wrapper::BranchWrapper) = wrapper.bus_ix_to
get_ix_range(wrapper::BranchWrapper) = wrapper.ix_range
get_ode_ouput_range(wrapper::BranchWrapper) = wrapper.ode_range
get_global_index(wrapper::BranchWrapper) = wrapper.global_index

get_system_base_power(wrapper::BranchWrapper) = wrapper.system_base_power
get_system_base_frequency(wrapper::BranchWrapper) = wrapper.system_base_frequency

PSY.get_name(wrapper::BranchWrapper) = PSY.get_name(wrapper.branch)
PSY.get_available(wrapper::BranchWrapper) = PSY.get_available(wrapper.branch)
PSY.get_active_power_flow(wrapper::BranchWrapper) = PSY.get_active_power(wrapper.branch)
PSY.get_reactive_power_flow(wrapper::BranchWrapper) = PSY.get_reactive_power(wrapper.branch)
PSY.get_arc(wrapper::BranchWrapper) = PSY.get_arc(wrapper.branch)
PSY.get_r(wrapper::BranchWrapper) = PSY.get_r(wrapper.branch)
PSY.get_x(wrapper::BranchWrapper) = PSY.get_x(wrapper.branch)
PSY.get_b(wrapper::BranchWrapper) = PSY.get_b(wrapper.branch)
PSY.get_rate(wrapper::BranchWrapper) = PSY.get_rate(wrapper.branch)
PSY.get_angle_limits(wrapper::BranchWrapper) = PSY.get_angle_limits(wrapper.branch)
PSY.get_ext(wrapper::BranchWrapper) = PSY.get_ext(wrapper.branch)
PSY.get_n_states(wrapper::BranchWrapper) = PSY.get_n_states(wrapper.branch)

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
