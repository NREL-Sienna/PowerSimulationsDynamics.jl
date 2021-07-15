"""
Wraps DynamicBranch devices from PowerSystems to handle changes in controls and connection
status, and allocate the required indexes of the state space.
"""
struct BranchWrapper
    branch::PSY.DynamicBranch
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
    )
        branch_states = PSY.get_states(branch)
        new(
            branch,
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
get_ode_range(wrapper::BranchWrapper) = wrapper.ode_range
get_global_index(wrapper::BranchWrapper) = wrapper.global_index

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
    inputs,
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
