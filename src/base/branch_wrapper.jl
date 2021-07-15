"""
Wraps DynamicBranch devices from PowerSystems to handle changes in controls and connection
status, and allocate the required indexes of the state space.
"""
struct BranchWrapper
    device::PSY.DynamicBranch
    connection_status::Base.RefValue{Bool}
    bus_ix_from::Int
    bus_ix_to::Int
    ix_range::Vector{Int}
    ode_range::Vector{Int}
    global_index::Base.ImmutableDict{Symbol, Int}
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
