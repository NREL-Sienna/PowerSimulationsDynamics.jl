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
