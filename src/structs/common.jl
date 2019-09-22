mutable struct Ports <: DynDevice
    state::Vector{Symbol}
    inner::Vector
end

function _make_state_mapping(components::Vector{DT}, states::Vector{Symbol}) where DT <: DynamicComponent

    local_state_mapping = Dict{DT, Vector{Int64}}()
    local_state_space_ix = 0
    for c in components
        component_state_index= Vector{Int64}(undef, length(c.states))
        for (ix,i) in enumerate(c.states)
            local_state_space_ix += 1
            component_state_index[ix] = local_state_space_ix
        end
        local_state_mapping[c] = component_state_index
    end

    return local_state_mapping

end

function _make_port_mapping(components::Vector{DT}, states::Vector{Symbol}) where DT <: DynamicComponent

    input_port_mapping = Dict{DT, Vector{Int64}}()

    for c in components
        index_component_in = Vector{Int64}()
        for i in c.ports.state
            tmp = [(ix,var) for (ix,var) in enumerate(states) if var == i]
            isempty(tmp) && continue
            push!(index_component_in, tmp[1][1])
        end
        input_port_mapping[c] = index_component_in
    end

    return input_port_mapping

end
