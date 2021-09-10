struct SimulationResults
    global_index::MAPPING_DICT
    voltage_buses::Vector{Int}
    current_buses::Vector{Int}
    solution::SciMLBase.AbstractODESolution
end
