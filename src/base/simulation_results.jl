struct SimulationResults
    global_index::MAPPING_DICT
    system::PSY.System
    solution::SciMLBase.AbstractODESolution
    function SimulationResults(inputs::SimulationResults, system::PSY.System, solution)
        new(make_global_state_map(inputs), system, solution)
    end
end
