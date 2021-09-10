struct SimulationResults
    global_index::MAPPING_DICT
    bus_lookup::Dict{Int, Int}
    system::PSY.System
    solution::SciMLBase.AbstractODESolution
    function SimulationResults(inputs::SimulationInputs, system::PSY.System, solution)
        new(make_global_state_map(inputs), get_lookup(inputs), system, solution)
    end
end

get_global_index(res::SimulationResults) = res.global_index
get_bus_count(res::SimulationResults) = get_n_buses(res.system)
get_bus_lookup(res::SimulationResults) = res.bus_lookup
get_system(res::SimulationResults) = res.system
