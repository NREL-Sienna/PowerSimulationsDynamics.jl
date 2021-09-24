struct SimulationResults
    global_index::MAPPING_DICT
    bus_lookup::Dict{Int, Int}
    system::PSY.System
    time_log::Dict{Symbol, Any}
    solution::SciMLBase.AbstractODESolution
    function SimulationResults(
        inputs::SimulationInputs,
        system::PSY.System,
        time_log,
        solution,
    )
        new(make_global_state_map(inputs), get_lookup(inputs), system, time_log, solution)
    end
end

get_global_index(res::SimulationResults) = res.global_index
get_bus_count(res::SimulationResults) = get_n_buses(res.system)
get_bus_lookup(res::SimulationResults) = res.bus_lookup
get_system(res::SimulationResults) = res.system
get_solution(res::SimulationResults) = res.solution
