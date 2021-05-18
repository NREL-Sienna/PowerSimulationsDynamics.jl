function _pp_state_series(solution, ix::Int64)    
    numel = length(solution.t)::Int64
    state = zeros(Float64, numel)
    for i in 1:numel
        state[i] = solution.u[i][ix]
    end
    return state
end

function pp_state_series(sim::Simulation, ref::Tuple{String, Symbol})
    global_state_index = get_global_index(sim.simulation_inputs)
    ix = get(global_state_index[ref[1]], ref[2], 0)
    return _pp_state_series(sim.solution, ix)
end

"""
Function to obtain the active power output time series of a Dynamic Injection series out of the DAE Solution. It receives the solution and the
string name of the Dynamic Injection device.

"""
function get_activepower_series(sim::Simulation, name::String)
    sim_inputs = sim.simulation_inputs
    n_buses = get_bus_count(sim_inputs)
    solution = sim.solution
    device = PSY.get_component(PSY.StaticInjection, sim_inputs.sys, name)
    bus_ix = get(get_lookup(sim_inputs), PSY.get_number(PSY.get_bus(device)), nothing)
    V_R = Float64[value[bus_ix] for value in solution.u]
    V_I = Float64[value[bus_ix + n_buses] for value in solution.u]
    dyn_device = PSY.get_dynamic_injector(device)
    I_R, I_I = compute_output_current(sim, dyn_device, V_R, V_I)
    return solution.t, V_R .* I_R + V_I .* I_I
end
