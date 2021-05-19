"""
Internal function to obtain as a Vector of Float64 of a specific state. It receives the solution and the
global index for a state.

"""
function _post_proc_state_series(solution, ix::Int64)
    numel = length(solution.t)::Int64
    state = zeros(Float64, numel)
    for i in 1:numel
        state[i] = solution.u[i][ix]
    end
    return state
end

"""
Function to obtain the state time series of a specific state. It receives the simulation, and a tuple
containing the name of the Dynamic Device and the symbol of the state.

"""
function post_proc_state_series(sim::Simulation, ref::Tuple{String, Symbol})
    global_state_index = get_global_index(sim.simulation_inputs)
    ix = get(global_state_index[ref[1]], ref[2], 0)
    return _post_proc_state_series(sim.solution, ix)
end

"""
Function to obtain voltage and output currents for a dynamic device. It receives the simulation, and the name
of the Dynamic Device.

"""
function post_proc_voltage_current_series(
    sim::Simulation,
    name::String,
)::NTuple{4, Vector{Float64}}
    #Note: Type annotation since get_dynamic_injector is type unstable and solution is Union{Nothing, DAESol}
    sim_inputs = sim.simulation_inputs
    n_buses = get_bus_count(sim_inputs)
    solution = sim.solution
    device = PSY.get_component(PSY.StaticInjection, sim_inputs.sys, name)
    bus_ix = get(get_lookup(sim_inputs), PSY.get_number(PSY.get_bus(device)), -1)
    V_R, V_I = post_proc_voltage_series(solution, bus_ix, n_buses)
    dyn_device = PSY.get_dynamic_injector(device)
    I_R, I_I = compute_output_current(sim, dyn_device, V_R, V_I)
    return V_R, V_I, I_R, I_I
end

"""
Function to obtain voltage using the bus index (and not the bus number). It receives the solution, the bus index and 
the total number of buses.

"""
function post_proc_voltage_series(solution, bus_ix::Int, n_buses::Int)
    bus_ix < 0 && error("Bus number $(bus_number) not found.")
    V_R = Float64[value[bus_ix] for value in solution.u]
    V_I = Float64[value[bus_ix + n_buses] for value in solution.u]
    return V_R, V_I
end

"""
Function to compute the active power output time series of a Dynamic Injection series out of the DAE Solution. It receives the solution and the
string name of the Dynamic Injection device.

"""
function post_proc_activepower_series(sim::Simulation, name::String)
    V_R, V_I, I_R, I_I = post_proc_voltage_current_series(sim, name)
    return V_R .* I_R + V_I .* I_I
end

"""
Function to compute the active power output time series of a Dynamic Injection series out of the DAE Solution. It receives the solution and the
string name of the Dynamic Injection device.

"""
function post_proc_reactivepower_series(sim::Simulation, name::String)
    V_R, V_I, I_R, I_I = post_proc_voltage_current_series(sim, name)
    return V_I .* I_R - V_R .* I_I
end
