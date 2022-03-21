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

"""
Internal function to obtain as a Vector of Float64 of a specific state. It receives the solution and the
global index for a state.

"""
function _post_proc_state_series(solution, ix::Int, dt::Union{Nothing, Float64})
    if dt === nothing
        ix_t = unique(i -> solution.t[i], eachindex(solution.t))
        ts = solution.t[ix_t]
        state = solution[ix, ix_t]
    else
        ts = range(0, stop = solution.t[end], step = dt)
        state = solution(collect(ts); idxs = ix)
    end
    return ts, state
end

"""
Function to obtain the state time series of a specific state. It receives the simulation, and a tuple
containing the name of the Dynamic Device and the symbol of the state.
"""
function post_proc_state_series(
    res::SimulationResults,
    ref::Tuple{String, Symbol},
    dt::Union{Nothing, Float64},
)
    global_state_index = get_global_index(res)
    ix = get(global_state_index[ref[1]], ref[2], 0)
    return _post_proc_state_series(get_solution(res), ix, dt)
end

"""
Function to obtain voltage and output currents for a dynamic device. It receives the simulation, and the name
of the Dynamic Device.
"""
function post_proc_voltage_current_series(
    res::SimulationResults,
    name::String,
    dt::Union{Nothing, Float64},
)::NTuple{5, Vector{Float64}}
    #Note: Type annotation since get_dynamic_injector is type unstable and solution is Union{Nothing, DAESol}
    system = get_system(res)
    bus_lookup = get_bus_lookup(res)
    n_buses = length(bus_lookup)
    solution = res.solution
    device = PSY.get_component(PSY.StaticInjection, system, name)
    if isnothing(device)
        error("Device $(name) not found in the system")
    end
    bus_ix = get(bus_lookup, PSY.get_number(PSY.get_bus(device)), -1)
    ts, V_R, V_I = post_proc_voltage_series(solution, bus_ix, n_buses, dt)
    dyn_device = PSY.get_dynamic_injector(device)
    _, I_R, I_I = compute_output_current(res, dyn_device, V_R, V_I, dt)
    return ts, V_R, V_I, I_R, I_I
end

"""
Function to obtain voltage using the bus index (and not the bus number). It receives the solution, the bus index and
the total number of buses.

"""
function post_proc_voltage_series(
    solution,
    bus_ix::Int,
    n_buses::Int,
    dt::Union{Nothing, Float64},
)
    bus_ix < 0 && error("Bus number $(bus_number) not found.")
    ts, V_R = _post_proc_state_series(solution, bus_ix, dt)
    _, V_I = _post_proc_state_series(solution, bus_ix + n_buses, dt)
    return collect(ts), collect(V_R), collect(V_I)
end

"""
Function to compute the real current output time series of a Dynamic Injection series out of the DAE Solution. It receives the solution and the
string name of the Dynamic Injection device.

"""
function post_proc_real_current_series(
    res::SimulationResults,
    name::String,
    dt::Union{Nothing, Float64},
)
    ts, _, _, I_R, _ = post_proc_voltage_current_series(res, name, dt)
    return ts, I_R
end
"""
Function to compute the imaginary current output time series of a Dynamic Injection series out of the DAE Solution. It receives the solution and the
string name of the Dynamic Injection device.

"""
function post_proc_imaginary_current_series(
    res::SimulationResults,
    name::String,
    dt::Union{Nothing, Float64},
)
    ts, _, _, _, I_I = post_proc_voltage_current_series(res, name, dt)
    return ts, I_I
end

"""
Function to compute the active power output time series of a Dynamic Injection series out of the DAE Solution. It receives the solution and the
string name of the Dynamic Injection device.

"""
function post_proc_activepower_series(
    res::SimulationResults,
    name::String,
    dt::Union{Nothing, Float64},
)
    ts, V_R, V_I, I_R, I_I = post_proc_voltage_current_series(res, name, dt)
    return ts, V_R .* I_R + V_I .* I_I
end

"""
Function to compute the active power output time series of a Dynamic Injection series out of the DAE Solution. It receives the solution and the
string name of the Dynamic Injection device.

"""
function post_proc_reactivepower_series(
    res::SimulationResults,
    name::String,
    dt::Union{Nothing, Float64},
)
    ts, V_R, V_I, I_R, I_I = post_proc_voltage_current_series(res, name, dt)
    return ts, V_I .* I_R - V_R .* I_I
end

"""
    get_state_series(
        res::SimulationResults,
        ref::Tuple{String, Symbol};
        dt::Union{Nothing, Float64} = nothing
    )
    end

Function to obtain series of states out of DAE Solution.

# Arguments

- `res::SimulationResults` : Simulation Results object that contains the solution
- `ref:Tuple{String, Symbol}` : Tuple used to identify the dynamic device, via its name, as a `String`, and the associated state as a `Symbol`.
"""
function get_state_series(res::SimulationResults, ref::Tuple{String, Symbol}; dt = nothing)
    return post_proc_state_series(res, ref, dt)
end

"""
    get_voltage_magnitude_series(
        res::SimulationResults,
        bus_number::Int
    )

Function to obtain the voltage magnitude series out of the DAE Solution.

# Arguments:

- `res::SimulationResults` : Simulation Results object that contains the solution
- `bus_number::Int` : Bus number identifier
"""
function get_voltage_magnitude_series(res::SimulationResults, bus_number::Int; dt = nothing)
    n_buses = get_bus_count(res)
    bus_ix = get(get_bus_lookup(res), bus_number, 0)
    ts, V_R, V_I = post_proc_voltage_series(res.solution, bus_ix, n_buses, dt)
    return ts, sqrt.(V_R .^ 2 .+ V_I .^ 2)
end

"""
    get_voltage_angle_series(
        res::SimulationResults,
        bus_number::Int
    )

Function to obtain the voltage angle series out of the DAE Solution.

# Arguments

- `res::SimulationResults` : Simulation Results object that contains the solution
- `bus_number::Int` : Bus number identifier
"""
function get_voltage_angle_series(res::SimulationResults, bus_number::Int; dt = nothing)
    n_buses = get_bus_count(res)
    bus_ix = get(get_bus_lookup(res), bus_number, 0)
    ts, V_R, V_I = post_proc_voltage_series(res.solution, bus_ix, n_buses, dt)
    return ts, atan.(V_I ./ V_R)
end

"""
    get_real_current_series(
            res::SimulationResults,
            name::String,
    )

Function to obtain the real current time series of a Dynamic Injection series out of the DAE Solution.

# Arguments

- `res::SimulationResults` : Simulation Results object that contains the solution
- `name::String` : Name to identify the specified device
"""
function get_real_current_series(res::SimulationResults, name::String; dt = nothing)
    return post_proc_real_current_series(res, name, dt)
end

"""
    get_imaginary_current_series(
            res::SimulationResults,
            name::String,
    )

Function to obtain the imaginary current time series of a Dynamic Injection series out of the DAE Solution.

# Arguments

- `res::SimulationResults` : Simulation Results object that contains the solution
- `name::String` : Name to identify the specified device
"""
function get_imaginary_current_series(res::SimulationResults, name::String; dt = nothing)
    return post_proc_imaginary_current_series(res, name, dt)
end

"""
    get_activepower_series(
            res::SimulationResults,
            name::String,
    )

Function to obtain the active power output time series of a Dynamic Injection series out of the DAE Solution.

# Arguments

- `res::SimulationResults` : Simulation Results object that contains the solution
- `name::String` : Name to identify the specified device
"""
function get_activepower_series(res::SimulationResults, name::String; dt = nothing)
    return post_proc_activepower_series(res, name, dt)
end

"""
    get_reactivepower_series(
            res::SimulationResults,
            name::String,
    )

Function to obtain the reactive power output time series of a Dynamic Injection series out of the DAE Solution.

# Arguments

- `res::SimulationResults` : Simulation Results object that contains the solution
- `name::String` : Name to identify the specified device
"""
function get_reactivepower_series(res::SimulationResults, name::String; dt = nothing)
    return post_proc_reactivepower_series(res, name, dt)
end

"""
    show_states_initial_value(res::SimulationResults)

Function to print initial states.

# Arguments

- `res::SimulationResults` : Simulation Results object that contains the solution
"""
function show_states_initial_value(res::SimulationResults)
    bus_size = get_bus_count(res)
    system = get_system(res)
    x0_init = res.solution.u[1]
    println("Voltage Variables")
    println("====================")
    buses_sorted =
        sort(collect(PSY.get_components(PSY.Bus, system)); by = x -> PSY.get_number(x))
    for bus in buses_sorted
        name = PSY.get_name(bus)
        println(name)
        println("====================")
        bus_n = PSY.get_number(bus)
        bus_ix = get_bus_lookup(res)[bus_n]
        V_R = x0_init[bus_ix]
        V_I = x0_init[bus_ix + bus_size]
        Vm = sqrt(V_R^2 + V_I^2)
        θ = angle(V_R + V_I * 1im)
        print("Vm ", round(Vm, digits = 4), "\n")
        print("θ ", round(θ, digits = 4), "\n")
        println("====================")
    end
    println("====================")
    for device in PSY.get_components(PSY.DynamicInjection, system)
        states = PSY.get_states(device)
        name = PSY.get_name(device)
        println("Differential States")
        println(name)
        println("====================")
        global_index = get_global_index(res)[name]
        for s in states
            print(s, " ", round(x0_init[global_index[s]], digits = 4), "\n")
        end
        println("====================")
    end
    dyn_branches = PSY.get_components(PSY.DynamicBranch, system)
    if !isempty(dyn_branches)
        println("====================")
        println("Line Current States")
        for br in dyn_branches
            states = PSY.get_states(br)
            name = PSY.get_name(br)
            printed_name = "Line " * name
            println("====================")
            println(printed_name)
            global_index = get_global_index(res)[name]
            for (i, s) in enumerate(states)
                print(s, " ", round(x0_init[global_index[s]], digits = 5), "\n")
            end
            println("====================")
        end
    end
    return
end
