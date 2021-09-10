
"""
Function to obtain series of states out of DAE Solution. It receives the solution, the dynamical system
and a tuple containing the string name of the Dynamic Injection device and the symbol of the state.

"""
function get_state_series(sim::SimulationResults, ref::Tuple{String, Symbol})
    if allunique(sim.solution.t)
        return sim.solution.t, post_proc_state_series(sim, ref)
    else
        @debug "found repeated time steps removing repetitions"
        ix = unique(i -> sim.solution.t[i], eachindex(sim.solution.t))
        return sim.solution.t[ix], post_proc_state_series(sim, ref)[ix]
    end
    return
end

"""
Function to obtain the voltage magnitude series out of the DAE Solution. It receives the solution, the dynamical system and the bus number.

"""
function get_voltage_magnitude_series(sim::SimulationResults, bus_number::Int)
    n_buses = get_bus_count(sim.inputs)
    bus_ix = get(get_lookup(sim.inputs), bus_number, 0)
    V_R, V_I = post_proc_voltage_series(sim.solution, bus_ix, n_buses)
    if allunique(sim.solution.t)
        return sim.solution.t, sqrt.(V_R .^ 2 + V_I .^ 2)
    else
        @debug "found repeated time steps removing repetitions"
        ix = unique(i -> sim.solution.t[i], eachindex(sim.solution.t))
        return sim.solution.t[ix], sqrt.(V_R[ix] .^ 2 + V_I[ix] .^ 2)
    end
end

"""
Function to obtain the voltage angle series out of the DAE Solution. It receives the solution, the dynamical system and the bus number.

"""
function get_voltage_angle_series(sim::SimulationResults, bus_number::Int)
    n_buses = get_bus_count(sim.inputs)
    bus_ix = get(get_lookup(sim.inputs), bus_number, 0)
    V_R, V_I = post_proc_voltage_series(sim.solution, bus_ix, n_buses)
    if allunique(sim.solution.t)
        return sim.solution.t, atan.(V_I ./ V_R)
    else
        @debug "found repeated time steps removing repetitions"
        ix = unique(i -> sim.solution.t[i], eachindex(sim.solution.t))
        return sim.solution.t[ix], atan.(V_I[ix] ./ V_R[ix])
    end
end

"""
Function to obtain the real current time series of a Dynamic Injection series out of the DAE Solution. It receives the solution and the
string name of the Dynamic Injection device.

"""
function get_real_current_series(sim::SimulationResults, name::String)
    if allunique(sim.solution.t)
        return sim.solution.t, post_proc_real_current_series(sim, name)
    else
        @debug "found repeated time steps removing repetitions"
        ix = unique(i -> sim.solution.t[i], eachindex(sim.solution.t))
        return sim.solution.t[ix], post_proc_real_current_series(sim, name)[ix]
    end
end
"""
Function to obtain the imaginary current time series of a Dynamic Injection series out of the DAE Solution. It receives the solution and the
string name of the Dynamic Injection device.

"""
function get_imaginary_current_series(sim::SimulationResults, name::String)
    if allunique(sim.solution.t)
        return sim.solution.t, post_proc_imaginary_current_series(sim, name)
    else
        @debug "found repeated time steps removing repetitions"
        ix = unique(i -> sim.solution.t[i], eachindex(sim.solution.t))
        return sim.solution.t[ix], post_proc_imaginary_current_series(sim, name)[ix]
    end
end

"""
Function to obtain the active power output time series of a Dynamic Injection series out of the DAE Solution. It receives the solution and the
string name of the Dynamic Injection device.

"""
function get_activepower_series(sim::SimulationResults, name::String)
    if allunique(sim.solution.t)
        return sim.solution.t, post_proc_activepower_series(sim, name)
    else
        @debug "found repeated time steps removing repetitions"
        ix = unique(i -> sim.solution.t[i], eachindex(sim.solution.t))
        return sim.solution.t[ix], post_proc_activepower_series(sim, name)[ix]
    end
end

"""
Function to obtain the reactive power output time series of a Dynamic Injection series out of the DAE Solution. It receives the solution and the
string name of the Dynamic Injection device.

"""
function get_reactivepower_series(sim::SimulationResults, name::String)
    if allunique(sim.solution.t)
        return sim.solution.t, post_proc_reactivepower_series(sim, name)
    else
        @debug "found repeated time steps removing repetitions"
        ix = unique(i -> sim.solution.t[i], eachindex(sim.solution.t))
        return sim.solution.t[ix], post_proc_reactivepower_series(sim, name)[ix]
    end
end

"""
Function to print initial states. It receives the vector of initial states and the dynamical system.
"""
function print_device_states(sim::SimulationResults)
    bus_size = get_bus_count(sim.inputs)
    system = get_system(sim)
    println("Voltage Variables")
    println("====================")
    buses_sorted =
        sort(collect(PSY.get_components(PSY.Bus, system)); by = x -> PSY.get_number(x))
    for bus in buses_sorted
        name = PSY.get_name(bus)
        println(name)
        println("====================")
        bus_n = PSY.get_number(bus)
        bus_ix = get_lookup(sim.inputs)[bus_n]
        V_R = sim.x0_init[bus_ix]
        V_I = sim.x0_init[bus_ix + bus_size]
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
        global_index = get_global_index(sim.inputs)[name]
        for s in states
            print(s, " ", round(sim.x0_init[global_index[s]], digits = 4), "\n")
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
            global_index = get_global_index(sim.inputs)[name]
            x0_br = Dict{Symbol, Float64}()
            for (i, s) in enumerate(states)
                print(s, " ", round(sim.x0_init[global_index[s]], digits = 5), "\n")
            end
            println("====================")
        end
    end
    return
end
