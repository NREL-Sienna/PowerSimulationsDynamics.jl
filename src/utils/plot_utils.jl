"""
Function to obtain series of states out of DAE Solution. It receives the solution, the dynamical system
and a tuple containing the symbol name of the Dynamic Injection device and the symbol of the state.

"""

function get_state_series(sim::Simulation, ref::Tuple{String, Symbol})
    global_state_index = get_global_index(sim.simulation_inputs)
    ix = get(global_state_index[ref[1]], ref[2], nothing)
    return sim.solution.t, [value[ix] for value in sim.solution.u]
end

"""
Function to obtain the voltage magnitude series out of the DAE Solution. It receives the solution, the dynamical system and the bus number.

"""
function get_voltagemag_series(sim::Simulation, bus_number::Int64)
    n_buses = get_bus_count(sim.simulation_inputs)
    bus_ix = get(get_lookup(sim.simulation_inputs), bus_number, nothing)
    if isnothing(bus_ix)
        @error("Bus number $(bus_number) not found.")
    else
        return sim.solution.t,
        [sqrt(value[bus_ix]^2 + value[bus_ix + n_buses]^2) for value in sim.solution.u]
    end
end

"""
Function to print initial states. It receives the vector of initial states and the dynamical system.
"""
function print_init_states(sim::Simulation)
    for (ix, val_sys) in get_global_index(sim.simulation_inputs)
        println("Differential States")
        println(ix)
        println("====================")
        for (k, val) in val_sys
            print(k, " ", sim.x0_init[val], "\n")
        end
        println("====================")
    end

    # println("Algebraic States") # TODO: Print Buses Voltages
    return
end

"""
Function to print initial states. It receives the vector of initial states and the dynamical system.
"""
function print_device_states(sim::Simulation)
    bus_size = get_bus_count(sim.simulation_inputs)
    system = get_system(sim.simulation_inputs)
    println("Voltage Variables")
    println("====================")
    buses_sorted =
        sort(collect(PSY.get_components(PSY.Bus, system)); by = x -> PSY.get_number(x))
    for bus in buses_sorted
        name = PSY.get_name(bus)
        println(name)
        println("====================")
        bus_n = PSY.get_number(bus)
        bus_ix = get_lookup(sim.simulation_inputs)[bus_n]
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
        global_index = get_global_index(sim.simulation_inputs)[name]
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
            global_index = get_global_index(sim.simulation_inputs)[name]
            x0_br = Dict{Symbol, Float64}()
            for (i, s) in enumerate(states)
                print(s, " ", round(sim.x0_init[global_index[s]], digits = 5), "\n")
            end
            println("====================")
        end
    end
    return
end
