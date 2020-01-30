"""
Function to obtain series of states out of DAE Solution. It receives the solution, the dynamical system
and a tuple containing the symbol name of the Dynamic Injection device and the symbol of the state.

"""

function get_state_series(sim::Simulation, ref::Tuple{String, Symbol})
    global_state_index = PSY.get_ext(sim.system)[GLOBAL_INDEX]
    ix = get(global_state_index[ref[1]], ref[2], nothing)
    return sim.solution.t, [value[ix] for value in sim.solution.u]
end

"""
Function to obtain the voltage magnitude series out of the DAE Solution. It receives the solution, the dynamical system and the bus number.

"""
function get_voltagemag_series(sim::Simulation, bus_number::Int64)
    n_buses = length(PSY.get_components(PSY.Bus, sim.system))
    return sim.solution.t, [sqrt(value[bus_number]^2 + value[bus_number+n_buses]^2) for value in sim.solution.u]
end


"""
Function to print initial states. It receives the vector of initial states and the dynamical system.
"""
function print_init_states(sim::Simulation)
    for (ix, val_sys) in PSY.get_ext(sim.system)[GLOBAL_INDEX]
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
