"""
Function to obtain series of states out of DAE Solution. It receives the solution, the dynamical system
and a tuple containing the symbol name of the Dynamic Injection device and the symbol of the state.

"""

function get_state_series(sol::DiffEqBase.DAESolution, dyn_system::DynamicSystem, ref::NTuple{2, Symbol})
    ix = get(dyn_system.global_state_index[ref[1]], ref[2], nothing)
    return sol.t, [value[ix] for value in sol.u]
end


"""
Function to obtain the voltage magnitude series out of the DAE Solution. It receives the solution, the dynamical system and the bus number.

"""

function get_voltagemag_series(sol::DiffEqBase.DAESolution, dyn_system::DynamicSystem, bus_number::Int64)
    bus_size = length(dyn_system.buses)
    return sol.t, [sqrt(value[bus_number]^2 + value[bus_number+bus_size]^2) for value in sol.u]
end


"""
Function to print initial states. It receives the vector of initial states and the dynamical system.
"""

function print_init_states(x0::Vector{Float64}, dyn_system::DynamicSystem)
    for (ix, val_sys) in dyn_system.global_state_index
        println(ix)
        println("====================")
        for (k, val) in val_sys
            print(k, " ", x0[val], "\n")
        end
        println("====================")
    end
end
