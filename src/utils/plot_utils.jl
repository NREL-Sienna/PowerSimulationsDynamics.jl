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