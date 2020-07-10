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
    bus_ix = get(PSY.get_ext(sim.system)[LOOKUP], bus_number, nothing)
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

"""
Function to print initial states. It receives the vector of initial states and the dynamical system.
"""
function print_device_states(sim::Simulation)
    for (ix, val_sys) in PSY.get_ext(sim.system)[GLOBAL_INDEX]
        ix_dyn_injector =
            PSY.get_dynamic_injector(PSY.get_component(PSY.StaticInjection, sim.system, ix))
        if !isnothing(ix_dyn_injector)
            println("Differential States")
            println(ix)
            println("====================")
            ix_states = PSY.get_states(ix_dyn_injector)
            for k in ix_states
                print(k, " ", sim.x0_init[val_sys[k]], "\n")
            end
            println("====================")
        end
    end

    #println("Algebraic States") # TODO: Print Buses Voltages
    return
end

function get_dict_init_states(sim::Simulation)
    bus_size = get_bus_count(sim.system)
    V_R = Vector{Float64}(undef, bus_size)
    V_I = Vector{Float64}(undef, bus_size)
    Vm = Vector{Float64}(undef, bus_size)
    θ = Vector{Float64}(undef, bus_size)
    for bus in PSY.get_components(PSY.Bus, sim.system)
        bus_n = PSY.get_number(bus)
        bus_ix = PSY.get_ext(sim.system)[LOOKUP][bus_n]
        V_R[bus_ix] = sim.x0_init[bus_ix]
        V_I[bus_ix] = sim.x0_init[bus_ix + bus_size]
        Vm[bus_ix] = sqrt(V_R[bus_ix]^2 + V_I[bus_ix]^2)
        θ[bus_ix] = angle(V_R[bus_ix] + V_I[bus_ix] * 1im)
    end
    results =
        Dict{String, Vector{Float64}}("V_R" => V_R, "V_I" => V_I, "Vm" => Vm, "θ" => θ)
    for (ix, val_sys) in PSY.get_ext(sim.system)[GLOBAL_INDEX]
        ix_dyn_injector =
            PSY.get_dynamic_injector(PSY.get_component(PSY.StaticInjection, sim.system, ix))
        if !isnothing(ix_dyn_injector)
            ix_states = PSY.get_states(ix_dyn_injector)
            x0_device = Vector{Float64}(undef, length(ix_states))
            for (i, state) in enumerate(ix_states)
                x0_device[i] = sim.x0_init[val_sys[state]]
            end
            results[ix] = x0_device
        end
    end
    return results
end
