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
    bus_size = get_bus_count(sim.system)
    println("Voltage Variables")
    println("====================")
    buses = PSY.get_components(PSY.Bus, sim.system)
    bus_ax = [PSY.get_number(bus) for bus in buses]
    sort!(bus_ax)
    look_up = PSY._make_ax_ref(bus_ax)
    buses_sorted = Vector{Any}(undef, bus_size)
    for b in buses
        buses_sorted[look_up[PSY.get_number(b)]] = b
    end
    for bus in buses_sorted
        name = PSY.get_name(bus)
        println(name)
        println("====================")
        bus_n = PSY.get_number(bus)
        bus_ix = PSY.get_ext(sim.system)[LOOKUP][bus_n]
        V_R = sim.x0_init[bus_ix]
        V_I = sim.x0_init[bus_ix + bus_size]
        Vm = sqrt(V_R^2 + V_I^2)
        θ = angle(V_R + V_I * 1im)
        print("Vm ", round(Vm, digits = 4), "\n")
        print("θ ", round(θ, digits = 4), "\n")
        println("====================")
    end
    println("====================")
    for device in PSY.get_components(PSY.DynamicInjection, sim.system)
        states = PSY.get_states(device)
        name = PSY.get_name(device)
        println("Differential States")
        println(name)
        println("====================")
        global_index = PSY.get_ext(sim.system)[GLOBAL_INDEX][name]
        for s in states
            print(s, " ", round(sim.x0_init[global_index[s]], digits = 4), "\n")
        end
        println("====================")
    end
    dyn_branches = PSY.get_components(PSY.DynamicBranch, sim.system)
    if !isempty(dyn_branches)
        println("====================")
        println("Line Current States")
        for br in dyn_branches
            states = PSY.get_states(br)
            name = PSY.get_name(br)
            printed_name = "Line " * name
            println("====================")
            println(printed_name)
            global_index = PSY.get_ext(sim.system)[GLOBAL_INDEX][name]
            x0_br = Dict{Symbol, Float64}()
            for (i, s) in enumerate(states)
                print(s, " ", round(sim.x0_init[global_index[s]], digits = 5), "\n")
            end
            println("====================")
        end
    end
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
    for device in PSY.get_components(PSY.DynamicInjection, sim.system)
        states = PSY.get_states(device)
        name = PSY.get_name(device)
        global_index = PSY.get_ext(sim.system)[GLOBAL_INDEX][name]
        x0_device = Vector{Float64}(undef, length(states))
        for (i, s) in enumerate(states)
            x0_device[i] = sim.x0_init[global_index[s]]
        end
        results[name] = x0_device
    end
    dyn_branches = PSY.get_components(PSY.DynamicBranch, sim.system)
    if !isempty(dyn_branches)
        for br in dyn_branches
            states = PSY.get_states(br)
            name = PSY.get_name(br)
            global_index = PSY.get_ext(sim.system)[GLOBAL_INDEX][name]
            x0_br = Vector{Float64}(undef, length(states))
            for (i, s) in enumerate(states)
                x0_br[i] = sim.x0_init[global_index[s]]
            end
            printed_name = "Line " * name
            results[printed_name] = x0_br
        end
    end
    return results
end

function get_dict_init_states_named(sim::Simulation)
    bus_size = get_bus_count(sim.system)
    V_R = Dict{Int64, Float64}()
    V_I = Dict{Int64, Float64}()
    Vm = Dict{Int64, Float64}()
    θ = Dict{Int64, Float64}()
    for bus in PSY.get_components(PSY.Bus, sim.system)
        bus_n = PSY.get_number(bus)
        bus_ix = PSY.get_ext(sim.system)[LOOKUP][bus_n]
        V_R[bus_n] = sim.x0_init[bus_ix]
        V_I[bus_n] = sim.x0_init[bus_ix + bus_size]
        Vm[bus_n] = sqrt(sim.x0_init[bus_ix]^2 + sim.x0_init[bus_ix + bus_size]^2)
        θ[bus_n] = angle(sim.x0_init[bus_ix] + sim.x0_init[bus_ix + bus_size] * 1im)
    end
    results = Dict{String, Any}("V_R" => V_R, "V_I" => V_I, "Vm" => Vm, "θ" => θ)
    for device in PSY.get_components(PSY.DynamicInjection, sim.system)
        states = PSY.get_states(device)
        name = PSY.get_name(device)
        global_index = PSY.get_ext(sim.system)[GLOBAL_INDEX][name]
        x0_device = Dict{Symbol, Float64}()
        for (i, s) in enumerate(states)
            x0_device[s] = sim.x0_init[global_index[s]]
        end
        results[name] = x0_device
    end
    dyn_branches = PSY.get_components(PSY.DynamicBranch, sim.system)
    if !isempty(dyn_branches)
        for br in dyn_branches
            states = PSY.get_states(br)
            name = PSY.get_name(br)
            global_index = PSY.get_ext(sim.system)[GLOBAL_INDEX][name]
            x0_br = Dict{Symbol, Float64}()
            for (i, s) in enumerate(states)
                x0_br[s] = sim.x0_init[global_index[s]]
            end
            printed_name = "Line " * name
            results[printed_name] = x0_br
        end
    end
    return results
end
