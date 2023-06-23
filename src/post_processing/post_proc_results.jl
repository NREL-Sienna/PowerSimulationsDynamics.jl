"""
    show_states_initial_value(sim::Simulation)

Function to print initial states.

# Arguments

- `sim::Simulation` : Simulation object that contains the initial condition
"""
function show_states_initial_value(sim::Simulation)
    inputs = sim.inputs
    bus_size = get_bus_count(inputs)
    system = get_system(sim)
    global_state_map = make_global_state_map(inputs)
    println("Voltage Variables")
    println("====================")
    buses_sorted =
        sort(collect(PSY.get_components(PSY.Bus, system)); by = x -> PSY.get_number(x))
    for bus in buses_sorted
        name = PSY.get_name(bus)
        println(name)
        println("====================")
        bus_n = PSY.get_number(bus)
        bus_ix = PSID.get_lookup(sim.inputs)[bus_n]
        V_R = sim.x0_init[bus_ix]
        V_I = sim.x0_init[bus_ix + bus_size]
        Vm = sqrt(V_R^2 + V_I^2)
        θ = angle(V_R + V_I * 1im)
        print("Vm ", round(Vm; digits = 4), "\n")
        print("θ ", round(θ; digits = 4), "\n")
        println("====================")
    end
    println("====================")
    for device in get_dynamic_injectors(inputs)
        states = PSY.get_states(device)
        name = PSY.get_name(device)
        println("Differential States")
        println(name)
        println("====================")
        global_index = global_state_map[name]
        for s in states
            print(s, " ", round(sim.x0_init[global_index[s]]; digits = 4), "\n")
        end
        println("====================")
    end
    dyn_branches = get_dynamic_branches(inputs)
    if !isempty(dyn_branches)
        println("====================")
        println("Line Current States")
        for br in dyn_branches
            states = PSY.get_states(br)
            name = PSY.get_name(br)
            printed_name = "Line " * name
            println("====================")
            println(printed_name)
            global_index = global_state_map[name]
            x0_br = Dict{Symbol, Float64}()
            for (i, s) in enumerate(states)
                print(s, " ", round(sim.x0_init[global_index[s]]; digits = 5), "\n")
            end
            println("====================")
        end
    end
    return
end
