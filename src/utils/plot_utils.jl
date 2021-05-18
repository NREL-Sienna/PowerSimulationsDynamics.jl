"""
Function to obtain series of states out of DAE Solution. It receives the solution, the dynamical system
and a tuple containing the string name of the Dynamic Injection device and the symbol of the state.

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
        error("Bus number $(bus_number) not found.")
    else
        return sim.solution.t,
        [sqrt(value[bus_ix]^2 + value[bus_ix + n_buses]^2) for value in sim.solution.u]
    end
end

"""
Function to obtain the active power output time series of a Dynamic Injection series out of the DAE Solution. It receives the solution and the
string name of the Dynamic Injection device.

"""
function get_activepower_series(sim::Simulation, name::String)
    sim_inputs = sim.simulation_inputs
    n_buses = get_bus_count(sim_inputs)
    solution = sim.solution
    device = PSY.get_component(PSY.StaticInjection, sim_inputs.sys, name)
    bus_ix = get(get_lookup(sim_inputs), PSY.get_number(PSY.get_bus(device)), nothing)
    V_R = [value[bus_ix] for value in solution.u]
    V_I = [value[bus_ix + n_buses] for value in solution.u]
    dyn_device = PSY.get_dynamic_injector(device)
    if typeof(dyn_device) <: PSY.DynamicGenerator
        I_R, I_I = compute_output_current(sim, dyn_device)
        return solution.t, V_R .* I_R + V_I .* I_I
    end
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

"""
Function to obtain the output current time series of a Base Machine model out of the DAE Solution. It receives the simulation inputs and
the dynamic device. It is dispatched for device type to compute the specific current.

"""
function compute_output_current(
    sim::Simulation,
    dynamic_device::PSY.DynamicGenerator{PSY.BaseMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain Data
    sys = sim.simulation_inputs.sys
    static_device =
        PSY.get_component(PSY.StaticInjection, sys, PSY.get_name(dynamic_device))
    n_buses = get_bus_count(sim.simulation_inputs)
    bus_ix = get(
        get_lookup(sim.simulation_inputs),
        PSY.get_number(PSY.get_bus(static_device)),
        nothing,
    )
    _, δ = get_state_series(sim, (PSY.get_name(dynamic_device), :δ))
    sim_length = length(δ)
    V_R = [value[bus_ix] for value in sim.solution.u]
    V_I = [value[bus_ix + n_buses] for value in sim.solution.u]

    #Get parameters
    machine = PSY.get_machine(dynamic_device)
    R = PSY.get_R(machine)
    Xd_p = PSY.get_Xd_p(machine)
    eq_p = PSY.get_eq_p(machine)
    basepower = PSY.get_base_power(dynamic_device)
    Sbase = PSY.get_base_power(sys)

    #RI to dq transformation
    V_dq = reduce(hcat, [ri_dq.(δ[i]) * [V_R[i]; V_I[i]] for i in 1:sim_length])
    V_d = @view V_dq[1, :]
    V_q = @view V_dq[2, :]

    #Obtain electric current
    i_d = (1.0 / (R^2 + Xd_p^2)) * (Xd_p * (eq_p .- V_q) - R * V_d)  #15.36
    i_q = (1.0 / (R^2 + Xd_p^2)) * (Xd_p * V_d + R * (eq_p .- V_q)) #15.36

    #Obtain output current
    I_RI = reduce(
        hcat,
        [(basepower / Sbase) * dq_ri(δ[i]) * [i_d[i]; i_q[i]] for i in 1:sim_length],
    )
    return I_RI[1, :], I_RI[2, :]
end

"""
Function to obtain the output current time series of a One-D-One-Q model out of the DAE Solution. It receives the simulation inputs and
the dynamic device. It is dispatched for device type to compute the specific current.

"""
function compute_output_current(
    sim::Simulation,
    dynamic_device::PSY.DynamicGenerator{PSY.OneDOneQMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain Data
    sys = sim.simulation_inputs.sys
    static_device =
        PSY.get_component(PSY.StaticInjection, sys, PSY.get_name(dynamic_device))
    n_buses = get_bus_count(sim.simulation_inputs)
    bus_ix = get(
        get_lookup(sim.simulation_inputs),
        PSY.get_number(PSY.get_bus(static_device)),
        nothing,
    )
    _, δ = get_state_series(sim, (PSY.get_name(dynamic_device), :δ))
    _, eq_p = get_state_series(sim, (PSY.get_name(dynamic_device), :eq_p))
    _, ed_p = get_state_series(sim, (PSY.get_name(dynamic_device), :ed_p))
    sim_length = length(δ)
    V_R = [value[bus_ix] for value in sim.solution.u]
    V_I = [value[bus_ix + n_buses] for value in sim.solution.u]

    #Get parameters
    #Get parameters
    machine = PSY.get_machine(dynamic_device)
    R = PSY.get_R(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    basepower = PSY.get_base_power(dynamic_device)
    Sbase = PSY.get_base_power(sys)

    #RI to dq transformation
    V_dq = reduce(hcat, [ri_dq.(δ[i]) * [V_R[i]; V_I[i]] for i in 1:sim_length])
    V_d = @view V_dq[1, :]
    V_q = @view V_dq[2, :]

    #Obtain electric current
    i_d = (1.0 / (R^2 + Xd_p * Xq_p)) * (Xq_p * (eq_p - V_q) + R * (ed_p - V_d))  #15.32
    i_q = (1.0 / (R^2 + Xd_p * Xq_p)) * (-Xd_p * (ed_p - V_d) + R * (eq_p - V_q))  #15.32

    #Obtain output current
    I_RI = reduce(
        hcat,
        [(basepower / Sbase) * dq_ri(δ[i]) * [i_d[i]; i_q[i]] for i in 1:sim_length],
    )
    return I_RI[1, :], I_RI[2, :]
end
