function calculate_initial_conditions!(sim::Simulation, inputs::SimulationInputs)
    @debug "Start state intialization routine"
    sys = get_system(inputs)
    initial_guess = sim.x0_init
    res = PSY.solve_powerflow!(sys)
    if !res
        error("Power Flow fail")
    end
    var_count = get_variable_count(inputs)

    @debug "Setting up initilization indexing"
    bus_size = get_bus_count(inputs)
    bus_vars_count = 2 * bus_size
    bus_range = 1:bus_vars_count
    injection_start = get_injection_pointer(inputs)
    injection_count = 1
    branches_start = get_branches_pointer(inputs)
    branches_count = 1

    @debug "Updating Voltage guess"
    for bus in PSY.get_components(PSY.Bus, sys)
        @debug PSY.get_name(bus)
        @debug "V_r" PSY.get_magnitude(bus) * cos(PSY.get_angle(bus))
        @debug "V_i" PSY.get_magnitude(bus) * sin(PSY.get_angle(bus))
        bus_n = PSY.get_number(bus)
        bus_ix = get_lookup(inputs)[bus_n]
        initial_guess[bus_ix] = PSY.get_magnitude(bus) * cos(PSY.get_angle(bus))
        initial_guess[bus_ix + bus_size] = PSY.get_magnitude(bus) * sin(PSY.get_angle(bus))
    end

    @debug "Updating Source internal voltages"
    sources = PSY.get_components(PSY.Source, sys)
    if !isempty(sources)
        for s in sources
            initialize_device(s)
        end
    end

    @debug "Updating Dynamic Injection Component Initial Guess"
    for d in PSY.get_components(PSY.DynamicInjection, sys)
        @debug PSY.get_name(d) typeof(d)
        bus = PSY.get_bus(d)
        bus_n = PSY.get_number(PSY.get_bus(d))
        bus_ix = get_lookup(inputs)[bus_n]
        n_states = PSY.get_n_states(d)
        ix_range = range(injection_start, length = n_states)
        injection_start = injection_start + n_states
        x0_device = initialize_device(d)
        @assert length(x0_device) == n_states
        initial_guess[ix_range] = x0_device
    end

    @debug "Updating Component Initial Guess"
    dyn_branches = PSY.get_components(PSY.DynamicBranch, sys)
    if !isempty(dyn_branches)
        for br in dyn_branches
            @debug PSY.get_name(br) typeof(br)
            arc = PSY.get_arc(br)
            n_states = PSY.get_n_states(br)
            from_bus_number = PSY.get_number(arc.from)
            to_bus_number = PSY.get_number(arc.to)
            bus_ix_from = get_lookup(inputs)[from_bus_number]
            bus_ix_to = get_lookup(inputs)[to_bus_number]
            ix_range = range(branches_start, length = n_states)
            branches_start = branches_start + n_states
            x0_branch = initialize_device(br)
            @assert length(x0_branch) == n_states
            initial_guess[ix_range] = x0_branch
        end
    else
        @debug "No Dynamic Branches in the system"
    end

    dx0 = zeros(var_count) #Define a vector of zeros for the derivative
    inif! = (out, x) -> system!(
        out,    #output of the function
        dx0,    #derivatives equal to zero
        x,      #states
        inputs,    #Parameters
        0.0,    #time equals to zero.
    )
    #Refine initial solution
    @debug "Start NLSolve Run"
    sys_solve = NLsolve.nlsolve(
        inif!,
        initial_guess,
        xtol = STRICT_NL_SOLVE_TOLERANCE,
        ftol = STRICT_NL_SOLVE_TOLERANCE,
        method = :trust_region,
    ) #Solve using initial guess x0
    if !NLsolve.converged(sys_solve)
        @warn("Initialization failed, initial conditions do not meet conditions for an stable equilibrium.\nTrying to solve again reducing numeric tolerance from $(STRICT_NL_SOLVE_TOLERANCE):")
        sys_solve = NLsolve.nlsolve(
            inif!,
            initial_guess,
            xtol = RELAXED_NL_SOLVE_TOLERANCE,
            ftol = RELAXED_NL_SOLVE_TOLERANCE,
            method = :trust_region,
        ) #Solve using initial guess x0
        if NLsolve.converged(sys_solve)
            @info("Initialization succeeded with a relaxed tolerance of $(RELAXED_NL_SOLVE_TOLERANCE). Saving solution")
        else
            @warn("Initialization failed again. Initial conditions do not meet conditions for an stable equilibrium\nSaving best result.")
        end
    end
    @debug "Write result to initial guess vector"
    initial_guess .= sys_solve.zero
    return NLsolve.converged(sys_solve)

end

"""
Returns a Dictionary with the resulting initial conditions of the simulation
"""
function get_initial_conditions(sim::Simulation)
    system = get_system(sim.simulation_inputs)
    bus_size = get_bus_count(sim.simulation_inputs)
    V_R = Dict{Int64, Float64}()
    V_I = Dict{Int64, Float64}()
    Vm = Dict{Int64, Float64}()
    θ = Dict{Int64, Float64}()
    for bus in PSY.get_components(PSY.Bus, system)
        bus_n = PSY.get_number(bus)
        bus_ix = get_lookup(sim.simulation_inputs)[bus_n]
        V_R[bus_n] = sim.x0_init[bus_ix]
        V_I[bus_n] = sim.x0_init[bus_ix + bus_size]
        Vm[bus_n] = sqrt(sim.x0_init[bus_ix]^2 + sim.x0_init[bus_ix + bus_size]^2)
        θ[bus_n] = angle(sim.x0_init[bus_ix] + sim.x0_init[bus_ix + bus_size] * 1im)
    end
    results = Dict{String, Any}("V_R" => V_R, "V_I" => V_I, "Vm" => Vm, "θ" => θ)
    for device in PSY.get_components(PSY.DynamicInjection, system)
        states = PSY.get_states(device)
        name = PSY.get_name(device)
        global_index = get_global_index(sim.simulation_inputs)[name]
        x0_device = Dict{Symbol, Float64}()
        for (i, s) in enumerate(states)
            x0_device[s] = sim.x0_init[global_index[s]]
        end
        results[name] = x0_device
    end
    dyn_branches = PSY.get_components(PSY.DynamicBranch, system)
    if !isempty(dyn_branches)
        for br in dyn_branches
            states = PSY.get_states(br)
            name = PSY.get_name(br)
            global_index = get_global_index(sim.simulation_inputs)[name]
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
