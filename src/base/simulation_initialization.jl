function calculate_initial_conditions!(sim::Simulation, inputs::SimulationInputs)
    @debug "Start state intialization routine"
    sys = get_system(sim)
    initial_guess = sim.x0_init
    res = PSY.solve_powerflow!(sys)
    if !res
        sim.status = BUILD_FAILED
        @error("PowerFlow failed to solve")
        return false
    end
    var_count = get_variable_count(inputs)

    @debug "Setting up initilization indexing"
    bus_size = get_bus_count(inputs)
    injection_start = get_injection_pointer(inputs)
    branches_start = get_branches_pointer(inputs)

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
            initialize_device!(s)
        end
    end

    @debug "Updating Dynamic Injection Component Initial Guess"
    for d in get_injectors_data(inputs)
        dynamic_device = PSY.get_dynamic_injector(d)
        @debug PSY.get_name(d) typeof(d) typeof(dynamic_device)
        n_states = PSY.get_n_states(dynamic_device)
        ix_range = range(injection_start, length = n_states)
        injection_start = injection_start + n_states
        x0_device = initialize_device!(d)
        @assert length(x0_device) == n_states
        initial_guess[ix_range] = x0_device
    end

    @debug "Updating Component Initial Guess"
    dyn_branches = PSY.get_components(PSY.DynamicBranch, sys)
    if !isempty(dyn_branches)
        for br in dyn_branches
            @debug PSY.get_name(br) typeof(br)
            n_states = PSY.get_n_states(br)
            ix_range = range(branches_start, length = n_states)
            branches_start = branches_start + n_states
            x0_branch = initialize_device!(br)
            @assert length(x0_branch) == n_states
            initial_guess[ix_range] = x0_branch
        end
    else
        @debug "No Dynamic Branches in the system"
    end

    if any(!isfinite, initial_guess)
        i = findall(!isfinite, initial_guess)
        invalid_initial_guess = String[]
        for (device, states) in get_global_index(sim)
            for state in states
                if state.second ∈ i
                    push!(invalid_initial_guess, "$device - $(p.first)")
                end
            end
        end
        error("Invalid initial guess values $invalid_initial_guess")
        sim.status = BUILD_FAILED
        return false
    end

    dx0 = zeros(var_count) #Define a vector of zeros for the derivative
    inif! = (out, x) -> system_implicit!(
        out,    #output of the function
        dx0,    #derivatives equal to zero
        x,      #states
        inputs,    #Parameters
        0.0,    #time equals to zero.
    )
    #Refine initial solution
    @debug "Start NLSolve System Run"
    @debug initial_guess
    try
        sys_solve = NLsolve.nlsolve(
            inif!,
            initial_guess,
            xtol = STRICT_NL_SOLVE_TOLERANCE,
            ftol = STRICT_NL_SOLVE_TOLERANCE,
            method = :trust_region,
        ) #Solve using initial guess x0
    catch e
        bt = catch_backtrace()
        @error "NLsolve failed to solve" exception = e, bt
        sim.status = BUILD_FAILED
        return false
    end
    if !NLsolve.converged(sys_solve)
        @warn(
            "Initialization failed, initial conditions do not meet conditions for an stable equilibrium.\nTrying to solve again reducing numeric tolerance from $(STRICT_NL_SOLVE_TOLERANCE):"
        )
        sys_solve = NLsolve.nlsolve(
            inif!,
            initial_guess,
            xtol = RELAXED_NL_SOLVE_TOLERANCE,
            ftol = RELAXED_NL_SOLVE_TOLERANCE,
            method = :trust_region,
        ) #Solve using initial guess x0
        if NLsolve.converged(sys_solve)
            @info(
                "Initialization succeeded with a relaxed tolerance of $(RELAXED_NL_SOLVE_TOLERANCE). Saving solution"
            )
        else
            @warn(
                "Initialization failed again. Initial conditions do not meet conditions for an stable equilibrium\nSaving best result."
            )
        end
    end
    @debug "Write result to initial guess vector"
    initial_guess .= sys_solve.zero
    return NLsolve.converged(sys_solve)
end
