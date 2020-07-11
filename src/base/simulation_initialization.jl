function calculate_initial_conditions!(sys::PSY.System, initial_guess::Vector{Float64})

    res = PSY.solve_powerflow!(sys)
    if !res
        error("Power Flow fail")
    end
    var_count = get_variable_count(sys)

    #Index Setup
    bus_size = get_bus_count(sys)
    bus_vars_count = 2 * bus_size
    bus_range = 1:bus_vars_count
    injection_start = get_injection_pointer(sys)
    injection_count = 1
    branches_start = get_branches_pointer(sys)
    branches_count = 1

    for bus in PSY.get_components(PSY.Bus, sys)
        #Write voltage initial guess
        bus_n = PSY.get_number(bus)
        bus_ix = PSY.get_ext(sys)[LOOKUP][bus_n]
        initial_guess[bus_ix] = PSY.get_voltage(bus) * cos(PSY.get_angle(bus))
        initial_guess[bus_ix + bus_size] = PSY.get_voltage(bus) * sin(PSY.get_angle(bus))
    end

    for d in PSY.get_components(PSY.DynamicInjection, sys)
        bus = PSY.get_bus(d)
        bus_n = PSY.get_number(PSY.get_bus(d))
        bus_ix = PSY.get_ext(sys)[LOOKUP][bus_n]
        n_states = PSY.get_n_states(d)
        ix_range = range(injection_start, length = n_states)
        injection_start = injection_start + n_states
        x0_device = initialize_device(d)
        @assert length(x0_device) == n_states
        initial_guess[ix_range] = x0_device
    end

    dx0 = zeros(var_count) #Define a vector of zeros for the derivative
    inif! = (out, x) -> system!(
        out,    #output of the function
        dx0,    #derivatives equal to zero
        x,      #states
        sys,    #Parameters
        0.0,    #time equals to zero.
    )
    #Refine initial solution
    sys_solve = NLsolve.nlsolve(
        inif!,
        initial_guess,
        xtol = :1e-9,
        ftol = :1e-9,
        method = :trust_region,
    ) #Solve using initial guess x0
    if !NLsolve.converged(sys_solve)
        @warn("Initialization failed, initial conditions do not meet conditions for an stable equilibrium.\nTrying to solve again reducing numeric tolerance:")
        sys_solve = NLsolve.nlsolve(
            inif!,
            initial_guess,
            xtol = :1e-6,
            ftol = :1e-6,
            method = :trust_region,
        ) #Solve using initial guess x0
        if NLsolve.converged(sys_solve)
            @info("Initialization succeeded with less tolerance. Saving solution")
        else
            @warn("Initialization failed again. Initial conditions do not meet conditions for an stable equilibrium\nSaving best result.")
        end
    end
    initial_guess .= sys_solve.zero
    return NLsolve.converged(sys_solve)

end
