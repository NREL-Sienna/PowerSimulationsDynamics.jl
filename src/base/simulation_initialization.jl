function calculate_initial_conditions!(sys::PSY.System, initial_guess::Vector{Float64})

    res = PSY.solve_powerflow!(sys)
    if !res
        error("Power Flow fail")
    end
    var_count = get_variable_count(sys)
   #=
    initial_guess = zeros(var_count)

    #Index Setup
    bus_size = get_bus_count(sys)
    bus_vars_count = 2 * bus_size
    bus_range = 1:bus_vars_count
    injection_start = get_injection_pointer(sys)
    injection_count = 1
    branches_start = get_branches_pointer(sys)
    branches_count = 1

    for d in PSY.get_components(PSY.DynamicInjection, sys)
        bus = PSY.get_bus(d)
        bus_n = PSY.get_number(PSY.get_bus(d)) # TODO: This requires that the bus numbers are indexed 1-N
        n_states = PSY.get_n_states(d)
        ix_range = range(injection_start, length = n_states)
        #injection_start = injection_start + n_states
        #Write Buses Voltages Initial Guess
        initial_guess[bus_n] = PSY.get_voltage(bus)*cos(PSY.get_angle(bus))
        initial_guess[bus_n + bus_size] = PSY.get_voltage(bus)*sin(PSY.get_angle(bus))

        dev = (out, x) ->
        LITS.initialize_device!(
            out,
            x,
            initial_guess[bus_n],
            initial_guess[bus_n + bus_size],
            d,
            sys,
        )
        dev_solve = NLsolve.nlsolve(dev,
                   [1.0, #eq_p
                    0.47, #ed_p
                    0.6, #δ
                    1.0, #ω
                    2.1, #Vf
                    0.28, #Vr1
                    -0.39, #Vr2,
                    1.0, #Vm
                    ],
    )
    end


=#
    dx0 = zeros(var_count) #Define a vector of zeros for the derivative

    inif! = (out, x) -> system!(
        out,    #output of the function
        dx0,    #derivatives equal to zero
        x,      #states
        sys,    #Parameters
        0.0,    #time equals to zero.
    )




    sys_solve = NLsolve.nlsolve(
        inif!,
        initial_guess,
        xtol = :1e-9,
        ftol = :1e-9,
        method = :trust_region,
    ) #Solve using initial guess x0
    if !NLsolve.converged(sys_solve)
        @warn("Initialization failed, initial conditions do not meet conditions for an stable equilibrium")
    end
    initial_guess .= sys_solve.zero
    return NLsolve.converged(sys_solve)

end
