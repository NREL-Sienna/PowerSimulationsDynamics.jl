function calculate_initial_conditions!(sys::PSY.System, initial_guess::Vector{Float64})

    bus_count = length(PSY.get_components(PSY.Bus, sys))
    var_count = get_variable_count(sys)
    dx0 = zeros(var_count) #Define a vector of zeros for the derivative

    #These lines stay
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
