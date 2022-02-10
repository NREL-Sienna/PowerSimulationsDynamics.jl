struct NLsolveWrapper
    zero::Vector{Float64}
    converged::Bool
    failed::Bool
end

NLsolveWrapper() = NLsolveWrapper(Vector{Float64}(), false, true)
converged(sol::NLsolveWrapper) = sol.converged
failed(sol::NLsolveWrapper) = sol.failed

function _get_model_closure(model::SystemModel{MassMatrixModel}, ::Vector{Float64})
    return (residual, x) -> model(residual, x, nothing, 0.0)
end

function _get_model_closure(model::SystemModel{ResidualModel}, x0::Vector{Float64})
    dx0 = zeros(length(x0))
    return (residual, x) -> model(residual, dx0, x, nothing, 0.0)
end

function _nlsolve_call(
    initial_guess::Vector{Float64},
    f_eval::Function,
    jacobian::JacobianFunctionWrapper,
    f_tolerance::Float64,
    solver::Symbol,
)
    df = NLsolve.OnceDifferentiable(
        f_eval,
        jacobian,
        initial_guess,
        similar(initial_guess),
        jacobian.Jv,
    )
    sys_solve = NLsolve.nlsolve(
        df,
        initial_guess;
        xtol = NLSOLVE_X_TOLERANCE,
        ftol = f_tolerance,
        method = solver,
    ) # Solve using initial guess x0
    return NLsolveWrapper(sys_solve.zero, NLsolve.converged(sys_solve), false)
end

function _nlsolve_call(
    initial_guess::Vector{Float64},
    f_eval::Function,
    f_tolerance::Float64,
    solver::Symbol,
)
    sys_solve = NLsolve.nlsolve(
        f_eval,
        initial_guess;
        xtol = NLSOLVE_X_TOLERANCE,
        ftol = f_tolerance,
        method = solver,
    ) # Solve using initial guess x0
    return NLsolveWrapper(sys_solve.zero, NLsolve.converged(sys_solve), false)
end

function _convergence_check(sys_solve::NLsolveWrapper, tol::Float64, solv::Symbol)
    if converged(sys_solve)
        @info(
            "Initialization succeeded with a tolerance of $(tol) using solver $(solv). Saving solution."
        )
    else
        @warn(
            "Initialization convergence failed, initial conditions do not meet conditions for an stable equilibrium.\nTrying to solve again reducing numeric tolerance or using another solver"
        )
    end
    return converged(sys_solve)
end

function refine_initial_condition!(
    sim::Simulation,
    model::SystemModel,
    jacobian::JacobianFunctionWrapper,
)
    @assert sim.status != BUILD_INCOMPLETE

    if sim.status == SIMULATION_INITIALIZED
        @info "Simulation already initialized. Refinement not executed"
        return
    end
    converged = false
    initial_guess = get_initial_conditions(sim)
    inputs = get_simulation_inputs(sim)
    bus_range = get_bus_range(inputs)
    powerflow_solution = deepcopy(initial_guess[bus_range])
    # @debug "NLsolve initial guess $initial_guess"
    f! = _get_model_closure(model, initial_guess)
    ini_res = similar(initial_guess)
    f!(ini_res, initial_guess)
    val, ix = findmax(ini_res)
    @debug "NLsolve initial residual: max = $(val) at $ix, total = $(sum(ini_res))"
    if sum(ini_res) > MAX_INIT_RESIDUAL
        state_map = make_global_state_map(inputs)
        gen_name = ""
        state = ""
        for (gen, states) in state_map
            for (state_name, index) in states
                if index == ix
                    gen_name = gen
                    state = state_name
                end
            end
        end
        @warn "The initial residual in $ix of the NLsolve function has a value of $val.
               Generator = $gen_name, state = $state
               NLsolve might not converge."
    end
    for tol in [STRICT_NLSOLVE_F_TOLERANCE, RELAXED_NLSOLVE_F_TOLERANCE]
        if converged
            break
        end
        for solv in [:trust_region, :newton]
            @debug "Start NLSolve System Run with $(solv) and F_tol = $tol"
            sys_solve = _nlsolve_call(initial_guess, f!, jacobian, tol, solv)
            #sys_solve = _nlsolve_call(initial_guess, f!, tol, solv)
            failed(sys_solve) && return BUILD_FAILED
            converged = _convergence_check(sys_solve, tol, solv)
            @debug "Write initial guess vector using $solv with tol = $tol convergence = $converged"
            initial_guess .= sys_solve.zero
            if converged
                break
            end
        end
    end
    sim.status = check_valid_values(initial_guess, inputs)
    if sim.status == BUILD_FAILED
        error("Initial conditions refinement failed to find a valid initial condition")
    end

    pf_diff = abs.(powerflow_solution .- initial_guess[bus_range])
    if maximum(pf_diff) > MINIMAL_ACCEPTABLE_NLSOLVE_F_TOLERANCE
        @warn "The resulting voltages in the initial conditions differ from the power flow results"
    end

    if !converged
        @warn("Initialization didn't converged to desired tolerances.\\
              Initial conditions do not meet conditions for an stable equilibrium. \\
              Simulation might fail")
    end
    return
end
