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
    return (residual, x) -> model(residual, x, dx0, nothing, 0.0)
end

function _nlsolve_call(
    initial_guess::Vector{Float64},
    model::SystemModel,
    jacobian::JacobianFunctionWrapper,
    tolerance::Float64,
    solver::Symbol,
)
    df = NLsolve.OnceDifferentiable(
        _get_model_closure(model, initial_guess),
        jacobian.Jf,
        initial_guess,
        similar(initial_guess),
        jacobian.Jv,
    )
    return try
        sys_solve = NLsolve.nlsolve(
            df,
            initial_guess;
            xtol = tolerance,
            ftol = tolerance,
            method = solver,
        ) #Solve using initial guess x0
        NLsolveWrapper(sys_solve.zero, NLsolve.converged(sys_solve), false)
    catch e
        bt = catch_backtrace()
        @error "NLsolve failed to solve" exception = e, bt
        NLsolveWrapper()
    end
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
    simulation::Simulation,
    model::SystemModel,
    jacobian::JacobianFunctionWrapper,
)
    @debug "Start NLSolve System Run"
    converged = false
    initial_guess = get_initial_conditions(simulation)
    @debug initial_guess
    for tol in [STRICT_NL_SOLVE_TOLERANCE, RELAXED_NL_SOLVE_TOLERANCE]
        if converged
            break
        end
        for solv in [:trust_region, :newton]
            sys_solve = _nlsolve_call(initial_guess, model, jacobian, tol, solv)
            failed(sys_solve) && return BUILD_FAILED
            converged = _convergence_check(sys_solve, tol, solv)
            @debug "Write result to initial guess vector under condition converged = $(converged)"
            initial_guess .= sys_solve.zero
            if converged
                break
            end
        end
    end
    if !converged
        @warn(
            "Initialization never converged to desired tolerances. Initial conditions do not meet conditions for an stable equilibrium. Simulation might diverge"
        )
        simulation.status = check_valid_values(initial_guess, inputs)
    else
        simulation.status = SIMULATION_INITIALIZED
    end
    return
end
