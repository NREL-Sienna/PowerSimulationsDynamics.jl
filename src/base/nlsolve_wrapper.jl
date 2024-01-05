struct NLsolveWrapper
    zero::Vector{Float64}
    converged::Bool
    failed::Bool
end

NLsolveWrapper() = NLsolveWrapper(Vector{Float64}(), false, true)
converged(sol::NLsolveWrapper) = sol.converged
failed(sol::NLsolveWrapper) = sol.failed

function _get_model_closure(
    model::SystemModel{MassMatrixModel, NoDelays},
    ::Vector{Float64},
)
    return (residual, x) -> model(residual, x, nothing, 0.0)
end

function _get_model_closure(
    model::SystemModel{MassMatrixModel, HasDelays},
    x0::Vector{Float64},
)
    h(p, t; idxs = nothing) = typeof(idxs) <: Number ? x0[idxs] : x0
    return (residual, x) -> model(residual, x, h, nothing, 0.0)
end

function _get_model_closure(
    model::SystemModel{ResidualModel, NoDelays},
    x0::Vector{Float64},
)
    dx0 = zeros(length(x0))
    return (residual, x) -> model(residual, dx0, x, nothing, 0.0)
end

function _nlsolve_call(
    initial_guess::Vector{Float64},
    f_eval::Function,
    jacobian::JacobianFunctionWrapper,
    f_tolerance::Float64,
    solver::Symbol,
    show_trace::Bool,
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
        iterations = MAX_NLSOLVE_INTERATIONS,
        show_trace = show_trace,
    ) # Solve using initial guess x0
    return NLsolveWrapper(sys_solve.zero, NLsolve.converged(sys_solve), false)
end

function _nlsolve_call(
    initial_guess::Vector{Float64},
    f_eval::Function,
    f_tolerance::Float64,
    solver::Symbol,
    show_trace::Bool,
)
    sys_solve = NLsolve.nlsolve(
        f_eval,
        initial_guess;
        xtol = NLSOLVE_X_TOLERANCE,
        ftol = f_tolerance,
        iterations = MAX_NLSOLVE_INTERATIONS,
        method = solver,
        show_trace = show_trace,
    ) # Solve using initial guess x0
    return NLsolveWrapper(sys_solve.zero, NLsolve.converged(sys_solve), false)
end

function _convergence_check(sys_solve::NLsolveWrapper, tol::Float64, solv::Symbol)
    if converged(sys_solve)
        @info(
            "Initialization non-linear solve succeeded with a tolerance of $(tol) using solver $(solv). Saving solution."
        )
    else
        @warn(
            "Initialization non-linear solve convergence failed, initial conditions do not meet conditions for an stable equilibrium.\nAttempting again with reduced numeric tolerance and using another solver"
        )
    end
    return converged(sys_solve)
end

function _sorted_residuals(residual::Vector{Float64})
    if isapprox(sum(abs.(residual)), 0.0; atol = STRICT_NLSOLVE_F_TOLERANCE)
        @debug "Residual is zero with tolerance $(STRICT_NLSOLVE_F_TOLERANCE)"
        return
    end
    ix_sorted = sortperm(abs.(residual); rev = true)
    show_residual = min(10, length(residual))
    for i in 1:show_residual
        ix = ix_sorted[i]
        @debug ix abs(residual[ix])
    end
    return
end

function _check_residual(
    residual::Vector{Float64},
    inputs::SimulationInputs,
    tolerance::Float64,
)
    @debug _sorted_residuals(residual)
    val, ix = findmax(residual)
    sum_residual = sum(abs.(residual))
    @info "Residual from initial guess: max = $(val) at $ix, total = $sum_residual"
    if sum_residual > tolerance
        state_map = make_global_state_map(inputs)
        for (k, val) in state_map
            inputs.global_state_map[k] = val
        end
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
        if gen_name != ""
            error("The initial residual in the state located at $ix has a value of $val.
                Generator = $gen_name, state = $state.
               Residual error is too large to continue")
        else
            bus_no = ix > bus_count ? ix - bus_count : ix
            component = ix > bus_count ? "imag" : "real"
            error("The initial residual in the state located at $ix has a value of $val.
                Voltage at bus = $bus_no, component = $component.
                Error is too large to continue")
        end
    end
    return
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
    f! = _get_model_closure(model, initial_guess)
    residual = similar(initial_guess)
    f!(residual, initial_guess)
    _check_residual(residual, inputs, MAX_INIT_RESIDUAL)
    for tol in [STRICT_NLSOLVE_F_TOLERANCE, RELAXED_NLSOLVE_F_TOLERANCE]
        if converged
            break
        end
        for solv in [:trust_region, :newton]
            @debug "Start NLSolve System Run with $(solv) and F_tol = $tol"
            show_trace = sim.console_level <= Logging.Info
            sys_solve = _nlsolve_call(initial_guess, f!, jacobian, tol, solv, show_trace)
            #sys_solve = _nlsolve_call(initial_guess, f!, tol, solv, show_trace)
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
        error(
            "Initial conditions refinement failed to find a valid initial condition. Run show_states_initial_value on your simulation",
        )
    end

    f!(residual, initial_guess)
    if !converged || (sum(residual) > MINIMAL_ACCEPTABLE_NLSOLVE_F_TOLERANCE)
        _check_residual(residual, inputs, MINIMAL_ACCEPTABLE_NLSOLVE_F_TOLERANCE)
        @warn("Initialization didn't found a solution to desired tolerances.\\
              Initial conditions do not meet conditions for an stable equilibrium. \\
              Simulation might fail")
    end

    pf_diff = abs.(powerflow_solution .- initial_guess[bus_range])
    if maximum(pf_diff) > MINIMAL_ACCEPTABLE_NLSOLVE_F_TOLERANCE
        @warn "The resulting voltages in the initial conditions differ from the power flow results"
    end

    return
end
