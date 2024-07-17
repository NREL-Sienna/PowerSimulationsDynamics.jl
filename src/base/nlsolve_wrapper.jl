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
    p::AbstractArray{Float64},
)
    return (residual, x, p) -> model(residual, x, p, 0.0)
end

function _get_model_closure(
    model::SystemModel{MassMatrixModel, HasDelays},
    x0::Vector{Float64},
    p::AbstractArray{Float64},
)
    h(p, t; idxs = nothing) = typeof(idxs) <: Number ? x0[idxs] : x0
    return (residual, x, p) -> model(residual, x, h, p, 0.0)
end

function _get_model_closure(
    model::SystemModel{ResidualModel, NoDelays},
    x0::Vector{Float64},
    p::AbstractArray{Float64},
)
    dx0 = zeros(length(x0))
    return (residual, x, p) -> model(residual, dx0, x, p, 0.0)
end

function _nlsolve_call(
    initial_guess::Vector{Float64},
    p::AbstractArray,
    f_eval::Function,
    jacobian::JacobianFunctionWrapper,
    f_tolerance::Float64,
    solver::NonlinearSolve.AbstractNonlinearSolveAlgorithm,
    show_trace::Bool,
)
    f = SciMLBase.NonlinearFunction(f_eval; jac = jacobian)
    prob = NonlinearSolve.NonlinearProblem{true}(f, initial_guess, p)
    sol = NonlinearSolve.solve(
        prob,
        solver;
        sensealg = SciMLSensitivity.SteadyStateAdjoint(),
        abstol = f_tolerance,
        reltol = f_tolerance,
        maxiters = MAX_NLSOLVE_INTERATIONS,
        show_trace = Val(show_trace),
    )
    return NLsolveWrapper(sol.u, SciMLBase.successful_retcode(sol), false)
end

function _convergence_check(
    sys_solve::NLsolveWrapper,
    tol::Float64,
    solv::Symbol,
)
    if converged(sys_solve)
        @warn(
            "Initialization non-linear solve succeeded with a tolerance of $(tol) using solver $(solv). Saving solution."
        )
    else
        @warn(
            "Initialization non-linear solve convergence failed with a tolerance of $(tol) using solver $(solv), initial conditions do not meet conditions for an stable equilibrium.\nAttempting again with reduced numeric tolerance and using another solver"
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
            get_global_state_map(inputs)[k] = val
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
            bus_count = get_bus_count(inputs)
            bus_no = ix > bus_count ? ix - bus_count : ix
            component = ix > bus_count ? "imag" : "real"
            error("The initial residual in the state located at $ix has a value of $val.
                Voltage at bus = $bus_no, component = $component.
                Error is too large to continue")
        end
    end
    return
end

function _refine_initial_condition!(x0, p, prob)
    h(p, t; idxs = nothing) = typeof(idxs) <: Number ? x0[idxs] : x0
    function ff(u, x0, p)
        if typeof(prob) <: SciMLBase.ODEProblem
            prob.f.f(u, x0, p, 0.0)
        elseif typeof(prob) <: SciMLBase.DDEProblem
            prob.f.f(u, x0, h, p, 0.0)
        end
        return
    end
    #residual = similar(x0)
    #ff(residual, x0, p) #Error: ERROR: AssertionError: length(getcolptr(S)) == size(S, 2) + 1 && (getcolptr(S))[end] - 1 == length(rowvals(S)) == length(nonzeros(S))
    #_check_residual(residual, inputs, MAX_INIT_RESIDUAL)
    solver = NonlinearSolve.TrustRegion()
    probnl = NonlinearSolve.NonlinearProblem{true}(ff, x0, p)
    #for tol in [STRICT_NLSOLVE_F_TOLERANCE, RELAXED_NLSOLVE_F_TOLERANCE] #ERROR: Enzyme execution failed., Enzyme: Non-constant keyword argument found for Tuple{UInt64, typeof(Core.kwcall), Duplicated{@NamedTuple{reltol::Float64, abstol::Float64}},
    for solver in [NonlinearSolve.TrustRegion(), NonlinearSolve.NewtonRaphson()]
        sol = NonlinearSolve.solve(
            probnl,
            solver;
            sensealg = SciMLSensitivity.SteadyStateAdjoint(),
            reltol = STRICT_NLSOLVE_F_TOLERANCE,
            abstol = STRICT_NLSOLVE_F_TOLERANCE,
            maxiters = MAX_NLSOLVE_INTERATIONS,
        )
        converged = SciMLBase.successful_retcode(sol)
        x0 .= sol.u
        if converged
            break
        end
    end
    #end 
    return nothing
end

function refine_initial_condition!(
    sim::Simulation,
    model::SystemModel,
    jacobian::JacobianFunctionWrapper,
)
    @assert sim.status != BUILD_INCOMPLETE
    converged = false
    initial_guess = get_x0(sim)
    inputs = get_simulation_inputs(sim)
    parameters = get_parameters(inputs)
    bus_range = get_bus_range(inputs)
    powerflow_solution = deepcopy(initial_guess[bus_range])
    f! = _get_model_closure(model, initial_guess, parameters)
    residual = similar(initial_guess)
    f!(residual, initial_guess, parameters)
    _check_residual(residual, inputs, MAX_INIT_RESIDUAL)
    for tol in [STRICT_NLSOLVE_F_TOLERANCE, RELAXED_NLSOLVE_F_TOLERANCE]
        if converged
            break
        end
        for solv in [:trust_region, :newton]
            @debug "Start NLSolve System Run with $(solv) and F_tol = $tol"
            show_trace = sim.console_level <= Logging.Info
            sys_solve = _nlsolve_call(
                initial_guess,
                parameters,
                f!,
                jacobian,
                tol,
                NonlinearSolve.NLsolveJL(; method = solv),
                show_trace,
            )
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

    f!(residual, initial_guess, parameters)
    if !converged || (sum(residual) > MINIMAL_ACCEPTABLE_NLSOLVE_F_TOLERANCE)
        _check_residual(residual, inputs, MINIMAL_ACCEPTABLE_NLSOLVE_F_TOLERANCE)
        @warn(
            "Initialization didn't found a solution to desired tolerances.\\
            Initial conditions do not meet conditions for an stable equilibrium. \\
            Simulation might fail"
        )
    end

    pf_diff = abs.(powerflow_solution .- initial_guess[bus_range])
    if maximum(pf_diff) > MINIMAL_ACCEPTABLE_NLSOLVE_F_TOLERANCE
        @warn "The resulting voltages in the initial conditions differ from the power flow results"
    end
    return
end
