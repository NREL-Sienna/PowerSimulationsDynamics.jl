function get_indices_in_parameter_vector(p, device_param_pairs)
    indices = Int[]
    for tuple in device_param_pairs
        label = join((tuple[1], "params", tuple[2:end]...), ".")
        ix = ComponentArrays.label2index(p, label)[1]
        if ix === nothing
            @error "Index not found for entry $tuple"
            return nothing
        end
        push!(indices, ix)
    end
    return indices
end

function get_required_initialization_level(sys, device_param_pairs)
    init_level = INITIALIZED
    #TODO - check parameters and return appropriate init_level 
    return init_level
end

function get_indices_in_state_vector(sim, state_data)
    @assert sim.results !== nothing
    res = sim.results
    global_state_index = get_global_index(res)
    state_ixs = Vector{Int64}(undef, length(state_data))
    for (ix, ref) in enumerate(state_data)
        if !haskey(global_state_index, ref[1])
            @error "$(keys(global_state_index))"
            error("State $(ref[2]) device $(ref[1]) not found in the system. ")
        end
        state_ixs[ix] = get(global_state_index[ref[1]], ref[2], 0)
    end
    return state_ixs
end

function get_parameter_values(sim, device_param_pairs)
    p = sim.inputs.parameters
    ixs = get_indices_in_parameter_vector(p, device_param_pairs)
    if ixs === nothing
        return nothing
    else
        return p[ixs]
    end
end

#TODO - think more carefully about how data should be included so it works with with Optimization API. 
#TODO - avoid code repetition with _execute! by defining functions appropriately.
#TODO - extend for initialization - first for H...
#TODO - extend for parameters that require initialization.
function get_sensitivity_functions(sim, param_data, state_data, solver, f_loss; kwargs...)
    init_level = get_required_initialization_level(sim.sys, param_data)
    p = sim.inputs.parameters
    param_ixs = get_indices_in_parameter_vector(p, param_data)
    state_ixs = get_indices_in_state_vector(sim, state_data)
    init_level = get_required_initialization_level(sim.sys, param_data)
    sim_inputs = deepcopy(sim.inputs_init)
    if get(kwargs, :auto_abstol, false)
        cb = AutoAbstol(true, get(kwargs, :abstol, 1e-9))
        callbacks =
            SciMLBase.CallbackSet((), tuple(push!(sim.callbacks, cb)...))
    else
        callbacks = SciMLBase.CallbackSet((), tuple(sim.callbacks...))
    end
    tstops = if !isempty(sim.tstops)
        [sim.tstops[1] ÷ 2, sim.tstops...]  #Note: Don't need to make a tuple because it is in ODEproblem. not a kwarg
    else
        []
    end

    prob_old = sim.problem
    f_old = prob_old.f
    f_new = SciMLBase.ODEFunction{true}(
        f_old.f;
        mass_matrix = f_old.mass_matrix,
    )

    prob_new = SciMLBase.remake(
        prob_old;
        f = f_new,
        tstops = tstops,
        advance_to_tstop = !isempty(tstops),
        initializealg = SciMLBase.NoInit(),
        callback = callbacks,
        kwargs...,
    )
    if param_ixs === nothing
        return nothing
    else
        x0 = deepcopy(sim.x0_init)
        sys = deepcopy(sim.sys)
        function f_enzyme(p, x0, sys, sim_inputs, prob, data, init_level)   #Make separate f_enzymes depending on init_level? 
            p_new = sim_inputs.parameters
            p_new[param_ixs] .= p
            if init_level == POWERFLOW_AND_DEVICES
                @error "POWERFLOW AND DEVICES -- not yet supported"
                #_initialize_powerflow_and_devices!(x0, inputs, sys)
            elseif init_level == DEVICES_ONLY
                @error "DEVICES ONLY"
            elseif init_level == INITIALIZED
                _initialize_devices_only!(x0, sim_inputs)
            end
            prob_new = SciMLBase.remake(prob; p = p_new, u0 = x0)
            sol = SciMLBase.solve(prob_new, solver)

            @assert length(state_ixs) == 1  #Hardcode for single state for now 
            ix = state_ixs[1]
            ix_t = unique(i -> sol.t[i], eachindex(sol.t))
            state = sol[ix, ix_t]

            return f_loss(state, data)
        end
        function f_forward(p, data)
            f_enzyme(
                p,
                x0,
                sys,
                deepcopy(sim.inputs_init),
                prob_new,
                data,
                init_level,
            )
        end
        function f_grad(p, data)
            dp = Enzyme.make_zero(p)
            dx0 = Enzyme.make_zero(x0)
            dsys = Enzyme.make_zero(sys)
            sim_inputs = deepcopy(sim.inputs_init)
            dsim_inputs = Enzyme.make_zero(sim_inputs)
            dprob_new = Enzyme.make_zero(prob_new)
            ddata = Enzyme.make_zero(data)
            Enzyme.autodiff(
                Enzyme.Reverse,
                f_enzyme,
                Enzyme.Active,
                Enzyme.Duplicated(p, dp),
                Enzyme.Duplicated(x0, dx0),
                Enzyme.Duplicated(sys, dsys),
                Enzyme.Duplicated(sim_inputs, dsim_inputs),
                Enzyme.Duplicated(prob_new, dprob_new),
                Enzyme.Duplicated(data, ddata),
                Enzyme.Const(init_level),
            )
            return dp
        end
        f_forward, f_grad
    end
end

#Inactive Rules
#Enzyme.EnzymeRules.inactive(::typeof(SimulationResults), args...) = nothing
#Enzyme.EnzymeRules.inactive(::typeof(get_activepower_branch_flow), args...) = nothing
