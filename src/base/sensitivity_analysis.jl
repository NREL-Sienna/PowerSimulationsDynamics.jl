function get_indices_in_parameter_vector(p, device_param_pairs)
    indices = Int[]
    for tuple in device_param_pairs
        label = join((tuple[1], "params", tuple[2:end]...), ".")
        ix = ComponentArrays.label2index(p, label)
        if ix === nothing
            @error "Index not found for entry $tuple"
            return nothing
        end
        if isa(ix, AbstractArray)
            indices = vcat(indices, ix)
        else
            push!(indices, ix)
        end
    end
    return indices
end

function get_required_initialization_level(p_metadata, param_ixs)
    #Check for invalid parameters 
    for metadata_entry in p_metadata[param_ixs]
        if metadata_entry.type === DEVICE_SETPOINT
            @error "The parameter given is unsupported because it is a device setpoint."
            return nothing
        end
        if metadata_entry.in_mass_matrix === DEVICE_SETPOINT
            @error "The parameter given is not yet supported because it appears in the mass matrix"
            return nothing
        end
    end
    #Check for parameters which change the power flow 
    for metadata_entry in p_metadata[param_ixs]
        if metadata_entry.type === NETWORK_PARAM || metadata_entry.type === NETWORK_SETPOINT
            return POWERFLOW_AND_DEVICES
        end
    end
    #Check for parameters which change device initialization
    for metadata_entry in p_metadata[param_ixs]
        if metadata_entry.impacts_ic == true
            return DEVICES_ONLY
        end
    end
    return INITIALIZED
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

function get_sensitivity_functions(sim, param_data, state_data, solver, f_loss; kwargs...)
    p_metadata = sim.inputs.parameters_metadata
    p = sim.inputs.parameters
    param_ixs = get_indices_in_parameter_vector(p, param_data)
    metadata_ixs = get_indices_in_parameter_vector(p_metadata, param_data)
    init_level = get_required_initialization_level(p_metadata, metadata_ixs)
    if init_level === nothing
        return nothing
    end
    state_ixs = get_indices_in_state_vector(sim, state_data)
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
    if typeof(prob_old) <: SciMLBase.ODEProblem
        f_new = SciMLBase.ODEFunction{true}(
            f_old.f;
            mass_matrix = f_old.mass_matrix,
            tgrad = (dT, u, p, t) -> dT .= false,
        )
    elseif typeof(prob_old) <: SciMLBase.DDEProblem
        f_new = SciMLBase.DDEFunction{true}(
            f_old.f;
            mass_matrix = f_old.mass_matrix,
        )
    else
        @error "Problem type not supported"
    end
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
                @info "Reinitializing devices only"
                _initialize_devices_only!(x0, sim_inputs)
                _refine_initial_condition!(x0, p_new, prob)
            elseif init_level == INITIALIZED
                @info "I.C.s not impacted by parameter change"
            end
            if typeof(prob) <: SciMLBase.AbstractODEProblem
                prob_new = SciMLBase.remake(prob; p = p_new, u0 = x0)
            elseif typeof(prob) <: SciMLBase.AbstractDDEProblem
                h(p, t; idxs = nothing) = begin
                    typeof(idxs) <: Number ? x0[idxs] : x0
                end
                prob_new = SciMLBase.remake(prob; h = h, p = p_new, u0 = x0)
            end
            sol = SciMLBase.solve(prob_new, solver)
            @assert length(state_ixs) == 1  #Hardcode for single state for now 
            ix = state_ixs[1]
            ix_t = unique(i -> sol.t[i], eachindex(sol.t))
            state = sol[ix, ix_t]
            return f_loss(state, data)
        end
        function f_Zygote(p, x0, sys, sim_inputs, prob, data, init_level)   #Make separate f_enzymes depending on init_level? 
            p_new = sim_inputs.parameters
            p_new_buff = Zygote.Buffer(p_new)
            for ix in eachindex(p_new)
                p_new_buff[ix] = p_new[ix]
            end
            for (i, ix) in enumerate(param_ixs)
                p_new_buff[ix] = p[i]
            end
            p_new = copy(p_new_buff)
            if init_level == POWERFLOW_AND_DEVICES
                @error "POWERFLOW AND DEVICES -- not yet supported"
                #_initialize_powerflow_and_devices!(x0, inputs, sys)
            elseif init_level == DEVICES_ONLY
                @error "Reinitializing not supported with Zygote"
                #_initialize_devices_only!(x0, sim_inputs)     #Mutation 
                #_refine_initial_condition!(x0, p_new, prob)   #Mutation 
            elseif init_level == INITIALIZED
                @info "I.C.s not impacted by parameter change"
            end
            if typeof(prob) <: SciMLBase.AbstractODEProblem
                prob_new = SciMLBase.remake(prob; p = p_new, u0 = x0)
            elseif typeof(prob) <: SciMLBase.AbstractDDEProblem
                h(p, t; idxs = nothing) = begin
                    typeof(idxs) <: Number ? x0[idxs] : x0
                end
                prob_new = SciMLBase.remake(prob; h = h, p = p_new, u0 = x0)
            end
            sol = SciMLBase.solve(prob_new, solver)
            @assert length(state_ixs) == 1  #Hardcode for single state for now 
            #Hack to avoid unique(i -> sol.t[i], eachindex(sol.t)) which mutates and is incompatible with Zygote: 
            ix_first = findfirst(x -> x == 1.0, sol.t)
            ix_last = findlast(x -> x == 1.0, sol.t)
            ix_t = vcat(1:ix_first, (ix_last + 1):(length(sol.t)))

            ix = state_ixs[1]
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
        function f_forward_zygote(p, data)
            f_Zygote(
                p,
                x0,
                sys,
                deepcopy(sim.inputs_init),
                prob_new,
                data,
                init_level,
            )
        end
        f_forward, f_grad, f_forward_zygote
    end
end

#Inactive Rules
#Enzyme.EnzymeRules.inactive(::typeof(SimulationResults), args...) = nothing
#Enzyme.EnzymeRules.inactive(::typeof(get_activepower_branch_flow), args...) = nothing
