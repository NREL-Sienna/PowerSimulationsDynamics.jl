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
        if metadata_entry.in_mass_matrix === true
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
        if haskey(global_state_index, ref[1])
            state_ixs[ix] = get(global_state_index[ref[1]], ref[2], 0)
        elseif ref[2] ∈ [:Vr, :Vi]
            bus_ix = get_lookup(get_simulation_inputs(sim))[ref[1]]
            if ref[2] == :Vr
                state_ixs[ix] = bus_ix
            else
                state_ixs[ix] = bus_ix + length(get_lookup(get_simulation_inputs(sim)))
            end
        else
            @error "Available devices: $(keys(global_state_index))"
            error("State $(ref[2]) device $(ref[1]) not found in the system. ")
        end
    end
    return state_ixs
end

function get_parameter_values(sim, device_param_pairs)
    p = sim.inputs.parameters
    p_metadata = sim.inputs.parameters_metadata
    ixs, _ = get_ix_and_level(p, p_metadata, device_param_pairs)
    if ixs === nothing
        return nothing
    else
        return p[ixs]
    end
end

function get_parameter_labels(sim, device_param_pairs)
    p = sim.inputs.parameters
    p_metadata = sim.inputs.parameters_metadata
    labels = ComponentArrays.labels(p)
    ixs, _ = get_ix_and_level(p, p_metadata, device_param_pairs)
    if ixs === nothing
        return nothing
    else
        return labels[ixs]
    end
end

function get_ix_and_level(p, p_metadata, param_data)
    param_ixs = get_indices_in_parameter_vector(p, param_data)
    metadata_ixs = get_indices_in_parameter_vector(p_metadata, param_data)
    init_level = get_required_initialization_level(p_metadata, metadata_ixs)
    return param_ixs, init_level
end

function get_ix_and_level(p, p_metadata, param_data::Symbol)
    @assert param_data == :All
    @warn "Device setpoints and network parameters not yet supported; returning all supported parameters."
    param_ixs = Int64[]
    @assert ComponentArrays.labels(p) == ComponentArrays.labels(p_metadata)
    for label in ComponentArrays.labels(p)
        ix_metadata = ComponentArrays.label2index(p_metadata, label)
        ix_param = ComponentArrays.label2index(p_metadata, label)[1]
        metadata_entry = p_metadata[ix_metadata][1]

        if (metadata_entry.type === DEVICE_PARAM) &&
           (metadata_entry.in_mass_matrix === false)
            push!(param_ixs, ix_param)
        end
    end
    return param_ixs, DEVICES_ONLY
end

function convert_perturbations_to_callbacks(sys, sim_inputs, perturbations)
    perturbations_count = length(perturbations)
    cb = Vector{SciMLBase.DiscreteCallback}(undef, perturbations_count)
    tstops = Float64[]
    for (ix, pert) in enumerate(perturbations)
        condition = (x, t, integrator) -> t in [pert.time]
        affect = get_affect(sim_inputs, sys, pert)
        cb[ix] = SciMLBase.DiscreteCallback(condition, affect)
        push!(tstops, pert.time)
    end
    callbacks = SciMLBase.CallbackSet((), tuple(cb...))
    return callbacks, tstops
end

function get_sensitivity_functions(
    sim,
    param_data,
    state_data,
    solver,
    f_loss,
    sys_reinit = false;
    kwargs...,
)
    p_metadata = sim.inputs.parameters_metadata
    p = sim.inputs.parameters
    param_ixs, init_level = get_ix_and_level(p, p_metadata, param_data)
    if init_level === nothing
        return nothing
    end
    state_ixs = get_indices_in_state_vector(sim, state_data)
    sim_inputs = deepcopy(sim.inputs_init)

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
        advance_to_tstop = true,
        initializealg = SciMLBase.NoInit(),
        kwargs...,
    )
    if param_ixs === nothing
        return nothing
    else
        x0 = deepcopy(sim.x0_init)
        sys = deepcopy(sim.sys)
        function f_enzyme(
            p,
            x0,
            sys,
            sim_inputs,
            prob,
            data,
            init_level,
            sys_reinit,
            callbacks,
            tstops,
        )
            p_new = sim_inputs.parameters
            p_new[param_ixs] .= p
            if init_level == POWERFLOW_AND_DEVICES
                @error "POWERFLOW AND DEVICES -- not yet supported"
                #_initialize_powerflow_and_devices!(x0, inputs, sys)
            elseif init_level == DEVICES_ONLY
                if sys_reinit
                    @info "Reinitializing inividual devices and full system"
                    _initialize_devices_only!(x0, sim_inputs)
                    _refine_initial_condition!(x0, p_new, prob)
                else
                    @info "Reinitializing devices only"
                    _initialize_devices_only!(x0, sim_inputs)
                end
            elseif init_level == INITIALIZED
                @info "I.C.s not impacted by parameter change"
            end
            if typeof(prob) <: SciMLBase.AbstractODEProblem
                prob_new = SciMLBase.remake(prob; p = p_new, u0 = x0, tstops = tstops)
            elseif typeof(prob) <: SciMLBase.AbstractDDEProblem
                h(p, t; idxs = nothing) = begin
                    typeof(idxs) <: Number ? x0[idxs] : x0
                end
                prob_new =
                    SciMLBase.remake(prob; h = h, p = p_new, u0 = x0, tstops = tstops)
            end
            sol = solve_with_callback(prob_new, callbacks, solver)
            ix_t = unique(i -> sol.t[i], eachindex(sol.t))
            states = [sol[ix, ix_t] for ix in state_ixs]
            return f_loss(states, data)
        end
        function f_Zygote(
            p,
            x0,
            sys,
            sim_inputs,
            prob,
            data,
            init_level,
            sys_reinit,
            perts, 
        )
            p_new = sim_inputs.parameters
            p_new_buff = Zygote.Buffer(p_new)
            for ix in eachindex(p_new)
                p_new_buff[ix] = p_new[ix]
            end
            for (i, ix) in enumerate(param_ixs)
                p_new_buff[ix] = p[i]
            end
            p_new = copy(p_new_buff)
            #Converting perturbations to callbacks, restore after : https://github.com/EnzymeAD/Enzyme.jl/issues/1650
            @error "Zygote assumes single perturbation"
            cb = Vector{SciMLBase.DiscreteCallback}(undef, 1)
            pert = perts[1]
            condition = (x, t, integrator) -> t in [pert.time]
            affect = get_affect(sim_inputs, sys, pert)
            cb = [SciMLBase.DiscreteCallback(condition, affect)]
            callbacks = SciMLBase.CallbackSet((), tuple(cb...))
            tstops = [pert.time]

            if init_level == POWERFLOW_AND_DEVICES
                @error "POWERFLOW AND DEVICES -- not yet supported"
                #_initialize_powerflow_and_devices!(x0, inputs, sys)
            elseif init_level == DEVICES_ONLY
                @error "Reinitializing not supported with Zygote"
                _initialize_devices_only!(x0, sim_inputs)     #Mutation 
            elseif init_level == INITIALIZED
                @info "I.C.s not impacted by parameter change"
            end
            if typeof(prob) <: SciMLBase.AbstractODEProblem
                prob_new = SciMLBase.remake(prob; p = p_new, u0 = x0, tstops = tstops)
            elseif typeof(prob) <: SciMLBase.AbstractDDEProblem
                h(p, t; idxs = nothing) = begin
                    typeof(idxs) <: Number ? x0[idxs] : x0
                end
                prob_new =
                    SciMLBase.remake(prob; h = h, p = p_new, u0 = x0, tstops = tstops)
            end
            sol = SciMLBase.solve(prob_new, solver; callback = callbacks)
            #Hack to avoid unique(i -> sol.t[i], eachindex(sol.t)) which mutates and is incompatible with Zygote: 
            ix_first = findfirst(x -> x == 1.0, sol.t)
            ix_last = findlast(x -> x == 1.0, sol.t)
            ix_t = vcat(1:ix_first, (ix_last + 1):(length(sol.t)))
            states = [sol[ix, ix_t] for ix in state_ixs]
            return f_loss(states, data)
        end
        function f_forward(p, perts, data)
            callbacks, tstops = convert_perturbations_to_callbacks(sys, sim_inputs, perts) 
            f_enzyme(
                p,
                x0,
                sys,
                deepcopy(sim.inputs_init),
                prob_new,
                data,
                init_level,
                sys_reinit,
                callbacks,
                tstops,
            )
        end
        function f_grad(p, perts, data)
            callbacks, tstops = convert_perturbations_to_callbacks(sys, sim_inputs, perts)       
            dp = Enzyme.make_zero(p)
            dx0 = Enzyme.make_zero(x0)
            dsys = Enzyme.make_zero(sys)
            sim_inputs = deepcopy(sim.inputs_init)
            dsim_inputs = Enzyme.make_zero(sim_inputs)
            dprob_new = Enzyme.make_zero(prob_new)
            ddata = Enzyme.make_zero(data)
            dcallbacks = Enzyme.make_zero(callbacks)
            dtstops = Enzyme.make_zero(tstops)
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
                Enzyme.Const(sys_reinit),
                Enzyme.Duplicated(callbacks, dcallbacks),
                Enzyme.Duplicated(tstops, dtstops),
            )
            return dp
        end
        function f_forward_zygote(p, perts, data)
            f_Zygote(
                p,
                x0,
                sys,
                deepcopy(sim.inputs_init),
                prob_new,
                data,
                init_level,
                sys_reinit,
                perts,
            )
        end
        f_forward, f_grad, f_forward_zygote
    end
end

#Inactive Rules; potential path to improving API
#Enzyme.EnzymeRules.inactive(::typeof(SimulationResults), args...) = nothing
#Enzyme.EnzymeRules.inactive(::typeof(get_activepower_branch_flow), args...) = nothing
