"""
Instantiate an ResidualModel for ForwardDiff calculations
"""
function ResidualModel(
    inputs::SimulationInputs,
    x0_init::Vector{T},
    ::Type{Ctype},
) where {T <: Float64, Ctype <: JacobianCache}
    U = ForwardDiff.Dual{
        typeof(ForwardDiff.Tag(system_residual!, T)),
        T,
        ForwardDiff.pickchunksize(length(x0_init)),
    }
    if isempty(inputs.delays)
        return SystemModel{ResidualModel, NoDelays}(
            inputs,
            Ctype{U}(system_residual!, inputs),
        )
    else
        error(
            "Cannot use ResidualModel for a system model with delays. Remove delays or use MassMatrixModel",
        )
    end
end

"""
Instantiate an ResidualModel for ODE inputs.
"""
function ResidualModel(inputs, ::Vector{Float64}, ::Type{Ctype}) where {Ctype <: SimCache}
    if isempty(inputs.delays)
        return SystemModel{ResidualModel, NoDelays}(inputs, Ctype(system_residual!, inputs))
    else
        error(
            "Cannot use ResidualModel for a system model with delays. Remove delays or use MassMatrixModel",
        )
    end
end

function (m::SystemModel{ResidualModel, NoDelays, C})(
    out::AbstractArray{T},
    du::AbstractArray{U},
    u::AbstractArray{V},
    p,
    t,
) where {
    C <: Cache,
    T <: ACCEPTED_REAL_TYPES,
    U <: ACCEPTED_REAL_TYPES,
    V <: ACCEPTED_REAL_TYPES,
}
    system_residual!(out, du, u, m.inputs, m.cache, t)
    return
end

function system_residual!(
    out::AbstractVector{T},
    dx::AbstractVector{U},
    x::AbstractVector{V},
    inputs::SimulationInputs,
    cache::Cache,
    t::Float64,
) where {T <: ACCEPTED_REAL_TYPES, U <: ACCEPTED_REAL_TYPES, V <: ACCEPTED_REAL_TYPES}
    update_global_vars!(cache, inputs, x)
    M = get_mass_matrix(inputs)

    # Global quantities
    system_ode_output = get_ode_output(cache, V)
    current_balance = get_current_balance(cache, V)
    fill!(current_balance, zero(V))
    bus_counts = get_bus_count(inputs)
    bus_range = get_bus_range(inputs)
    voltage_r = @view x[1:bus_counts]
    voltage_i = @view x[(bus_counts + 1):bus_range[end]]
    current_r = @view current_balance[1:bus_counts]
    current_i = @view current_balance[(bus_counts + 1):bus_range[end]]
    global_vars = get_global_vars(cache, V)
    inner_vars = get_inner_vars(cache, V)

    for dynamic_device in get_dynamic_injectors(inputs)
        ix_range = get_ix_range(dynamic_device)
        device_ode_output = @view system_ode_output[get_ode_ouput_range(dynamic_device)]
        device_inner_vars = @view inner_vars[get_inner_vars_index(dynamic_device)]
        device_states = @view x[ix_range]
        bus_ix = get_bus_ix(dynamic_device)

        device!(
            device_states,
            device_ode_output,
            voltage_r[bus_ix],
            voltage_i[bus_ix],
            view(current_r, bus_ix),
            view(current_i, bus_ix),
            global_vars,
            device_inner_vars,
            dynamic_device,
            nothing,
            t,
        )
        M_ = @view M[ix_range, ix_range]
        out[ix_range] .= device_ode_output .- M_ * dx[ix_range]
    end

    for static_load in get_static_loads(inputs)
        bus_ix = get_bus_ix(static_load)
        device!(
            voltage_r[bus_ix],
            voltage_i[bus_ix],
            view(current_r, bus_ix),
            view(current_i, bus_ix),
            global_vars,
            inner_vars,
            static_load,
            t,
        )
    end

    for static_device in get_static_injectors(inputs)
        bus_ix = get_bus_ix(static_device)
        device!(
            voltage_r[bus_ix],
            voltage_i[bus_ix],
            view(current_r, bus_ix),
            view(current_i, bus_ix),
            global_vars,
            inner_vars,
            static_device,
            t,
        )
    end

    if has_dyn_lines(inputs)
        branches_ode = get_branches_ode(cache, T)
        for dynamic_branch in get_dynamic_branches(inputs)
            dyn_branch = get_branch(dynamic_branch) # DynamicBranch
            branch = PSY.get_branch(dyn_branch) # Line or Transformer2W
            ix_range = get_ix_range(dynamic_branch)
            branch_output_ode = @view branches_ode[get_ode_ouput_range(dynamic_branch)]
            branch_states = @view x[ix_range]
            bus_ix_from = get_bus_ix_from(dynamic_branch)
            bus_ix_to = get_bus_ix_to(dynamic_branch)
            branch!(
                branch_states,
                branch_output_ode,
                #Get Voltage data
                voltage_r[bus_ix_from],
                voltage_i[bus_ix_from],
                voltage_r[bus_ix_to],
                voltage_i[bus_ix_to],
                #Get Current data
                view(current_r, bus_ix_from),
                view(current_i, bus_ix_from),
                view(current_r, bus_ix_to),
                view(current_i, bus_ix_to),
                dynamic_branch,
                branch,
            )
            M_ = @view M[ix_range, ix_range]
            out[ix_range] .= branch_output_ode .- M_ * dx[ix_range]
        end
    end
    voltages = @view x[bus_range]
    M_ = @view M[bus_range, bus_range]
    out[bus_range] .= network_model(inputs, current_balance, voltages) .- M_ * dx[bus_range]
    return
end

"""
Instantiate a MassMatrixModel for ODE inputs.
"""
function MassMatrixModel(inputs, ::Vector{Float64}, ::Type{Ctype}) where {Ctype <: SimCache}
    if isempty(inputs.delays)
        return SystemModel{MassMatrixModel, NoDelays}(
            inputs,
            Ctype(system_mass_matrix!, inputs),
        )
    else
        return SystemModel{MassMatrixModel, HasDelays}(
            inputs,
            Ctype(system_mass_matrix!, inputs),
        )
    end
end

"""
Instantiate a MassMatrixModel for ForwardDiff calculations
"""
function MassMatrixModel(
    inputs::SimulationInputs,
    x0_init::Vector{T},
    ::Type{Ctype},
) where {T <: Float64, Ctype <: JacobianCache}
    U = ForwardDiff.Dual{
        typeof(ForwardDiff.Tag(system_mass_matrix!, T)),
        T,
        ForwardDiff.pickchunksize(length(x0_init)),
    }
    if isempty(inputs.delays)
        return SystemModel{MassMatrixModel, NoDelays}(
            inputs,
            Ctype{U}(system_mass_matrix!, inputs),
        )
    else
        return SystemModel{MassMatrixModel, HasDelays}(
            inputs,
            Ctype{U}(system_mass_matrix!, inputs),
        )
    end
end

function (m::SystemModel{MassMatrixModel, NoDelays, C})(
    du::AbstractArray{T},
    u::AbstractArray{U},
    p,
    t,
) where {C <: Cache, U <: ACCEPTED_REAL_TYPES, T <: ACCEPTED_REAL_TYPES}
    system_mass_matrix!(du, u, nothing, m.inputs, m.cache, t)
end

function (m::SystemModel{MassMatrixModel, HasDelays, C})(
    du::AbstractArray{T},
    u::AbstractArray{U},
    h,
    p,
    t,
) where {C <: Cache, U <: ACCEPTED_REAL_TYPES, T <: ACCEPTED_REAL_TYPES}
    system_mass_matrix!(du, u, h, m.inputs, m.cache, t)
end

function system_mass_matrix!(
    dx::Vector{T},
    x::Vector{V},
    h,
    inputs::SimulationInputs,
    cache::Cache,
    t,
) where {T <: ACCEPTED_REAL_TYPES, V <: ACCEPTED_REAL_TYPES}
    update_global_vars!(cache, inputs, x)

    # Global quantities
    system_ode_output = get_ode_output(cache, V)
    current_balance = get_current_balance(cache, V)
    fill!(current_balance, zero(V))
    bus_counts = get_bus_count(inputs)
    bus_range = get_bus_range(inputs)
    voltage_r = @view x[1:bus_counts]
    voltage_i = @view x[(bus_counts + 1):bus_range[end]]
    current_r = @view current_balance[1:bus_counts]
    current_i = @view current_balance[(bus_counts + 1):bus_range[end]]
    global_vars = get_global_vars(cache, V)
    inner_vars = get_inner_vars(cache, V)

    for dynamic_device in get_dynamic_injectors(inputs)
        ix_range = get_ix_range(dynamic_device)
        device_ode_output = @view system_ode_output[get_ode_ouput_range(dynamic_device)]
        device_inner_vars = @view inner_vars[get_inner_vars_index(dynamic_device)]
        device_states = @view x[ix_range]
        bus_ix = get_bus_ix(dynamic_device)

        device!(
            device_states,
            device_ode_output,
            voltage_r[bus_ix],
            voltage_i[bus_ix],
            view(current_r, bus_ix),
            view(current_i, bus_ix),
            global_vars,
            device_inner_vars,
            dynamic_device,
            h,
            t,
        )
        dx[ix_range] .= device_ode_output
    end

    for static_load in get_static_loads(inputs)
        bus_ix = get_bus_ix(static_load)
        device!(
            voltage_r[bus_ix],
            voltage_i[bus_ix],
            view(current_r, bus_ix),
            view(current_i, bus_ix),
            global_vars,
            inner_vars,
            static_load,
            t,
        )
    end

    for static_device in get_static_injectors(inputs)
        bus_ix = get_bus_ix(static_device)
        device!(
            voltage_r[bus_ix],
            voltage_i[bus_ix],
            view(current_r, bus_ix),
            view(current_i, bus_ix),
            global_vars,
            inner_vars,
            static_device,
            t,
        )
    end

    if has_dyn_lines(inputs)
        branches_ode = get_branches_ode(cache, T)
        for dynamic_branch in get_dynamic_branches(inputs)
            dyn_branch = get_branch(dynamic_branch) # DynamicBranch
            branch = PSY.get_branch(dyn_branch) # Line or Transformer2W
            ix_range = get_ix_range(dynamic_branch)
            branch_output_ode = @view branches_ode[get_ode_ouput_range(dynamic_branch)]
            branch_states = @view x[ix_range]
            bus_ix_from = get_bus_ix_from(dynamic_branch)
            bus_ix_to = get_bus_ix_to(dynamic_branch)
            branch!(
                branch_states,
                branch_output_ode,
                #Get Voltage data
                voltage_r[bus_ix_from],
                voltage_i[bus_ix_from],
                voltage_r[bus_ix_to],
                voltage_i[bus_ix_to],
                #Get Current data
                view(current_r, bus_ix_from),
                view(current_i, bus_ix_from),
                view(current_r, bus_ix_to),
                view(current_i, bus_ix_to),
                dynamic_branch,
                branch,
            )
            dx[ix_range] .= branch_output_ode
        end
    end
    voltages = @view x[bus_range]
    dx[bus_range] .= network_model(inputs, current_balance, voltages)
    return
end
