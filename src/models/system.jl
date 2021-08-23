struct SystemModel{T<:PSID.SimulationModel}
    inputs::SimulationInputs
    cache::Cache
end

"""
    Instantiate an ImplicitModel for ForwardDiff calculations
"""
function ImplicitModel(inputs, x0_init::Vector{T}, ::Type{Ctype}) where {T <: Number, Ctype <: JacobianCache}
    U = ForwardDiff.Dual{typeof(ForwardDiff.Tag(system_implicit!, T)), T, Val{ForwardDiff.pickchunksize(length(x0_init))}}
    return SystemModel{ImplicitModel}(inputs, Ctype{T,U}(inputs))
end

"""
Instantiate an ImplicitModel for ODE inputs.
"""
function ImplicitModel(inputs, ::Vector{T}, ::Type{Ctype}) where {Ctype <: SimCache{T}} where {T <: Real}
    return SystemModel{ImplicitModel}(inputs, Ctype(inputs))
end

function (m::SystemModel{ImplicitModel})(out::AbstractArray{T}, du::AbstractArray{T}, u::AbstractArray{T}, p, t) where {T <: Real}
    return system_implicit!(out, du, u, f.inputs, f.cache, t)
end

function system_implicit!(
    out::Vector{T},
    dx::Vector{U},
    x::Vector{V},
    inputs::SimulationInputs,
    cache::Cache,
    t::Float64,
) where {T, U, V <: Real}
    I_injections_r = get_current_injections_r(cache, T)
    I_injections_i = get_current_injections_i(cache, T)
    injection_ode = get_injection_ode(cache, T)
    branches_ode = get_branches_ode(cache, T)
    M = get_mass_matrix(inputs)
    update_global_vars!(cache, inputs, x)
    fill!(I_injections_r, 0.0)
    fill!(I_injections_i, 0.0)

    #Network quantities
    bus_counts = get_bus_count(inputs)
    bus_range = get_bus_range(inputs)
    voltages = @view x[bus_range]
    V_r = @view x[1:bus_counts]
    V_i = @view x[bus_counts+1:bus_range[end]]

    for dynamic_device in get_dynamic_injectors_data(inputs)
        ix_range = get_ix_range(dynamic_device)
        ode_range = get_ode_range(dynamic_device)
        device!(
            x,
            injection_ode,
            V_r,
            V_i,
            I_injections_r,
            I_injections_i,
            dynamic_device,
            inputs,
            cache,
            t,
        )
        M_ = @view M[ix_range, ix_range]
        out[ix_range] .= injection_ode[ode_range] .- M_ * dx[ix_range]
    end

    for static_device in get_static_injectiors_data(inputs)
        device!(
            V_r,
            V_i,
            I_injections_r,
            I_injections_i,
            static_device,
            inputs,
            cache,
            t,
        )
    end

    if has_dyn_lines(inputs)
        for br in get_dynamic_branches(inputs)
            ix_range = get_ix_range(dynamic_branch)
            ode_range = get_ode_range(dynamic_branch)
            branch!(
                x,
                dx,
                branches_ode,
                #Get Voltage data
                view(V_r, bus_ix_from),
                view(V_i, bus_ix_from),
                view(V_r, bus_ix_to),
                view(V_i, bus_ix_to),
                #Get Current data
                view(I_injections_r, bus_ix_from),
                view(I_injections_i, bus_ix_from),
                view(I_injections_r, bus_ix_to),
                view(I_injections_i, bus_ix_to),
                ix_range,
                ode_range,
                br,
                inputs,
            )
            M_ = @view M[ix_range, ix_range]
            out[ix_range] .= branches_ode[ode_range] .- M_ * dx[ix_range]
        end
    end

    out[bus_range] .= network_model(inputs, cache, voltages) .- M[bus_range, bus_range] * dx[bus_range]
end

"""
Instantiate a MassMatrixModel for ODE inputs.
"""
function MassMatrixModel(inputs, ::Vector{T}, ::Type{Ctype}) where {Ctype <: SimCache{T}} where T <: Number
    return SystemModel{MassMatrixModel}(inputs, Ctype(inputs))
end

function MassMatrixModel(inputs, x0_init::Vector{T}, ::Type{Ctype}) where {T <: Number, Ctype <: JacobianCache}
    U = ForwardDiff.Dual{typeof(ForwardDiff.Tag(system_mass_matrix!, T)), T, Val{ForwardDiff.pickchunksize(length(x0_init))}}
    return SystemModel{MassMatrixModel}(inputs, Ctype{T,U}(inputs))
end

function (m::SystemModel{MassMatrixModel})(du::AbstractArray{T}, u::AbstractArray{T}, p, t) where {T <: Real}
    system_mass_matrix!(du, u, f.inputs, f.cache, t)
end

function system_mass_matrix!(
    dx::Vector{T},
    x::Vector{T},
    inputs::SimulationInputs,
    cache::Cache,
    t,
) where {T <: Real}
    I_injections_r = get_current_injections_r(cache, T)
    I_injections_i = get_current_injections_i(cache, T)
    injection_ode = get_injection_ode(cache, T)
    branches_ode = get_branches_ode(cache, T)
    update_global_vars!(cache, inputs, x)
    fill!(I_injections_r, 0.0)
    fill!(I_injections_i, 0.0)

    #Network quantities
    bus_counts = get_bus_count(inputs)
    bus_range = get_bus_range(inputs)
    voltages = @view x[bus_range]
    V_r = @view x[1:bus_counts]
    V_i = @view x[bus_counts+1:bus_range[end]]

    for dynamic_device in get_injectors_data(inputs)
        ix_range = get_ix_range(dynamic_device)
        ode_range = get_ode_range(dynamic_device)
        device!(
            x,
            injection_ode,
            view(V_r, bus_ix),
            view(V_i, bus_ix),
            view(I_injections_r, bus_ix),
            view(I_injections_i, bus_ix),
            ix_range,
            ode_range,
            dynamic_device,
            inputs,
            t,
        )
        dx[ix_range] .= injection_ode[ode_range]
    end

    for static_device in get_static_injectiors_data(inputs)
        device!(
            V_r,
            V_i,
            I_injections_r,
            I_injections_i,
            static_device,
            inputs,
            cache,
            t,
        )
    end

    if has_dyn_lines(inputs)
        for dynamic_branch in get_dynamic_branches(inputs)
            ix_range = get_ix_range(dynamic_branch)
            ode_range = get_ode_range(dynamic_branch)
            branch!(
                x,
                dx,
                branches_ode,
                #Get Voltage data
                view(V_r, bus_ix_from),
                view(V_i, bus_ix_from),
                view(V_r, bus_ix_to),
                view(V_i, bus_ix_to),
                #Get Current data
                view(I_injections_r, bus_ix_from),
                view(I_injections_i, bus_ix_from),
                view(I_injections_r, bus_ix_to),
                view(I_injections_i, bus_ix_to),
                ix_range,
                ode_range,
                br,
                inputs,
            )
            dx[ix_range] .= branches_ode[ode_range]
        end
    end

    dx[bus_range] .= network_model(inputs, cache, voltages)
end
