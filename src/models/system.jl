struct System{Ctype <: Cache}
    inputs::SimulationInputs
    cache::Ctype
    function System(inputs, x0_init::Vector{T}) where T <: Number
        U = ForwardDiff.Dual{typeof(ForwardDiff.Tag(system_implicit!, T)), T, Val{ForwardDiff.pickchunksize(length(x0_init))}}
        cache = JacobianCache{T, U}
        new{typeof(cache)}(inputs, cache)
    end
end

function update_global_vars!(cache::Cache, inputs::SimulationInputs, x::AbstractArray{U}) where {U <: Real}
    index = get_global_vars_update_pointers(inputs)[GLOBAL_VAR_SYS_FREQ_INDEX]
    index == 0 && return
    # get_global_vars(cache)[:ω_sys] = x[index]
    return
end

function system_implicit!(
    out::Vector{T},
    dx::Vector{T},
    x::Vector{T},
    inputs::SimulationInputs,
    cache::Cache,
    t::Float64,
) where {T <: Real}
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
    V_r = @view x[1:bus_counts]
    V_i = @view x[bus_counts+1:get_bus_range(inputs)[end]]

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

    for d in get_static_injections_data(inputs)
        bus_n = PSY.get_number(PSY.get_bus(d))
        bus_ix = get_lookup(inputs)[bus_n]
        device!(
            view(V_r, bus_ix),
            view(V_i, bus_ix),
            view(I_injections_r, bus_ix),
            view(I_injections_i, bus_ix),
            d,
            inputs,
            t,
        )
    end

    if has_dyn_lines(inputs)
        for br in get_dynamic_branches(inputs)
            arc = PSY.get_arc(br)
            n_states = PSY.get_n_states(br)
            from_bus_number = PSY.get_number(arc.from)
            to_bus_number = PSY.get_number(arc.to)
            bus_ix_from = get_lookup(inputs)[from_bus_number]
            bus_ix_to = get_lookup(inputs)[to_bus_number]
            ix_range = range(branches_start, length = n_states)
            ode_range = range(branches_count, length = n_states)
            branches_count = branches_count + n_states
            branches_start += n_states
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

    out[bus_range] .=
        Ybus_current_kirchoff(inputs, V_r, V_i, I_injections_r, I_injections_i) .-
        M[bus_range, bus_range] * dx[bus_range]
end

function system_mass_matrix!(
    dx::AbstractArray{U},
    x::AbstractArray{U},
    inputs::SimulationInputs,
    cache::Cache,
    t,
) where {U <: Real}
    I_injections_r = get_current_injections_r(cache)
    I_injections_i = get_current_injections_i(cache)
    injection_ode = get_injection_ode(cache, T)
    branches_ode = get_branches_ode(cache, T)

    #Index Setup
    bus_size = get_bus_count(inputs)
    bus_vars_count = 2 * bus_size
    bus_range = 1:bus_vars_count
    injection_start = get_injection_pointer(inputs)
    injection_count = 1
    branches_start = get_branches_pointer(inputs)
    branches_count = 1
    update_global_vars!(inputs, x)

    #Network quantities
    V_r = @view x[1:bus_size]
    V_i = @view x[(bus_size + 1):bus_vars_count]
    fill!(I_injections_r, 0.0)
    fill!(I_injections_i, 0.0)

    for d in get_injectors_data(inputs)
        dynamic_device = PSY.get_dynamic_injector(d)
        bus_n = PSY.get_number(PSY.get_bus(d))
        bus_ix = get_lookup(inputs)[bus_n]
        n_states = PSY.get_n_states(dynamic_device)
        ix_range = range(injection_start, length = n_states)
        ode_range = range(injection_count, length = n_states)
        injection_count = injection_count + n_states
        injection_start = injection_start + n_states
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

    for d in get_static_injections_data(inputs)
        bus_n = PSY.get_number(PSY.get_bus(d))
        bus_ix = get_lookup(inputs)[bus_n]
        device!(
            view(V_r, bus_ix),
            view(V_i, bus_ix),
            view(I_injections_r, bus_ix),
            view(I_injections_i, bus_ix),
            d,
            inputs,
            t,
        )
    end

    if has_dyn_lines(inputs)
        for br in get_dynamic_branches(inputs)
            arc = PSY.get_arc(br)
            n_states = PSY.get_n_states(br)
            from_bus_number = PSY.get_number(arc.from)
            to_bus_number = PSY.get_number(arc.to)
            bus_ix_from = get_lookup(inputs)[from_bus_number]
            bus_ix_to = get_lookup(inputs)[to_bus_number]
            ix_range = range(branches_start, length = n_states)
            ode_range = range(branches_count, length = n_states)
            branches_count = branches_count + n_states
            branches_start += n_states
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

    dx[bus_range] .= Ybus_current_kirchoff(inputs, V_r, V_i, I_injections_r, I_injections_i)
end
