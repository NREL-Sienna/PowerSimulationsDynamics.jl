function update_global_vars!(sys::PSY.System, x::AbstractArray)
    index = PSY.get_ext(sys)[GLOBAL_VARS][:ω_sys_index]
    index == 0 && return
    PSY.get_ext(sys)[GLOBAL_VARS][:ω_sys] = x[index]
    return
end

function system!(out::Vector{<:Real}, dx, x, sys::PSY.System, t::Float64)

    I_injections_r = PSY.get_ext(sys)[AUX_ARRAYS][1]
    I_injections_i = PSY.get_ext(sys)[AUX_ARRAYS][2]
    injection_ode = PSY.get_ext(sys)[AUX_ARRAYS][3]
    branches_ode = PSY.get_ext(sys)[AUX_ARRAYS][4]

    #Index Setup
    bus_size = get_bus_count(sys)
    bus_vars_count = 2 * bus_size
    bus_range = 1:bus_vars_count
    injection_start = get_injection_pointer(sys)
    injection_count = 1
    branches_start = get_branches_pointer(sys)
    branches_count = 1
    update_global_vars!(sys, x)

    #Network quantities
    V_r = @view x[1:bus_size]
    V_i = @view x[(bus_size + 1):bus_vars_count]
    Sbase = PSY.get_basepower(sys)
    fill!(I_injections_r, 0.0)
    fill!(I_injections_i, 0.0)

    for d in PSY.get_components(PSY.DynamicInjection, sys)
        bus_n = PSY.get_number(PSY.get_bus(d)) # TODO: This requires that the bus numbers are indexed 1-N
        n_states = PSY.get_n_states(d)
        ix_range = range(injection_start, length = n_states)
        ode_range = range(injection_count, length = n_states)
        injection_count = injection_count + n_states
        injection_start = injection_start + n_states
        device!(
            x,
            injection_ode,
            view(V_r, bus_n),
            view(V_i, bus_n),
            view(I_injections_r, bus_n),
            view(I_injections_i, bus_n),
            ix_range,
            ode_range,
            d,
            sys,
        )
        out[ix_range] = injection_ode[ode_range] - dx[ix_range]
    end

    for d in PSY.get_components(PSY.StaticInjection, sys)
        isa(d, PSY.Generator) && continue
        bus_n = PSY.get_number(PSY.get_bus(d))

        device!(
            view(V_r, bus_n),
            view(V_i, bus_n),
            view(I_injections_r, bus_n),
            view(I_injections_i, bus_n),
            d,
            sys,
        )
    end

    dyn_branches = PSY.get_components(DynamicLine, sys)
    if !(isempty(dyn_branches))
        for br in dyn_branches
            arc = PSY.get_arc(br)
            n_states = PSY.get_n_states(br)
            from_bus_number = PSY.get_number(arc.from)
            to_bus_number = PSY.get_number(arc.to)
            ix_range = range(branches_start, length = n_states)
            ode_range = range(branches_count, length = n_states)
            branches_count = branches_count + n_states
            branch!(
                x,
                dx,
                branches_ode,
                #Get Voltage data
                view(V_r, from_bus_number),
                view(V_i, from_bus_number),
                view(V_r, to_bus_number),
                view(V_i, to_bus_number),
                #Get Current data
                view(I_injections_r, from_bus_number),
                view(I_injections_i, from_bus_number),
                view(I_injections_r, to_bus_number),
                view(I_injections_i, to_bus_number),
                ix_range,
                ode_range,
                br,
                sys,
            )
            out[ix_range] = branches_ode[ode_range] - dx[ix_range]
        end
    end

    kirchoff_laws!(sys, V_r, V_i, I_injections_r, I_injections_i, dx)
    out[bus_range] = PSY.get_ext(sys)[AUX_ARRAYS][6]

end
