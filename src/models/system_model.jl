function system_model!(out::Vector{T}, dx, x, sys, t) where {T<:Real}

    #Index Setup
    bus_size = length(PSY.get_components(PSY.Bus, sys))
    bus_range = 1:(2*bus_size)
    bus_vars_count = length(bus_range)
    injection_start = get_injection_pointer(sys)
    injection_count = 1
    branches_start = get_branches_pointer(sys)
    branches_count = 1

    #Network quantities
    V_r = @view x[1:bus_size]
    V_i = @view x[(bus_size+1):bus_vars_count]
    Sbase = PSY.get_basepower(sys)
    I_injections_r = zeros(T, bus_size)
    I_injections_i = zeros(T, bus_size)
    injection_ode = zeros(T, get_n_injection_states(sys))
    branches_ode = zeros(T, get_n_branches_states(sys))

    for d in PSY.get_components(PSY.DynamicInjection, sys)
        bus_n = PSY.get_number(PSY.get_bus(d)) # TODO: This requires that the bus numbers are indexed 1-N
        ix_range = range(injection_start, length = PSY.get_n_states(d))
        ode_range = range(injection_count, length = PSY.get_n_states(d))
        injection_count = injection_count + PSY.get_n_states(d)
        injection_start = injection_start + PSY.get_n_states(d)
        device_model!(
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
    #=
    if !isnothing(sys.dyn_branch)
        for br in sys.dyn_branch
            arc = br.arc
            n_buses = get_bus_size(sys)
            from_bus_number = PSY.get_number(arc.from)
            to_bus_number = PSY.get_number(arc.to)
            ix_dx = [from_bus_number,
                    from_bus_number+n_buses,
                    to_bus_number,
                    to_bus_number+n_buses]
            ix_range = range(branches_start, length=br.n_states)
            ode_range = range(branches_count, length=br.n_states)
            branches_count = branches_count + br.n_states
            branch_model!(x,
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
                        ix_dx,
                        ode_range,
                        br,
                        sys)

            out[ix_range] = branches_ode[ode_range] - dx[ix_range]
        end
    end
    =#

    for d in PSY.get_components(PSY.StaticInjection, sys)
        bus_n = PSY.get_number(PSY.get_bus(d))

        device_model!(
            view(V_r, bus_n),
            view(V_i, bus_n),
            view(I_injections_r, bus_n),
            view(I_injections_i, bus_n),
            d,
            sys,
        )
    end

    out[bus_range] = kcl(PSY.get_ext(sys)[YBUS], V_r, V_i, I_injections_r, I_injections_i)

end
