function system_model!(out, dx, x, (controls, sys), t)

    #Index Setup
    bus_size = get_bus_size(sys)
    bus_range = get_bus_range(sys)
    bus_vars_count = length(bus_range)
    injection_start = get_injection_pointer(sys)
    injection_count = 1
    branches_start = get_branches_pointer(sys)
    branches_count = 1

    #Network quantities
    V_r = @view x[1:bus_size]
    V_i = @view x[bus_size+1:bus_vars_count]
    Sbase = get_sys_base(sys)
    I_injections_r = zeros(bus_size)
    I_injections_i = zeros(bus_size)
    injection_ode = zeros(get_n_injection_states(sys))
    branches_ode = zeros(get_n_branches_states(sys))

    for d in sys.dyn_injections
        bus_n = get_bus_number(d)
        ix_range = range(injection_start, length=d.n_states)
        ode_range = range(injection_count, length=d.n_states)
        injection_count = injection_count + d.n_states
        injection_start = injection_start + d.n_states
        device_model!(x,
                      injection_ode,
                      view(V_r, bus_n),
                      view(V_i, bus_n),
                      view(I_injections_r, bus_n),
                      view(I_injections_i, bus_n),
                      ix_range,
                      ode_range,
                      controls,
                      d,
                      sys)
        out[ix_range] = injection_ode[ode_range] - dx[ix_range]
    end
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
                        controls,
                        br,
                        sys)

            out[ix_range] = branches_ode[ode_range] - dx[ix_range]
        end
    end

    for d in sys.injections
        bus_n = get_bus_number(d)

        device_model!( view(V_r, bus_n),
                       view(V_i, bus_n),
                       view(I_injections_r, bus_n),
                       view(I_injections_i, bus_n),
                       d,
                       sys)
    end

    out[bus_range] = kcl(sys.Ybus, V_r, V_i,
                               I_injections_r, I_injections_i)


end
