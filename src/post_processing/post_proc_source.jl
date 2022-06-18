"""
Function to obtain voltage and output currents for a source. It receives the simulation resutls and 
    an optional argument of the time step of the results.
"""

function post_proc_source_voltage_current_series(
    res::SimulationResults;
    dt = nothing,
)::NTuple{5, Vector{Float64}}
    system = get_system(res)
    bus_lookup = get_bus_lookup(res)
    n_buses = length(bus_lookup)
    solution = res.solution
    source = collect(PSY.get_components(PSY.Source, system))[1]
    device = PSY.get_component(PSY.StaticInjection, system, source.name)
    if isnothing(device)
        error("Device $(name) not found in the system")
    end

    Vt_sim = get_state_series(res, (source.name, :Vt), dt = dt)
    θt_sim = get_state_series(res, (source.name, :θt), dt = dt)

    Vs_R = Vt_sim[2] .* cos.(θt_sim[2])
    Vs_I = Vt_sim[2] .* sin.(θt_sim[2])

    R_th = PSY.get_R_th(source)
    X_th = PSY.get_X_th(source)
    Z_sq = R_th^2 + X_th^2

    bus_ix = get(bus_lookup, PSY.get_number(PSY.get_bus(device)), -1)

    ts, Vb_R, Vb_I = post_proc_voltage_series(solution, bus_ix, n_buses, dt)

    I_R = R_th * (Vs_R - Vb_R) / Z_sq + X_th * (Vs_I - Vb_I) / Z_sq
    I_I = R_th * (Vs_I - Vb_I) / Z_sq - X_th * (Vs_R - Vb_R) / Z_sq

    return ts, Vs_R, Vs_I, I_R, I_I
end