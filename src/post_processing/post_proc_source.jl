"""
Function to obtain voltage and output currents for a source. It receives the simulation resutls and 
    an optional argument of the time step of the results.
"""
function post_proc_source_voltage_current_series(
    res::SimulationResults,
    name::String,
    dt = nothing,
)
    system = get_system(res)
    bus_lookup = get_bus_lookup(res)
    n_buses = length(bus_lookup)
    solution = res.solution
    device = PSY.get_component(PSY.StaticInjection, system, name)
    if isnothing(device)
        error("Device $(name) not found in the system")
    end
    if typeof(device) != PSY.Source
        error("Device $(name) is not a Source.")
    end

    Vt_sim = PSY.get_internal_voltage(device)
    θt_sim = PSY.get_internal_angle(device)

    Vs_R = Vt_sim * cos.(θt_sim)
    Vs_I = Vt_sim * sin.(θt_sim)

    R_th = PSY.get_R_th(device)
    X_th = PSY.get_X_th(device)
    Z_sq = R_th^2 + X_th^2

    bus_ix = get(bus_lookup, PSY.get_number(PSY.get_bus(device)), -1)

    ts, Vb_R, Vb_I = post_proc_voltage_series(solution, bus_ix, n_buses, dt)

    I_R = R_th * (Vs_R .- Vb_R) / Z_sq + X_th * (Vs_I .- Vb_I) / Z_sq
    I_I = R_th * (Vs_I .- Vb_I) / Z_sq - X_th * (Vs_R .- Vb_R) / Z_sq

    return ts, Vs_R, Vs_I, I_R, I_I
end

"""
Function to obtain output real current for a source. It receives the simulation results,
the Source name and an optional argument of the time step of the results.

"""
function get_source_real_current_series(res::SimulationResults, name::String, dt = nothing)
    ts, _, _, I_R, _ = post_proc_source_voltage_current_series(res, name, dt)
    return ts, I_R
end

"""
Function to obtain output imaginary current for a source. It receives the simulation results,
the Source name and an optional argument of the time step of the results.

"""
function get_source_imaginary_current_series(
    res::SimulationResults,
    name::String,
    dt = nothing,
)
    ts, _, _, _, I_I = post_proc_source_voltage_current_series(res, name, dt)
    return ts, I_I
end
