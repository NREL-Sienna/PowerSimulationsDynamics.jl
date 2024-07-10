"""
Function to obtain voltage and output currents for a source. It receives the simulation resutls and 
    an optional argument of the time step of the results.
"""
function post_proc_source_voltage_current_series(
    res::SimulationResults,
    name::String,
    dt = nothing,
    unique_timestamps = true,
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

    ts, Vb_R, Vb_I =
        post_proc_voltage_series(solution, bus_ix, n_buses, dt, unique_timestamps)

    I_R = R_th * (Vs_R .- Vb_R) / Z_sq + X_th * (Vs_I .- Vb_I) / Z_sq
    I_I = R_th * (Vs_I .- Vb_I) / Z_sq - X_th * (Vs_R .- Vb_R) / Z_sq

    return ts, Vs_R, Vs_I, I_R, I_I
end

"""
Function to obtain output real current for a source. It receives the simulation results,
the Source name and an optional argument of the time step of the results.

"""
function get_source_real_current_series(
    res::SimulationResults,
    name::String,
    dt = nothing,
    unique_timestamps = true,
)
    ts, _, _, I_R, _ =
        post_proc_source_voltage_current_series(res, name, dt, unique_timestamps)
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
    unique_timestamps = true,
)
    ts, _, _, _, I_I =
        post_proc_source_voltage_current_series(res, name, dt, unique_timestamps)
    return ts, I_I
end

"""
Function to obtain the output current time series of a PeriodicVariableSource model out of the DAE Solution. It receives the simulation inputs,
the dynamic device and bus voltage. It is dispatched for device type to compute the specific current.
compute_output_current(::SimulationResults, ::PeriodicVariableSource, ::Vector{Float64}, ::Vector{Float64}, ::Nothing)
"""
function compute_output_current(
    res::SimulationResults,
    dynamic_device::PSY.PeriodicVariableSource,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    dt::Union{Nothing, Float64, Vector{Float64}},
    unique_timestamps::Bool,
)
    name = PSY.get_name(dynamic_device)
    ts, Vt_internal = post_proc_state_series(res, (name, :Vt), dt, unique_timestamps)
    _, θt_internal = post_proc_state_series(res, (name, :θt), dt, unique_timestamps)

    Vr_internal = Vt_internal .* cos.(θt_internal)
    Vi_internal = Vt_internal .* sin.(θt_internal)

    R_th = PSY.get_R_th(dynamic_device)
    X_th = PSY.get_X_th(dynamic_device)
    Z_sq = R_th^2 + X_th^2

    I_R = R_th * (Vr_internal .- V_R) / Z_sq + X_th * (Vi_internal .- V_I) / Z_sq
    I_I = R_th * (Vi_internal .- V_I) / Z_sq - X_th * (Vr_internal .- V_R) / Z_sq

    return ts, I_R, I_I
end
