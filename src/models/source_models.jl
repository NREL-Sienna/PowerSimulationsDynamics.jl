function mdl_source!(voltage_r, voltage_i, current_r, current_i, device::PSY.Source)
    #Load device parameters
    V_R =
        device.ext[CONTROL_REFS][V_source_index] *
        cos(device.ext[CONTROL_REFS][θ_source_index])
    V_I =
        device.ext[CONTROL_REFS][V_source_index] *
        sin(device.ext[CONTROL_REFS][θ_source_index])
    R_th = PSY.get_R_th(device)
    X_th = PSY.get_X_th(device)
    Zmag = R_th^2 + X_th^2

    #update current
    current_r[1] += R_th * (V_R - voltage_r[1]) / Zmag + X_th * (V_I - voltage_i[1]) / Zmag #in system pu flowing out
    current_i[1] += R_th * (V_I - voltage_i[1]) / Zmag - X_th * (V_R - voltage_r[1]) / Zmag #in system pu flowing out

    return
end


function device!(
    x,
    output_ode::Vector{T},
    voltage_r,
    voltage_i,
    current_r,
    current_i,
    ix_range::UnitRange{Int},
    ode_range::UnitRange{Int},
    dynamic_device::PSY.PeriodicVariableSource,
    inputs::SimulationInputs,
    t,
) where {T <: Real}
    #Obtain local device states
    internal_states = @view x[ix_range]

    ω_θ = PSY.get_angle_frequencies(dynamic_device)
    ω_V = PSY.get_angle_frequencies(dynamic_device)
    output_ode[ode_range][ix_range[1]] = dV = 0
    for (ix, A) in PSY.get_internal_voltage_coefficients(dynamic_device)
        t <= 0 && continue
        dV += ω_V[ix]*(A[1]*cos(ω_V[ix]* t) - A[2]*sin(ω_V[ix]* t))
    end

    output_ode[ode_range][ix_range[2]] = dθ = 0
    for (ix, A) in PSY.get_internal_angle_coefficients(dynamic_device)
        t <= 0 && continue
        dθ += ω_θ[ix]*(A[1]*cos(ω_θV[ix]* t) - A[2]*sin(ω_θ[ix]* t))
    end

    # Internal Voltage states
    V_R = internal_states[1]*cos(internal_states[2])
    V_I = internal_states[1]*sin(internal_states[2])

    R_th = PSY.get_R_th(device)
    X_th = PSY.get_X_th(device)
    Zmag = R_th^2 + X_th^2

    #update current
    current_r[1] += R_th * (V_R - voltage_r[1]) / Zmag + X_th * (V_I - voltage_i[1]) / Zmag #in system pu flowing out
    current_i[1] += R_th * (V_I - voltage_i[1]) / Zmag - X_th * (V_R - voltage_r[1]) / Zmag #in system pu flowing out

    return
end
