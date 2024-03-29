function mdl_source!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    static_device::StaticWrapper{PSY.Source},
) where {T <: ACCEPTED_REAL_TYPES}
    #Load device parameters
    V_R = get_V_ref(static_device) * cos(get_θ_ref(static_device))
    V_I = get_V_ref(static_device) * sin(get_θ_ref(static_device))
    # TODO: field forwarding
    R_th = PSY.get_R_th(static_device.device)
    X_th = PSY.get_X_th(static_device.device)
    Zmag = R_th^2 + X_th^2

    #update current
    current_r[1] += R_th * (V_R - voltage_r) / Zmag + X_th * (V_I - voltage_i) / Zmag #in system pu flowing out
    current_i[1] += R_th * (V_I - voltage_i) / Zmag - X_th * (V_R - voltage_r) / Zmag #in system pu flowing out

    return
end
