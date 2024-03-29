function mdl_source!(
    device_parameters,
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    static_device::StaticWrapper{PSY.Source},
) where {T <: ACCEPTED_REAL_TYPES}
    #Load device parameters
    R_th, X_th, internal_voltage, internal_angle = device_parameters
    V_R = internal_voltage * cos(internal_angle)
    V_I = internal_voltage * sin(internal_angle)
    Zmag = R_th^2 + X_th^2

    #update current
    current_r[1] += R_th * (V_R - voltage_r) / Zmag + X_th * (V_I - voltage_i) / Zmag #in system pu flowing out
    current_i[1] += R_th * (V_I - voltage_i) / Zmag - X_th * (V_R - voltage_r) / Zmag #in system pu flowing out

    return
end
