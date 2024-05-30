function mdl_source!(
    p,
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    static_device::StaticWrapper{PSY.Source},
) where {T <: ACCEPTED_REAL_TYPES}
    #Load device parameters
    R_th = p[:params][:R_th]
    X_th = p[:params][:X_th]
    #Load device references
    V_ref = p[:refs][:V_ref]
    θ_ref = p[:refs][:θ_ref]
    #@error V_ref
    V_R = V_ref * cos(θ_ref)
    V_I = V_ref * sin(θ_ref)
    Zmag = R_th^2 + X_th^2

    #update current
    current_r[1] += R_th * (V_R - voltage_r) / Zmag + X_th * (V_I - voltage_i) / Zmag #in system pu flowing out
    current_i[1] += R_th * (V_I - voltage_i) / Zmag - X_th * (V_R - voltage_r) / Zmag #in system pu flowing out

    return
end
