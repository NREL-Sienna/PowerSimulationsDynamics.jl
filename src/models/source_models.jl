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
