function mdl_source!(voltage_r, voltage_i, current_r, current_i, device::StaticWrapper{PSY.Source})
    #Load device parameters
    bus_ix = get_bus_ix(device)
    V_R = get_V_ref(device) * cos(get_θ_ref(device))
    V_I = get_V_ref(device) * sin(get_θ_ref(device))
    # TODO: field forwarding
    R_th = PSY.get_R_th(device.device)
    X_th = PSY.get_X_th(device.device)
    Zmag = R_th^2 + X_th^2

    #update current
    current_r[bus_ix] += R_th * (V_R - voltage_r[bus_ix]) / Zmag + X_th * (V_I - voltage_i[bus_ix]) / Zmag #in system pu flowing out
    current_i[bus_ix] += R_th * (V_I - voltage_i[bus_ix]) / Zmag - X_th * (V_R - voltage_r[bus_ix]) / Zmag #in system pu flowing out

    return
end
