function mdl_source!(
    voltage_r,
    voltage_i,
    current_r,
    current_i,
    device::PSY.Source,
    sys::PSY.System,
)

    #Load device parameters
    bus = PSY.get_bus(device)
    #V_R = PSY.get_magnitude(bus) * cos(PSY.get_angle(bus))
    #V_I = PSY.get_magnitude(bus) * sin(PSY.get_angle(bus))
    V_R = PSY.get_internal_voltage(device) * cos(PSY.get_internal_angle(device))
    V_I = PSY.get_internal_voltage(device) * sin(PSY.get_internal_angle(device))
    R_th = PSY.get_R_th(device)
    X_th = PSY.get_X_th(device)
    Zmag = R_th^2 + X_th^2

    #update current
    current_r[1] += R_th * (V_R - voltage_r[1]) / Zmag + X_th * (V_I - voltage_i[1]) / Zmag #in system pu flowing out
    current_i[1] += R_th * (V_I - voltage_i[1]) / Zmag - X_th * (V_R - voltage_r[1]) / Zmag #in system pu flowing out

    return
end
