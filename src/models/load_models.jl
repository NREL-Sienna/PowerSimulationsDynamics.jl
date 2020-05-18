function mdl_Zload!(
    voltage_r,
    voltage_i,
    current_r,
    current_i,
    device::PSY.PowerLoad,
    sys::PSY.System,
)

    #Load squared voltage magnitude at steady state
    bus = PSY.get_bus(device)
    Vmag_sq = PSY.get_voltage(bus)^2
    #Vmag_sq = 1.0

    #Load device parameters
    P = PSY.get_activepower(device)
    Q = PSY.get_reactivepower(device)

    #Compute impedance and load: Z = |V|^2/conj(S)
    #Z_load =  Vmag_sq/(P-Q*1im) #in pu
    #I = -(voltage_r[1] + voltage_i[1]*1im)/Z_load #in pu flowing out
    #current_r[1] += real(I)
    #current_i[1] += imag(I)

    #Update current
    #This model creates an equivalent RL/RC circuit based on steady state voltage
    current_r[1] += -(1.0 / Vmag_sq) * (voltage_r[1] * P + voltage_i[1] * Q) #in system pu flowing out
    current_i[1] += -(1.0 / Vmag_sq) * (voltage_i[1] * P - voltage_r[1] * Q) #in system pu flowing out

    return
end
