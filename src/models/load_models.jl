function mdl_Zload!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    static_device::PSY.PowerLoad,
) where {T <: Real}
    #Load squared voltage magnitude at steady state
    Vmag_sq = PSY.get_magnitude(PSY.get_bus(static_device))^2
    #Vmag_sq = 1.0

    #Load device parameters
    P = PSY.get_active_power(static_device)
    Q = PSY.get_reactive_power(static_device)

    #Compute impedance and load: Z = |V|^2/conj(S)
    #Z_load =  Vmag_sq/(P-Q*1im) #in pu
    #I = -(voltage_r[1] + voltage_i[1]*1im)/Z_load #in pu flowing out
    #current_r[1] += real(I)
    #current_i[1] += imag(I)

    #Update current
    #This model creates an equivalent RL/RC circuit based on steady state voltage
    current_r[1] += -(1.0 / Vmag_sq) * (voltage_r * P + voltage_i * Q) #in system pu flowing out
    current_i[1] += -(1.0 / Vmag_sq) * (voltage_i * P - voltage_r * Q) #in system pu flowing out

    return
end
