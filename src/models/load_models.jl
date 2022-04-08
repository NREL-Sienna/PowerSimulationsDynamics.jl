"""
Model for ZIP Load model given by:

P_zip = P_power + P_current * (V / V0) + P_impedance * (V / V0)^2
Q_zip = Q_power + Q_current * (V / V0) + Q_impedance * (V / V0)^2

with V = sqrt(V_r^2 + V_i^2) and V0 the voltage magnitude from the power flow solution

The current taken for the load is computed as:
I_zip = (P_zip + j Q_zip) / (V_r + j V_i)

For constant impedance it is obtained:
Iz_re = (1 / V0)^2 * (V_r * P_impedance + V_i * Q_impedance)
Iz_im = (1 / V0)^2 * (V_i * P_impedance - V_r * Q_impedance)

For constant current it is obtained:
Ii_re = (1 / V0) * ( (V_r * P_current + V_i * Q_current) / V )
Ii_im = (1 / V0) * ( (V_i * P_current - V_r * Q_current) / V )

For constant power it is obtained:
Ip_re =  (V_r * P_power + V_i * Q_power) / V^2
Ip_im =  (V_i * P_power - V_r * Q_power) / V^2

"""
function mdl_zip_load!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    wrapper::ZIPLoadWrapper,
) where {T <: ACCEPTED_REAL_TYPES}
    # Read power flow voltages
    #V0_mag_inv = 1.0 / get_V_ref(wrapper)
    V0_mag_inv = 1.0 / PSY.get_magnitude(PSY.get_bus(wrapper))
    V0_mag_sq_inv = V0_mag_inv^2

    V_mag_inv = 1.0 / sqrt(voltage_r^2 + voltage_i^2)
    V_mag_sq_inv = V_mag_inv^2

    #Load device parameters
    P_power = get_P_power(wrapper)
    P_current = get_P_current(wrapper)
    P_impedance = get_P_impedance(wrapper)
    Q_power = get_Q_power(wrapper)
    Q_current = get_Q_current(wrapper)
    Q_impedance = get_Q_impedance(wrapper)

    Iz_re = V0_mag_sq_inv * (voltage_r * P_impedance + voltage_i * Q_impedance)
    Iz_im = V0_mag_sq_inv * (voltage_i * P_impedance - voltage_r * Q_impedance)

    Ii_re = V0_mag_inv * V_mag_inv * (voltage_r * P_current + voltage_i * Q_current)
    Ii_im = V0_mag_inv * V_mag_inv * (voltage_i * P_current - voltage_r * Q_current)

    Ip_re = V_mag_sq_inv * (voltage_r * P_power + voltage_i * Q_power)
    Ip_im = V_mag_sq_inv * (voltage_i * P_power - voltage_r * Q_power)

    #Update current
    current_r[1] += -(Iz_re + Ii_re + Ip_re) #in system pu flowing out
    current_i[1] += -(Iz_im + Ii_im + Ip_im) #in system pu flowing out

    return
end
