"""
Model for ZIP Load model given by:

P_zip = P_power + P_current * (V / V0) + P_impedance * (V / V0)^2
Q_zip = Q_power + Q_current * (V / V0) + Q_impedance * (V / V0)^2

with V = sqrt(V_r^2 + V_i^2) and V0 the voltage magnitude from the power flow solution

The current taken for the load is computed as:
I_zip = (P_zip + j Q_zip)^* / (V_r + j V_i)^*
I_zip = (P_zip - j Q_zip) / (V_r - j V_i)

For constant impedance it is obtained:
Iz_re = (1 / V0)^2 * (V_r * P_impedance + V_i * Q_impedance)
Iz_im = (1 / V0)^2 * (V_i * P_impedance - V_r * Q_impedance)

For constant current it is obtained:
Ii_re = (1 / V0) * ( (V_r * P_current + V_i * Q_current) / V )
Ii_im = (1 / V0) * ( (V_i * P_current - V_r * Q_current) / V )

For constant power it is obtained:
Ip_re =  (V_r * P_power + V_i * Q_power) / V^2
Ip_im =  (V_i * P_power - V_r * Q_power) / V^2

Model for Exponential Load model given by:

P_exp = P0 * (V / V0)^α
Q_exp = Q0 * (V / V0)^β

The current taken for the load is computed as:
I_exp = (P_exp + j Q_exp)^* / (V_r + j V_i)^*
I_exp = (P_exp - j Q_exp) / (V_r - j V_i)

It results:
Ir_exp = V_r * P0 * (V^(α - 2) / V0^α) + V_i * Q0 * (V^(β - 2)/ V0^β)
Ii_im  = V_i * P0 * (V^(α - 2) / V0^α) - V_r * Q0 * (V^(β - 2)/ V0^β)

"""
function mdl_zip_load!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    wrapper::StaticLoadWrapper,
) where {T <: ACCEPTED_REAL_TYPES}
    # Read power flow voltages
    #V0_mag_inv = 1.0 / get_V_ref(wrapper)
    V0_mag_inv = 1.0 / PSY.get_magnitude(PSY.get_bus(wrapper))
    V0_mag_sq_inv = V0_mag_inv^2

    V_mag = hypot(voltage_r, voltage_i)
    V_mag_inv = 1.0 / V_mag
    V_mag_sq_inv = V_mag_inv^2

    # Load device parameters
    P_power = get_P_power(wrapper)
    P_current = get_P_current(wrapper)
    P_impedance = get_P_impedance(wrapper)
    Q_power = get_Q_power(wrapper)
    Q_current = get_Q_current(wrapper)
    Q_impedance = get_Q_impedance(wrapper)
    exp_params = get_exp_params(wrapper)
    Ir_exp = zero(T)
    Ii_exp = zero(T)

    # Compute ZIP currents
    Iz_re = V0_mag_sq_inv * (voltage_r * P_impedance + voltage_i * Q_impedance)
    Iz_im = V0_mag_sq_inv * (voltage_i * P_impedance - voltage_r * Q_impedance)

    Ii_re = V0_mag_inv * V_mag_inv * (voltage_r * P_current + voltage_i * Q_current)
    Ii_im = V0_mag_inv * V_mag_inv * (voltage_i * P_current - voltage_r * Q_current)

    Ip_re = V_mag_sq_inv * (voltage_r * P_power + voltage_i * Q_power)
    Ip_im = V_mag_sq_inv * (voltage_i * P_power - voltage_r * Q_power)

    # Compute Exponential Currents
    if !isempty(exp_params)
        for tuple in exp_params
            P0 = tuple.P_exp
            α = tuple.P_coeff
            Q0 = tuple.Q_exp
            β = tuple.Q_coeff

            # Compute Auxiliary terms
            V_coeff_active = V_mag^(α - 2.0) * V0_mag_inv^α
            V_coeff_reactive = V_mag^(β - 2.0) * V0_mag_inv^β

            # Compute currents
            Ir_exp += voltage_r * P0 * V_coeff_active + voltage_i * Q0 * V_coeff_reactive
            Ii_exp += voltage_i * P0 * V_coeff_active - voltage_r * Q0 * V_coeff_reactive
        end
    end

    # Update current
    current_r[1] += -(Iz_re + Ii_re + Ip_re + Ir_exp) #in system pu flowing out
    current_i[1] += -(Iz_im + Ii_im + Ip_im + Ii_exp) #in system pu flowing out

    return
end
