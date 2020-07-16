function Ybus_current_kirchoff(sys, V_r, V_i, I_injections_r, I_injections_i)
    Ybus = PSY.get_ext(sys)[YBUS]
    # Note: BLAS doesn't work because the the type of Vr and Vi is not Matrix of Complex
    I_bus = PSY.get_ext(sys)[AUX_ARRAYS][5]
    I_balance = PSY.get_ext(sys)[AUX_ARRAYS][6]
    LinearAlgebra.mul!(I_bus, Ybus, (V_r + V_i .* 1im))
    bus_count = get_bus_count(sys)
    for n in 1:bus_count
        I_balance[n] = real(I_bus[n]) - I_injections_r[n]
        I_balance[n + bus_count] = imag(I_bus[n]) - I_injections_i[n]
    end

    return
end

function voltage_kirchoff(sys, V_r, V_i, dx)
    voltage_buses = get_voltage_bus_no(sys)
    isempty(voltage_buses) && return
    I_bus = PSY.get_ext(sys)[AUX_ARRAYS][5]
    I_balance = PSY.get_ext(sys)[AUX_ARRAYS][6]
    sys_f = PSY.get_frequency(sys)
    ω_b = 2.0 * π * sys_f
    n_buses = length(PSY.get_components(PSY.Bus, sys))
    shunts = get_total_shunts(sys)
    for bus_no in voltage_buses
        shunt_multiplier = shunts[bus_no]
        I_balance[bus_no] =
            -ω_b * I_balance[bus_no] * shunt_multiplier + ω_b * V_i[bus_no] - dx[bus_no]
        I_balance[bus_no + n_buses] =
            -ω_b * I_balance[bus_no + n_buses] * shunt_multiplier - ω_b * V_r[bus_no] -
            dx[bus_no + n_buses]
    end
end

"""
    Kirchoff law for buses.
    I_gen_r[j]: Real current injection from all generators connected at bus j.
                It is zero if no generator is connected at bus j.
    I_gen_i[j]: Real current injection from all generators connected at bus j.
                It is zero if no generator is connected at bus j.
"""
function kirchoff_laws!(sys, V_r, V_i, I_injections_r, I_injections_i, dx)
    Ybus_current_kirchoff(sys, V_r, V_i, I_injections_r, I_injections_i)
    voltage_kirchoff(sys, V_r, V_i, dx)
    return
end
