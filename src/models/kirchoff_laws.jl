function Ybus_current_kirchoff(Ybus, V_r, V_i, I_injections_r, I_injections_i, I_bus, I_balance)
    # Note: BLAS doesn't work because the the type of Vr and Vi is not Matrix of Complex
    LinearAlgebra.mul!(I_bus, Ybus, (V_r + V_i .* 1im))
    bus_count = get_bus_count(sys)
    for n in eachindex(I_bus)
        I_balance[n] = real(I_bus[n]) - I_injections_r[n]
        I_balance[n + bus_count] = imag(I_bus[n]) - I_injections_i[n]
    end

    return
end

function voltage_kirchoff(voltage_buses, V_r, V_i, dx)
    isempty(voltage_buses) && return

    sys_f = PSY.get_frequency(sys)
    ω_b = 2.0 * π * sys_f
    n_buses = length(PSY.get_components(PSY.Bus, sys))
    shunts = get_total_shunts(sys)
    for bus_ix in voltage_buses
        shunt_multiplier = shunts[bus_ix]
        I_balance[bus_ix] =
            -ω_b * I_balance[bus_ix] * shunt_multiplier + ω_b * V_i[bus_ix] - dx[bus_ix]
        I_balance[bus_ix + n_buses] =
            -ω_b * I_balance[bus_ix + n_buses] * shunt_multiplier - ω_b * V_r[bus_ix] -
            dx[bus_ix + n_buses]
    end
end

"""
    Kirchoff law for buses.
    I_gen_r[j]: Real current injection from all generators connected at bus j.
                It is zero if no generator is connected at bus j.
    I_gen_i[j]: Real current injection from all generators connected at bus j.
                It is zero if no generator is connected at bus j.
"""
function kirchoff_laws!(inputs::SimulationInputs, V_r, V_i, I_injections_r, I_injections_i, dx)
    I_bus = get_aux_arrays(inputs)[5]
    I_balance = get_aux_arrays(inputs)[6]
    Ybus = get_Ybus(inputs)
    Ybus_current_kirchoff(Ybus, V_r, V_i, I_injections_r, I_injections_i, I_bus, I_balance)
    voltage_kirchoff(get_voltage_buses_ix(inputs), V_r, V_i, dx)
    return
end
