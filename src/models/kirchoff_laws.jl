"""
    Kirchoff law for buses.
    I_gen_r[j]: Real current injection from all generators connected at bus j.
                It is zero if no generator is connected at bus j.
    I_gen_i[j]: Real current injection from all generators connected at bus j.
                It is zero if no generator is connected at bus j.

"""
#TODO: Improve performance of this function
function kirchoff_laws(sys, V_r, V_i, I_injections_r, I_injections_i, dx)
    Ybus = PSY.get_ext(sys)[YBUS]
    I_bus = Ybus * (V_r + V_i .* 1im)
    I_balance = [real(I_bus) - I_injections_r; imag(I_bus) - I_injections_i]

    voltage_buses = PSY.get_ext(sys)["voltage_buses_ix"]
    isempty(voltage_buses) && return I_balance

    sys_f = PSY.get_frequency(sys)
    ω_b = 2.0 * π * sys_f
    n_buses = length(PSY.get_components(PSY.Bus, sys))

    for bus_no in voltage_buses
        shunt_multiplier = PSY.get_ext(sys)["total_shunts"][bus_no]
        I_balance[bus_no] = -ω_b*I_balance[bus_no]*shunt_multiplier + ω_b*V_i[bus_no] - dx[bus_no]
        I_balance[bus_no + n_buses] = -ω_b*I_balance[bus_no+n_buses]*shunt_multiplier - ω_b*V_r[bus_no] - dx[bus_no+n_buses]
    end

    return I_balance
end
#=
            ix_dx = [
                from_bus_number,
                from_bus_number + bus_size,
                to_bus_number,
                to_bus_number + bus_size,
            ]
            =#
